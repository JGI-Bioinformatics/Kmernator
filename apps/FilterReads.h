//
// Kmernator/apps/FilterReads.h
//
// Author: Rob Egan
//
// Copyright 2010 The Regents of the University of California.
// All rights reserved.
//
// The United States Government has rights in this work pursuant
// to contracts DE-AC03-76SF00098, W-7405-ENG-36 and/or
// W-7405-ENG-48 between the United States Department of Energy
// and the University of California.
//
// Redistribution and use in source and binary forms are permitted
// provided that: (1) source distributions retain this entire
// copyright notice and comment, and (2) distributions including
// binaries display the following acknowledgement:  "This product
// includes software developed by the University of California,
// JGI-PSF and its contributors" in the documentation or other
// materials provided with the distribution and in all advertising
// materials mentioning features or use of this software.  Neither the
// name of the University nor the names of its contributors may be
// used to endorse or promote products derived from this software
// without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE.
//

#ifndef FILTER_READS_H_
#define FILTER_READS_H_


#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "Options.h"
#include "ReadSet.h"
#include "FilterKnownOddities.h"
#include "DuplicateFragmentFilter.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "ReadSelector.h"
#include "KmerTrackingData.h"
#include "Utils.h"
#include "Log.h"

#include <boost/lexical_cast.hpp>

using namespace std;

// TODO add outputformat of fasta
class _FilterReadsBaseOptions : public OptionsBaseInterface {
public:
	_FilterReadsBaseOptions() :  sizeHistoryFile("") {}
	virtual ~_FilterReadsBaseOptions() {}

	std::string &getSizeHistoryFile() {
		return sizeHistoryFile;
	}
	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		KmerBaseOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();
		ReadSelectorOptions::_resetDefaults();
		FilterKnownOdditiesOptions::_resetDefaults();
		DuplicateFragmentFilterOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		// set options specific to this program
		p.add("kmer-size", 1);
		p.add("input-file", -1);
		po::options_description opts("General Filtering Options");
		opts.add_options()

				("size-history-file", po::value<std::string>()->default_value(sizeHistoryFile), "if set, a text file with accumulated kmer counts will be generated (for EstimateSize.R)");

		desc.add(opts);
		GeneralOptions::_setOptions(desc, p);
		KmerBaseOptions::_setOptions(desc, p);
		KmerSpectrumOptions::_setOptions(desc, p);
		ReadSelectorOptions::_setOptions(desc, p);
		FilterKnownOdditiesOptions::_setOptions(desc,p);
		DuplicateFragmentFilterOptions::_setOptions(desc, p);
	}
	bool _parseOptions( po::variables_map &vm) {

		bool ret = true;

		ret &= GeneralOptions::_parseOptions(vm);
		ret &= KmerBaseOptions::_parseOptions(vm);
		ret &= KmerSpectrumOptions::_parseOptions(vm);
		ret &= ReadSelectorOptions::_parseOptions(vm);
		ret &= FilterKnownOdditiesOptions::_parseOptions(vm);
		ret &= DuplicateFragmentFilterOptions::_parseOptions(vm);

		setOpt("size-history-file", sizeHistoryFile);

		if (Options::getOptions().getOutputFile().empty() && Logger::isMaster())
		{
			LOG_WARN(1, "no output file specified... This is a dry run!");
		}

		if (Options::getOptions().getInputFiles().empty() && Logger::isMaster()) {
			setOptionsErrorMsg("Please specify at least one input file");
			ret = false;
		}

		return ret ;
	}

protected:
	std::string sizeHistoryFile;
};
typedef OptionsBaseTemplate< _FilterReadsBaseOptions > FilterReadsBaseOptions;

template<typename _ReadSelector>
long selectReads(unsigned int minDepth, ReadSet &reads, _ReadSelector &selector, std::string outputFilename)
{
	typedef typename _ReadSelector::OFM OFM;
	LOG_VERBOSE_OPTIONAL(1, true, "selectReads with minDepth " << minDepth << ", minLength " << ReadSelectorOptions::getOptions().getMinReadLength() << ": " << reads.getSize() << " reads");
	LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());

	long oldPicked = 0;
	long picked = 0;

	int maximumKmerDepth = ReadSelectorOptions::getOptions().getMaxKmerDepth();


	if (maximumKmerDepth > 0) {
		outputFilename += "-MaxDepth" + boost::lexical_cast<std::string>(maximumKmerDepth);
		OFM ofmap = selector.getOFM(outputFilename);
		std::string normalizationMethod = ReadSelectorOptions::getOptions().getNormalizationMethod();
		if (normalizationMethod == "RANDOM") {
			picked += selector.pickCoverageNormalizedSubset(maximumKmerDepth, minDepth, ReadSelectorOptions::getOptions().getMinReadLength(), reads.hasPairs(), ReadSelectorOptions::getOptions().getBothPairs());
		} else if (normalizationMethod == "OPTIMAL") {
			for (int depth = 1; depth <= maximumKmerDepth; depth++) {

				LOG_VERBOSE_OPTIONAL(2, true, "Picking depth " << depth << " layer of reads");
				if (reads.hasPairs())
					picked += selector.pickBestCoveringSubsetPairs(depth,
							minDepth, ReadSelectorOptions::getOptions().getMinReadLength(), ReadSelectorOptions::getOptions().getBothPairs());
				else
					picked += selector.pickBestCoveringSubsetReads(depth,
							minDepth, ReadSelectorOptions::getOptions().getMinReadLength());
				LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());
			}
		} else {
			LOG_WARN(1, "INVALID normalization-method: " << normalizationMethod);
		}

		if (picked > 0 && !outputFilename.empty()) {
			LOG_VERBOSE_OPTIONAL(1, true, "Writing " << picked << " reads to output file(s)");
			selector.writePicks(ofmap, oldPicked);
		}
		LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());
		oldPicked += picked;


	} else {

		int maxDepth = ReadSelectorOptions::getOptions().getPartitionByDepth();
		bool isPartitioned = (maxDepth > 0);
		if (!isPartitioned) {
			maxDepth = minDepth;
		}
		// record potentially modified options
		int oldMinReadLength = ReadSelectorOptions::getOptions().getMinReadLength();
		int oldMinPassingInPair = ReadSelectorOptions::getOptions().getMinPassingInPair();
		bool hasRemainderTrim = false;

		for (unsigned int depth = maxDepth; depth >= minDepth; depth /= 2) {

			float tmpMinDepth = std::max(minDepth, depth);
			if (KmerBaseOptions::getOptions().getKmerSize() == 0) {
				tmpMinDepth = 0;
				depth = 0;
			}
			string ofname = outputFilename;
			if (hasRemainderTrim) {
				ofname += "-Remainder";
			} else if (isPartitioned && tmpMinDepth > 0) {
				ofname += "-PartitionDepth" + boost::lexical_cast< string >( tmpMinDepth );
			}
			OFM ofmap = selector.getOFM(ofname);
			LOG_VERBOSE(1, "Selecting reads over depth: " << tmpMinDepth);

			if (reads.hasPairs()) {
				LOG_DEBUG(3, "getBothPairs: " << ReadSelectorOptions::getOptions().getBothPairs() << " " << ReadSelectorOptions::getOptions().getMinPassingInPair());
				picked = selector.pickAllPassingPairs(tmpMinDepth,
						ReadSelectorOptions::getOptions().getMinReadLength(),
						ReadSelectorOptions::getOptions().getBothPairs());
			} else {
				picked = selector.pickAllPassingReads(tmpMinDepth,
						ReadSelectorOptions::getOptions().getMinReadLength());
			}
			LOG_VERBOSE(2, "At or above coverage: " << tmpMinDepth << " Picked " << picked
					<< " / " << reads.getSize() << " reads");
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			if (!outputFilename.empty()) {
				LOG_VERBOSE(1, "Writing " << picked << " reads to output files");
				selector.writePicks(ofmap, oldPicked);
			}
			oldPicked += picked;

			if (depth == minDepth) {
				if ((!hasRemainderTrim)
					&& isPartitioned
					&& ReadSelectorOptions::getOptions().getRemainderTrim() > 0
					&& (ReadSelectorOptions::getOptions().getMinPassingInPair() != 1
						|| ((int) ReadSelectorOptions::getOptions().getMinReadLength()) != ReadSelectorOptions::getOptions().getRemainderTrim()
						)) {
					// modify options
					ReadSelectorOptions::getOptions().getMinPassingInPair() = 1;
					ReadSelectorOptions::getOptions().getMinReadLength() = ReadSelectorOptions::getOptions().getRemainderTrim();
					hasRemainderTrim=true;
					depth *= 2;
				} else {
					break;
				}
			}
		}
		if (hasRemainderTrim) {
			// reset modified options
			ReadSelectorOptions::getOptions().getMinReadLength() = oldMinReadLength;
			ReadSelectorOptions::getOptions().getMinPassingInPair() = oldMinPassingInPair;
		}
	}
	LOG_VERBOSE(1, "Done.  Cleaning up. " << MemoryUtils::getMemoryUsage());

	return oldPicked;
};



#endif /* FILTERREADS_H_ */
