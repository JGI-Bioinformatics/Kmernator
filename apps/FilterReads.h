//
// Kmernator/apps/FilterReads.h
//
// Author: Rob Egan, Craig Furman
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
	_FilterReadsBaseOptions() : maxKmerDepth(-1), partitionByDepth(-1), bothPairs(1), remainderTrim(-1), sizeHistoryFile("") {}
	virtual ~_FilterReadsBaseOptions() {}

	int &getMaxKmerDepth() {
		return maxKmerDepth;
	}
	int &getPartitionByDepth() {
		return partitionByDepth;
	}
	int &getMinPassingInPair() {
		return bothPairs;
	}
	int &getRemainderTrim() {
		return remainderTrim;
	}
	bool getBothPairs() {
		return bothPairs == 2;
	}
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

				("max-kmer-output-depth", po::value<int>()->default_value(maxKmerDepth), "maximum number of times a kmer will be output among the selected reads (mutually exclusive with partition-by-depth).  This is not a criteria on the kmer spectrum, just a way to reduce the redundancy of the output")

				("partition-by-depth", po::value<int>()->default_value(partitionByDepth), "partition filtered reads by powers-of-two coverage depth (mutually exclusive with max-kmer-depth)")

				("min-passing-in-pair", po::value<int>()->default_value(bothPairs), "1 or 2 reads in a pair must pass filters")

				("remainder-trim", po::value<int>()->default_value(remainderTrim), "if set, a final round letting single reads and lesser trimmed reads will be selected")

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

		setOpt<int>("max-kmer-output-depth", maxKmerDepth);
		setOpt<int>("partition-by-depth", partitionByDepth);
		setOpt<int>("min-passing-in-pair", bothPairs);
		setOpt<int>("remainder-trim", remainderTrim);
		setOpt<std::string>("size-history-file", sizeHistoryFile);

		// verify mutually exclusive options are not set
		if ( (getMaxKmerDepth() > 0 && getPartitionByDepth() >  0) )
		{
			LOG_ERROR(1, "You can not specify both max-kmer-depth and partition-by-depth");
			ret = false;
		}
		if (Options::getOptions().getOutputFile().empty() && Logger::isMaster())
		{
			LOG_WARN(1, "no output file specified... This is a dry run!");
		}

		if (Options::getOptions().getInputFiles().empty() && Logger::isMaster()) {
			LOG_ERROR(1, "Please specify at least one input file");
			ret = false;
		}

		return ret ;
	}

protected:
	int maxKmerDepth, partitionByDepth, bothPairs, remainderTrim;
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

	int maximumKmerDepth = FilterReadsBaseOptions::getOptions().getMaxKmerDepth();


	if (maximumKmerDepth > 0) {
		OFM ofmap = selector.getOFM(outputFilename);
		for (int depth = 1; depth <= maximumKmerDepth; depth++) {
			LOG_VERBOSE_OPTIONAL(2, true, "Picking depth " << depth << " layer of reads");
			if (reads.hasPairs())
				picked += selector.pickBestCoveringSubsetPairs(depth,
						minDepth, ReadSelectorOptions::getOptions().getMinReadLength(), FilterReadsBaseOptions::getOptions().getBothPairs());
			else
				picked += selector.pickBestCoveringSubsetReads(depth,
						minDepth, ReadSelectorOptions::getOptions().getMinReadLength());
			LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());
		}

		if (picked > 0 && !outputFilename.empty()) {
			LOG_VERBOSE_OPTIONAL(1, true, "Writing " << picked << " reads to output file(s)");
			selector.writePicks(ofmap, oldPicked);
		}
		LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());
		oldPicked += picked;


	} else {

		int maxDepth = FilterReadsBaseOptions::getOptions().getPartitionByDepth();
		bool isPartitioned = (maxDepth > 0);
		if (!isPartitioned) {
			maxDepth = minDepth;
		}
		// record potentially modified options
		int oldMinReadLength = ReadSelectorOptions::getOptions().getMinReadLength();
		int oldMinPassingInPair = FilterReadsBaseOptions::getOptions().getMinPassingInPair();
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
				LOG_DEBUG(3, "getBothPairs: " << FilterReadsBaseOptions::getOptions().getBothPairs() << " " << FilterReadsBaseOptions::getOptions().getMinPassingInPair());
				picked = selector.pickAllPassingPairs(tmpMinDepth,
						ReadSelectorOptions::getOptions().getMinReadLength(),
						FilterReadsBaseOptions::getOptions().getBothPairs());
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
					&& FilterReadsBaseOptions::getOptions().getRemainderTrim() > 0
					&& (FilterReadsBaseOptions::getOptions().getMinPassingInPair() != 1
						|| ((int) ReadSelectorOptions::getOptions().getMinReadLength()) != FilterReadsBaseOptions::getOptions().getRemainderTrim()
						)) {
					// modify options
					FilterReadsBaseOptions::getOptions().getMinPassingInPair() = 1;
					ReadSelectorOptions::getOptions().getMinReadLength() = FilterReadsBaseOptions::getOptions().getRemainderTrim();
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
			FilterReadsBaseOptions::getOptions().getMinPassingInPair() = oldMinPassingInPair;
		}
	}
	LOG_VERBOSE(1, "Done.  Cleaning up. " << MemoryUtils::getMemoryUsage());

	return oldPicked;
};



#endif /* FILTERREADS_H_ */
