//
// Kmernator/apps/FilterReads.h
//
// Author: Rob Egan
//
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/

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
		po::options_description opts("FilterReads <options> [[kmer-size] [input-file ...]]\n\tNote: --kmer-size and --input-file can either be specified as positional argumens at the end or within <options>\n\nGeneral Filtering Options");
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

	string suffix;
	if (ReadSelectorOptions::getOptions().getSeparateOutputs()) {
		if (KmerBaseOptions::getOptions().getKmerSize() > 0) {
			outputFilename += "-MinDepth" + boost::lexical_cast<std::string>(minDepth);
		}
		suffix = FormatOutput::getDefaultSuffix();
	}

	if (maximumKmerDepth > 0) {
		if (ReadSelectorOptions::getOptions().getSeparateOutputs())
			outputFilename += "-MaxDepth" + boost::lexical_cast<std::string>(maximumKmerDepth);
		OFM ofmap = selector.getOFM(outputFilename, suffix);
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
		float oldMinReadLength = ReadSelectorOptions::getOptions().getMinReadLength();
		int oldMinPassingInPair = ReadSelectorOptions::getOptions().getMinPassingInPair();
		bool hasRemainderTrim = false;

		for (unsigned int depth = maxDepth; depth >= minDepth; depth /= 2) {

			float tmpMinDepth = std::max(minDepth, depth);
			if (KmerBaseOptions::getOptions().getKmerSize() == 0) {
				tmpMinDepth = 0;
				depth = 0;
			}
			string ofname = outputFilename;
			if (hasRemainderTrim && ReadSelectorOptions::getOptions().getSeparateOutputs()) {
				ofname += "-Remainder";
			} else if (isPartitioned && tmpMinDepth > 0 && ReadSelectorOptions::getOptions().getSeparateOutputs()) {
				ofname += "-PartitionDepth" + boost::lexical_cast< string >( tmpMinDepth );
			}
			OFM ofmap = selector.getOFM(ofname, suffix);
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
					&& ReadSelectorOptions::getOptions().getRemainderTrim() > 0.0
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
