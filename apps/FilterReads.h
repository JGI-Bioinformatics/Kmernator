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
	static int getMaxKmerDepth() {
		return getVarMap()["max-kmer-output-depth"].as<int> ();
	}
	static int getPartitionByDepth() {
		return getVarMap()["partition-by-depth"].as<int> ();
	}
	static bool getBothPairs() {
		return getVarMap()["min-passing-in-pair"].as<int>() == 2;
	}
	void _resetDefaults() {
		KmerOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		// set options specific to this program
		p.add("kmer-size", 1);
		p.add("input-file", -1);
		po::options_description opts("General Filtering Options");
		opts.add_options()

		("max-kmer-output-depth", po::value<int>()->default_value(-1),
				"maximum number of times a kmer will be output among the selected reads (mutually exclusive with partition-by-depth).  This is not a criteria on the kmer spectrum, just a way to reduce the redundancy of the output")

		("partition-by-depth", po::value<int>()->default_value(-1),
				"partition filtered reads by powers-of-two coverage depth (mutually exclusive with max-kmer-depth)")

		("min-passing-in-pair", po::value<int>()->default_value(1),
				"1 or 2 reads in a pair must pass filters");

		desc.add(opts);
		KmerOptions::_setOptions(desc, p);
	}
	bool _parseOptions( po::variables_map &vm) {

		bool ret = true;

		// verify mutually exclusive options are not set
		if ( (getMaxKmerDepth() > 0 && getPartitionByDepth() >  0) )
		{
			throw std::invalid_argument("You can not specify both max-kmer-depth and partition-by-depth");
		}
		if (Options::getOptions().getOutputFile().empty() && Logger::isMaster())
		{
			LOG_WARN(1, "no output file specified... This is a dry run!");
		}

		if (Options::getOptions().getInputFiles().empty() && Logger::isMaster()) {
			LOG_ERROR(1, "Please specify at least one input file");
			ret = false;
		}

		ret &= KmerOptions::_parseOptions(vm);
		return ret ;
	}
};
typedef OptionsBaseTemplate< _FilterReadsBaseOptions > FilterReadsBaseOptions;

template<typename _ReadSelector>
long selectReads(unsigned int minDepth, ReadSet &reads, _ReadSelector &selector, std::string outputFilename)
{
	typedef typename _ReadSelector::OFM OFM;
	LOG_VERBOSE_OPTIONAL(1, true, "selectReads with minDepth " << minDepth << ", minLength " << Options::getOptions().getMinReadLength() << ": " << reads.getSize() << " reads");
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
						minDepth, Options::getOptions().getMinReadLength(), FilterReadsBaseOptions::getOptions().getBothPairs());
			else
				picked += selector.pickBestCoveringSubsetReads(depth,
						minDepth, Options::getOptions().getMinReadLength());
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
		if (maxDepth < 0) {
			maxDepth = 1;
		}

		for (unsigned int depth = maxDepth; depth >= 1; depth /= 2) {

			string ofname = outputFilename;
			if (maxDepth > 1) {
				ofname += "-PartitionDepth" + boost::lexical_cast< string >( depth );
			}
			OFM ofmap = selector.getOFM(ofname);
			float tmpMinDepth = std::max(minDepth, depth);
			if (KmerOptions::getOptions().getKmerSize() == 0) {
				tmpMinDepth = 0;
				depth = 0;
			}
			LOG_VERBOSE(1, "Selecting reads over depth: " << depth << " (" << tmpMinDepth << ") ");

			if (reads.hasPairs()) {
				picked = selector.pickAllPassingPairs(tmpMinDepth,
						Options::getOptions().getMinReadLength(),
						FilterReadsBaseOptions::getOptions().getBothPairs());
			} else {
				picked = selector.pickAllPassingReads(tmpMinDepth,
						Options::getOptions().getMinReadLength());
			}
			LOG_VERBOSE(2, "At or above coverage: " << depth << " Picked " << picked
			<< " / " << reads.getSize() << " reads");
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			if (picked > 0 && !outputFilename.empty()) {
				LOG_VERBOSE(1, "Writing " << picked << " reads  to output files");
				selector.writePicks(ofmap, oldPicked);
			}
			oldPicked += picked;

			if (minDepth > depth) {
				break;
			}

		}
	}
	LOG_VERBOSE(1, "Done.  Cleaning up. " << MemoryUtils::getMemoryUsage());

	return oldPicked;
};



#endif /* FILTERREADS_H_ */
