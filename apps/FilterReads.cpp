//
// Kmernator/apps/FilterReads.cpp
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

#include "FilterReads.h"

typedef TrackingDataWithDirection DataType;
typedef KmerSpectrum<DataType, DataType> KS;
typedef ReadSelector<DataType> RS;
class _FilterReadsOptions : public _FilterReadsBaseOptions {
public:
	void _resetDefaults() {
		_FilterReadsBaseOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
		FilterKnownOdditiesOptions::_resetDefaults();
		DuplicateFragmentFilterOptions::_resetDefaults();

		KmerOptions::getOptions().getSaveKmerMmap() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		_FilterReadsBaseOptions::_setOptions(desc, p);
		GeneralOptions::_setOptions(desc, p);
		KmerOptions::_setOptions(desc, p);
		FilterKnownOdditiesOptions::_setOptions(desc, p);
		DuplicateFragmentFilterOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= KmerOptions::_parseOptions(vm);
		ret &= FilterKnownOdditiesOptions::_parseOptions(vm);
		ret &= DuplicateFragmentFilterOptions::_parseOptions(vm);

		ret &= _FilterReadsBaseOptions::_parseOptions(vm);
		return ret;
	}
};
typedef OptionsBaseTemplate< _FilterReadsOptions > FilterReadsOptions;

int main(int argc, char *argv[]) {

	try {
		if (!FilterReadsOptions::parseOpts(argc, argv))
			throw invalid_argument("Please fix the command line arguments");
	} catch (...) {
		std::cerr << FilterReadsOptions::getDesc() << std::endl << std::endl;
		std::cerr << "Please fix the command line arguments and/or OpenMP environment" << std::endl;
		exit(1);
	}

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	ReadSet reads;

	try {
		OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();
		LOG_VERBOSE(1, "Reading Input Files");
		reads.appendAllFiles(inputs);
		LOG_VERBOSE(1, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
				<< " Bases ");
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		LOG_VERBOSE(1, "Identifying Pairs: ");
		long numPairs = reads.identifyPairs();
		LOG_VERBOSE(1, "Pairs + single = " << numPairs);
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		if (FilterKnownOdditiesOptions::getOptions().getSkipArtifactFilter() == 0) {

			LOG_VERBOSE(1, "Preparing artifact filter: ");
			FilterKnownOddities filter;
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			LOG_VERBOSE(2, "Applying sequence artifact filter to Input Files");
			unsigned long filtered = filter.applyFilter(reads);
			LOG_VERBOSE(1, "filter affected (trimmed/removed) " << filtered << " Reads ");;
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		}
		if (DuplicateFragmentFilterOptions::getOptions().getDeDupMode() > 0 && DuplicateFragmentFilterOptions::getOptions().getDeDupEditDistance() >= 0) {
			LOG_VERBOSE(2, "Applying DuplicateFragmentPair Filter to Input Files");
			unsigned long duplicateFragments = DuplicateFragmentFilter::filterDuplicateFragments(reads);
			LOG_VERBOSE(1, "filter removed duplicate fragment pair reads: " << duplicateFragments);
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
		}

		KS spectrum(0);

		Kmernator::MmapFileVector spectrumMmaps;
		if (KmerOptions::getOptions().getKmerSize() > 0 && !KmerOptions::getOptions().getLoadKmerMmap().empty()) {
			spectrum.restoreMmap(KmerOptions::getOptions().getLoadKmerMmap());
		} else if (KmerOptions::getOptions().getKmerSize() > 0) {

			long numBuckets = KS::estimateWeakKmerBucketSize(reads);
			LOG_DEBUG(1, "targeting " << numBuckets << " buckets for reads ");

			spectrum = KS(numBuckets);
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			spectrumMmaps = spectrum.buildKmerSpectrumInParts(reads, KmerOptions::getOptions().getBuildPartitions(), outputFilename.empty() ? "" : outputFilename + "-mmap");
			spectrum.optimize();
			spectrum.trackSpectrum(true);
			std::string sizeHistoryFile = FilterReadsOptions::getOptions().getSizeHistoryFile();
			if (!sizeHistoryFile.empty()) {
				LOG_VERBOSE(1, "Writing size history file to: " << sizeHistoryFile);
				OfstreamMap ofm(sizeHistoryFile, "");
				ofm.getOfstream("") << spectrum.getSizeTracker().toString();
			} else {
				LOG_VERBOSE(1, "Kmer Size History:" << std::endl << spectrum.getSizeTracker().toString());
			}

			if (Log::isVerbose(1))
				spectrum.printHistograms(Log::Verbose("Kmer Histogram"));

			if (Options::getOptions().getVariantSigmas() > 0.0) {
				spectrum.purgeVariants();
				if (Log::isVerbose(1)) {
					spectrum.printHistograms(Log::Verbose("Variant-Removed Kmer Histogram"));
				}
			}
		}

		if (KmerOptions::getOptions().getKmerSize() > 0) {
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			if (Options::getOptions().getGCHeatMap() && ! outputFilename.empty()) {
				LOG_VERBOSE(1, "Creating GC Heat Map ");
				LOG_DEBUG(1,  MemoryUtils::getMemoryUsage());
				OfstreamMap ofmap(outputFilename + "-GC", ".txt");
				spectrum.printGC(ofmap.getOfstream(""));
			}

			if (KmerOptions::getOptions().getMinDepth() > 1) {
				LOG_DEBUG(1, "Clearing singletons from memory");
				spectrum.singleton.clear();
				LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
			} else {
				spectrum.optimize(true);
			}
		}


		unsigned int minDepth = KmerOptions::getOptions().getMinDepth();
		unsigned int depthRange = Options::getOptions().getDepthRange();
		unsigned int depthStep = 2;
		if (depthRange < minDepth) {
			depthRange = minDepth;
		}

		if (!outputFilename.empty()) {

			for(unsigned int thisDepth = depthRange ; thisDepth >= minDepth; thisDepth /= depthStep) {
				std::string pickOutputFilename = outputFilename;
				if (KmerOptions::getOptions().getKmerSize() > 0) {
					pickOutputFilename += "-MinDepth" + boost::lexical_cast<std::string>(thisDepth);
					LOG_VERBOSE(1, "Trimming reads with minDepth: " << thisDepth);
				} else {
					LOG_VERBOSE(1, "Trimming reads that pass Artifact Filter with length: " << Options::getOptions().getMinReadLength());
				}

				RS selector(reads, spectrum.weak);
				selector.scoreAndTrimReads(minDepth);

				selectReads(thisDepth, reads, selector, pickOutputFilename);
			}
		}
		LOG_DEBUG(1, "Clearing spectrum");
		spectrum.reset();

	} catch (...) {
		LOG_ERROR(1, "caught an error!" << StackTrace::getStackTrace());
	}


	LOG_VERBOSE(1, "Finished");

	return 0;
}


