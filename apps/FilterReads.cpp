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

int main(int argc, char *argv[]) {

	Options::getSaveKmerMmap() = 1;
	try {
		if (!FilterReadsOptions::parseOpts(argc, argv))
			throw invalid_argument("Please fix the command line arguments");
	} catch (...) {
		std::cerr << FilterReadsOptions::getDesc() << std::endl << std::endl;
		std::cerr << "Please fix the command line arguments and/or OpenMP environment" << std::endl;
		exit(1);
	}

	MemoryUtils::getMemoryUsage();
    std::string outputFilename = Options::getOutputFile();

	ReadSet reads;
	KmerSizer::set(Options::getKmerSize());

	Options::FileListType inputs = Options::getInputFiles();
	LOG_VERBOSE(1, "Reading Input Files");
	reads.appendAllFiles(inputs);
	LOG_VERBOSE(1, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
			<< " Bases ");
	LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

	LOG_VERBOSE(1, "Identifying Pairs: ");
	long numPairs = reads.identifyPairs();
	LOG_VERBOSE(1, "Pairs + single = " << numPairs);
	LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

	if (Options::getSkipArtifactFilter() == 0) {

	  LOG_VERBOSE(1, "Preparing artifact filter: ");
      FilterKnownOddities filter;
      LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

	  LOG_VERBOSE(2, "Applying sequence artifact filter to Input Files");
	  unsigned long filtered = filter.applyFilter(reads);
	  LOG_VERBOSE(1, "filter affected (trimmed/removed) " << filtered << " Reads ");;
	  LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

	}
	if (Options::getDeDupMode() > 0 && Options::getDeDupEditDistance() >= 0) {
	  LOG_VERBOSE(2, "Applying DuplicateFragmentPair Filter to Input Files");
	  unsigned long duplicateFragments = DuplicateFragmentFilter::filterDuplicateFragments(reads);
	  LOG_VERBOSE(1, "filter removed duplicate fragment pair reads: " << duplicateFragments);
	  LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
	}

	KS spectrum(0);

	Kmernator::MmapFileVector spectrumMmaps;

	if (Options::getKmerSize() > 0) {

	  long numBuckets = KS::estimateWeakKmerBucketSize(reads, 64);
	  LOG_DEBUG(1, "targeting " << numBuckets << " buckets for reads ");

	  spectrum = KS(numBuckets);
	  LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

	  TrackingData::minimumWeight = Options::getMinKmerQuality();

	  spectrumMmaps = spectrum.buildKmerSpectrumInParts(reads, Options::getBuildPartitions());
	  if (Options::getVariantSigmas() > 0.0) {
		  spectrum.purgeVariants();
		  if (Log::isVerbose(1)) {
			  spectrum.printHistograms(Log::Verbose("Variant-Removed Histogram"));
		  }
	  }

	  LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

	  if (Options::getGCHeatMap() && ! outputFilename.empty()) {
		  LOG_VERBOSE(1, "Creating GC Heat Map ");
		  LOG_DEBUG(1,  MemoryUtils::getMemoryUsage());
		  OfstreamMap ofmap(outputFilename + "-GC", ".txt");
		  spectrum.printGC(ofmap.getOfstream(""));
	  }

	  if (Options::getMinDepth() > 1) {
		  LOG_DEBUG(1, "Clearing singletons from memory");
		  spectrum.singleton.clear();
		  LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
	  }
	}


	unsigned int minDepth = Options::getMinDepth();
	unsigned int depthRange = Options::getDepthRange();
	unsigned int depthStep = 2;
	if (depthRange < minDepth) {
		depthRange = minDepth;
	}

	if (!outputFilename.empty()) {

		for(unsigned int thisDepth = depthRange ; thisDepth >= minDepth; thisDepth /= depthStep) {
			std::string pickOutputFilename = outputFilename;
			if (Options::getKmerSize() > 0) {
				pickOutputFilename += "-MinDepth" + boost::lexical_cast<std::string>(thisDepth);
				LOG_VERBOSE(1, "Trimming reads with minDepth: " << thisDepth);
			} else {
				LOG_VERBOSE(1, "Trimming reads that pass Artifact Filter with length: " << Options::getMinReadLength());
			}

			RS selector(reads, spectrum.weak);
			selector.scoreAndTrimReads(minDepth);

			selectReads(thisDepth, reads, selector, pickOutputFilename);
		}
	}

	LOG_DEBUG(1, "Clearing spectrum");
	spectrum.reset();

	LOG_VERBOSE(1, "Finished");

	return 0;
}


// $Log: FilterReads.cpp,v $
// Revision 1.22  2010-05-24 21:48:50  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.21.2.3  2010-05-20 18:35:25  regan
// bugfix in output naming
//
// Revision 1.21.2.2  2010-05-20 18:25:58  regan
// fixed to not output unless asked for
//
// Revision 1.21.2.1  2010-05-19 21:36:57  regan
// refactored duplicate fragment filter code
// added duplicate fragment on single ended reads
//
// Revision 1.21  2010-05-18 20:50:18  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.20.2.4  2010-05-18 16:43:49  regan
// added GC heatmap output .. still refining
//
// Revision 1.20.2.3  2010-05-12 20:47:42  regan
// minor refactor.
// adjusted output file names
// support of option to output a range of min-depths
//
// Revision 1.20.2.2  2010-05-12 18:24:25  regan
// bugfix
//
// Revision 1.20.2.1  2010-05-12 17:57:11  regan
// help destructor ordering
//
// Revision 1.20  2010-05-06 21:46:57  regan
// merged changes from PerformanceTuning-20100501
//
//
