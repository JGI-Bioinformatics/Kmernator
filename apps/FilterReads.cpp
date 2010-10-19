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

typedef TrackingDataMinimal4f DataType;
typedef KmerSpectrum<DataType, DataType> KS;
typedef ReadSelector<DataType> RS;

#include "FilterReadsTemplate.h"


int main(int argc, char *argv[]) {

	if (!FilterReadsOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	MemoryUtils::getMemoryUsage();
    std::string outputFilename = Options::getOutputFile();

	ReadSet reads;
	KmerSizer::set(Options::getKmerSize());

	Options::FileListType inputs = Options::getInputFiles();
	LOG_VERBOSE(1, "Reading Input Files");
	reads.appendAllFiles(inputs);
	LOG_VERBOSE(1, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
			<< " Bases ");
	LOG_VERBOSE(1, MemoryUtils::getMemoryUsage());

	LOG_VERBOSE(1, "Identifying Pairs: ");
	long numPairs = reads.identifyPairs();
	LOG_VERBOSE(1, "Pairs + single = " << numPairs);
	LOG_VERBOSE(1, MemoryUtils::getMemoryUsage());

	if (Options::getSkipArtifactFilter() == 0) {

	  LOG_VERBOSE(1, "Preparing artifact filter: ");
      FilterKnownOddities filter;
      LOG_VERBOSE(1, MemoryUtils::getMemoryUsage());

	  LOG_VERBOSE(1, "Applying sequence artifact filter to Input Files");
	  unsigned long filtered = filter.applyFilter(reads);
	  LOG_VERBOSE(1, "filter affected (trimmed/removed) " << filtered << " Reads ");;
	  LOG_VERBOSE(1, MemoryUtils::getMemoryUsage());

	  LOG_VERBOSE(1, "Applying DuplicateFragmentPair Filter to Input Files");
	  unsigned long duplicateFragments = filter.filterDuplicateFragments(reads);
	  LOG_VERBOSE(1, "filter affected  (removed) " << duplicateFragments);
	  LOG_VERBOSE(1, MemoryUtils::getMemoryUsage());
	}

	KS spectrum(0);

	Kmernator::MmapFileVector spectrumMmaps;

	if (Options::getKmerSize() > 0) {

	  long numBuckets = KS::estimateWeakKmerBucketSize(reads, 64);
	  LOG_VERBOSE(1, "targeting " << numBuckets << " buckets for reads ");

	  spectrum = KS(numBuckets);
	  LOG_VERBOSE(1, MemoryUtils::getMemoryUsage());

	  TrackingData::minimumWeight = Options::getMinKmerQuality();

	  spectrumMmaps = spectrum.buildKmerSpectrumInParts(reads, Options::getBuildPartitions());
	  LOG_VERBOSE(1, MemoryUtils::getMemoryUsage());

	  if (Options::getGCHeatMap() && ! outputFilename.empty()) {
		  LOG_VERBOSE(1, "Creating GC Heat Map " <<  MemoryUtils::getMemoryUsage());
		  OfstreamMap ofmap(outputFilename + "-GC", ".txt");
		  spectrum.printGC(ofmap.getOfstream(""));
	  }

	  if (Options::getMinDepth() > 1) {
        LOG_VERBOSE(1, "Clearing singletons from memory");
        spectrum.singleton.clear();
	    LOG_VERBOSE(1, MemoryUtils::getMemoryUsage());
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
			}
			LOG_VERBOSE(1, "Trimming reads with minDepth: " << thisDepth);
			RS selector(reads, spectrum.weak);
			selector.scoreAndTrimReads(minDepth);

			selectReads(thisDepth, reads, spectrum, selector, pickOutputFilename);
		}
	}

	spectrum.reset();

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
