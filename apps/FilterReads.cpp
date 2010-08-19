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

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "Options.h"
#include "ReadSet.h"
#include "FilterKnownOddities.h"
#include "KmerSpectrum.h"
#include "ReadSelector.h"
#include "Utils.h"

#include <boost/lexical_cast.hpp>

using namespace std;

typedef TrackingDataMinimal4f DataType;
typedef KmerSpectrum<DataType, DataType> KS;
typedef ReadSelector<DataType> RS;

// TODO add outputformat of fasta
class FilterReadsOptions : Options {
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
	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("kmer-size", 1);
		getPosDesc().add("input-file", -1);

		getDesc().add_options()

		("max-kmer-output-depth", po::value<int>()->default_value(-1),
				"maximum number of times a kmer will be output among the selected reads (mutually exclusive with partition-by-depth).  This is not a criteria on the kmer spectrum, just a way to reduce the redundancy of the output")

		("partition-by-depth", po::value<int>()->default_value(-1),
				"partition filtered reads by powers-of-two coverage depth (mutually exclusive with max-kmer-depth)")

		("min-passing-in-pair", po::value<int>()->default_value(1),
				"1 or 2 reads in a pair must pass filters");

		bool ret = Options::parseOpts(argc, argv);

		if (ret) {
			// verify mutually exclusive options are not set
			if ( (getMaxKmerDepth() > 0 && getPartitionByDepth() >  0) )
			{
				throw std::invalid_argument("You can not specify both max-kmer-depth and partition-by-depth");
			}
			if (Options::getOutputFile().empty())
			{
				std::cerr << "WARNING: no output file specified... This is a dry run!" << std::endl;
			}
		}
		return ret;
	}
};

long selectReads(unsigned int minDepth, ReadSet &reads, KS &spectrum, std::string outputFileName);

int main(int argc, char *argv[]) {
	if (!FilterReadsOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	MemoryUtils::getMemoryUsage();
    std::string outputFilename = Options::getOutputFile();

	ReadSet reads;
	KmerSizer::set(Options::getKmerSize());

	Options::FileListType inputs = Options::getInputFiles();
	cerr << "Reading Input Files" << endl;
	reads.appendAllFiles(inputs);
	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
			<< " Bases " << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	cerr << "Identifying Pairs: ";
	long numPairs = reads.identifyPairs();
	cerr << "Pairs + single = " << numPairs << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	if (Options::getSkipArtifactFilter() == 0) {

	  cerr << "Preparing artifact filter: ";
      FilterKnownOddities filter;
      cerr << MemoryUtils::getMemoryUsage() << endl;

	  cerr << "Applying sequence artifact filter to Input Files" << endl;
	  unsigned long filtered = filter.applyFilter(reads);
	  cerr << "filter affected (trimmed/removed) " << filtered << " Reads " << endl;;
	  cerr << MemoryUtils::getMemoryUsage() << endl;

	  cerr << "Applying DuplicateFragmentPair Filter to Input Files" << endl;
	  unsigned long duplicateFragments = filter.filterDuplicateFragments(reads);
	  cerr << "filter affected  (removed) " << duplicateFragments << endl;
	  cerr << MemoryUtils::getMemoryUsage() << endl;
	}

	KS spectrum(0);

	KoMer::MmapFileVector spectrumMmaps;

	if (Options::getKmerSize() > 0) {

	  long numBuckets = KS::estimateWeakKmerBucketSize(reads, 64);
	  cerr << "targeting " << numBuckets << " buckets for reads " << endl;

	  spectrum = KS(numBuckets);
	  cerr << MemoryUtils::getMemoryUsage() << endl;

	  TrackingData::minimumWeight = Options::getMinKmerQuality();

	  spectrumMmaps = spectrum.buildKmerSpectrumInParts(reads, Options::getBuildPartitions());
	  cerr << MemoryUtils::getMemoryUsage() << endl;

	  if (Options::getGCHeatMap() && ! outputFilename.empty()) {
		  cerr << "Creating GC Heat Map " <<  MemoryUtils::getMemoryUsage() << endl;
		  OfstreamMap ofmap(outputFilename + "-GC", ".txt");
		  spectrum.printGC(ofmap.getOfstream(""));
	  }

	  if (Options::getMinDepth() > 1) {
        cerr << "Clearing singletons from memory" << endl;
        spectrum.singleton.clear();
	    cerr << MemoryUtils::getMemoryUsage() << endl;
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
			selectReads(thisDepth, reads, spectrum, pickOutputFilename);
		}
	}

	spectrum.reset();

	return 0;
}

long selectReads(unsigned int minDepth, ReadSet &reads, KS &spectrum, std::string outputFilename)
{
	cerr << "Trimming reads: ";
	RS selector(reads, spectrum.weak, minDepth);
	cerr << MemoryUtils::getMemoryUsage() << endl;
	cerr << "Picking reads: " << endl;

	long oldPicked = 0;
	long picked = 0;

	int maximumKmerDepth = FilterReadsOptions::getMaxKmerDepth();

	OfstreamMap ofmap(outputFilename, ".fastq");

	if (maximumKmerDepth > 0) {
		for (int depth = 1; depth < maximumKmerDepth; depth++) {
			cerr << "Picking depth " << depth << " layer of reads" << endl;
			if (reads.hasPairs())
				picked += selector.pickBestCoveringSubsetPairs(depth,
						minDepth, Options::getMinReadLength(), FilterReadsOptions::getBothPairs());
			else
				picked += selector.pickBestCoveringSubsetReads(depth,
						minDepth, Options::getMinReadLength());
			cerr << MemoryUtils::getMemoryUsage() << endl;
		}

		if (picked > 0 && !outputFilename.empty()) {
			cerr << "Writing " << picked << " reads to output file(s)" << endl;
			selector.writePicks(ofmap, oldPicked);
		}
		cerr << MemoryUtils::getMemoryUsage() << endl;
		oldPicked += picked;


	} else {

		int maxDepth = FilterReadsOptions::getPartitionByDepth();
		if (maxDepth < 0) {
			maxDepth = 1;
		}

		for (unsigned int depth = maxDepth; depth >= 1; depth /= 2) {

			string ofname = outputFilename;
			if (maxDepth > 1) {
				ofname += "-PartitionDepth" + boost::lexical_cast< string >( depth );
			}
			ofmap = OfstreamMap(ofname, ".fastq");
			float tmpMinDepth = std::max(minDepth, depth);
			if (Options::getKmerSize() == 0) {
				tmpMinDepth = 0;
				depth = 0;
			}
			cerr << "Selecting reads over depth: " << depth << " (" << tmpMinDepth << ") " << endl;

			if (reads.hasPairs()) {
				picked = selector.pickAllPassingPairs(tmpMinDepth,
						Options::getMinReadLength(),
						FilterReadsOptions::getBothPairs());
			} else {
				picked = selector.pickAllPassingReads(tmpMinDepth,
						Options::getMinReadLength());
			}
			cerr << "At or above coverage: " << depth << " Picked " << picked
			<< " / " << reads.getSize() << " reads" << endl;
			cerr << MemoryUtils::getMemoryUsage() << endl;

			if (picked > 0 && !outputFilename.empty()) {
				cerr << "Writing " << picked << " reads  to output files" << endl;
				selector.writePicks(ofmap, oldPicked);
			}
			oldPicked += picked;

			if (minDepth > depth) {
				break;
			}
		}
	}
	ofmap.clear();
	cerr << "Done.  Cleaning up" << MemoryUtils::getMemoryUsage() << endl;
	selector.clear();

	return oldPicked;
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
