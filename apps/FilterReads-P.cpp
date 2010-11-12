//
// Kmernator/apps/FilterReads-P.cpp
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
#include "DistributedFunctions.h"

typedef TrackingDataWithDirection DataType;
typedef DistributedKmerSpectrum<DataType, DataType> KS;
typedef DistributedReadSelector<DataType> RS;

int main(int argc, char *argv[]) {

	// assign defaults
	Options::getMmapInput() = 0;
	Options::getVerbosity() = 2;

	int threadSupport = MPI::Init_thread(MPI_THREAD_MULTIPLE);
	mpi::environment env(argc, argv);
	mpi::communicator world;

	try {

		validateMPIWorld(world, threadSupport);

		if (!FilterReadsOptions::parseOpts(argc, argv))
			throw std::invalid_argument("Please fix the command line arguments");

		if (FilterReadsOptions::getMaxKmerDepth() > 0 && world.size() > 1)
			throw std::invalid_argument("Distributed version does not support max-kmer-output-depth option");

	} catch (...) {
		std::cerr << std::endl << "Please fix the options and/or MPI environment" << std::endl;
		exit(1);
	}
	world.barrier();

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOutputFile();

	ReadSet reads;
	KmerSizer::set(Options::getKmerSize());

	Options::FileListType inputs = Options::getInputFiles();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files");

	reads.appendAllFiles(inputs, world.rank(), world.size());

	LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Identifying Pairs: ");

	unsigned long counts[3], totalCounts[3];
	unsigned long &readCount = counts[0] = reads.getSize();
	unsigned long &numPairs  = counts[1] = reads.identifyPairs();
	unsigned long &baseCount = counts[2] = reads.getBaseCount();
	LOG_VERBOSE(2, "loaded " << readCount << " Reads, " << baseCount << " Bases ");
	LOG_VERBOSE(2, "Pairs + single = " << numPairs);

	all_reduce(world, (unsigned long*) counts, 3, (unsigned long*) totalCounts, std::plus<unsigned long>());
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Loaded " << totalCounts[0] << " distributed reads, " << totalCounts[1] << " distributed pairs, " << totalCounts[2] << " distributed bases");

	setGlobalReadSetOffset(world, reads);

	if (Options::getSkipArtifactFilter() == 0) {

		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Preparing artifact filter: ");

		FilterKnownOddities filter;
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		LOG_VERBOSE_OPTIONAL(2, world.rank() == 0, "Applying sequence artifact filter to Input Files");

		unsigned long filtered = filter.applyFilter(reads);

		LOG_VERBOSE(2, "local filter affected (trimmed/removed) " << filtered << " Reads ");
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		unsigned long allFiltered;
		reduce(world, filtered, allFiltered, std::plus<unsigned long>(), 0);
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "distributed filter (trimmed/removed) " << allFiltered << " Reads ");

	}

	if ( Options::getDeDupMode() > 0 && Options::getDeDupEditDistance() >= 0) {
		if (world.size() == 1) {
			LOG_VERBOSE(2, "Applying DuplicateFragmentPair Filter to Input Files");
			unsigned long duplicateFragments = DuplicateFragmentFilter::filterDuplicateFragments(reads);

			LOG_VERBOSE(2, "filter removed duplicate fragment pair reads: " << duplicateFragments);
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			unsigned long allDuplicateFragments;
			reduce(world, duplicateFragments, allDuplicateFragments, std::plus<unsigned long>(), 0);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "distributed removed duplicate fragment pair reads: " << allDuplicateFragments);
		} else {
			if (world.rank() == 0)
				LOG_WARN(1, "Distributed DuplicateFragmentPair Filter is not supported (yet)." << std::endl
					<< "If you want this feature please run the non-MPI FilterReads");
		}

	}

	long numBuckets = 0;
	if (Options::getKmerSize() > 0) {

		numBuckets = KS::estimateWeakKmerBucketSize(reads, 64);

		numBuckets = all_reduce(world, numBuckets, mpi::maximum<int>());
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "targeting " << numBuckets << " buckets for reads");
	}
	KS spectrum(world, numBuckets);
	Kmernator::MmapFileVector spectrumMmaps;

	if (Options::getKmerSize() > 0) {
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		TrackingData::minimumWeight = Options::getMinKmerQuality();

		spectrum.buildKmerSpectrum(reads);
		if (Options::getVariantSigmas() > 0.0) {
			long purgedVariants = spectrum.purgeVariants();
			long totalPurgedVariants = all_reduce(world, purgedVariants, std::plus<long>());
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Distributed Purged " << totalPurgedVariants << " kmer variants");

			std::string hist = spectrum.getHistogram(false);

			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Variant Purged Histogram\n" << hist);
			world.barrier();

		}

		spectrumMmaps = spectrum.writeKmerMaps(Options::getOutputFile() + "-mmap");
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

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
			}
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Trimming reads with minDepth: " << thisDepth);
			RS selector(world, reads, spectrum.weak);
			selector.scoreAndTrimReads(minDepth);

			// TODO implement a more efficient algorithm to output data in order

			// rank 0 will overwrite, all others will append
			if (world.rank() != 0)
				OfstreamMap::getAppend() = true;

			// let only one rank at a time write to the files
			LOG_VERBOSE(1, "Writing Files");
			int rank = 0;
			while (rank < world.size()) {
				if (rank == world.rank()) {
					LOG_VERBOSE_OPTIONAL(1, true, "Writing files part " << (rank+1) << " of " << world.size());
					selectReads(thisDepth, reads, selector, pickOutputFilename);
				}
				world.barrier();
				rank++;
			}
		}
	}

	LOG_DEBUG(2, "Clearing spectrum");
	spectrum.reset();
	LOG_DEBUG(1, "Finished, waiting for rest of collective");

	world.barrier();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");

	MPI::Finalize();

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
