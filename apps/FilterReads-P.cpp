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

class _MPIFilterReadsOptions : public _FilterReadsBaseOptions, public _MPIOptions {
public:
	void _resetDefaults() {
		_FilterReadsBaseOptions::_resetDefaults();
		_MPIOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
		// assign defaults
		GeneralOptions::getOptions().getMmapInput() = 0;
		GeneralOptions::getOptions().getVerbose() = 2;
		KmerOptions::getOptions().getSaveKmerMmap() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		_FilterReadsBaseOptions::_setOptions(desc, p);
		_MPIOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= _MPIOptions::_parseOptions(vm);
		ret &= _FilterReadsBaseOptions::_parseOptions(vm);

		return ret;
	}
};
typedef OptionsBaseTemplate< _MPIFilterReadsOptions > MPIFilterReadsOptions;

int main(int argc, char *argv[]) {

	mpi::communicator world = initializeWorldAndOptions< MPIFilterReadsOptions >(argc, argv);

	if (MPIFilterReadsOptions::getOptions().getMaxKmerDepth() > 0 && world.size() > 1)
		LOG_THROW("Distributed version does not support max-kmer-output-depth option");

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	ReadSet reads;

	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();
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

	setGlobalReadSetOffsets(world, reads);

	if (Options::getOptions().getSkipArtifactFilter() == 0) {

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

	if ( Options::getOptions().getDeDupMode() > 0 && Options::getOptions().getDeDupEditDistance() >= 0) {
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
	if (KmerOptions::getOptions().getKmerSize() > 0) {

		numBuckets = KS::estimateWeakKmerBucketSize(reads);

		numBuckets = all_reduce(world, numBuckets, mpi::maximum<int>());
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "targeting " << numBuckets << " buckets for reads");
	}
	KS spectrum(world, numBuckets);
	Kmernator::MmapFileVector spectrumMmaps;
	if (KmerOptions::getOptions().getKmerSize() > 0 && !KmerOptions::getOptions().getLoadKmerMmap().empty()) {
		spectrum.restoreMmap(KmerOptions::getOptions().getLoadKmerMmap());
	} else if (KmerOptions::getOptions().getKmerSize() > 0) {
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		spectrum.buildKmerSpectrum(reads);

		std::string sizeHistoryFile = MPIFilterReadsOptions::getOptions().getSizeHistoryFile();
		if (!sizeHistoryFile.empty()) {
			spectrum.trackSpectrum(true);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Writing size history file to: " << sizeHistoryFile);
			KS::SizeTracker reducedSizeTracker = spectrum.reduceSizeTracker(world);
			if (world.rank() == 0) {
				OfstreamMap ofm(sizeHistoryFile, "");
				ofm.getOfstream("") << reducedSizeTracker.toString();
			}
		}

		if (Log::isVerbose(1)) {
			std::string hist = spectrum.getHistogram(false);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Kmer Histogram\n" << hist);
		}
	}
	if (KmerOptions::getOptions().getKmerSize() > 0) {

		if (Options::getOptions().getVariantSigmas() > 0.0) {
			long purgedVariants = spectrum.purgeVariants();
			long totalPurgedVariants = all_reduce(world, purgedVariants, std::plus<long>());
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Distributed Purged " << totalPurgedVariants << " kmer variants");

			std::string hist = spectrum.getHistogram(false);

			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Variant Purged Kmer Histogram\n" << hist);
			world.barrier();

		}

		if (!outputFilename.empty() && KmerOptions::getOptions().getSaveKmerMmap() > 0) {
			spectrumMmaps = spectrum.writeKmerMaps(outputFilename + "-mmap");
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
        }

		if (KmerOptions::getOptions().getMinDepth() > 1) {
			LOG_DEBUG(1, "Clearing singletons from memory");
			spectrum.singleton.clear();
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
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
				LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Trimming reads with minDepth: " << thisDepth);
			} else {
				LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Trimming reads that pass Artifact Filter with length: " << Options::getOptions().getMinReadLength());
			}
			RS selector(world, reads, spectrum.weak);
			selector.scoreAndTrimReads(minDepth);

			// TODO implement a more efficient algorithm to output data in order

			// rank 0 will overwrite, all others will append
			if (world.rank() != 0)
				OfstreamMap::getDefaultAppend() = true;

			// let only one rank at a time write to the files
			LOG_VERBOSE(1, "Writing Files");

			selectReads(thisDepth, reads, selector, pickOutputFilename);

		}
	}

	LOG_DEBUG(2, "Clearing spectrum");
	spectrum.reset();
	LOG_DEBUG(1, "Finished, waiting for rest of collective");

	world.barrier();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");

	MPI_Finalize();

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
