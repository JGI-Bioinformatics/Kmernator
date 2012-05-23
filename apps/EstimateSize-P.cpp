//
// Kmernator/apps/EstimateSize-P.cpp
//
// Author: Rob Egan
//
// Copyright 2012 The Regents of the University of California.
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
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "KmerTrackingData.h"
#include "Utils.h"
#include "Log.h"

#include <boost/lexical_cast.hpp>

using namespace std;

#include "DistributedFunctions.h"

typedef TrackingDataWithDirection DataType;
typedef DistributedKmerSpectrum<DataType, DataType> KS;
typedef DistributedReadSelector<DataType> RS;

class _MPIEstimateSizeOptions: public OptionsBaseInterface {
public:
	_MPIEstimateSizeOptions() :
		samplePartitions(50),
		maxSampleFraction(0.05) {
	}
	~_MPIEstimateSizeOptions() {
	}

	long &getSamplePartitions() {
		return samplePartitions;
	}
	double &getMaxSampleFraction() {
		return maxSampleFraction;
	}
	void _resetDefaults() {
		MPIOptions::_resetDefaults();
		KmerBaseOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();
		// assign defaults
		GeneralOptions::getOptions().getMmapInput() = 0;
		GeneralOptions::getOptions().getVerbose() = 1;
		KmerSpectrumOptions::getOptions().getMinDepth() = 1;
		KmerSpectrumOptions::getOptions().getSaveKmerMmap() = 0;
		KmerSpectrumOptions::getOptions().getKmerSubsample() = 1000;
	}

	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {
		p.add("kmer-size", 1);
		p.add("input-file", -1);

		po::options_description opts("EstimateSize Options");
		opts.add_options()

						("sample-partitions",po::value<long>()->default_value(samplePartitions),"The number of partitions to break up the file")

						("max-sample-fraction", po::value<double>()->default_value(maxSampleFraction), "The maximum amount of data to read")


						;
		desc.add(opts);
		KmerSpectrumOptions::_setOptions(desc, p);

		MPIOptions::_setOptions(desc, p);
		GeneralOptions::_setOptions(desc, p);
		KmerBaseOptions::_setOptions(desc, p);
		KmerSpectrumOptions::_setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);
		ret &= KmerBaseOptions::_parseOptions(vm);
		ret &= KmerSpectrumOptions::_parseOptions(vm);

		setOpt<long> ("sample-partitions", samplePartitions);
		setOpt<double> ("max-sample-fraction", maxSampleFraction);

		return ret;
	}
protected:
	long samplePartitions;
	double maxSampleFraction;

};
typedef OptionsBaseTemplate< _MPIEstimateSizeOptions > MPIEstimateSizeOptions;


int main(int argc, char *argv[]) {

	mpi::communicator world = initializeWorldAndOptions< MPIEstimateSizeOptions >(argc, argv);

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files");

	long partitions = MPIEstimateSizeOptions::getOptions().getSamplePartitions();
	double maxFraction = MPIEstimateSizeOptions::getOptions().getMaxSampleFraction();
	assert(maxFraction < 1.0);
	double fraction = 0.0;
	long totalPartitions = (long) partitions / maxFraction;

	unsigned long totalReads = 0;
	unsigned long totalBases = 0;
	long numBuckets = 0;
	KS spectrum(world, numBuckets);
	for (long iter = 0 ; iter < partitions && fraction < maxFraction; iter++) {
		fraction += (double) 1. / (double) totalPartitions;

		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0,  "Starting iteration " << iter << " at " << fraction*100 << "%");

		ReadSet reads;
		reads.appendAllFiles(inputs, world.rank()*totalPartitions + iter, world.size()*totalPartitions);
		setGlobalReadSetOffsets(world, reads);

		unsigned long counts[3], totalCounts[3];
		unsigned long &readCount = counts[0] = reads.getSize();
		unsigned long &baseCount = counts[2] = reads.getBaseCount();
		LOG_VERBOSE(2, "loaded " << readCount << " Reads, " << baseCount << " Bases ");

		all_reduce(world, (unsigned long*) counts, 3, (unsigned long*) totalCounts, std::plus<unsigned long>());
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Loaded " << totalCounts[0] << " distributed reads, " << totalCounts[2] << " distributed bases");
		totalReads += totalCounts[0];
		totalBases += totalCounts[2];

		if (KmerBaseOptions::getOptions().getKmerSize() > 0) {

			if (numBuckets == 0) {

				numBuckets = KS::estimateWeakKmerBucketSize(reads) * partitions / KS::getKmerSubsample();

				numBuckets = all_reduce(world, numBuckets, mpi::maximum<int>());

				LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "targeting " << numBuckets << " buckets for reads");
				spectrum = KS(world, numBuckets);
			}

			spectrum.buildKmerSpectrum(reads);
			spectrum.trackSpectrum(true);
			LOG_DEBUG_OPTIONAL(1, true, "SizeTracker: " << spectrum.getSizeTracker().toString());

			if (Log::isDebug(1)) {
				KS::MPIHistogram h = spectrum._getHistogram(false);
				std::string hist = h.toString();
				LOG_DEBUG_OPTIONAL(1, world.rank() == 0, "Collective Kmer Histogram\n" << hist);
			}
		}

	}

	std::string hist = spectrum.getHistogram(false);
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Kmer Histogram\n" << hist);

	KS::SizeTracker reducedSizeTracker = spectrum.reduceSizeTracker(world);
	std::string reducedSizeTrackerFile = outputFilename;
	if (reducedSizeTrackerFile.empty()) {
		reducedSizeTrackerFile = UniqueName::generateUniqueName("tmp-estimateSize");
	}
	float errorRate = TrackingData::getErrorRate();
	LOG_DEBUG_OPTIONAL(1, true, "Kmer error rate: " << errorRate);
	float commonErrorRate = 0.0;
	MPI_Reduce(&errorRate, &commonErrorRate, 1, MPI_FLOAT, MPI_SUM, 0, world);
	commonErrorRate /= (float) world.size();
	if (world.rank() == 0) {
		{
			LOG_VERBOSE_OPTIONAL(1, true, "Writing size tracking file to:" << reducedSizeTrackerFile);
			LOG_DEBUG_OPTIONAL(1, true, "SizeTracker:\n" << reducedSizeTracker.toString());
			OfstreamMap ofm(reducedSizeTrackerFile, "");
			ofm.getOfstream("") << reducedSizeTracker.toString();
		}
		std::string basePath = FileUtils::getBasePath(argv[0]);
		std::stringstream cmdss;
		cmdss << "Rscript " << basePath << "/EstimateSize.R " << reducedSizeTrackerFile; // << " " << (commonErrorRate*1.25);
		std::string command = cmdss.str();

		if (!FileUtils::getBasePath("Rscript").empty() || basePath.empty()) {

			LOG_DEBUG_OPTIONAL(1, true, "Executing: " << command);
			IPipestream ipipe(command);
			double errorRate = 0.0, genomeSize = 0.0;
			bool readValues = false;
			while (ipipe.good() && !ipipe.eof()) {
				std::string line;
				getline(ipipe, line);
				LOG_DEBUG_OPTIONAL(2, true, "Read: " << line);
				if (line.find("errorRate") != std::string::npos) {
					LOG_DEBUG_OPTIONAL(2, true, "Found headers in: " << line);
					readValues = true;
					continue;
				}
				if (readValues) {
					readValues = false;
					LOG_DEBUG_OPTIONAL(2, true, "Reading errorRate and GenomeSize from " << line);
					std::stringstream ss;
					ss << line;
					ss >> errorRate;
					ss >> genomeSize;
				}
			}
			LOG_VERBOSE_OPTIONAL(1, true, "Estimated Kmer-quality errorRate: " << commonErrorRate);
			LOG_VERBOSE_OPTIONAL(1, true, "Distributed readCount: " << totalReads);
			LOG_VERBOSE_OPTIONAL(1, true, "Estimated fractionRead: " << fraction);
			LOG_VERBOSE_OPTIONAL(1, true, "Estimated errorRate: " << errorRate);
			LOG_VERBOSE_OPTIONAL(1, true, "Estimated genomeSize: " << genomeSize);

			ipipe.close();

			double totalRawKmers = reducedSizeTracker.getLastElement().rawKmers / fraction;
			double estimatedUniqueKmers = totalRawKmers * errorRate + genomeSize;
			LOG_VERBOSE_OPTIONAL(1, true, "Estimated totalRawKmers: " << totalRawKmers);
			LOG_VERBOSE_OPTIONAL(1, true, "Estimated totalUniqueKmers: " << estimatedUniqueKmers);
			if (reducedSizeTrackerFile.compare(outputFilename) != 0) {
				LOG_DEBUG_OPTIONAL(1, true, "Removing temporary size tracking file: " << reducedSizeTrackerFile);
				unlink(reducedSizeTrackerFile.c_str());
			}
		} else {
			LOG_ERROR(1, "Could not find Rscript or EstimateSize.R.  To get the coefficients, please run '" << command << "'");
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

