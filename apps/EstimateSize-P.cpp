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

class _MPIEstimateSizeOptions : public OptionsBaseInterface {
public:
        long getNumPoints() {
            return 250;
        } 
        long getTotalPoints() {
            return 1000;
        }
	void _resetDefaults() {
		MPIOptions::_resetDefaults();
                KmerOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
		// assign defaults
		GeneralOptions::getOptions().getMmapInput() = 0;
		GeneralOptions::getOptions().getVerbose() = 1;
                KmerOptions::getOptions().getMinDepth() = 1;
		KmerOptions::getOptions().getSaveKmerMmap() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
                p.add("kmer-size", 1);
                p.add("input-file", -1);

		MPIOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc, p);
                KmerOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);
		ret &= KmerOptions::_parseOptions(vm);

		return ret;
	}
};
typedef OptionsBaseTemplate< _MPIEstimateSizeOptions > MPIEstimateSizeOptions;

int main(int argc, char *argv[]) {

	mpi::communicator world = initializeWorldAndOptions< MPIEstimateSizeOptions >(argc, argv);

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files");
  
        long numPoints = MPIEstimateSizeOptions::getOptions().getNumPoints();
        long totalPoints = MPIEstimateSizeOptions::getOptions().getTotalPoints();

        if (world.rank() == 0) {
            std::cout << "totalBases\tuniqueKmers\trawKmers\n0\t0\t0\n" << KmerSizer::getSequenceLength() << "\t1\t1\n";
        }

        unsigned long totalBases = 0;
        long numBuckets = 0;
        KS spectrum(world, numBuckets);
        for (long iter = 0 ; iter < numPoints; iter++) {
            LOG_VERBOSE(1, "Starting iteration " << iter << " of " << numPoints);

	    ReadSet reads;
	    reads.appendAllFiles(inputs, world.rank()*totalPoints + iter, world.size()*totalPoints);

	    unsigned long counts[3], totalCounts[3];
	    unsigned long &readCount = counts[0] = reads.getSize();
	    unsigned long &numPairs  = counts[1] = reads.identifyPairs();
	    unsigned long &baseCount = counts[2] = reads.getBaseCount();
	    LOG_VERBOSE(2, "loaded " << readCount << " Reads, " << baseCount << " Bases ");

	    all_reduce(world, (unsigned long*) counts, 3, (unsigned long*) totalCounts, std::plus<unsigned long>());
	    LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Loaded " << totalCounts[0] << " distributed reads, " << totalCounts[1] << " distributed pairs, " << totalCounts[2] << " distributed bases");
            totalBases += totalCounts[2];

	    if (numBuckets == 0 && KmerOptions::getOptions().getKmerSize() > 0) {

		numBuckets = KS::estimateWeakKmerBucketSize(reads);

		numBuckets = all_reduce(world, numBuckets, mpi::maximum<int>()) * numPoints;
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "targeting " << numBuckets << " buckets for reads");
	        spectrum = KS(world, numBuckets);
	    }
	    if (KmerOptions::getOptions().getKmerSize() > 0) {

		spectrum.buildKmerSpectrum(reads);
                KS::MPIHistogram h = spectrum._getHistogram(false);
                if (world.rank() == 0) {
		  h.finish();
                  std:: cout << totalBases << "\t" << h.getCount() << "\t" << h.getTotalCount() << std::endl;
                }
		if (Log::isDebug(1)) {
			std::string hist = h.toString();
			LOG_DEBUG_OPTIONAL(1, world.rank() == 0, "Collective Kmer Histogram\n" << hist);
		}
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

