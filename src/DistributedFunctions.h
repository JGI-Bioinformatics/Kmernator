//
// Kmernator/src/DistributedFunctions.h
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

#ifndef DISTRIBUTED_FUNCTIONS_H_
#define DISTRIBUTED_FUNCTIONS_H_

#include "config.h"
#include "Options.h"
#include "Kmer.h"
#include "KmerSpectrum.h"

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

void reduceOMPThreads(mpi::communicator &world) {
	int numThreads = omp_get_max_threads();
	numThreads = all_reduce(world, numThreads, mpi::minimum<int>());
	omp_set_num_threads(numThreads);
	if (world.rank() == 0)
		LOG_DEBUG_MT(1, world.rank() << ": set OpenMP threads to " << numThreads);
}

template<typename So, typename We, typename Si = TrackingDataSingleton>
class DistributedKmerSpectrum : public KmerSpectrum<So, We, Si>
{
	typedef KmerSpectrum<So, We, Si> KS;
	typedef Kmer::NumberType NumberType;
	typedef Kmer::IndexType IndexType;
	typedef Kmernator::MmapFile MmapFile;
	typedef typename KS::Histogram Histogram;

public:
	DistributedKmerSpectrum(unsigned long buckets = 0, bool separateSingletons = true)
	: KS(buckets, separateSingletons) {}
	~DistributedKmerSpectrum() {}
	DistributedKmerSpectrum &operator=(const KS &other) {
		*((KS*) this) = other;
		return *this;
	}

	template<typename D>
	MmapFile writeKmerMap(mpi::communicator &world, D &kmerMap, std::string filepath) {
		NumberType numBuckets = kmerMap.getNumBuckets();
		NumberType *mySizeCounts = new NumberType[ numBuckets ];
		for(IndexType i = 0 ; i < numBuckets; i++)
			mySizeCounts[i] = kmerMap.getBucketByIdx(i).size();
		NumberType *ourSizeCounts = new NumberType[ numBuckets ];
		all_reduce(world, mySizeCounts, numBuckets, ourSizeCounts, std::plus<NumberType>());

		// rename variable for clarity
		NumberType *offsetArray = mySizeCounts;
		NumberType offset = sizeof(NumberType) * (2+numBuckets);
		NumberType totalSize = 0;
		for(IndexType i = 0; i < numBuckets; i++) {
			totalSize += ourSizeCounts[i];
			offsetArray[i] = offset;
			offset += D::BucketType::sizeToStore( ourSizeCounts[i] );
		}

		delete [] ourSizeCounts;

		MmapFile mmap;
		NumberType totalMmapSize = kmerMap.getSizeToStoreCountsAndIndexes() + totalSize * D::BucketType::getElementByteSize();
		if (world.rank() == 0) {
			mmap = MmapTempFile::buildNewMmap(totalMmapSize , filepath);
			NumberType *numbers = (NumberType *) mmap.data();
			*(numbers++) = numBuckets;
			*(numbers++) = kmerMap.getBucketMask();
			// store offsetArray in mmap
			for(IndexType i = 0; i < numBuckets; i++) {
				*(numbers++) = offsetArray[i];
			}

			world.barrier();
		} else {
			world.barrier();
			mmap = MmapFile(filepath, std::ios_base::in | std::ios_base::out, totalMmapSize);
		}

		// store our part of mmap (interleaved DMP)
		for(IndexType i = 0 ; i < numBuckets; i++) {
			if (kmerMap.getBucketByIdx(i).size() > 0)
				kmerMap.getBucketByIdx(i).store(mmap.data() + offsetArray[i]);
		}

		delete [] mySizeCounts;	// aka offsetArray
		world.barrier();

		// swap memory maps
		const D full = D::restore(mmap.data());
		kmerMap.swap( const_cast<D&>(full) );
		madvise(const_cast<char*>(mmap.data()), mmap.size(), MADV_RANDOM);

		return mmap;
	};

	void _buildKmerSpectrum( mpi::communicator &world, ReadSet &store, bool isSolid = false ) {
//TODO
	}
	void printHistograms(mpi::communicator &world, std::ostream &os, bool printSolidOnly = false) {
		Histogram histogram(127);

		if (!printSolidOnly) {
			this->setWeakHistogram(histogram);
			this->setSingletonHistogram(histogram);
		}
		this->setSolidHistogram(histogram);


		if (world.rank() == 0)
			os << histogram.toString();
	}

	Kmernator::MmapFileVector buildKmerSpectrumInParts(mpi::communicator &world, ReadSet &store ) {
		bool isSolid = false; // not supported for references...
		int numParts = world.size();
		if (numParts == 1) {
			this->buildKmerSpectrum(store);
			return Kmernator::MmapFileVector();
		}

		assert((numParts & (numParts-1)) == 0); // numParts must be a power of 2

		_buildKmerSpectrum(world, store, isSolid);

		// purge low counts
		if (Options::getMinDepth() > 1) {
			LOG_VERBOSE_MT(1, world.rank() << ": Clearing memory from singletons: " << this->singleton.size() << std::endl << MemoryUtils::getMemoryUsage());
			this->singleton.clear();
		}
		if (Options::getMinDepth() > 2) {
			LOG_VERBOSE_MT(1, world.rank() << ": Purging low count kmers (< " << Options::getMinDepth() << ")" << std::endl << MemoryUtils::getMemoryUsage());
			this->purgeMinDepth(Options::getMinDepth());
		}

		// communicate sizes and allocate permanent file
		LOG_VERBOSE(1, "Merging partial spectrums" << std::endl << MemoryUtils::getMemoryUsage() );
		Kmernator::MmapFileVector ourSpectrum(2);
		ourSpectrum[0] = writeKmerMap(world, this->weak, Options::getTmpDir() + "/weak-kmer-mmap");

		if (Options::getMinDepth() <= 1) {
			ourSpectrum[1] = writeKmerMap(world, this->singleton, Options::getTmpDir() + "/singleton-kmer-mmap");
		}

		LOG_VERBOSE(1, "Finished merging partial spectrums" << std::endl << MemoryUtils::getMemoryUsage());

		if (Log::isVerbose(1)) {
			this->printStats(Log::Verbose("Final Stats"), store.getSize(), isSolid, true);
			if (!isSolid) {
				// TODO printHistograms(world, Log::Verbose("Final Histogram"));
			}
		}

		return ourSpectrum;

	}
};

#endif /* DISTRIBUTED_FUNCTIONS_H_ */
