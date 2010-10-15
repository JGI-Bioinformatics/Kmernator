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
#include "MPIBuffer.h"

#include "boost/optional.hpp"
#include <boost/thread/thread.hpp>

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

void reduceOMPThreads(mpi::communicator &world) {
	int numThreads = omp_get_max_threads();
	numThreads = all_reduce(world, numThreads, mpi::minimum<int>());
	omp_set_num_threads(numThreads);
	if (world.rank() == 0)
		LOG_DEBUG(1, "set OpenMP threads to " << numThreads);
}

template<typename So, typename We, typename Si = TrackingDataSingleton>
class DistributedKmerSpectrum : public KmerSpectrum<So, We, Si>
{
public:
	typedef KmerSpectrum<So, We, Si> KS;
	typedef Kmer::NumberType NumberType;
	typedef Kmer::IndexType IndexType;
	typedef Kmernator::MmapFile MmapFile;
	typedef typename KS::Histogram Histogram;
	typedef typename KS::ReadSetSizeType ReadSetSizeType;
	typedef typename KS::PositionType PositionType;
	typedef typename KS::WeightType WeightType;
	typedef typename KS::DataPointers DataPointers;

public:
	DistributedKmerSpectrum(unsigned long buckets = 0, bool separateSingletons = true)
	: KS(buckets, separateSingletons) {
	}
	~DistributedKmerSpectrum() {
	}
	DistributedKmerSpectrum &operator=(const KS &other) {
		*((KS*) this) = other;
		return *this;
	}

	template<typename D>
	MmapFile writeKmerMapMPI(mpi::communicator &world, D &kmerMap, std::string filepath) {
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

	/*
	 * each rank & thread listens for worldSize IRecv messages (rank, data tag)
	 * each rank scans ReadSet, each thread scheduled dynamically
	 * each thread builds message buffer of kmers for that rank & thread (numThreads * worldSize * numThreads) buffers, sends when full.
	 * send first receives & processes all pending messages, initiates new irecv if msg was non-zero, if zero, checkpoint increments
	 * done when all requests are MPI_REQUEST_NULL
	 *
	 * message (for a specific rank & thread):
	 * readIdx, readPos, weight (rc if neg) + kmer
	 *
	 */
	class StoreKmerMessageHeader;

	typedef MPIMessageBufferBase< StoreKmerMessageHeader > StoreKmerMessageBuffersBase;
	typedef MPIRecvMessageBuffer< StoreKmerMessageHeader > RecvStoreKmerMessageBufferBase;
	typedef MPISendMessageBuffer< StoreKmerMessageHeader > SendStoreKmerMessageBufferBase;

	class RecvStoreKmerMessageBuffer : public RecvStoreKmerMessageBufferBase
	{
		DistributedKmerSpectrum &_spectrum;
		DataPointers _pointers;

	public:
		RecvStoreKmerMessageBuffer(mpi::communicator &world, DistributedKmerSpectrum &spectrum, int messageSize, int tag)
			: RecvStoreKmerMessageBufferBase(world, messageSize, tag), _spectrum(spectrum), _pointers(spectrum) {
		}
		~RecvStoreKmerMessageBuffer() {}
		inline DataPointers &getDataPointer() {
			return _pointers;
		}
		inline DistributedKmerSpectrum &getSpectrum() {
			return _spectrum;
		}
	};
	class SendStoreKmerMessageBuffer : public SendStoreKmerMessageBufferBase
	{
		DistributedKmerSpectrum &_spectrum;

	public:
		SendStoreKmerMessageBuffer(mpi::communicator &world, DistributedKmerSpectrum &spectrum, int messageSize)
			: SendStoreKmerMessageBufferBase(world, messageSize), _spectrum(spectrum) {
		}
		~SendStoreKmerMessageBuffer() {}
		inline DistributedKmerSpectrum &getSpectrum() {
			return _spectrum;
		}
	};

	class StoreKmerMessageHeader {
	public:
		ReadSetSizeType readIdx;
		PositionType readPos;
		WeightType weight; // weight is negative if kmer is rc of observed direction
		// Kmer is next bytes, dynamically determined by KmerSizer::getTwoBitLength()
		// kmer is least complement

		// THIS IS DANGEROUS unless allocated an extra Kmer!
		Kmer *getKmer() {
			return (Kmer*) (((char*)this)+sizeof(*this));
		}
		void set(ReadSetSizeType _readIdx, PositionType _readPos, WeightType _weight, const Kmer &_kmer) {
			readIdx = _readIdx;
			readPos = _readPos;
			weight = _weight;
			*(getKmer()) = _kmer;
		}
		void process(RecvStoreKmerMessageBufferBase *bufferCallback) {
			RecvStoreKmerMessageBuffer *kbufferCallback = (RecvStoreKmerMessageBuffer*) bufferCallback;
			LOG_DEBUG(4, bufferCallback->getWorld().rank() << ": " << omp_get_thread_num() << ": message: " << readIdx << " " << readPos << " " << weight << " " << getKmer()->toFasta());
			kbufferCallback->getSpectrum().append(kbufferCallback->getDataPointer(), *getKmer(), weight, readIdx, readPos);
		}
	};


	void _buildKmerSpectrumMPI( mpi::communicator &world, ReadSet &store, bool isSolid = false ) {
		int numThreads = omp_get_max_threads();
		int rank = world.rank();
		int messageSize = sizeof(StoreKmerMessageHeader) + KmerSizer::getTwoBitLength();

		NumberType distributedThreadMask = world.size() - 1;

		LOG_VERBOSE(1, "starting _buildSpectrum");

		SendStoreKmerMessageBuffer *sendBuffers[numThreads][numThreads];
		RecvStoreKmerMessageBuffer *recvBuffers[numThreads];

		// share ReadSet sizes for globally unique readIdx calculations
		long readSetSize = store.getSize();
		long readSetSizesInput[ world.size() ];
		long readSetSizes[ world.size() ];
		for(int i = 0; i < world.size(); i++)
			readSetSizesInput[i] = i == world.rank() ? readSetSize : 0;
		mpi::all_reduce(world, (long*) readSetSizesInput, world.size(), (long*) readSetSizes, mpi::maximum<long>());
		long globalReadSetOffset = 0;
		for(int i = 0; i < world.rank(); i++)
			globalReadSetOffset += readSetSizes[i];
		LOG_DEBUG(2, "reduced readSetSizes " << globalReadSetOffset);

		LOG_DEBUG(2, "building spectrum using " << numThreads << " threads (" << omp_get_max_threads() << ")");

		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();
			LOG_DEBUG(3, "allocating buffers for thread");
			recvBuffers[threadId] = new RecvStoreKmerMessageBuffer(world, *this, messageSize, threadId);
			for(int recvThread = 0 ; recvThread < numThreads; recvThread++) {
				sendBuffers[threadId][recvThread] = new SendStoreKmerMessageBuffer(world, *this, messageSize);
				sendBuffers[threadId][recvThread]->addCallback( *recvBuffers[threadId] );
			}


			for(long readIdx = threadId ; readIdx < readSetSize; readIdx+=numThreads)
			{

				const Read &read = store.getRead( readIdx );

				if (read.isDiscarded())
					continue;

				KmerWeights kmers = KmerReadUtils::buildWeightedKmers(read, true, true);
				ReadSetSizeType globalReadIdx = readIdx + globalReadSetOffset;
				LOG_DEBUG(2, "Read " << readIdx << " (" << globalReadIdx << ") " << kmers.size() );

				for (PositionType readPos = 0 ; readPos < kmers.size(); readPos++) {
					int rankDest, threadDest;
					WeightType weight = kmers.valueAt(readPos);
					this->solid.getThreadIds(kmers[readPos], threadDest, numThreads, rankDest, distributedThreadMask);

					if (rankDest == rank && threadDest == threadId) {
						this->append(recvBuffers[ threadId ]->getDataPointer(), kmers[readPos], weight, globalReadIdx, readPos);
					} else {
						sendBuffers[ threadId ][ threadDest ]->bufferMessage(rankDest, threadDest)->set(globalReadIdx, readPos, weight, kmers[readPos]);
					}
				}

				if (threadId == 0 && readIdx > 0 && readIdx % 100000 == 0 )
					LOG_VERBOSE(1, "processed " << readIdx << " reads");

			}

			LOG_DEBUG(2, "finished generating kmers from reads");

			// receiving any new messages
			recvBuffers[threadId]->receiveAllIncomingMessages();

			LOG_DEBUG(2, "sending final flush");
			// send all pending buffers
			for(int destThread = 0; destThread < numThreads; destThread++) {
				sendBuffers[threadId][destThread]->flushAllMessageBuffers(destThread);
			}
			LOG_DEBUG(2, "sending final messages")
			for(int destThread = 0; destThread < numThreads; destThread++) {
				sendBuffers[threadId][destThread]->finalize(destThread);
				delete sendBuffers[threadId][destThread];
			}

			LOG_DEBUG(2, "receiving final messages");
			recvBuffers[threadId]->finalize(numThreads);
			delete recvBuffers[threadId];

		}
		world.barrier();
		LOG_DEBUG(1, "finished _buildKmerSpectrum");
	}

	void printHistogramsMPI(mpi::communicator &world, std::ostream &os, bool printSolidOnly = false) {
		Histogram histogram(127);

		if (!printSolidOnly) {
			this->setWeakHistogram(histogram);
			this->setSingletonHistogram(histogram);
		}
		this->setSolidHistogram(histogram);

		Histogram dest(127);
		mpi::reduce(world, histogram, dest, std::plus<Histogram>(), 0);

		if (world.rank() == 0)
			os << dest.toString();
	}

	Kmernator::MmapFileVector buildKmerSpectrumMPI(mpi::communicator &world, ReadSet &store ) {
		bool isSolid = false; // not supported for references...
		int numParts = world.size();

		assert((numParts & (numParts-1)) == 0); // numParts must be a power of 2

		_buildKmerSpectrumMPI(world, store, isSolid);

		// purge low counts
		if (Options::getMinDepth() > 1) {
			LOG_VERBOSE(1, "Clearing memory from singletons: " << this->singleton.size() << std::endl << MemoryUtils::getMemoryUsage());
			this->singleton.clear();
		}
		if (Options::getMinDepth() > 2) {
			LOG_VERBOSE(1, "Purging low count kmers (< " << Options::getMinDepth() << ")" << std::endl << MemoryUtils::getMemoryUsage());
			this->purgeMinDepth(Options::getMinDepth());
		}

		// communicate sizes and allocate permanent file
		LOG_VERBOSE(1, "Merging partial spectrums" << std::endl << MemoryUtils::getMemoryUsage() );
		Kmernator::MmapFileVector ourSpectrum(2);
		ourSpectrum[0] = writeKmerMapMPI(world, this->weak, Options::getTmpDir() + "/weak-kmer-mmap");

		if (Options::getMinDepth() <= 1) {
			ourSpectrum[1] = writeKmerMapMPI(world, this->singleton, Options::getTmpDir() + "/singleton-kmer-mmap");
		}

		LOG_VERBOSE(1, "Finished merging partial spectrums" << std::endl << MemoryUtils::getMemoryUsage());

		if (Log::isVerbose(1)) {
			this->printStats(Log::Verbose("Final Stats"), store.getSize(), isSolid, true);
			if (!isSolid) {
				//printHistogramsMPI(world, Log::Verbose("Final Histogram"));
			}
		}

		return ourSpectrum;

	}

	/*
	 * ReadTrim
	 * each rank IRecv message for each file to output from previous rank (except rank that first opens file)

message:
receiveMessage[worldSize-1] { numReads, trimLengths[numReads], twoBitLengths[numReads], TwoBitReadBytes[ numReads ] }, numReads, trimLengths[numReads]
    last is reduction from worldSize cycles ago, apply to input files and output
reduce in-place all trimLengths for message
sendMessage[worldSize-1] { numReads, trimLengths[numReads], twoBitLengths[numRead], twoBitReadBytes[numReads] }, numReads, trimLengths[numReads]
    last is reduction for (rank+1)%worldSize
    pop first message part (buff+offset), push new message part + popped reduction

done when empty cycle is received


	 *
	 */


	/*
	 * ReadSet
	 *
	 * stream reads (no mmap or Read in-memory storage)
	 * open/write files by (inputFileIdx + rank) % inputFileSize
keep track of globalReadIds, or at least read-file boundaries
	 *
	 */

};

#endif /* DISTRIBUTED_FUNCTIONS_H_ */
