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

#include "boost/optional.hpp"
#include <boost/thread/thread.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

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
	class StoreKmerMessageHeader {
	public:
		ReadSetSizeType readIdx;
		PositionType readPos;
		WeightType weight; // weight is negative if kmer is rc of observed direction
		// Kmer is next bytes, dynamically determined by KmerSizer::getTwoBitLength()
		// kmer is least complement
	};

	class StoreKmerMessageBuffers {
		static const int STORE_KMER_MESSAGE_BUFFER_SIZE = 16 * 1024;
	public:

		mpi::communicator &_world;
		DistributedKmerSpectrum &_spectrum;
		int _numCheckpoints;
		int _messageSize;
		int _maxNumMessages;
		int _numThreads;
		char **_messageInBuffers;
		char **_messageOutBuffers;
		int  **_offsets;
		mpi::request **_requests;
		DataPointers **_pointers;
		StoreKmerMessageBuffers(mpi::communicator &world, DistributedKmerSpectrum &spectrum) : _world(world), _spectrum(spectrum), _numCheckpoints(0) {
			_messageSize = sizeof(StoreKmerMessageHeader) + KmerSizer::getTwoBitLength();
			_maxNumMessages = STORE_KMER_MESSAGE_BUFFER_SIZE / _messageSize;
			_numThreads = omp_get_max_threads();
			_messageInBuffers = new char * [ _numThreads ];
			_messageOutBuffers = new char * [ _numThreads ];
			_offsets = new int * [ _numThreads ];
			_requests = new mpi::request * [ _numThreads ];
			_pointers = new DataPointers * [ _numThreads ];
			#pragma omp parallel num_threads(_numThreads)
			{
				int threadId = omp_get_thread_num();
				_messageInBuffers[threadId] = new char[ STORE_KMER_MESSAGE_BUFFER_SIZE * _world.size() ];
				_messageOutBuffers[threadId] = new char[ STORE_KMER_MESSAGE_BUFFER_SIZE * _world.size() * _numThreads ];
				_offsets[threadId] = new int[ _world.size() * _numThreads ];
				_requests[threadId] = new mpi::request[ _world.size() ];
				for(int destRank = 0; destRank < _world.size(); destRank++) {
					for(int destThreadId = 0; destThreadId < _numThreads; destThreadId++)
						_offsets[destRank*_numThreads+destThreadId] = 0;
					_requests[threadId][destRank] = _world.irecv(destRank, threadId, _messageInBuffers[threadId] + destRank*STORE_KMER_MESSAGE_BUFFER_SIZE, _messageSize);
				}
				_pointers[threadId] = new DataPointers(_spectrum);
			}
		}
		~StoreKmerMessageBuffers() {
			#pragma omp parallel num_threads(_numThreads)
			{
				int threadId = omp_get_thread_num();
				delete [] _messageInBuffers[threadId];
				delete [] _messageOutBuffers[threadId];
				delete [] _offsets[threadId];
				for(int destRank = 0; destRank < _world.size(); destRank++)
					_requests[threadId][destRank].cancel();
				delete [] _requests[threadId];
				delete _pointers[threadId];
			}
			delete [] _messageInBuffers;
			delete [] _messageOutBuffers;
			delete [] _offsets;
			delete [] _requests;
			delete [] _pointers;
		}
		inline DataPointers &getDataPointer(int threadId) {
			return *_pointers[threadId];
		}
		StoreKmerMessageHeader *_processStoreKmerMessage(StoreKmerMessageHeader *msg) {
			Kmer *kmer = (Kmer*) (msg+1);
			_spectrum.append(getDataPointer( omp_get_thread_num() ), *kmer, msg->weight, msg->readIdx, msg->readPos);
			return (StoreKmerMessageHeader*) ((char*) msg) + _messageSize;
		}
		void bufferStoreKmerMessage(int rankDest, int threadDest, ReadSetSizeType &readIdx, PositionType &readPos, WeightType &weight, Kmer &kmer) {
			int threadId = omp_get_thread_num();
			int tDest = _numThreads*rankDest + threadDest;
			char *buffStart = _messageOutBuffers[threadId] + STORE_KMER_MESSAGE_BUFFER_SIZE * tDest;
			int &offset = _offsets[threadId][ tDest ];
			if (offset + _messageSize > STORE_KMER_MESSAGE_BUFFER_SIZE)
				flushStoreKmerMessageBuffer(rankDest, threadDest);
			StoreKmerMessageHeader *msg = (StoreKmerMessageHeader*) (buffStart + offset);
			msg->readIdx = readIdx;
			msg->readPos = readPos;
			msg->weight = weight;
			*((Kmer*)(++msg)) = kmer;
			offset += _messageSize;
		}
		void flushStoreKmerMessageBuffer(int rankDest, int threadDest, bool checkReceive = true, bool sendZeroMessage = false) {
			int threadId = omp_get_thread_num();
			if (checkReceive)
				receiveAllIncomingMessages(threadId);
			int tDest = _numThreads*rankDest + threadDest;
			int &offset = _offsets[threadId][ tDest ];
			if (offset > 0 || sendZeroMessage) {
				char *buffStart = _messageOutBuffers[threadId] + STORE_KMER_MESSAGE_BUFFER_SIZE * tDest;
				_world.send(rankDest, threadDest, buffStart, offset);
			}
			offset = 0;
		}
		void flushAllStoreKmerMessageBuffers(bool checkReceive = true, bool sendZeroMessage = false) {
			int threadId = omp_get_thread_num();
			if (checkReceive)
				receiveAllIncomingMessages(threadId);
			for(int i = 0 ; i < _world.size(); i++)
				for(int j = 0; j < _numThreads; j++)
					flushStoreKmerMessageBuffer(i, j, checkReceive, sendZeroMessage);
		}
		int receiveAllIncomingMessages(int t) {
			int messages = 0;
			boost::optional< mpi::status > status;
			while (true) {
				status = _world.iprobe(mpi::any_source, t);
				if (!status)
					break;
				messages++;
				int source = status.get().source();
				int size = status.get().count<char>().get();
				mpi::request &request = _requests[t][source];
				assert( !!request.test() );
				assert( request.test().get().source() == source );
				assert( request.test().get().tag() == status.get().tag() );

				char *buffer = _messageInBuffers[ t ] + ( source * STORE_KMER_MESSAGE_BUFFER_SIZE );
				StoreKmerMessageHeader *msg = (StoreKmerMessageHeader*) buffer;
				StoreKmerMessageHeader *last = (StoreKmerMessageHeader*) (buffer+size);
				if (size == 0) {
					checkpoint();
				} else {
					while (msg != last)
						msg = _processStoreKmerMessage(msg);
					request = _world.irecv(source, t, buffer, _messageSize);
				}
			}
			return messages;
		}
		void checkpoint() {
			#pragma omp atomic
			_numCheckpoints++;
		}
		inline int getNumCheckpoints() const {
			return _numCheckpoints;
		}

	};

	void _buildKmerSpectrum( mpi::communicator &world, ReadSet &store, bool isSolid = false ) {
		int numThreads = omp_get_max_threads();
		int rank = world.rank();

		NumberType distributedThreadMask = world.size() - 1;

		StoreKmerMessageBuffers buffers(world, *this);

		// share ReadSet sizes for globally unique readIdx calculations
		long readSetSize = store.getSize();
		long readSetSizes[ world.size() ];
		for(int i = 0; i < world.size(); i++)
			readSetSizes[i] = i == world.rank() ? readSetSize : 0;
		mpi::all_reduce(world, (long*) readSetSizes, world.size(), (long*) readSetSizes, mpi::maximum<long>());
		long globalReadSetOffset = 0;
		for(int i = 0; i < world.rank(); i++)
			globalReadSetOffset += readSetSizes[i];

		#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
		for(long readIdx = 0 ; readIdx < readSetSize; readIdx++)
		{
			const Read &read = store.getRead( readIdx );
			int threadId = omp_get_thread_num();

			if (read.isDiscarded())
				continue;
			KmerWeights kmers = KmerReadUtils::buildWeightedKmers(read, true, true);
			ReadSetSizeType globalReadIdx = readIdx + globalReadSetOffset;
			for (PositionType readPos = 0 ; readPos < kmers.size(); readPos++) {
				int rankDest, threadDest;
				WeightType weight = kmers.valueAt(readPos);
				this->solid.getThreadIds(kmers[readPos], threadDest, numThreads, rankDest, distributedThreadMask);
				if (rankDest == rank && threadDest == threadId)
					this->append(buffers.getDataPointer( threadId ), kmers[readPos], weight, globalReadIdx, readPos);
				else
					buffers.bufferStoreKmerMessage(rankDest, threadDest, globalReadIdx, readPos, weight, kmers[readPos]);
			}
		}

		// flush buffers and wait
		#pragma omp parallel num_threads(numThreads)
		{
			// send all pending buffers
			buffers.flushAllStoreKmerMessageBuffers(true);
			// send zero message buffer as checkpoint signal to stop
			buffers.flushAllStoreKmerMessageBuffers(true, true);
			buffers.checkpoint();
			while (buffers.getNumCheckpoints() < numThreads * world.size() ) {
				buffers.receiveAllIncomingMessages(omp_get_thread_num());
				boost::this_thread::sleep( boost::posix_time::milliseconds(2) );
			}
		}
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
