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
	class StoreKmerMessageHeader {
	public:
		ReadSetSizeType readIdx;
		PositionType readPos;
		WeightType weight; // weight is negative if kmer is rc of observed direction
		// Kmer is next bytes, dynamically determined by KmerSizer::getTwoBitLength()
		// kmer is least complement
	};

	class StoreKmerMessageBuffers {
		static const int STORE_KMER_MESSAGE_BUFFER_SIZE = 1 * 1024;
	public:
		typedef boost::optional< mpi::request > OptionalRequest;
		typedef boost::optional< mpi::status > OptionalStatus;
		mpi::communicator &_world;
		DistributedKmerSpectrum &_spectrum;
		int _numCheckpoints;
		int _messageSize;
		int _maxNumMessages;
		int _numThreads;
		char **_messageInBuffers;
		char **_messageOutBuffers;
		int  **_offsets;
		OptionalRequest **_requests;
		DataPointers **_pointers;
		StoreKmerMessageBuffers(mpi::communicator &world, DistributedKmerSpectrum &spectrum) : _world(world), _spectrum(spectrum), _numCheckpoints(0) {
			_messageSize = sizeof(StoreKmerMessageHeader) + KmerSizer::getTwoBitLength();
			_maxNumMessages = STORE_KMER_MESSAGE_BUFFER_SIZE / _messageSize;
			_numThreads = omp_get_max_threads();
			_messageInBuffers = new char * [ _numThreads ];
			_messageOutBuffers = new char * [ _numThreads ];
			_offsets = new int * [ _numThreads ];
			_requests = new OptionalRequest * [ _numThreads ];
			_pointers = new DataPointers * [ _numThreads ];

			#pragma omp parallel num_threads(_numThreads)
			{
				int threadId = omp_get_thread_num();
				_messageInBuffers[threadId] = new char[ STORE_KMER_MESSAGE_BUFFER_SIZE * _world.size() ];
				_messageOutBuffers[threadId] = new char[ STORE_KMER_MESSAGE_BUFFER_SIZE * _world.size() * _numThreads ];
				_offsets[threadId] = new int[ _world.size() * _numThreads ];

				LOG_DEBUG_MT(3, _world.rank() << ": " << threadId << ": allocated buffers " << (void*) _messageInBuffers[threadId] << " " << (void*) _messageOutBuffers[threadId] << " " << (void*) _offsets[threadId]);

				_requests[threadId] = new OptionalRequest[ _world.size() ];
				for(int destRank = 0; destRank < _world.size(); destRank++) {
					for(int destThreadId = 0; destThreadId < _numThreads; destThreadId++)
						_offsets[threadId][destRank*_numThreads+destThreadId] = 0;
					_requests[threadId][destRank] = _world.irecv(destRank, threadId, _messageInBuffers[threadId] + destRank*STORE_KMER_MESSAGE_BUFFER_SIZE, STORE_KMER_MESSAGE_BUFFER_SIZE);
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
		int processStoreKmerMessages(char *start, int size) {
			StoreKmerMessageHeader *msg, *end;
			msg = (StoreKmerMessageHeader*) start;
			end = (StoreKmerMessageHeader*) (start+size);
			int count = 0;
			while (msg != end) {
				msg = _processStoreKmerMessage(msg);
			}
			return count;
		}
		StoreKmerMessageHeader *_processStoreKmerMessage(StoreKmerMessageHeader *msg) {
			Kmer *kmer = (Kmer*) (((char*)msg)+sizeof(StoreKmerMessageHeader));
			LOG_DEBUG(4, _world.rank() << ": " << omp_get_thread_num() << ": message: " << msg->readIdx << " " << msg->readPos << " " << msg->weight << " " << kmer->toFasta());
			_spectrum.append(getDataPointer( omp_get_thread_num() ), *kmer, msg->weight, msg->readIdx, msg->readPos);
			return (StoreKmerMessageHeader*) (((char*) msg) + _messageSize);
		}
		void bufferStoreKmerMessage(int rankDest, int threadDest, ReadSetSizeType &readIdx, PositionType &readPos, WeightType &weight, Kmer &kmer) {
			int threadId = omp_get_thread_num();
			int tDest = _numThreads*rankDest + threadDest;

			char *buffStart = _messageOutBuffers[threadId] + STORE_KMER_MESSAGE_BUFFER_SIZE * tDest;
			LOG_DEBUG(4, _world.rank() << ": " << threadId << ": buffering message to " << rankDest << ", " << threadDest << " @ " << tDest << " " << (void*) (buffStart));
			int &offset = _offsets[threadId][tDest];
			if (offset + _messageSize > STORE_KMER_MESSAGE_BUFFER_SIZE)
				flushStoreKmerMessageBuffer(rankDest, threadDest);

			StoreKmerMessageHeader *msg = (StoreKmerMessageHeader*) (buffStart + offset);
			LOG_DEBUG(4, _world.rank() << ": " << threadId << ": buffering message to " << rankDest << ", " << threadDest << " @ " << tDest << " " << (void*) msg
					<< " : " << readIdx << " " << readPos << " " << weight << " " << kmer.toFasta());

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
			int &offset = _offsets[threadId][tDest];
			LOG_DEBUG(3, _world.rank() << ": " << threadId << ": flushing outgoing message buffer to " << rankDest << ", " << threadDest << " @ " << tDest << " size: " << offset);

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
		int receiveAllIncomingMessages(int threadDest) {
			int threadId = omp_get_thread_num();
			int messages = 0;
			assert( threadDest == threadId );


			LOG_DEBUG_MT(3, _world.rank() << ": " << threadId << ": receiving all messages for thread " << threadDest);
			while (true) {
				LOG_DEBUG_MT(4, _world.rank() << ": " << threadId << ": calling iprobe");

				int source, tag, size;
				OptionalStatus optionalStatus;
				OptionalRequest *optionalRequestPtr;
				for(int testSource = 0 ; testSource < _world.size(); testSource++) {
					optionalRequestPtr = _requests[threadDest] + testSource;
					if (!*optionalRequestPtr)
						continue;
					optionalStatus = (*optionalRequestPtr)->test();
					if (!!optionalStatus) {
						source = optionalStatus.get().source();
						tag = optionalStatus.get().tag();
						size = optionalStatus.get().count<char>().get();
						LOG_DEBUG_MT(2,  _world.rank() << ": " << threadId << ": received message from " << source << ", " << tag << " size: " << size << " request " << testSource << ", " << threadId);
						assert( testSource == source );
						assert( threadId == tag );
						break;
					}
				}
				if (!optionalStatus)
					break;

				messages++;

				LOG_DEBUG_MT(4, _world.rank() << ": " << threadId << ": receiving message from " << source << " size " << size << " t=" << threadId);

				char *buffer = _messageInBuffers[ threadDest ] + ( source * STORE_KMER_MESSAGE_BUFFER_SIZE );
				if (size == 0) {
					checkpoint();
					LOG_DEBUG_MT(2,_world.rank() << ": " << threadId << ": got checkpoint from " << source);
					*optionalRequestPtr = OptionalRequest();
				} else {
					processStoreKmerMessages(buffer, size);
					*optionalRequestPtr = _world.irecv(source, threadDest, buffer, STORE_KMER_MESSAGE_BUFFER_SIZE);
				}
			}
			LOG_DEBUG_MT(3, _world.rank() << ": " << threadId << ": processed messages: " << messages);
			return messages;
		}
		void checkpoint() {
			#pragma omp atomic
			_numCheckpoints++;
			LOG_DEBUG_MT(2, _world.rank() << ": " << omp_get_thread_num() << ": checkpoint received:" << _numCheckpoints);
		}
		inline int getNumCheckpoints() const {
			return _numCheckpoints;
		}

	};

	void _buildKmerSpectrumMPI( mpi::communicator &world, ReadSet &store, bool isSolid = false ) {
		int numThreads = omp_get_max_threads();
		int rank = world.rank();

		NumberType distributedThreadMask = world.size() - 1;

		LOG_VERBOSE_MT(1, world.rank() << ": starting _buildSpectrum");
		StoreKmerMessageBuffers buffers(world, *this);

		LOG_DEBUG_MT(2, world.rank() << ": allocated buffers");

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
		LOG_DEBUG_MT(2, world.rank() << ": reduced readSetSizes " << globalReadSetOffset);

		long myCount = 0;
		#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
		for(long readIdx = 0 ; readIdx < readSetSize; readIdx++)
		{
			int threadId = threadId;
			const Read &read = store.getRead( readIdx );

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
			if (++myCount & 1024*128 == 0)
				LOG_DEBUG_MT(1, world.rank() << ": " << threadId << ": processed " << myCount << " reads");
		}

		LOG_DEBUG_MT(1, world.rank() << ": finished generating kmers from reads");

		// flush buffers and wait
		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();
			// send all pending buffers
			buffers.flushAllStoreKmerMessageBuffers(true);
			// send zero message buffer as checkpoint signal to stop
			LOG_DEBUG_MT(2, world.rank() << ": " << threadId << ": sending stop message");
			buffers.flushAllStoreKmerMessageBuffers(true, true);
			buffers.checkpoint();
			LOG_DEBUG_MT(2, world.rank() << ": " << threadId << ": sent checkpoints: " << buffers.getNumCheckpoints());
			while (buffers.getNumCheckpoints() < numThreads * world.size() ) {
				buffers.receiveAllIncomingMessages(omp_get_thread_num());
				boost::this_thread::sleep( boost::posix_time::milliseconds(2) );
			}
		}
		LOG_DEBUG_MT(1, world.rank() << ": finished _buildKmerSpectrum");
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
		if (numParts == 1) {
			this->buildKmerSpectrum(store);
			return Kmernator::MmapFileVector();
		}

		assert((numParts & (numParts-1)) == 0); // numParts must be a power of 2

		_buildKmerSpectrumMPI(world, store, isSolid);

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
