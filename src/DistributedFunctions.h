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
#include "ReadSet.h"
#include "ReadSelector.h"

#include "boost/optional.hpp"
#include <boost/thread/thread.hpp>
#include <vector>

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

ReadSet::ReadSetSizeType setGlobalReadSetOffset(mpi::communicator &world, ReadSet &store) {
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
	store.setGlobalOffset( globalReadSetOffset );
	return globalReadSetOffset;
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

protected:
	mpi::communicator world;

public:
	DistributedKmerSpectrum(mpi::communicator &_world, unsigned long buckets = 0, bool separateSingletons = true)
	: KS(buckets, separateSingletons), world(_world) {
	}
	~DistributedKmerSpectrum() {
	}
	DistributedKmerSpectrum &operator=(const DistributedKmerSpectrum &other) {
		*((KS*) this) = other;
		world = other.world;
		return *this;
	}

	template<typename D>
	MmapFile writeKmerMapMPI(D &kmerMap, std::string filepath) {
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
		// store the arrays in mmap blocked by mpi rank, to minimized mapped footprint
		for(int rank = 0 ; rank < world.size(); rank++) {
			for(IndexType i = rank; i < numBuckets; i+= world.size()) {
				totalSize += ourSizeCounts[i];
				offsetArray[i] = offset;
				offset += D::BucketType::sizeToStore( ourSizeCounts[i] );
			}
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
	 * StoreKmer / buildKmerSpectrumMPI
	 *
	 * each rank & thread listens for worldSize IRecv messages (rank, thread-tag)
	 * each rank scans ReadSet, each thread scheduled dynamically
	 * each thread builds own message buffer of kmers for each rank & thread numThreads * ( worldSize * numThreads ) buffers, sends when full.
	 * every send first receives & processes all pending messages, initiates new irecv, if zero sized, checkpoint increments
	 * done when all checkpoints have been received
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
			LOG_DEBUG(4, "StoreKmerMessage: " << readIdx << " " << readPos << " " << weight << " " << getKmer()->toFasta());
			kbufferCallback->getSpectrum().append(kbufferCallback->getDataPointer(), *getKmer(), weight, readIdx, readPos);
		}
	};

	void _buildKmerSpectrumMPI( mpi::communicator &world, ReadSet &store, bool isSolid = false ) {
		int numThreads = omp_get_max_threads();
		int rank = world.rank();
		int numParts = world.size();
		assert((numParts & (numParts-1)) == 0); // numParts must be a power of 2

		int messageSize = sizeof(StoreKmerMessageHeader) + KmerSizer::getTwoBitLength();

		long readSetSize = store.getSize();

		NumberType distributedThreadMask = numParts - 1;

		LOG_VERBOSE(1, "starting _buildSpectrum");

		SendStoreKmerMessageBuffer *sendBuffers[numThreads][numThreads];
		RecvStoreKmerMessageBuffer *recvBuffers[numThreads];

		LOG_DEBUG(2, "building spectrum using " << numThreads << " threads (" << omp_get_max_threads() << ")");

		ReadSetSizeType globalReadSetOffset = store.getGlobalOffset();
		assert( world.rank() == 0 ? (globalReadSetOffset == 0) : (globalReadSetOffset > 0) );

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
				LOG_DEBUG(3, "Read " << readIdx << " (" << globalReadIdx << ") " << kmers.size() );

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

	class MPIHistogram : public Histogram {
	public:
		typedef typename Histogram::HistogramElement HistogramElement;
		typedef typename Histogram::BucketsType BucketsType;
		MPIHistogram(unsigned int _zoomMax, double _logBase = 2.0) : Histogram(_zoomMax, _logBase) {}

		void reduce(mpi::communicator &world) {

			BucketsType copy(this->buckets);
			unsigned long *inbuffer = new unsigned long[this->buckets.size()];
			unsigned long *outbuffer = new unsigned long[this->buckets.size()];
			double *inbuffer_d = new double[this->buckets.size()];
			double *outbuffer_d = new double[this->buckets.size()];

			for(unsigned int i = 0; i < copy.size(); i++) {
				inbuffer[i] = this->buckets[i].visits;
			}
			mpi::all_reduce(world, inbuffer, copy.size(), outbuffer, std::plus<unsigned long>());
			for(unsigned int i = 0; i < copy.size(); i++) {
				copy[i].visits = outbuffer[i];
				inbuffer[i] = copy[i].visitedCount;
			}
			mpi::all_reduce(world, inbuffer, copy.size(), outbuffer, std::plus<unsigned long>());
			for(unsigned int i = 0; i < copy.size(); i++) {
				copy[i].visitedCount = outbuffer[i];
				inbuffer_d[i] = copy[i].visitedWeight;
			}
			mpi::all_reduce(world, inbuffer_d, copy.size(), outbuffer_d, std::plus<double>());
			for(unsigned int i = 0; i < copy.size(); i++) {
				copy[i].visitedWeight = outbuffer_d[i];
			}

			delete [] inbuffer;
			delete [] outbuffer;
			delete [] inbuffer_d;
			delete [] outbuffer_d;

			this->buckets.swap(copy);
		}
	};

	Kmernator::MmapFileVector buildKmerSpectrumMPI(ReadSet &store ) {
		bool isSolid = false; // not supported for references...

		_buildKmerSpectrumMPI(world, store, isSolid);

		MPIHistogram histogram(127);
		histogram.set(*this, false);
		LOG_DEBUG(2, "Individual raw histogram\n" << histogram.toString());

		histogram.reduce(world);
		if (world.rank() == 0)
			LOG_VERBOSE(1, "Collective raw histogram\n" << histogram.toString());
		world.barrier();

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
		ourSpectrum[0] = writeKmerMapMPI(this->weak, Options::getTmpDir() + "/weak-kmer-mmap");

		if (Options::getMinDepth() <= 1) {
			ourSpectrum[1] = writeKmerMapMPI(this->singleton, Options::getTmpDir() + "/singleton-kmer-mmap");
		}

		LOG_VERBOSE(1, "Finished merging partial spectrums" << std::endl << MemoryUtils::getMemoryUsage());

		return ourSpectrum;

	}

}; // DistributedKmerSpectrum

template<typename M>
class DistributedReadSelector : public ReadSelector<M>
{
public:
	typedef ReadSelector<M> RS;
	typedef M DataType;
	typedef typename RS::KMType KMType;
	typedef typename RS::ScoreType ScoreType;
	typedef typename RS::ReadSetSizeType ReadSetSizeType;
	typedef typename RS::ReadTrimType ReadTrimType;
	typedef typename RS::KA KA;
	typedef typename RS::ElementType ElementType;

	typedef Kmer::NumberType NumberType;

	typedef std::vector<ScoreType> KmerValueVector;
	typedef typename KmerValueVector::const_iterator KmerValueVectorIterator;
	typedef std::vector<ReadSetSizeType> ReadIdxVector;
	typedef typename ReadIdxVector::const_iterator ReadIdxVectorIterator;

protected:
	mpi::communicator _world;

public:
	DistributedReadSelector(mpi::communicator &world, const ReadSet &reads, const KMType &map)
		: RS(reads, map), _world(world) {}

	/*
	 * ReadTrim
	 * each rank IRecv FileReadyMessage tag=2*numThreads for each file to output from previous rank (except rank that first opens file)
     * each rank & thread Irecv RequestMessage tag = threadId
     * each rank & thread Irecv ResponseMessage tag = threadId + numThreads
     *
     * each rank & thread maintains linear buffer of read kmer values
     * RequestMessage process( callback ) adds message to response
     * ResponseMessage process( callback ) populates data in kmer values linear buffer
     *
     * periodically (by batch) communication flushes, checkpoints, and read trims get updated, files get written.
     *
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

	class RequestKmerMessageHeader;
	class RespondKmerMessageHeader;

	typedef MPIMessageBufferBase< RequestKmerMessageHeader > RequestKmerMessageBuffersBase;
	typedef MPIRecvMessageBuffer< RequestKmerMessageHeader > RecvRequestKmerMessageBufferBase;
	typedef MPISendMessageBuffer< RequestKmerMessageHeader > SendRequestKmerMessageBufferBase;

	typedef MPIMessageBufferBase< RespondKmerMessageHeader > RespondKmerMessageBuffersBase;
	typedef MPIRecvMessageBuffer< RespondKmerMessageHeader > RecvRespondKmerMessageBufferBase;
	typedef MPISendMessageBuffer< RespondKmerMessageHeader > SendRespondKmerMessageBufferBase;

	class RecvRespondKmerMessageBuffer : public RecvRespondKmerMessageBufferBase
	{
	public:
		KmerValueVector &_kmerValues;
		RecvRespondKmerMessageBuffer(mpi::communicator &world, int messageSize, int tag, KmerValueVector &kmerValues)
			: RecvRespondKmerMessageBufferBase(world, messageSize, tag), _kmerValues(kmerValues) {
		}
	};
	class SendRespondKmerMessageBuffer : public SendRespondKmerMessageBufferBase
	{
	public:
		SendRespondKmerMessageBuffer(mpi::communicator &world, int messageSize)
			: SendRespondKmerMessageBufferBase(world, messageSize) {
		}
	};

	class RecvRequestKmerMessageBuffer : public RecvRequestKmerMessageBufferBase
	{
	public:
		SendRespondKmerMessageBuffer &_sendResponse;
		RS &_readSelector;
		int _numThreads;
		RecvRequestKmerMessageBuffer(mpi::communicator &world, int messageSize, int tag, SendRespondKmerMessageBuffer &sendResponse, RS &readSelector, int numThreads)
			: RecvRequestKmerMessageBufferBase(world, messageSize, tag), _sendResponse(sendResponse), _readSelector(readSelector), _numThreads(numThreads) {
		}
	};
	class SendRequestKmerMessageBuffer : public SendRequestKmerMessageBufferBase
	{
	public:
		SendRequestKmerMessageBuffer(mpi::communicator &world,int messageSize)
			: SendRequestKmerMessageBufferBase(world, messageSize) {
		}
	};

	class RespondKmerMessageHeader {
	public:
		long requestId;
		ScoreType score;

		void set(long _requestId, ScoreType _score) {
			requestId = _requestId;
			score = _score;
		}
		// store response in kmer value vector
		void process(RecvRespondKmerMessageBufferBase *bufferCallback) {
			RecvRespondKmerMessageBuffer *kbufferCallback = (RecvRespondKmerMessageBuffer*) bufferCallback;
			LOG_DEBUG(4, "RespondKmerMessage: " << requestId << " " << score);
			kbufferCallback->_kmerValues[requestId] = score;
		}
	};

	class RequestKmerMessageHeader {
	public:
		long requestId;
		// Kmer is next bytes, dynamically determined by KmerSizer::getTwoBitLength()
		// kmer is least complement

		// THIS IS DANGEROUS unless allocated an extra Kmer!
		Kmer *getKmer() {
			return (Kmer*) (((char*)this)+sizeof(*this));
		}
		void set(long _requestId, const Kmer &_kmer) {
			requestId = _requestId;
			*(getKmer()) = _kmer;
		}
		// lookup kmer in map and build response message
		void process(RecvRequestKmerMessageBufferBase *bufferCallback) {
			RecvRequestKmerMessageBuffer *kbufferCallback = (RecvRequestKmerMessageBuffer*) bufferCallback;
			LOG_DEBUG(4, "RequestKmerMessage: " << requestId << " " << getKmer()->toFasta());

			ScoreType score = kbufferCallback->_readSelector.getValue( *getKmer() );
			kbufferCallback->_sendResponse.bufferMessage(kbufferCallback->getSource(), kbufferCallback->getTag() + kbufferCallback->_numThreads)->set( requestId, score );
		}
	};

	void _batchKmerLookup(const Read &read, SequenceLengthType trimLength, ReadIdxVector &readOffsetVector, KmerValueVector &batchBuffer, SendRequestKmerMessageBuffer &sendReq, int &thisThreadId, int &numThreads, NumberType &distributedThreadBitMask) {
		KA kmers = this->getKmersForRead(read);
		SequenceLengthType numKmers = 0;
		if (trimLength > KmerSizer::getSequenceLength()) {
			numKmers = trimLength - KmerSizer::getSequenceLength();
		}

		ReadSetSizeType offset = batchBuffer.size();
		readOffsetVector.push_back( offset );
		batchBuffer.insert(batchBuffer.end(), numKmers, ScoreType(0));
		for(SequenceLengthType kmerIdx = 0; kmerIdx < numKmers; kmerIdx++) {
			int localThreadId, distributedThreadId;
			this->_map.getThreadIds(kmers[kmerIdx], localThreadId, numThreads, distributedThreadId, distributedThreadBitMask);
			if (distributedThreadId == _world.rank()) {
				// handle this directly
				batchBuffer[offset + kmerIdx] =  this->getValue(kmers[kmerIdx]);
			} else {
				sendReq.bufferMessage( distributedThreadId, thisThreadId )->set( offset + kmerIdx, kmers[kmerIdx]) ;
			}
		}
	}
	void scoreAndTrimReadsMPI(ScoreType minimumKmerScore, int correctionAttempts = 0) {
		this->_trims.resize(this->_reads.getSize());
		bool useKmers = Options::getKmerSize() != 0;

		ReadSetSizeType readsSize = this->_reads.getSize();
		ReadSetSizeType batchSize = 10;
		ReadSetSizeType batchReadIdx = 0;
		int respondMessageSize = sizeof(RespondKmerMessageHeader);
		int requestMessageSize = sizeof(RequestKmerMessageHeader) + KmerSizer::getTwoBitLength();

		int numThreads = omp_get_max_threads();
		NumberType distributedThreadBitMask = _world.size() - 1;

		KmerValueVector batchBuffer[ numThreads ];
		ReadIdxVector readIndexBuffer[ numThreads ];
		ReadIdxVector readOffsetBuffer[ numThreads ];

		RecvRequestKmerMessageBuffer *recvReq[numThreads];
		SendRequestKmerMessageBuffer *sendReq[numThreads];

		RecvRespondKmerMessageBuffer *recvResp[numThreads];
		SendRespondKmerMessageBuffer *sendResp[numThreads];

		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();
			// initialize read/kmer buffers
			batchBuffer[threadId].reserve((batchSize * 200 / numThreads) + 1);
			readIndexBuffer[threadId].reserve((batchSize / numThreads) + 1);
			readOffsetBuffer[threadId].reserve((batchSize / numThreads) + 1);

			// initialize message buffers

			recvResp[threadId] = new RecvRespondKmerMessageBuffer(_world, respondMessageSize, threadId + numThreads, batchBuffer[threadId]);
			sendResp[threadId] = new SendRespondKmerMessageBuffer(_world, respondMessageSize);
			sendResp[threadId]->addCallback( *recvResp[threadId] );

			recvReq[threadId] = new RecvRequestKmerMessageBuffer(_world, requestMessageSize, threadId, *sendResp[threadId], *this, numThreads);
			sendReq[threadId] = new SendRequestKmerMessageBuffer(_world, requestMessageSize);
			sendReq[threadId]->addCallback( *recvReq[threadId] );

		}

		#pragma omp parallel num_threads(numThreads) firstprivate(batchReadIdx)
		while (batchReadIdx < readsSize) {
			int threadId = omp_get_thread_num();
			batchBuffer[threadId].resize(0);
			readIndexBuffer[threadId].resize(0);
			readOffsetBuffer[threadId].resize(0);

			LOG_DEBUG(1, "Starting batch for kmer lookups: " << batchReadIdx);

			for(ReadSetSizeType i = threadId ; i < batchSize ; i+=numThreads) {

				ReadSetSizeType readIdx = batchReadIdx + i;
				if (readIdx >= readsSize)
					continue;
				const Read &read = this->_reads.getRead(readIdx);
				if (read.isDiscarded()) {
					continue;
				}
				readIndexBuffer[threadId].push_back(readIdx);

				ReadTrimType &trim = this->_trims[readIdx];

				Sequence::BaseLocationVectorType markups = read.getMarkups();
				SequenceLengthType markupLength = TwoBitSequence::firstMarkupNorX(markups);
				if ( markupLength == 0 ) {
					trim.trimLength = read.getLength();
				} else {
					// trim at first N or X markup
					trim.trimLength = markupLength - 1;
				}
				if (useKmers) {
					_batchKmerLookup(read, trim.trimLength, readOffsetBuffer[threadId], batchBuffer[threadId], *sendReq[threadId], threadId, numThreads, distributedThreadBitMask);
				}

			}

			LOG_DEBUG(1, "kmer lookups finished, flushing communications");
			sendReq[threadId]->flushAllMessageBuffers(threadId);
			sendReq[threadId]->finalize(threadId);
			recvReq[threadId]->finalize();
			sendResp[threadId]->flushAllMessageBuffers(threadId+numThreads);
			sendResp[threadId]->finalize(threadId+numThreads);
			recvResp[threadId]->finalize();

			LOG_DEBUG(1, "assigning trim values");
			for(ReadSetSizeType i = 0; i < readIndexBuffer[threadId].size() ; i++ ) {
				ReadSetSizeType &readIdx = readIndexBuffer[threadId][i];

				KmerValueVectorIterator buffBegin = (batchBuffer[threadId].begin() + readOffsetBuffer[threadId][i] );
				KmerValueVectorIterator buffEnd = ( (i+1) < readIndexBuffer[threadId].size() ? (batchBuffer[threadId].begin() + readOffsetBuffer[threadId][i+1]) : batchBuffer[threadId].end() );
				ReadTrimType &trim = this->_trims[readIdx];

				trim.score = 0;
				trim.trimLength = 0;
				while (buffBegin != buffEnd) {
					ScoreType score = *(buffBegin++);
					if (score >= minimumKmerScore) {
						trim.trimLength++;
						trim.score += score;
					 } else
						 break;
				}
				double reportScore;
				if (trim.trimLength > 0) {
					// calculate average score (before adding kmer length)
					reportScore= trim.score /= (ScoreType) trim.trimLength;
					if (useKmers) {
						trim.trimLength += KmerSizer::getSequenceLength() - 1;
					}
				} else {
					// keep available so that pairs will be selected together
					trim.score = -1.0;
					reportScore = 0.0;
				}
				if (!trim.label.empty())
					trim.label += " ";
				trim.label += "Trim:" + boost::lexical_cast<std::string>( trim.trimLength ) + " Score:" + boost::lexical_cast<std::string>( reportScore );
			}

			LOG_DEBUG(1, "Finished assigning trim values: " << batchReadIdx);
			batchReadIdx += batchSize;

			// local & world threads are okay to start without sync
		}

		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();
			// initialize message buffers

			delete recvResp[threadId];
			delete sendResp[threadId];

			delete recvReq[threadId];
			delete sendReq[threadId];
		}

	}



};


/*
 * ReadSet
 *
 * stream reads (no mmap or Read in-memory storage)
 * open/write files by (inputFileIdx + rank) % inputFileSize
keep track of globalReadIds, or at least read-file boundaries
 *
 */



#endif /* DISTRIBUTED_FUNCTIONS_H_ */
