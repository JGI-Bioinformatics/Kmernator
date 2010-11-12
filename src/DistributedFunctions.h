//
// Kmernator/src/DistributedFunctions.h
//
// Author: Rob Egan
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
#include "Log.h"

#include "boost/optional.hpp"
#include <boost/thread/thread.hpp>
#include <vector>

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

void reduceOMPThreads(mpi::communicator &world) {
	Options::validateOMPThreads();
	int numThreads = Options::getMaxThreads();
	numThreads = all_reduce(world, numThreads, mpi::minimum<int>());
	omp_set_num_threads(numThreads);
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "set OpenMP threads to " << numThreads);
}

void validateMPIWorld(mpi::communicator &world, int threadSupport) {
	if (Options::getGatheredLogs())
		Logger::setWorld(&world, Options::getDebug() >= 2);
	if ((world.size() & (world.size()-1)) != 0) {
		throw std::invalid_argument(
				(std::string("The number of mpi processes must be a power-of-two.\nPlease adjust the number of processes. ")
		+ boost::lexical_cast<std::string>(world.size())).c_str());
	}
	if (threadSupport != MPI_THREAD_MULTIPLE) {
		LOG_WARN(1, "Your version of MPI does not support MPI_THREAD_MULTIPLE, reducing OpenMP threads to 1")
		omp_set_num_threads(1);
	}
	reduceOMPThreads(world);
}


ReadSet::ReadSetSizeType setGlobalReadSetOffset(mpi::communicator &world, ReadSet &store) {
	// share ReadSet sizes for globally unique readIdx calculations
	long readSetSize = store.getSize();
	long readSetSizesInput[ world.size() ];
	long readSetSizes[ world.size() ];
	for(int i = 0; i < world.size(); i++)
		readSetSizesInput[i] = i == world.rank() ? readSetSize : 0;
	LOG_DEBUG(2, "setGlobalReadSetOffset: all_reduce:" << readSetSize);
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
	typedef typename KS::WeakElementType WeakElementType;
	typedef typename KS::WeakBucketType WeakBucketType;

protected:
	mpi::communicator world;
	NumberType distributedThreadMask;

public:
	DistributedKmerSpectrum(mpi::communicator &_world, unsigned long buckets = 0, bool separateSingletons = true)
	: KS(buckets, separateSingletons), world(_world) {
		NumberType numParts = world.size();
		distributedThreadMask = numParts - 1;
		assert((numParts & (distributedThreadMask)) == 0); // numParts must be a power of 2
	}
	~DistributedKmerSpectrum() {
	}
	DistributedKmerSpectrum &operator=(const DistributedKmerSpectrum &other) {
		*((KS*) this) = other;
		world = other.world;
		return *this;
	}

	template<typename D>
	MmapFile writeKmerMap(D &kmerMap, std::string filepath) {
		NumberType numBuckets = kmerMap.getNumBuckets();
		NumberType *mySizeCounts = new NumberType[ numBuckets ];
		for(IndexType i = 0 ; i < numBuckets; i++)
			mySizeCounts[i] = kmerMap.getBucketByIdx(i).size();
		NumberType *ourSizeCounts = new NumberType[ numBuckets ];
		LOG_DEBUG(2, "writeKmerMap(): all_reduce() " << numBuckets << ", " << ourSizeCounts[0] << " " << ourSizeCounts[1] << " " << ourSizeCounts[2] << " " << ourSizeCounts[3]);
		all_reduce(world, mySizeCounts, numBuckets, ourSizeCounts, std::plus<NumberType>());
		LOG_DEBUG(3, "Distributed Size Counts: " << ourSizeCounts[0] << " " << ourSizeCounts[1] << " " << ourSizeCounts[2] << " " << ourSizeCounts[3]);

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
			msync(mmap.data(), (char*) numbers - mmap.data(), MS_SYNC);
			LOG_DEBUG(2, "writeKmerMap(): barrier");
			world.barrier();
		} else {
			LOG_DEBUG(2, "writeKmerMap(): barrier");
			world.barrier();
			mmap = MmapFile(filepath, std::ios_base::in | std::ios_base::out, totalMmapSize);
			NumberType *numbers = (NumberType *) mmap.data();
			LOG_DEBUG_OPTIONAL(1, true, *numbers << " vs " << numBuckets);
			assert(*(numbers++) == numBuckets);
			assert(*(numbers++) == kmerMap.getBucketMask());

			for(IndexType i = 0 ; i < numBuckets; i++) {
				assert(*(numbers++) == offsetArray[i]);
			}
		}

		// store our part of mmap (interleaved DMP)
		for(IndexType i = 0 ; i < numBuckets; i++) {
			long size = kmerMap.getBucketByIdx(i).size();
			if (size > 0) {
				kmerMap.getBucketByIdx(i).store(mmap.data() + offsetArray[i]);
				msync(mmap.data() + offsetArray[i], size, MS_ASYNC);
			}
		}

		delete [] mySizeCounts;	// aka offsetArray
		LOG_DEBUG(2, "writeKmerMap(): barrier2");
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
		bool _isSolid;

	public:
		RecvStoreKmerMessageBuffer(mpi::communicator &world, DistributedKmerSpectrum &spectrum, int messageSize, int tag, bool isSolid = false)
			: RecvStoreKmerMessageBufferBase(world, messageSize, tag), _spectrum(spectrum), _pointers(spectrum), _isSolid(isSolid) {
		}
		inline DataPointers &getDataPointer() {
			return _pointers;
		}
		inline DistributedKmerSpectrum &getSpectrum() {
			return _spectrum;
		}
		inline bool isSolid() const {
			return _isSolid;
		}
	};
	class SendStoreKmerMessageBuffer : public SendStoreKmerMessageBufferBase
	{
		DistributedKmerSpectrum &_spectrum;

	public:
		SendStoreKmerMessageBuffer(mpi::communicator &world, DistributedKmerSpectrum &spectrum, int messageSize)
			: SendStoreKmerMessageBufferBase(world, messageSize), _spectrum(spectrum) {
		}
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
			kbufferCallback->getSpectrum().append(kbufferCallback->getDataPointer(), *getKmer(), weight, readIdx, readPos, kbufferCallback->isSolid());
		}
	};

	void _buildKmerSpectrumMPI(ReadSet &store, bool isSolid) {
		int numThreads = omp_get_max_threads();
		int rank = world.rank();

		int messageSize = sizeof(StoreKmerMessageHeader) + KmerSizer::getTwoBitLength();

		long readSetSize = store.getSize();


		LOG_VERBOSE(2, "starting _buildSpectrumMPI");

		SendStoreKmerMessageBuffer *sendBuffers[numThreads][numThreads];
		RecvStoreKmerMessageBuffer *recvBuffers[numThreads];

		LOG_DEBUG(2, "building spectrum using " << numThreads << " threads (" << omp_get_max_threads() << ")");

		ReadSetSizeType globalReadSetOffset = store.getGlobalOffset();
		assert( world.rank() == 0 ? (globalReadSetOffset == 0) : (globalReadSetOffset > 0) );

		std::stringstream ss;
		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();
			LOG_DEBUG(3, "allocating buffers for thread");
			recvBuffers[threadId] = new RecvStoreKmerMessageBuffer(world, *this, messageSize, threadId, isSolid);
			for(int recvThread = 0 ; recvThread < numThreads; recvThread++) {
				sendBuffers[threadId][recvThread] = new SendStoreKmerMessageBuffer(world, *this, messageSize);
				sendBuffers[threadId][recvThread]->addReceiveAllCallback( recvBuffers[threadId] );
			}

			#pragma omp master
			{
				LOG_DEBUG(2, "message buffers ready");
				world.barrier();
			}
			#pragma omp barrier

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
					if ( TrackingData::isDiscard( (weight<0.0) ? 0.0-weight : weight ) )  {
						LOG_DEBUG(4, "discarded kmer " << readIdx << "@" << readPos << " " << weight << " " << kmers[readPos].toFasta());
					} else {

						this->solid.getThreadIds(kmers[readPos], threadDest, numThreads, rankDest, distributedThreadMask);

						if (rankDest == rank && threadDest == threadId) {
							this->append(recvBuffers[ threadId ]->getDataPointer(), kmers[readPos], weight, globalReadIdx, readPos, isSolid);
						} else {
							sendBuffers[ threadId ][ threadDest ]->bufferMessage(rankDest, threadDest)->set(globalReadIdx, readPos, weight, kmers[readPos]);
						}
					}
				}

				if (threadId == 0 && readIdx % 1000000 == 0) {
					if (world.rank() == 0) {
						LOG_VERBOSE_OPTIONAL(1, true, "distributed processing " << (readIdx * world.size()) << " reads");
					} else {
						LOG_DEBUG(2, "local processed: " << readIdx << " reads");
					}

				}
			}

			LOG_DEBUG(2, "finished generating kmers from reads");

			// receiving any new messages
			recvBuffers[threadId]->receiveAllIncomingMessages();

			LOG_DEBUG(3, "sending final flush");
			// send all pending buffers
			for(int destThread = 0; destThread < numThreads; destThread++) {
				sendBuffers[threadId][destThread]->flushAllMessageBuffers(destThread);
			}
			recvBuffers[threadId]->receiveAllIncomingMessages();


			LOG_DEBUG(3, "sending final messages")

			for(int destThread = 0; destThread < numThreads; destThread++) {
				sendBuffers[threadId][destThread]->finalize(destThread);
			}
			LOG_DEBUG(3, "receiving final messages");
			recvBuffers[threadId]->finalize(numThreads);

			#pragma omp critical
			{
				if (Log::isDebug(2)) {
					for(int destThread = 0; destThread < numThreads; destThread++)
						ss << "R" << world.rank() << ": sendBuffers["<<threadId<<"]["<<destThread<<"] sent " << sendBuffers[threadId][destThread]->getNumDeliveries() << "/" << sendBuffers[threadId][destThread]->getNumMessages() << std::endl;
					ss  << "R" << world.rank() << ":recvBuffers["<<threadId<<"] received " << recvBuffers[threadId]->getNumDeliveries() << "/" << recvBuffers[threadId]->getNumMessages() << std::endl;
				}
			}
			for(int destThread = 0; destThread < numThreads; destThread++) {
				delete sendBuffers[threadId][destThread];
			}
			delete recvBuffers[threadId];

		} // omp parallel

		std::string s = ss.str();
		LOG_DEBUG(1, s);
		LOG_DEBUG(2, "_buildKmerSpectrumMPI() final barrier");
		world.barrier();
		LOG_DEBUG(1, "finished _buildKmerSpectrumMPI");
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

			LOG_DEBUG(2, "MPIHistogram::reduce() all_reduce");
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

	void buildKmerSpectrum( ReadSet &store ) {
		return this->buildKmerSpectrum(store, false);
	}

	void buildKmerSpectrum(ReadSet &store, bool isSolid) {

		_buildKmerSpectrumMPI(store, isSolid);

		std::string hist = getHistogram(isSolid);
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective raw histogram\n" << hist);
		LOG_DEBUG(2, "buildKmerSpectrum() barrier");
		world.barrier();

		// purge low counts
		if (Options::getMinDepth() > 1) {
			LOG_VERBOSE(2, "Clearing memory from singletons: " << this->singleton.size() );
			LOG_DEBUG(2, MemoryUtils::getMemoryUsage());
			this->singleton.clear();
		}
		if (Options::getMinDepth() > 2) {
			LOG_VERBOSE(2, "Purging low count kmers (< " << Options::getMinDepth() << ")");
			LOG_DEBUG(2, MemoryUtils::getMemoryUsage());
			this->purgeMinDepth(Options::getMinDepth());
		}
		LOG_DEBUG(3, this->weak.toString());

	}

	virtual std::string getHistogram(bool solidOnly = false) {
		MPIHistogram histogram(127);
		histogram.set(*this, solidOnly);
		LOG_DEBUG(2, "Individual histogram\n" << histogram.toString());

		histogram.reduce(world);

		return histogram.toString();
	}


	Kmernator::MmapFileVector writeKmerMaps(string fileprefix = Options::getOutputFile()) {
		// communicate sizes and allocate permanent file
		LOG_VERBOSE(2, "Merging partial spectrums" );
		LOG_DEBUG(3, MemoryUtils::getMemoryUsage() );
		Kmernator::MmapFileVector ourSpectrum(3);
		if (this->hasSolids){
			ourSpectrum[0] = this->writeKmerMap(this->solid, fileprefix + "-kmer-mmap");
		}
		ourSpectrum[1] = this->writeKmerMap(this->weak, fileprefix + "-kmer-mmap");

		if (Options::getMinDepth() <= 1 && this->hasSingletons) {
			ourSpectrum[2] = this->writeKmerMap(this->singleton, fileprefix + "-singleton-kmer-mmap");
		}

		LOG_DEBUG(1, "Finished merging partial spectrums" << std::endl << MemoryUtils::getMemoryUsage());

		return ourSpectrum;
	}

	/*
	 * PurgeVariantKmer / buildKmerSpectrumMPI
	 *
	 * message (for a specific rank & thread 0):
	 * threshold + kmer
	 *
	 */
	class PurgeVariantKmerMessageHeader;

	typedef MPIMessageBufferBase< PurgeVariantKmerMessageHeader > PurgeVariantKmerMessageBuffersBase;
	typedef MPIRecvMessageBuffer< PurgeVariantKmerMessageHeader > RecvPurgeVariantKmerMessageBufferBase;
	typedef MPISendMessageBuffer< PurgeVariantKmerMessageHeader > SendPurgeVariantKmerMessageBufferBase;

	class RecvPurgeVariantKmerMessageBuffer : public RecvPurgeVariantKmerMessageBufferBase
	{
		DistributedKmerSpectrum &_spectrum;
		DataPointers _pointers;
		double _variantSigmas, _minDepth;

	public:
		RecvPurgeVariantKmerMessageBuffer(mpi::communicator &world, DistributedKmerSpectrum &spectrum, int messageSize, int srcTag, double variantSigmas, double minDepth)
			: RecvPurgeVariantKmerMessageBufferBase(world, messageSize, srcTag), _spectrum(spectrum), _pointers(spectrum), _variantSigmas(variantSigmas), _minDepth(minDepth) {
		}
		inline DataPointers &getDataPointer() {
			return _pointers;
		}
		inline DistributedKmerSpectrum &getSpectrum() {
			return _spectrum;
		}
		inline double getVariantSigmas() {
			return _variantSigmas;
		}
		inline double getMinDepth() {
			return _minDepth;
		}
	};
	class SendPurgeVariantKmerMessageBuffer : public SendPurgeVariantKmerMessageBufferBase
	{
		DistributedKmerSpectrum &_spectrum;

	public:
		SendPurgeVariantKmerMessageBuffer(mpi::communicator &world, DistributedKmerSpectrum &spectrum, int messageSize)
			: SendPurgeVariantKmerMessageBufferBase(world, messageSize), _spectrum(spectrum) {
		}
		inline DistributedKmerSpectrum &getSpectrum() {
			return _spectrum;
		}
	};

	class PurgeVariantKmerMessageHeader {
	public:
		float threshold;
		// Kmer is next bytes, dynamically determined by KmerSizer::getTwoBitLength()
		// kmer is least complement

		// THIS IS DANGEROUS unless allocated an extra Kmer!
		Kmer *getKmer() {
			return (Kmer*) (((char*)this)+sizeof(*this));
		}
		void set(float _threshold, const Kmer &_kmer) {
			threshold = _threshold;
			*(getKmer()) = _kmer;
		}
		void process(RecvPurgeVariantKmerMessageBufferBase *bufferCallback) {
			RecvPurgeVariantKmerMessageBuffer *kbufferCallback = (RecvPurgeVariantKmerMessageBuffer*) bufferCallback;
			LOG_DEBUG(4, "PurgeVariantKmerMessage: " << threshold << " " << getKmer()->toFasta());
			DistributedKmerSpectrum &spectrum = kbufferCallback->getSpectrum();
			DataPointers &pointers = kbufferCallback->getDataPointer();
			double dummy;
			bool purged = spectrum._setPurgeVariant(pointers, *getKmer(), threshold, dummy);

			if (purged)
				spectrum.variantWasPurged();
		}
	};
	void variantWasPurged(long count = 1) {
		#pragma omp atomic
		_purgedVariants += count;
	}
private:
	std::vector<SendPurgeVariantKmerMessageBuffer*> sendPurgeVariant;
	std::vector<RecvPurgeVariantKmerMessageBuffer*> recvPurgeVariant;
	long _purgedVariants;
	double variantSigmas;
	double minDepth;
	int _variantNumThreads;

	void _preVariants(double variantSigmas, double minDepth) {
		_purgedVariants = 0;
		long messageSize = sizeof(PurgeVariantKmerMessageHeader) + KmerSizer::getByteSize();

		int &numThreads = _variantNumThreads = omp_get_max_threads();
		sendPurgeVariant.resize(numThreads*numThreads, NULL);
		recvPurgeVariant.resize(numThreads, NULL);
		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();
			recvPurgeVariant[threadId] = new RecvPurgeVariantKmerMessageBuffer(world, *this, messageSize, threadId, variantSigmas, minDepth);
			for (int t = 0 ; t < numThreads; t++) {
				sendPurgeVariant[threadId*numThreads+t] = new SendPurgeVariantKmerMessageBuffer(world, *this, messageSize);
				sendPurgeVariant[threadId*numThreads+t]->addReceiveAllCallback( recvPurgeVariant[threadId] );
			}
		}
		LOG_DEBUG(2, "_preVariants(): barrier");
		world.barrier();
	}

	long _postVariants(long purgedKmers) {
		_purgedVariants += this->KS::_postVariants(purgedKmers);
		int &numThreads = _variantNumThreads;
		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();

			for(int t = 0 ; t < numThreads; t++) {
				delete sendPurgeVariant[threadId*numThreads+t];
			}
			delete recvPurgeVariant[threadId];
		}
		sendPurgeVariant.clear();
		recvPurgeVariant.clear();

		long allPurged;
		LOG_DEBUG(2, "_postVariants(): all_reduce");
		mpi::all_reduce(world, _purgedVariants, allPurged, std::plus<long>());
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Distributed Purged " << allPurged << " kmer-variants");

		return allPurged;
	}
	void _variantThreadSync(long processed, long remaining, double maxDepth) {
		// call parent
		this->KS::_variantThreadSync(processed, remaining, maxDepth);

		int &numThreads = _variantNumThreads;
		int threadId = omp_get_thread_num();
		recvPurgeVariant[threadId]->receiveAllIncomingMessages();
		for (int t = 0 ; t < numThreads; t++) {
			sendPurgeVariant[threadId*numThreads+t]->flushAllMessageBuffers(t);
			LOG_DEBUG(3, "Flushed: " << t);
		}
		recvPurgeVariant[threadId]->receiveAllIncomingMessages();
		LOG_DEBUG(3, "Received all incoming");
		for (int t = 0 ; t < numThreads; t++) {
			sendPurgeVariant[threadId*numThreads+t]->finalize(t);
		}
		recvPurgeVariant[threadId]->finalize(numThreads);
		LOG_DEBUG(3, "_variantThreadSync() finished:" << maxDepth <<std::endl);
	}
	long _variantBatchSync(long remaining, long purgedKmers, double maxDepth, double threshold) {
		// call parent
		remaining = this->KS::_variantBatchSync(remaining, purgedKmers, maxDepth, threshold);

		long allRemaining;
		LOG_DEBUG(2, "_variantBatchSync(): all_reduce + barrier");
		mpi::all_reduce(world, remaining, allRemaining, std::plus<long>());
		if (threshold > 0.0) {
			long allPurged;
			mpi::reduce(world, purgedKmers+_purgedVariants, allPurged, std::plus<long>(), 0);
			LOG_VERBOSE_OPTIONAL(2, world.rank() == 0, "Distributed Purged " << allPurged << " variants below: " << maxDepth << " / " <<  threshold << ".  Remaining: " << allRemaining);
		}
		world.barrier();
		return allRemaining;
	}
	// recursively purge kmers within editdistance
	long _purgeVariants(DataPointers &pointers, const Kmer &kmer, WeakBucketType &variants, double threshold, short editDistance) {
		int rank = world.rank();
		int threadId = omp_get_thread_num();
		int &numThreads = _variantNumThreads;
		if (editDistance == 0)
			return 0;

		WeakBucketType::permuteBases(kmer, variants, editDistance, true);

		for(SequenceLengthType i = 0 ; i < variants.size(); i++) {
			int rankDest, threadDest;
			Kmer &varKmer = variants[i];
			this->solid.getThreadIds(varKmer, threadDest, 1, rankDest, distributedThreadMask);

			if (rankDest == rank) {
				double dummy;
				if (this->_setPurgeVariant(pointers, varKmer, threshold, dummy))
					variantWasPurged();
			} else {
				sendPurgeVariant[threadId*numThreads+threadDest]->bufferMessage(rankDest, threadDest)->set(threshold, varKmer);
			}
		}
		return 0;
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
		: RS(reads, map), _world(world) {
		LOG_DEBUG(3, this->_map.toString());
	}

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
		RS *_readSelector;
		int _numThreads;
		RecvRequestKmerMessageBuffer(mpi::communicator &world, int messageSize, int tag, SendRespondKmerMessageBuffer &sendResponse, RS &readSelector, int numThreads)
			: RecvRequestKmerMessageBufferBase(world, messageSize, tag), _sendResponse(sendResponse), _readSelector(&readSelector), _numThreads(numThreads) {
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
			int destSource = kbufferCallback->getRecvSource();
			int destTag = kbufferCallback->getRecvTag() + kbufferCallback->_numThreads;
			assert(destTag == omp_get_thread_num() + kbufferCallback->_numThreads);

			ScoreType score = kbufferCallback->_readSelector->getValue( *getKmer() );
			kbufferCallback->_sendResponse.bufferMessage(destSource, destTag)->set( requestId, score );
			LOG_DEBUG(4, "RequestKmerMessage: " << getKmer()->toFasta() << " " << requestId << " " << score);
		}
	};

	int _batchKmerLookup(const Read &read, SequenceLengthType markupLength, ReadSetSizeType offset, KmerValueVector &batchBuffer, SendRequestKmerMessageBuffer &sendReq, int &thisThreadId, int &numThreads, NumberType &distributedThreadBitMask) {
		KA kmers = this->getKmersForRead(read);
		SequenceLengthType numKmers = kmers.size();

		this->_setNumKmers(markupLength, numKmers);
		assert(numKmers <= kmers.size());
		assert(offset == batchBuffer.size());

		if (numKmers > 0)
			batchBuffer.insert(batchBuffer.end(), numKmers, ScoreType(0));

		for(SequenceLengthType kmerIdx = 0; kmerIdx < numKmers; kmerIdx++) {
			int localThreadId, distributedThreadId;
			this->_map.getThreadIds(kmers[kmerIdx], localThreadId, numThreads, distributedThreadId, distributedThreadBitMask);
			if (distributedThreadId == _world.rank()) {
				// handle this directly
				batchBuffer[offset + kmerIdx] = this->getValue(kmers[kmerIdx]);
			} else {
				sendReq.bufferMessage( distributedThreadId, thisThreadId )->set( offset + kmerIdx, kmers[kmerIdx]) ;
			}
		}
		return numKmers;
	}
	void scoreAndTrimReads(ScoreType minimumKmerScore, int correctionAttempts = 0) {
		this->_trims.resize(this->_reads.getSize());
		bool useKmers = Options::getKmerSize() != 0;

		int numThreads = omp_get_max_threads();
		NumberType distributedThreadBitMask = _world.size() - 1;

		ReadSetSizeType readsSize = this->_reads.getSize();
		ReadSetSizeType batchSize = Options::getBatchSize();
		ReadSetSizeType batchReadIdx = 0;
		int maxKmers = std::max(1, (int) this->_reads.getMaxSequenceLength() - (int) KmerSizer::getSequenceLength());
		int reserveBB = ((batchSize * maxKmers) / numThreads) + 1;
		int reserveOffsets = (batchSize / numThreads) + 1;

		int respondMessageSize = sizeof(RespondKmerMessageHeader);
		int requestMessageSize = sizeof(RequestKmerMessageHeader) + KmerSizer::getTwoBitLength();

		LOG_DEBUG(1, "Starting scoreAndTrimReadsMPI - trimming: " << readsSize << " using " << numThreads << " threads");

		KmerValueVector batchBuffer[ numThreads ];
		ReadIdxVector readIndexBuffer[ numThreads ];
		ReadIdxVector readOffsetBuffer[ numThreads ];

		RecvRequestKmerMessageBuffer *recvReq[numThreads];
		SendRequestKmerMessageBuffer *sendReq[numThreads];

		RecvRespondKmerMessageBuffer *recvResp[numThreads];
		SendRespondKmerMessageBuffer *sendResp[numThreads];

		LOG_DEBUG(2, "scoreAndTrimReads(): all_reduce: " << readsSize);
		ReadSetSizeType mostReads = mpi::all_reduce(_world, readsSize, mpi::maximum<ReadSetSizeType>());
		LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "Largest number of reads to batch: " << mostReads);

		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();

			// initialize message buffers

			recvResp[threadId] = new RecvRespondKmerMessageBuffer(_world, respondMessageSize, threadId + numThreads, batchBuffer[threadId]);
			sendResp[threadId] = new SendRespondKmerMessageBuffer(_world, respondMessageSize);
			sendResp[threadId]->addReceiveAllCallback( recvResp[threadId] );

			recvReq[threadId] = new RecvRequestKmerMessageBuffer(_world, requestMessageSize, threadId, *sendResp[threadId], *this, numThreads);
			sendReq[threadId] = new SendRequestKmerMessageBuffer(_world, requestMessageSize);
			sendReq[threadId]->addReceiveAllCallback( recvReq[threadId] );
			sendReq[threadId]->addReceiveAllCallback( recvResp[threadId] );

			recvReq[threadId]->addFlushAllCallback( sendResp[threadId], threadId + numThreads);
		}

		LOG_DEBUG(2, "scoreAndTrimReads(): barrier. message buffers ready");
		_world.barrier();

		#pragma omp parallel num_threads(numThreads) firstprivate(batchReadIdx)
		while (batchReadIdx < mostReads) {
			int threadId = omp_get_thread_num();

			// initialize read/kmer buffers
			batchBuffer[threadId].resize(0);
			readIndexBuffer[threadId].resize(0);
			readOffsetBuffer[threadId].resize(0);
			batchBuffer[threadId].reserve(reserveBB);
			readIndexBuffer[threadId].reserve(reserveOffsets);
			readOffsetBuffer[threadId].reserve(reserveOffsets);

			LOG_VERBOSE_OPTIONAL(1, _world.rank() == 0 && threadId == 0, "trimming batch: " << batchReadIdx * _world.size());

			LOG_DEBUG(3, "Starting batch for kmer lookups: " << batchReadIdx);

			for(ReadSetSizeType i = threadId ; i < batchSize ; i+=numThreads) {

				ReadSetSizeType readIdx = batchReadIdx + i;
				if (readIdx >= readsSize)
					continue;
				const Read &read = this->_reads.getRead(readIdx);
				if (read.isDiscarded()) {
					continue;
				}
				ReadSetSizeType offset = batchBuffer[threadId].size();
				readOffsetBuffer[threadId].push_back( offset );
				readIndexBuffer[threadId].push_back(readIdx);

				Sequence::BaseLocationVectorType markups = read.getMarkups();
				SequenceLengthType markupLength = TwoBitSequence::firstMarkupNorX(markups);

				if (useKmers) {
					_batchKmerLookup(read, markupLength, offset, batchBuffer[threadId], *sendReq[threadId], threadId, numThreads, distributedThreadBitMask);
				} else {
					this->trimReadByMarkupLength(read, this->_trims[readIdx], markupLength);
				}

			}
			assert(readOffsetBuffer[threadId].size() == readIndexBuffer[threadId].size());

			LOG_DEBUG(3, "Starting communication sync for kmer lookups: " << batchReadIdx);

			sendReq[threadId]->flushAllMessageBuffers(threadId);
			while (sendReq[threadId]->getNumMessages() != recvResp[threadId]->getNumMessages()) {
				recvReq[threadId]->receiveAllIncomingMessages();
				sendResp[threadId]->flushAllMessageBuffers(threadId+numThreads);
				recvResp[threadId]->receiveAllIncomingMessages();
			}

			LOG_DEBUG(3, "Finishing communication sync for kmer lookups: " << batchReadIdx);
			LOG_DEBUG(4, "readIndexBuffer: " << readIndexBuffer[threadId].size() << "/" << readIndexBuffer[threadId][readIndexBuffer[threadId].size()-1]);
			LOG_DEBUG(4, "batchBuffer: " << batchBuffer[threadId].size());
			LOG_DEBUG(4, "Waiting for Request buffers to finalize");

			sendReq[threadId]->finalize(threadId);
			recvReq[threadId]->finalize();

			LOG_DEBUG(4, "Waiting for Response buffers to finalize");

			sendResp[threadId]->finalize(threadId+numThreads);
			recvResp[threadId]->finalize();

			LOG_DEBUG(4, "Delivery Request sent: " << sendReq[threadId]->getNumDeliveries() << " received: " << recvReq[threadId]->getNumDeliveries() << " "
					 <<  "Response sent: " << sendResp[threadId]->getNumDeliveries() << " received: " << recvResp[threadId]->getNumDeliveries());
			LOG_DEBUG(4, "Messages request/response: " << sendReq[threadId]->getNumMessages() << "/" << recvResp[threadId]->getNumMessages() << " "
					 << recvReq[threadId]->getNumMessages() << "/" << sendResp[threadId]->getNumMessages());
			assert( sendReq[threadId]->getNumMessages() == recvResp[threadId]->getNumMessages() );
			assert( recvReq[threadId]->getNumMessages() == sendResp[threadId]->getNumMessages() );

			LOG_DEBUG(3, "Starting trim for kmer lookups: " << batchReadIdx);
			for(ReadSetSizeType i = 0; i < readIndexBuffer[threadId].size() ; i++ ) {
				ReadSetSizeType &readIdx = readIndexBuffer[threadId][i];

				KmerValueVectorIterator buffBegin = (batchBuffer[threadId].begin() + readOffsetBuffer[threadId][i] );
				KmerValueVectorIterator buffEnd = ( ((i+1) < readOffsetBuffer[threadId].size()) ? (batchBuffer[threadId].begin() + readOffsetBuffer[threadId][i+1]) : batchBuffer[threadId].end() );
				ReadTrimType &trim = this->_trims[readIdx];

				if (useKmers) {
					this->trimReadByMinimumKmerScore(minimumKmerScore, trim, buffBegin, buffEnd);
				}

				this->setTrimHeaders(trim, useKmers);
			}

			LOG_DEBUG(2, "Finished assigning trim values: " << batchReadIdx);
			batchReadIdx += batchSize;

			// local & world threads are okay to start without sync
		}

		LOG_DEBUG(2, "scoreAndTrimReads(): barrier.  Finished trimming, waiting for remote processes");
		_world.barrier();

		std::stringstream ss;
		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();
			// initialize message buffers

			LOG_DEBUG(2, "Releasing Request/Response message buffers");
			#pragma omp critical
			{
				ss		<< "recvResp["<<threadId<<"] received " << recvResp[threadId]->getNumDeliveries() << "/" << recvResp[threadId]->getNumMessages()
						<< "sendResp["<<threadId<<"] sent     " << sendResp[threadId]->getNumDeliveries() << "/" << sendResp[threadId]->getNumMessages()
						<< "recvReq["<<threadId<<"] received " << recvReq[threadId]->getNumDeliveries() << "/" << recvReq[threadId]->getNumMessages()
						<< "sendReq["<<threadId<<"] sent     " << sendReq[threadId]->getNumDeliveries() << "/" << sendReq[threadId]->getNumMessages() << std::endl;
			}
			delete recvResp[threadId];
			delete sendResp[threadId];

			delete recvReq[threadId];
			delete sendReq[threadId];
		}
		std::string s = ss.str();
		LOG_DEBUG(1, s);
		LOG_DEBUG(2, "scoreAndTrimReads(): barrier.  Finished scoreAndTrimReadsMPI");
		_world.barrier();
	}

	// TODO
	// rescoreByBestCoveringSubset*

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
