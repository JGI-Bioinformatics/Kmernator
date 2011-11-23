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
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <vector>

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

// collective
void reduceOMPThreads(mpi::communicator &world) {
	Options::getOptions().validateOMPThreads();
#ifdef _USE_OPENMP
	int numThreads = omp_get_max_threads();
	numThreads = all_reduce(world, numThreads, mpi::minimum<int>());
	omp_set_num_threads(numThreads);
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "set OpenMP threads to " << numThreads);
#endif
}

// collective
void validateMPIWorld(mpi::communicator &world) {
	int provided;
	MPI_Query_thread(&provided);
#ifdef _USE_OPENMP
	if (provided != MPI_THREAD_FUNNELED && omp_get_max_threads() > 1) {
		if (world.rank() == 0)
			LOG_WARN(1, "Your version of MPI does not support MPI_THREAD_FUNNELED (" << provided << "), reducing OpenMP threads to 1")
		omp_set_num_threads(1);
	}
#endif
	reduceOMPThreads(world);
}

template< typename OptionsTempl >
mpi::communicator initializeWorldAndOptions(int argc, char *argv[]) {
	int threadProvided;
	int threadRequest = omp_get_max_threads() == 1 ? MPI_THREAD_SINGLE : MPI_THREAD_FUNNELED;
	MPI_Init_thread(&argc, &argv, threadRequest, &threadProvided);
	mpi::environment env(argc, argv);
	mpi::communicator world;
	MPI_Comm_set_errhandler( world, MPI::ERRORS_THROW_EXCEPTIONS );

	try {
		Logger::setWorld(&world);

		if (!OptionsTempl::parseOpts(argc, argv))
			LOG_THROW("Please fix the command line arguments." << std::endl << OptionsTempl::getOptionsErrorMsg());

		if (GeneralOptions::getOptions().getGatheredLogs())
			Logger::setWorld(&world, Options::getOptions().getDebug() >= 2);

		validateMPIWorld(world);

	} catch (...) {
		MPI_Finalize();
		exit(1);
	}
	world.barrier();
	return world;
}

// collective
std::string getRankSubdir(mpi::communicator &world, std::string prefix) {
	int subRank = world.rank() / 256;
	std::stringstream ss;
	ss << prefix << "/rank-subdirs-" << world.size() << "-0x" << std::hex << subRank;
	std::string subDir = ss.str();
	if (world.rank() == subRank * 256) {
		LOG_VERBOSE_OPTIONAL(1, true, "Making rank-subdirs: " << subDir);
	}
	mkdir(subDir.c_str(), 0777); // all ranks must mkdir if writing to local disks
	std::string subSubDir = subDir + "/" + boost::lexical_cast<std::string>(world.rank()) + "of" + boost::lexical_cast<std::string>(world.size());
	LOG_DEBUG(2, "getRankSubdir(" << prefix << "): " << subSubDir);
	world.barrier();
	mkdir(subSubDir.c_str(), 0777);
	return subSubDir;
}

// collective
void niceBarrier(mpi::communicator &world, int waitMs = 1) {
	mpi::communicator tmpWorld(world, mpi::comm_duplicate);
	int rank = tmpWorld.rank();
	int size = tmpWorld.size();
	char buf, buf2;
	buf  = '\0';
	buf2 = '\0';
	mpi::request rreq, rreq2;
	mpi::request sreq, sreq2;

	LOG_DEBUG_OPTIONAL(2, true, "Entering niceBarrier");
	if (rank == 0)
		sreq = tmpWorld.isend((rank+1) % size, 0, buf);

	rreq = tmpWorld.irecv((rank+size-1) % size, 0, buf2);
	while (! rreq.test() ) {
		boost::this_thread::sleep( boost::posix_time::milliseconds(waitMs) );
	}

	if (rank != 0)
		sreq = tmpWorld.isend((rank+1) % size, 0, buf);

	while (! sreq.test() ) {
		boost::this_thread::sleep( boost::posix_time::milliseconds(waitMs) );
	}

	if (rank == 0)
		sreq2 = tmpWorld.isend((rank+1) % size, 1, buf);

	rreq2 = tmpWorld.irecv((rank+size-1) % size, 1, buf2);
	while (! rreq2.test() ) {
		boost::this_thread::sleep( boost::posix_time::milliseconds(waitMs) );
	}

	if (rank != 0)
		sreq2 = tmpWorld.isend((rank+1) % size, 1, buf);

	while (! sreq2.test() ) {
		boost::this_thread::sleep( boost::posix_time::milliseconds(waitMs) );
	}

	LOG_DEBUG_OPTIONAL(2, true, "Exiting niceBarrier");

}

// collective
void setGlobalReadSetOffsets(mpi::communicator &world, ReadSet &store) {
	// share ReadSet sizes for globally unique readIdx calculations
	ReadSet::ReadSetSizeType readSetSize = store.getSize();
	ReadSet::ReadSetSizeType readSetSizesInput[world.size()], readSetSizesOutput[world.size()];
	for(int i = 0; i < world.size(); i++)
		readSetSizesInput[i] = i == world.rank() ? readSetSize : 0;
	LOG_DEBUG(2, "setGlobalReadSetOffset: all_reduce:" << readSetSize);

	mpi::all_reduce(world, (ReadSet::ReadSetSizeType*) readSetSizesInput, world.size(),  (ReadSet::ReadSetSizeType*) readSetSizesOutput, mpi::maximum<ReadSet::ReadSetSizeType>());
	ReadSet::ReadIdxVector readSizes(readSetSizesOutput, readSetSizesOutput + world.size());
	store.setGlobalOffsets(world.rank(), readSizes);
	LOG_DEBUG(2, "globalOffset: " << store.getGlobalOffset(world.rank()) << " of " << store.getGlobalSize());
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

public:
	DistributedKmerSpectrum(mpi::communicator &_world, unsigned long buckets = 0, bool separateSingletons = true)
	: KS(buckets * _world.size(), separateSingletons), world(_world) {
	}
	~DistributedKmerSpectrum() {
	}
	DistributedKmerSpectrum &operator=(const DistributedKmerSpectrum &other) {
		*((KS*) this) = other;
		world = other.world;
		return *this;
	}
	mpi::communicator &getWorld() {
		return world;
	}
	template<typename D>
	MmapFile writeKmerMap(D &kmerMap, std::string filepath) {
		MmapFile mmap;
		NumberType alignment = mmap.alignment();

		NumberType totalMmapSize = sizeof(NumberType) * 2; // numBuckets + bucketMask
		NumberType numBuckets = kmerMap.getNumBuckets();
		totalMmapSize += numBuckets * sizeof(NumberType); // offsets for each bucket
		totalMmapSize += alignment - (totalMmapSize % alignment); // align

		LOG_DEBUG_OPTIONAL(1, true, "numBuckets: " << numBuckets);
		// Get the size of each bucket
		NumberType *mySizeCounts = new NumberType[ numBuckets ];
		for(IndexType i = 0 ; i < numBuckets; i++) {
			mySizeCounts[i] = kmerMap.getBucketByIdx(i).size();
			// DMP buckets should only be populated at rank blocks
			assert( mySizeCounts[i] == 0ul || (long) (world.size() * i / numBuckets) == world.rank());
		}
		// Globally share bucket sizes
		NumberType *ourSizeCounts = new NumberType[ numBuckets ];
		all_reduce(world, mySizeCounts, numBuckets, ourSizeCounts, std::plus<NumberType>());

		// rename variable for clarity
		NumberType *offsetArray = mySizeCounts;

		// store the arrays in mmap blocked and aligned by mpi rank
		NumberType myOffset = 0, mySize = 0;
		int lastRank = world.size();
		for(IndexType i = 0; i < numBuckets; i++) {
			int rank = world.size() * i / numBuckets;
			if (rank == world.rank() && lastRank != rank)
				myOffset = totalMmapSize;

			offsetArray[i] = totalMmapSize;
			totalMmapSize += D::BucketType::sizeToStore( ourSizeCounts[i] );

			if (world.size() * (i+1) / numBuckets != (IndexType) rank) {
				totalMmapSize += alignment - (totalMmapSize % alignment); // align for each rank
				if (rank == world.rank())
					mySize = totalMmapSize - myOffset;
			}
			LOG_DEBUG(1, "myOffset " << myOffset << " mySize " << mySize << " totalMmapSize " << totalMmapSize);
			lastRank = rank;
		}

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
			int rank = (world.size() * i / numBuckets);
			if (rank == world.rank()) {
				assert(myOffset <= offsetArray[i]);
				assert(myOffset + mySize >= offsetArray[i] +  kmerMap.getBucketByIdx(i).size());
				kmerMap.getBucketByIdx(i).store(mmap.data() + offsetArray[i]);
			} else {
				assert( kmerMap.getBucketByIdx(i).size() == 0);
			}
		}
		// share with other processes
		msync(mmap.data() + myOffset, mySize, MS_ASYNC);

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

	class StoreKmerMessageHeader {
	public:

		ReadSetSizeType readIdx;
		PositionType readPos;
		WeightType weight; // weight is negative if kmer is rc of observed direction
		ExtensionMessagePacket extensionMsgPacket;

		// Kmer is next bytes, dynamically determined by KmerSizer::getTwoBitLength()
		// kmer is least complement

		// THIS IS DANGEROUS unless allocated an extra Kmer!
		Kmer *getKmer() {
			return (Kmer*) (((char*)this)+sizeof(*this));
		}

		Extension getLeft() const {
			return extensionMsgPacket.getLeft();
		}
		Extension getRight() const {
			return extensionMsgPacket.getRight();
		}
		void set(ReadSetSizeType _readIdx, PositionType _readPos, const WeightedExtensionMessagePacket &wemsgPkt, const Kmer &_kmer) {
			readIdx = _readIdx;
			readPos = _readPos;
			weight = wemsgPkt.getWeight();
			extensionMsgPacket.setExtensions(wemsgPkt.getLeft(), wemsgPkt.getRight());
			*(getKmer()) = _kmer;
		}
	};

	class  StoreKmerMessageHeaderProcessor {
	public:
		DistributedKmerSpectrum &_spectrum;
		std::vector<DataPointers> _pointers;
		StoreKmerMessageHeaderProcessor(DistributedKmerSpectrum &spectrum, bool isSolid = false) : _spectrum(spectrum), _isSolid(isSolid) {
			assert(!omp_in_parallel());
			_pointers.resize(omp_get_max_threads(), DataPointers(spectrum));
		}
		inline DataPointers &getDataPointer() {
			return _pointers[omp_get_thread_num()];
		}
		inline DistributedKmerSpectrum &getSpectrum() {
			return _spectrum;
		}
		inline bool isSolid() const {
			return _isSolid;
		}
		bool _isSolid;
		int process(StoreKmerMessageHeader *msg, MessagePackage &msgPkg) {
			assert(msgPkg.tag == omp_get_thread_num());
			LOG_DEBUG(5, "StoreKmerMessage: " << msg->readIdx << " " << msg->readPos << " " << msg->weight << " " << msg->getKmer()->toFasta());
			getSpectrum().append(getDataPointer(), *msg->getKmer(), msg->weight, msg->readIdx, msg->readPos, isSolid(), msg->getLeft(), msg->getRight());
			return 0;
		}
	};

	typedef MPIAllToAllMessageBuffer< StoreKmerMessageHeader, StoreKmerMessageHeaderProcessor > StoreKmerMessageBuffer;

	void _buildKmerSpectrumMPI(const ReadSet &store, bool isSolid) {
		int numThreads = omp_get_max_threads();
		int rank = world.rank();
		int worldSize = world.size();

		int messageSize = sizeof(StoreKmerMessageHeader) + KmerSizer::getByteSize();

		long readSetSize = store.getSize();

		LOG_VERBOSE(2, "starting _buildSpectrumMPI");

		StoreKmerMessageBuffer *msgBuffers;

		LOG_DEBUG(2, "building spectrum using " << numThreads << " threads (" << omp_get_max_threads() << ")");

		ReadSetSizeType globalReadSetOffset = store.getGlobalOffset(world.rank());
		assert( world.rank() == 0 ? (globalReadSetOffset == 0) : (globalReadSetOffset > 0) );
		msgBuffers = new StoreKmerMessageBuffer(world, messageSize, StoreKmerMessageHeaderProcessor(*this,isSolid));

		std::stringstream ss;
		#pragma omp parallel num_threads(numThreads)
		{
			int threadId = omp_get_thread_num();

			#pragma omp master
			{
				LOG_DEBUG(2, "message buffers ready");
				world.barrier();
			}
			#pragma omp barrier

			// allow the master thread to only handle communications
			int loopThreadId = threadId, loopNumThreads = numThreads;
			if (numThreads > 1) {
				loopThreadId--; loopNumThreads--;
			}
			if (loopThreadId >= 0) {
			  for(long readIdx = loopThreadId ; readIdx < readSetSize; readIdx+=loopNumThreads)
			  {

				const Read &read = store.getRead( readIdx );

				if (read.isDiscarded())
					continue;

				DataPointers pointers(*this);
				KmerWeightedExtensions kmers = KmerReadUtils::buildWeightedKmers(read, true, true);
				ReadSetSizeType globalReadIdx = readIdx + globalReadSetOffset;
				LOG_DEBUG(3, "_buildKmerSpectrumMPI(): Read " << readIdx << " (" << globalReadIdx << ") " << kmers.size() );

				for (PositionType readPos = 0 ; readPos < kmers.size(); readPos++) {
					int rankDest, threadDest;
					const WeightedExtensionMessagePacket &v = kmers.valueAt(readPos);
					WeightType weight = v.getWeight();
					if ( TrackingData::isDiscard( (weight<0.0) ? 0.0-weight : weight ) )  {
						LOG_DEBUG(4, "discarded kmer " << readIdx << "@" << readPos << " " << weight << " " << kmers[readPos].toFasta());
					} else {

						this->getThreadIds(kmers[readPos], threadDest, numThreads, rankDest, worldSize, true);

						if (rankDest == rank && threadDest == threadId) {
							this->append(pointers, kmers[readPos], v.getWeight(), globalReadIdx, readPos, isSolid, v.getLeft(), v.getRight());
						} else {
							msgBuffers->bufferMessage(rankDest, threadDest)->set(globalReadIdx, readPos, v, kmers[readPos]);
						}
					}
				}

				if (loopThreadId == 0 && readIdx % 1000000 == 0) {
					if (world.rank() == 0) {
						LOG_VERBOSE_OPTIONAL(1, true, "distributed processing " << (readIdx * world.size()) << " reads");
					} else {
						LOG_DEBUG(2, "local processed: " << readIdx << " reads");
					}
				}
			  }
			}

			LOG_DEBUG(2, "finished generating kmers from reads");

			msgBuffers->finalize();

		} // omp parallel

		delete msgBuffers;

		LOG_DEBUG(3, "_buildKmerSpectrumMPI() final barrier");
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

	void buildKmerSpectrum(const ReadSet &store ) {
		return this->buildKmerSpectrum(store, false);
	}

	void buildKmerSpectrum(const ReadSet &store, bool isSolid) {

		_buildKmerSpectrumMPI(store, isSolid);

		std::string hist = getHistogram(isSolid);
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Raw Histogram\n" << hist);
		LOG_DEBUG(2, "buildKmerSpectrum() barrier");
		world.barrier();

		// purge low counts
		if (KmerOptions::getOptions().getMinDepth() > 1) {
			LOG_VERBOSE(2, "Clearing memory from singletons: " << this->singleton.size() );
			LOG_DEBUG(2, MemoryUtils::getMemoryUsage());
			this->singleton.clear();
		}
		if (KmerOptions::getOptions().getMinDepth() > 2) {
			LOG_VERBOSE(2, "Purging low count kmers (< " << KmerOptions::getOptions().getMinDepth() << ")");
			LOG_DEBUG(2, MemoryUtils::getMemoryUsage());
			this->purgeMinDepth(KmerOptions::getOptions().getMinDepth());
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


	Kmernator::MmapFileVector writeKmerMaps(string mmapFilename = Options::getOptions().getOutputFile()) {
		// communicate sizes and allocate permanent file
		LOG_VERBOSE(2, "Merging partial spectrums" );
		LOG_DEBUG(3, MemoryUtils::getMemoryUsage() );
		Kmernator::MmapFileVector ourSpectrum;
		if (this->hasSolids){
			ourSpectrum.push_back(this->writeKmerMap(this->solid, mmapFilename + "-solid"));
		}
		ourSpectrum.push_back(this->writeKmerMap(this->weak, mmapFilename));

		if (KmerOptions::getOptions().getMinDepth() <= 1 && this->hasSingletons) {
			ourSpectrum.push_back(this->writeKmerMap(this->singleton, mmapFilename + "-singleton"));
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

	};
	class PurgeVariantKmerMessageHeaderProcessor;

	class PurgeVariantKmerMessageHeaderProcessor {
	public:
		DistributedKmerSpectrum &_spectrum;
		std::vector< DataPointers > _pointers;
		double _variantSigmas, _minDepth;
		PurgeVariantKmerMessageHeaderProcessor(DistributedKmerSpectrum &spectrum,  double variantSigmas, double minDepth)
		: _spectrum(spectrum), _variantSigmas(variantSigmas), _minDepth(minDepth) {
			for(int i = 0; i < omp_get_max_threads(); i++)
				_pointers.push_back( DataPointers(_spectrum) );
		}
		inline DataPointers &getDataPointer() {
			return _pointers[omp_get_thread_num()];
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
		int process(PurgeVariantKmerMessageHeader *msg, MessagePackage &msgPkg) {
			LOG_DEBUG(5, "PurgeVariantKmerMessage: " << msg->threshold << " " << msg->getKmer()->toFasta());
			DistributedKmerSpectrum &spectrum = getSpectrum();
			DataPointers &pointers = getDataPointer();
			double dummy;
			bool purged = spectrum._setPurgeVariant(pointers, *msg->getKmer(), msg->threshold, dummy);

			if (purged)
				spectrum.variantWasPurged();
			return 0;
		}
	};

	typedef MPIAllToAllMessageBuffer< PurgeVariantKmerMessageHeader, PurgeVariantKmerMessageHeaderProcessor > PurgeVariantKmerMessageBuffer;

	void variantWasPurged(long count = 1) {
		#pragma omp atomic
		_purgedVariants += count;
	}
private:
	PurgeVariantKmerMessageBuffer *msgPurgeVariant;

	long _purgedVariants;
	double variantSigmas;
	double minDepth;
	int _variantNumThreads;

	void _preVariants(double variantSigmas, double minDepth) {
		_purgedVariants = 0;
		long messageSize = sizeof(PurgeVariantKmerMessageHeader) + KmerSizer::getByteSize();

		msgPurgeVariant = new PurgeVariantKmerMessageBuffer(world, messageSize, PurgeVariantKmerMessageHeaderProcessor(*this, variantSigmas, minDepth));

		LOG_DEBUG(2, "_preVariants(): barrier");
		world.barrier();
	}

	long _postVariants(long purgedKmers) {
		_purgedVariants += this->KS::_postVariants(purgedKmers);
		delete msgPurgeVariant;

		return _purgedVariants;
	}
	void _variantThreadSync(long processed, long remaining, double maxDepth) {
		// call parent
		this->KS::_variantThreadSync(processed, remaining, maxDepth);

		msgPurgeVariant->finalize();
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
		int worldSize = world.size();
		if (editDistance == 0)
			return 0;

		WeakBucketType::permuteBases(kmer, variants, editDistance, true);

		for(SequenceLengthType i = 0 ; i < variants.size(); i++) {
			int rankDest, threadDest;
			Kmer &varKmer = variants[i];
			this->solid.getThreadIds(varKmer, threadDest, 1, rankDest, worldSize);

			if (rankDest == rank) {
				double dummy;
				if (this->_setPurgeVariant(pointers, varKmer, threshold, dummy))
					variantWasPurged();
			} else {
				msgPurgeVariant->bufferMessage(rankDest, threadDest)->set(threshold, varKmer);
			}
		}
		msgPurgeVariant->sendReceive();
		return 0;
	}


}; // DistributedKmerSpectrum


/*
 * DistributedOfstreamMap
 */
class DistributedOfstreamMap : public OfstreamMap
{

private:
	mpi::communicator _world;
	std::string _tempPrefix;
	std::string _realOutputPrefix;

	static string getTempPath(std::string tempPath) {
		return tempPath + UniqueName::generateUniqueName("/.tmp-output");
	}
protected:
	virtual void close() {
		LOG_VERBOSE_OPTIONAL(2, _world.rank() == 0, "Concatenating all MPI rank files");
		KeySet keys = getGlobalKeySet();
		if (isBuildInMemory())
			writeGlobalFiles(keys);
		OfstreamMap::close();
		if (!isBuildInMemory())
			concatenateMPI(keys);
	}

public:
	DistributedOfstreamMap(mpi::communicator &world, std::string outputFilePathPrefix = Options::getOptions().getOutputFile(), std::string suffix = FormatOutput::getDefaultSuffix(), std::string tempPath = Options::getOptions().getTmpDir())
	 :  OfstreamMap(getTempPath(tempPath), suffix), _world(world), _tempPrefix(), _realOutputPrefix(outputFilePathPrefix) {
		_tempPrefix = OfstreamMap::getOutputPrefix();
		LOG_DEBUG(3, "DistributedOfstreamMap(world, " << outputFilePathPrefix << ", " << suffix << "," << tempPath <<")");
		setBuildInMemory(Options::getOptions().getBuildOutputInMemory());
	}

	~DistributedOfstreamMap() {
		LOG_DEBUG_OPTIONAL(2, _world.rank() == 0, "~DistributedOfstreamMap()");
		this->clear();
	}

	virtual std::string getRank() const {
		return std::string("--MPIRANK-") + boost::lexical_cast<std::string>(_world.rank());
	}
	virtual void clear() {
		LOG_DEBUG_OPTIONAL(2, true, "DistributedOfstreamMap::clear()");
		this->close();
		_clear();
	}
	virtual std::string getRealFilePath(std::string key) const {
		return _realOutputPrefix + key + getSuffix();
	}
	// gets global keys to rank0.  All other ranks may have partial set...
	KeySet getGlobalKeySet() {
		LOG_DEBUG_OPTIONAL(2, _world.rank()==0, "Calling DistributedOfstreamMap::getGlobalKeySet()");

		// Send all filenames (minus Rank) to master
		KeySet keys = getKeySet();

		if (_world.rank() != 0) {
			_world.send(0, 0, keys);
		} else {
			for(int i = 1 ; i < _world.size() ; i++) {
				KeySet newKeys;
				_world.recv(i, 0, newKeys);
				keys.insert(newKeys.begin(), newKeys.end());
			}
			LOG_DEBUG_OPTIONAL(2, true, "getGlobalKeySet(): Collectively writing " << keys.size() << " files");
			for(KeySet::iterator it = keys.begin(); it != keys.end(); it++)
				LOG_DEBUG(3, "File key: " << *it);
		}

		return keys;
	}
	void writeGlobalFiles(KeySet &keys) {
		LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "Calling DistributedOfstreamMap::writeGlobalFiles()");
		assert(isBuildInMemory());

		int size = _world.size();
		MPI_Comm world = _world;

		// synchronize all files
		int numFiles = keys.size();
		mpi::broadcast(_world, numFiles, 0);

		KeySet::iterator itF = keys.begin();
		for(int fileNum = 0; fileNum < numFiles; fileNum++) {
			std::string key;
			if (_world.rank() == 0) {
				assert(itF != keys.end());
				key = *(itF++);
			}
			mpi::broadcast(_world, key, 0);
			std::string fullPath = getRealFilePath(key);
			LOG_VERBOSE_OPTIONAL(1, _world.rank() == 0, "writeGlobalFiles(): Collectively writing: " << fullPath);

			std::string contents;

			Iterator it = this->_map->find(key);

			if (it == this->_map->end()) {
				LOG_WARN(1, "Could not find " << key << " in DistributedOfstreamMap");
			} else {
				assert(it->second.isStringStream());
				contents = it->second.getFinalString();
			}

			long long int mySize = contents.length();
			long long int sendPos[size], recvPos[size], totalSize = 0, myStart = 0;
			for(int i = 0; i < size; i++)
				sendPos[i] = _world.rank() == i ? mySize : 0;
			MPI_Allreduce(&sendPos, &recvPos, size, MPI_LONG_LONG_INT, MPI_SUM, _world);
			for(int i = 0; i < size; i++) {
				if (_world.rank() == i)
					myStart = totalSize;
				totalSize += recvPos[i];
			}

			LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "Opening " << fullPath);
			int err;
			MPI_File ourFile;
			err = MPI_File_open(world, const_cast<char*>(fullPath.c_str()), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &ourFile);
			if (err != MPI_SUCCESS) {
				LOG_ERROR(1, "Could not open " << fullPath << " collectively");
				throw;
			}
			err = MPI_File_set_size(ourFile, totalSize);
			if (err != MPI_SUCCESS) {
				LOG_ERROR(1, "Could not set the size for " << fullPath << " to " << totalSize);
				throw;
			}
			LOG_DEBUG(2, "Writing " << mySize << " at " << myStart << " to " << fullPath);
			char *data = const_cast<char*>(contents.data());
			long long int  offset = 0;
			long long int maxwrite = 0xf000000; // keep writes to less than max int size at a time to avoid MPI overflows
			while (offset < mySize) {
				long long int thisWriteSize = std::min(maxwrite, mySize - offset);
				MPI_File_write_at(ourFile, myStart+offset, data+offset, thisWriteSize, MPI_BYTE, MPI_STATUS_IGNORE);
				offset += thisWriteSize;
			}
			LOG_DEBUG_OPTIONAL(1, _world.rank()==0, "Closing " << fullPath);
			MPI_File_close(&ourFile);
		}

	}
	void concatenateMPI(KeySet &keys) {
		LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "Calling DistributedOfstreamMap::concatenateMPI()");

		// synchronize all files
		int numFiles = keys.size();
		mpi::broadcast(_world, numFiles, 0);

		KeySet::iterator itF = keys.begin();
		for(int fileNum = 0; fileNum < numFiles; fileNum++) {
			std::string key;
			if (_world.rank() == 0) {
				assert(itF != keys.end());
				key = *(itF++);
			}
			mpi::broadcast(_world, key, 0);
			std::string fullPath = getRealFilePath(key);
			LOG_VERBOSE_OPTIONAL(1, _world.rank() == 0, "concatenateMPI(): Collectively writing: " << fullPath);

			Iterator it = this->_map->find(key);
			std::string myFilePath;
			if (it == this->_map->end()) {
				LOG_WARN(1, "Could not find " << myFilePath << " in DistributedOfstreamMap");
			} else {
				myFilePath = getFilePath(key);
			}
			mergeFiles(_world, myFilePath, fullPath, true);
		}
	}

	static void mergeFiles(mpi::communicator &world, std::string rankFile, std::string globalFile, bool unlinkAfter = false) {
		long long int mySize = 0;
		char *buf[2];
		int bufSize = 1024*1024*32;
		buf[0] = new char[bufSize];
		buf[1] = new char[bufSize];
		int bufId = 0;
		MPI_Info info(MPI_INFO_NULL);

		int rank = world.rank();
		int size = world.size();
		int err;
		MPI_File myFile;
		err = MPI_File_open(MPI_COMM_SELF, const_cast<char*>(rankFile.c_str()), MPI_MODE_RDONLY, info, &myFile);
		if (err != MPI_SUCCESS) {
			LOG_WARN(1, "Could not open " << rankFile << " myself.  Merging 0 bytes");
		} else {
			err = MPI_File_get_size(myFile, &mySize);
		}

		long long int sendPos[size], recvPos[size], totalSize = 0, myStart = 0, myPos = 0;
		for(int i = 0; i < size; i++)
			sendPos[i] = rank == i ? mySize : 0;
		MPI_Allreduce(&sendPos, &recvPos, size, MPI_LONG_LONG_INT, MPI_SUM, world);
		for(int i = 0; i < size; i++) {
			if (rank == i)
				myStart = totalSize;
			totalSize += recvPos[i];
		}

		LOG_DEBUG_OPTIONAL(2, rank==0, "Writing to '" << globalFile << "' at " << myStart << " for " << totalSize << " bytes");

		MPI_File ourFile;
		err = MPI_File_open(world, const_cast<char*>(globalFile.c_str()), MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &ourFile);
		if (err != MPI_SUCCESS) {
			LOG_THROW("Could not open " << globalFile << " collectively");
		}
		err = MPI_File_set_size(ourFile, totalSize);
		if (err != MPI_SUCCESS) {
			LOG_THROW("Could not set the size for " << globalFile << " to " << totalSize);
		}

		MPI_Status status;
		MPI_Request writeRequest = MPI_REQUEST_NULL;
		myPos = myStart;
		while (myPos < myStart + mySize) {
			err = MPI_File_read(myFile, buf[bufId % 2], bufSize, MPI_BYTE, &status);
			if (err != MPI_SUCCESS) {
				LOG_THROW("Could not read from " << rankFile);
			}
			int sendBytes;
			err = MPI_Get_count(&status, MPI_BYTE, &sendBytes);
			if (sendBytes == 0 || err != MPI_SUCCESS)
				break;

			err = MPI_Wait(&writeRequest, MPI_STATUS_IGNORE);
			if (err != MPI_SUCCESS) {
				LOG_THROW("Could not wait for write of " << globalFile);
			}
			err = MPI_File_iwrite_at(ourFile, myPos, buf[bufId++ % 2], sendBytes, MPI_BYTE, &writeRequest);
			if (err != MPI_SUCCESS) {
				LOG_THROW("Could not write to " << globalFile);
			}
			myPos += sendBytes;
		}
		err = MPI_Wait(&writeRequest, MPI_STATUS_IGNORE);
		if (err != MPI_SUCCESS) throw;

		err = MPI_File_close(&ourFile);
		if (err != MPI_SUCCESS) throw;
		err = MPI_File_close(&myFile);
		if (err != MPI_SUCCESS) throw;

		delete [] buf[0];
		delete [] buf[1];

		if (unlinkAfter)
			unlink(rankFile.c_str());
	}

};


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
	typedef DistributedOfstreamMap OFM;
	typedef Kmer::NumberType NumberType;

	typedef std::vector<ScoreType> KmerValueVector;
	typedef typename KmerValueVector::const_iterator KmerValueVectorIterator;
	typedef std::vector<KmerValueVector> KmerValueVectorVector;
	typedef std::vector<ReadSetSizeType> ReadIdxVector;
	typedef typename ReadIdxVector::const_iterator ReadIdxVectorIterator;

protected:
	mpi::communicator _world;

public:
	DistributedReadSelector(mpi::communicator &world, const ReadSet &reads, const KMType &map)
		: RS(reads, map), _world(world) {
		LOG_DEBUG(3, this->_map.toString());
	}
	OFM getOFM(std::string outputFile, std::string suffix = FormatOutput::getDefaultSuffix()) {
		LOG_DEBUG(2, "DistributedReadSelector::getOFM(" << outputFile << ", " << suffix << ")");
		return OFM(_world, outputFile, suffix);
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

	class ReqRespKmerMessageHeader {
	public:
		long requestId;
		// either ScoreType or Kmer is the next part of the message
		// Kmer is next bytes, dynamically determined by KmerSizer::getTwoBitLength()
		// kmer is least complement

		// THIS IS DANGEROUS unless allocated an extra Kmer or ScoreType!
		Kmer *getKmer() {
			return (Kmer*) (((char*)this)+sizeof(*this));
		}
		ScoreType &getScore() {
			return *((ScoreType*) (((char*)this)+sizeof(*this)));
		}

		void set(long _requestId, const Kmer &_kmer) {
			requestId = _requestId;
			*(getKmer()) = _kmer;
		}

		void set(long _requestId, ScoreType _score) {
			requestId = _requestId;
			getScore() = _score;
		}
	};


	class ReqRespKmerMessageHeaderProcessor;
    typedef MPIAllToAllMessageBuffer< ReqRespKmerMessageHeader, ReqRespKmerMessageHeaderProcessor > ReqRespKmerMessageBuffer;

    class ReqRespKmerMessageHeaderProcessor {
    public:
		KmerValueVectorVector &_kmerValues;
		RS *_readSelector;
		int _numThreads;

		ReqRespKmerMessageHeaderProcessor(KmerValueVectorVector &kmerValues, RS &readSelector, int numThreads): _kmerValues(kmerValues), _readSelector(&readSelector), _numThreads(numThreads)  {}

 		// store response in kmer value vector
		int processRespond(ReqRespKmerMessageHeader *msg, MessagePackage &msgPkg) {
			LOG_DEBUG(5, "RespondKmerMessage: " << msg->requestId << " " << msg->getScore() << " recv Source: " << msgPkg.source << " recvTag: " << msgPkg.tag);
			assert( msgPkg.tag == omp_get_thread_num() + _numThreads);
			_kmerValues[omp_get_thread_num()][msg->requestId] = msg->getScore();
			return sizeof(ScoreType);
		}

		// lookup kmer in map and build response message
		int processRequest(ReqRespKmerMessageHeader *msg, MessagePackage &msgPkg) {
			int destSource = msgPkg.source;
			int destTag = msgPkg.tag + _numThreads;
			ScoreType score = _readSelector->getValue( *msg->getKmer() );
			LOG_DEBUG(5, "RequestKmerMessage: " << msg->getKmer()->toFasta() << " " << msg->requestId << " " << score << " destSource: " << destSource << " sdestTag: " << destTag);
			assert(msgPkg.tag == omp_get_thread_num());

			((ReqRespKmerMessageBuffer*)msgPkg.bufferCallback)->bufferMessage(destSource, destTag, sizeof(ScoreType))->set(msg->requestId, score);
			return KmerSizer::getByteSize();
		}

		int process(ReqRespKmerMessageHeader *msg, MessagePackage &msgPkg) {
			int recvTag = msgPkg.tag;
			if (recvTag >= _numThreads)
				return processRespond(msg, msgPkg);
			else
				return processRequest(msg, msgPkg);
		}
	};


	int _batchKmerLookup(const Read &read, SequenceLengthType markupLength, ReadSetSizeType offset, KmerValueVector &batchBuffer,
			ReqRespKmerMessageBuffer &sendReq, int &thisThreadId, int &numThreads, int &rank, int &worldSize) {
		KA kmers = this->getKmersForRead(read);
		SequenceLengthType numKmers = kmers.size();

		this->_setNumKmers(markupLength, numKmers);
		assert(numKmers <= kmers.size());
		assert(offset == batchBuffer.size());

		if (numKmers > 0)
			batchBuffer.insert(batchBuffer.end(), numKmers, ScoreType(-1));

		for(SequenceLengthType kmerIdx = 0; kmerIdx < numKmers; kmerIdx++) {
			int localThreadId, distributedThreadId;
			this->_map.getThreadIds(kmers[kmerIdx], localThreadId, numThreads, distributedThreadId, worldSize);
			ReadSetSizeType requestId = offset+kmerIdx;
			if (distributedThreadId == rank) {
				// handle this directly
				batchBuffer[requestId] = this->getValue(kmers[kmerIdx]);
			} else {
				LOG_DEBUG(5, "_batchKmerLookup() sendReq: destRank: " << distributedThreadId << " destTag: " << thisThreadId << " requestId: " << (requestId) << " kmer: " << kmers[kmerIdx].toFasta());
				sendReq.bufferMessage( distributedThreadId, thisThreadId, KmerSizer::getByteSize() )->set( requestId, kmers[kmerIdx]) ;
			}
		}
		return numKmers;
	}
	void scoreAndTrimReads(ScoreType minimumKmerScore, int correctionAttempts = 0) {
		this->_trims.resize(this->_reads.getSize());
		bool useKmers = KmerOptions::getOptions().getKmerSize() != 0;

		int numThreads = omp_get_max_threads();
		int rank = _world.rank();
		int worldSize = _world.size();

		ReadSetSizeType readsSize = this->_reads.getSize();
		ReadSetSizeType batchSize = Options::getOptions().getBatchSize();
		ReadSetSizeType batchReadIdx = 0;
		int maxKmers = std::max(1, (int) this->_reads.getMaxSequenceLength() - (int) KmerSizer::getSequenceLength());
		int reserveBB = ((batchSize * maxKmers) / numThreads) + 1;
		int reserveOffsets = (batchSize / numThreads) + 1;

		LOG_DEBUG(1, "Starting scoreAndTrimReadsMPI - trimming: " << readsSize << " using " << numThreads << " threads");

		KmerValueVectorVector batchBuffer; batchBuffer.resize(numThreads);
		ReadIdxVector readIndexBuffer[ numThreads ];
		ReadIdxVector readOffsetBuffer[ numThreads ];

		ReqRespKmerMessageBuffer *reqRespBuffer;

		LOG_DEBUG(2, "scoreAndTrimReads(): all_reduce: " << readsSize);
		ReadSetSizeType mostReads = mpi::all_reduce(_world, readsSize, mpi::maximum<ReadSetSizeType>());
		LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "Largest number of reads to batch: " << mostReads);

		reqRespBuffer = new ReqRespKmerMessageBuffer(_world, sizeof(ReqRespKmerMessageHeader),
				ReqRespKmerMessageHeaderProcessor(batchBuffer, *this, numThreads), 2);

		LOG_DEBUG(2, "scoreAndTrimReads(): barrier. message buffers ready");
		_world.barrier();

		#pragma omp parallel num_threads(numThreads) firstprivate(batchReadIdx)
		{
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

			// allow master thread to only handle communications
			int loopThreadId = threadId, loopNumThreads = numThreads;
			if (loopNumThreads > 1) {
				loopThreadId--; loopNumThreads--;
			}
			if (loopThreadId >= 0) {
			  for(ReadSetSizeType i = loopThreadId ; i < batchSize ; i+=loopNumThreads) {

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
					_batchKmerLookup(read, markupLength, offset, batchBuffer[threadId], *reqRespBuffer, threadId, numThreads, rank, worldSize);
				} else {
					this->trimReadByMarkupLength(read, this->_trims[readIdx], markupLength);
				}
			  }
			}
			assert(readOffsetBuffer[threadId].size() == readIndexBuffer[threadId].size());

			reqRespBuffer->sendReceive(); // flush/send all pending requests for this thread's batch
			reqRespBuffer->sendReceive();
			reqRespBuffer->sendReceive(); // receive all pending responses for this threads's batch
			reqRespBuffer->sendReceive();

			LOG_DEBUG(3, "Starting trim for kmer lookups: " << batchReadIdx);
			for(ReadSetSizeType i = 0; i < readIndexBuffer[threadId].size() ; i++ ) {
				ReadSetSizeType &readIdx = readIndexBuffer[threadId][i];

				KmerValueVectorIterator buffBegin = (batchBuffer[threadId].begin() + readOffsetBuffer[threadId][i] );
				KmerValueVectorIterator buffEnd = ( ((i+1) < readOffsetBuffer[threadId].size()) ? (batchBuffer[threadId].begin() + readOffsetBuffer[threadId][i+1]) : batchBuffer[threadId].end() );
				ReadTrimType &trim = this->_trims[readIdx];
				for(KmerValueVectorIterator it = buffBegin; it != buffEnd; it++)
					if (*it == -1)
						LOG_WARN(1, "readIdx: " << readIdx << " pos: " << (it - buffBegin) << " did not get updated!");

				if (useKmers) {
					this->trimReadByMinimumKmerScore(minimumKmerScore, trim, buffBegin, buffEnd);
				}

				this->setTrimHeaders(trim, useKmers);
			}

			//LOG_DEBUG(2, "Finished assigning trim values: " << batchReadIdx);
			batchReadIdx += batchSize;

			// local & world threads are okay to start without sync
		  }
		  reqRespBuffer->finalize();
		}

		LOG_DEBUG(2, "scoreAndTrimReads(): barrier.  Finished trimming, waiting for remote processes");
		_world.barrier();

		delete reqRespBuffer;
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
