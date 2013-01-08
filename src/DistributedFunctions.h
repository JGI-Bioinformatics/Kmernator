//
// Kmernator/src/DistributedFunctions.h
//
// Author: Rob Egan
//
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/

#ifndef DISTRIBUTED_FUNCTIONS_H_
#define DISTRIBUTED_FUNCTIONS_H_

#include "mpi.h"

#include "config.h"
#include "Options.h"
#include "Log.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "MPIBuffer.h"
#include "ReadSet.h"
#include "ReadSelector.h"
#include "MPIUtils.h"
#include "DistributedOfstreamMap.h"

#include "boost/optional.hpp"
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <vector>


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
	typedef typename KS::SizeTracker SizeTracker;
	typedef typename SizeTracker::Elements SizeTrackerElements;

protected:
	mpi::communicator world;

public:
	DistributedKmerSpectrum(mpi::communicator &_world, unsigned long buckets = 0, bool separateSingletons = true)
	: KS(buckets, separateSingletons), world(_world, mpi::comm_duplicate) {
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
		assert(store.isGlobal());
		if (store.getGlobalSize() == 0)
			return;
		int numThreads = omp_get_max_threads();
		int rank = world.rank();
		int worldSize = world.size();

		int messageSize = sizeof(StoreKmerMessageHeader) + KmerSizer::getByteSize();

		long readSetSize = store.getSize();

		LOG_VERBOSE(2, "starting _buildSpectrumMPI with " << omp_get_max_threads() << " threads");

		StoreKmerMessageBuffer *msgBuffers;

		LOG_DEBUG(2, "building spectrum using " << numThreads << " threads (" << omp_get_max_threads() << ")");

		ReadSetSizeType globalReadSetOffset = store.getGlobalOffset(world.rank());
		assert( world.rank() == 0 ? (globalReadSetOffset == 0) : (store.getSize() == 0 || globalReadSetOffset > 0) );
		msgBuffers = new StoreKmerMessageBuffer(world, messageSize, StoreKmerMessageHeaderProcessor(*this,isSolid));

		long kmerSubsample = KS::getKmerSubsample();

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
			bool isRunningInLoop = true;
			if (numThreads > 1) {
				if (loopThreadId == 0)
					isRunningInLoop = false;
				loopThreadId--; loopNumThreads--;
			}
			if (isRunningInLoop) {
				KmerReadUtils kru;
				long progressCount = 0, progressMark = 1000000 / world.size();
				for(long readIdx = loopThreadId ; readIdx < readSetSize; readIdx+=loopNumThreads)
				{

					const Read &read = store.getRead( readIdx );

					if (read.isDiscarded())
						continue;

					DataPointers pointers(*this);
					KmerWeightedExtensions &kmers = kru.buildWeightedKmers(read, true, true);
					ReadSetSizeType globalReadIdx = readIdx + globalReadSetOffset;
					LOG_DEBUG(3, "_buildKmerSpectrumMPI(): Read " << readIdx << " (" << globalReadIdx << ") " << kmers.size() );

					for (PositionType readPos = 0 ; readPos < kmers.size(); readPos++) {
						int rankDest, threadDest;
						if (kmerSubsample > 1 && kmers[readPos].hash() % kmerSubsample != 0) {
							continue;
						}
						const WeightedExtensionMessagePacket &v = kmers.valueAt(readPos);
						WeightType weight = v.getWeight();
						if ( TrackingData::isDiscard( (weight<0.0) ? 0.0-weight : weight ) )  {
							LOG_DEBUG(4, "discarded kmer " << readIdx << "@" << readPos << " " << weight << " " << kmers[readPos].toFasta());
						} else {

							this->getThreadIds(kmers[readPos], threadDest, loopNumThreads, rankDest, worldSize, true);

							if (rankDest == rank && threadDest == loopThreadId) {
								this->append(pointers, kmers[readPos], v.getWeight(), globalReadIdx, readPos, isSolid, v.getLeft(), v.getRight());
							} else {
								msgBuffers->bufferMessage(rankDest, numThreads == loopNumThreads ? threadDest : threadDest+1)->set(globalReadIdx, readPos, v, kmers[readPos]);
							}
						}
					}

					if (loopThreadId == 0 &&  progressCount++ % progressMark == 0) {
						LOG_VERBOSE_OPTIONAL(1, (progressCount-1) % world.size() == world.rank(), "distributed processing " << (readIdx * world.size()) << " reads. " << this->solid.size()* world.size() << "/" << this->weak.size()* world.size() << "/" << this->singleton.size()* world.size() << " kmers");
					}
				}
			}

			LOG_DEBUG(2, "finished generating kmers from reads");

			msgBuffers->finalize();

		} // omp parallel
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "distributed processed " << store.getGlobalSize() << " reads");

		delete msgBuffers;

		LOG_DEBUG(3, "_buildKmerSpectrumMPI() final barrier");
		world.barrier();
		LOG_DEBUG(1, "finished _buildKmerSpectrumMPI");
	}

	SizeTracker reduceSizeTracker(mpi::communicator &world) {
		SizeTracker s = this->getSizeTracker();
		SizeTrackerElements &elements = s.elements;
		int inct = elements.size();
		int ct, tmpct;
		LOG_DEBUG_OPTIONAL(1, true, "my sizeTracker size: " << inct);
		MPI_Allreduce(&inct, &ct, 1, MPI_INT, MPI_MAX, world);
		LOG_DEBUG_OPTIONAL(1, true, "our sizeTracker size: " << ct);

		SizeTracker newTracker;
		newTracker.resize(ct);
		SizeTrackerElements &newElements = newTracker.elements;

		long in[ct*4], out[ct*4];
		long *tmp = &in[0];
		for(int i = 0; i < ct; i++) {
			// repeat the last record, if shorter than the global max size
			tmpct = i;
			if (tmpct >= inct) tmpct = inct - 1;
			*(tmp++) = elements[tmpct].rawKmers;
			*(tmp++) = elements[tmpct].rawGoodKmers;
			*(tmp++) = elements[tmpct].uniqueKmers;
			*(tmp++) = elements[tmpct].singletonKmers;
		}
		MPI_Allreduce(in, out, ct*4, MPI_LONG_LONG_INT, MPI_SUM, world);
		tmp = &out[0];
		for(int i = 0; i < ct ; i++) {
			newElements[i].rawKmers = *(tmp++);
			newElements[i].rawGoodKmers = *(tmp++);
			newElements[i].uniqueKmers = *(tmp++);
			newElements[i].singletonKmers = *(tmp++);
		}
		return newTracker;
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
		this->buildKmerSpectrum(store, false);
	}

	void buildKmerSpectrum(const ReadSet &store, bool isSolid) {

		_buildKmerSpectrumMPI(store, isSolid);

		if(Log::isVerbose(2)) {
			std::string hist = getHistogram(isSolid);
			LOG_VERBOSE_OPTIONAL(2, world.rank() == 0, "Collective Raw Histogram\n" << hist);
		}
		LOG_DEBUG(2, "buildKmerSpectrum() barrier");
		world.barrier();

		// purge low counts
		if (KmerSpectrumOptions::getOptions().getMinDepth() > 1) {
			LOG_VERBOSE(2, "Clearing memory from singletons: " << this->singleton.size() );
			LOG_DEBUG(2, MemoryUtils::getMemoryUsage());
			this->singleton.clear();
		}
		if (KmerSpectrumOptions::getOptions().getMinDepth() > 2) {
			LOG_VERBOSE(2, "Purging low count kmers (< " << KmerSpectrumOptions::getOptions().getMinDepth() << ")");
			LOG_DEBUG(2, MemoryUtils::getMemoryUsage());
			this->purgeMinDepth(KmerSpectrumOptions::getOptions().getMinDepth());
		}
		LOG_DEBUG(3, this->weak.toString());

	}

	MPIHistogram _getHistogram(bool solidOnly = false) {
		MPIHistogram histogram(255);
		histogram.set(*this, solidOnly);
		LOG_DEBUG(2, "Individual histogram\n" << histogram.toString());
		histogram.reduce(world);
		return histogram;
	}

	virtual std::string getHistogram(bool solidOnly = false) {
		return _getHistogram(solidOnly).toString();
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

		if (KmerSpectrumOptions::getOptions().getMinDepth() <= 1 && this->hasSingletons) {
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
	typedef typename KmerValueVector::iterator KmerValueVectorIterator;
	typedef std::vector<KmerValueVector> KmerValueVectorVector;
	typedef std::vector<ReadSetSizeType> ReadIdxVector;
	typedef typename ReadIdxVector::const_iterator ReadIdxVectorIterator;

protected:
	mpi::communicator _world;

public:
	DistributedReadSelector(mpi::communicator &world, const ReadSet &reads, const KMType &map)
	: RS(reads, map), _world(world, mpi::comm_duplicate) {
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
			ReqRespKmerMessageBuffer &sendReq, int &thisThreadId, int &numThreads, int &rank, int &worldSize, KA &kmers) {
		this->getKmersForRead(read, kmers);
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
	void scoreAndTrimReads(ScoreType minimumKmerScore, enum RS::KmerScoringType scoringType = RS::_KS_MAX_SCORING) {
		if (scoringType == RS::_KS_MAX_SCORING)
			scoringType = this->_defaultScoringType;

		this->_trims.resize(this->_reads.getSize());
		bool useKmers = KmerBaseOptions::getOptions().getKmerSize() != 0;

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
		std::vector< KA > _kmers(numThreads, KA());

#pragma omp parallel num_threads(numThreads) firstprivate(batchReadIdx)
		{
			if (numThreads != omp_get_num_threads() && omp_get_thread_num() == 0) {
				LOG_WARN(1, "Using less threads than expected: " << omp_get_num_threads() << " not " << numThreads);
				numThreads = omp_get_num_threads();
			}
#pragma omp barrier
			while (batchReadIdx < mostReads) {
				int threadId = omp_get_thread_num();
				KA &kmers = _kmers[omp_get_thread_num()];
				// initialize read/kmer buffers
				batchBuffer[threadId].resize(0);
				readIndexBuffer[threadId].resize(0);
				readOffsetBuffer[threadId].resize(0);
				batchBuffer[threadId].reserve(reserveBB);
				readIndexBuffer[threadId].reserve(reserveOffsets);
				readOffsetBuffer[threadId].reserve(reserveOffsets);

				LOG_DEBUG(3, "Starting batch for kmer lookups: " << batchReadIdx);

				// allow master thread to only handle communications
				int loopThreadId = threadId, loopNumThreads = numThreads;
				bool isRunningInLoop = true;
				if (numThreads > 1) {
					if (loopThreadId == 0)
						isRunningInLoop = false;
					loopThreadId--; loopNumThreads--;
				}
				if (isRunningInLoop) {
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
							_batchKmerLookup(read, markupLength, offset, batchBuffer[threadId], *reqRespBuffer, threadId, numThreads, rank, worldSize, kmers);
						} else {
							this->trimReadByMarkupLength(read, this->_trims[readIdx], markupLength);
						}
					}
					LOG_VERBOSE_OPTIONAL(1, _world.rank() == 0 && loopThreadId == 0, "trimming batch: " << batchReadIdx * _world.size());
				
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

						for(KmerValueVectorIterator it = buffBegin; it != buffEnd; it++) {
							if (*it == -1)
								LOG_WARN(1, "readIdx: " << readIdx << " pos: " << (it - buffBegin) << " did not get updated!");
							if (*it < minimumKmerScore)
								*it = 0.0;
						}

						if (useKmers) {
							this->trimReadByMinimumKmerScore(minimumKmerScore, trim, buffBegin, buffEnd);
							this->scoreReadByScoringType(buffBegin + trim.trimOffset, buffBegin + trim.trimOffset + trim.trimLength, trim, scoringType);
						}

						this->setTrimHeaders(trim, useKmers);
					}

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
