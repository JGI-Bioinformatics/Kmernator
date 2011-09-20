/*
 * Meraculous.h
 *
 *  Created on: Sep 16, 2011
 *      Author: regan
 */

#ifndef MERACULOUS_H_
#define MERACULOUS_H_

#include "KmerTrackingData.h"
#include "DistributedFunctions.h"

typedef DistributedKmerSpectrum<ExtensionTrackingData, ExtensionTrackingData, ExtensionTrackingDataSingleton> _MeraculousDistributedKmerSpectrum;
class MeraculousDistributedKmerSpectrum : public _MeraculousDistributedKmerSpectrum {
public:

	typedef _MeraculousDistributedKmerSpectrum::DataPointers DataPointers;
	typedef _MeraculousDistributedKmerSpectrum::WeakElementType WeakElementType;

	MeraculousDistributedKmerSpectrum(mpi::communicator &_world, unsigned long buckets = 0, bool separateSingletons = true)
	: _MeraculousDistributedKmerSpectrum(_world, buckets, separateSingletons) {
	}
	~MeraculousDistributedKmerSpectrum() {
	}
	MeraculousDistributedKmerSpectrum &operator=(const MeraculousDistributedKmerSpectrum &other) {
		*((_MeraculousDistributedKmerSpectrum*) this) = (_MeraculousDistributedKmerSpectrum) other;
		return *this;
	}

	class StoreKmerExtensionMessageHeader : public ExtensionMessagePacket {
	public:
		// first 4 bytes are the ExtensionMessagePacket
		// Kmer is next bytes, dynamically determined by KmerSizer::getTwoBitLength()
		// kmer is least complement

		// THIS IS DANGEROUS unless allocated an extra Kmer!
		Kmer *getKmer() {
			return (Kmer*) (((char*)this)+sizeof(*this));
		}
		void set(Extension left, Extension right, const Kmer &_kmer) {
			setExtensions(left, right);
			*(getKmer()) = _kmer;
		}
	};

	class  StoreKmerExtensionMessageHeaderProcessor {
	public:
		MeraculousDistributedKmerSpectrum &_spectrum;
		std::vector<DataPointers> _pointers;
		StoreKmerExtensionMessageHeaderProcessor(MeraculousDistributedKmerSpectrum &spectrum) : _spectrum(spectrum) {
			assert(!omp_in_parallel());
			_pointers.resize(omp_get_max_threads(), DataPointers(spectrum));
		}
		inline DataPointers &getDataPointer() {
			return _pointers[omp_get_thread_num()];
		}
		inline DistributedKmerSpectrum &getSpectrum() {
			return _spectrum;
		}
		int process(StoreKmerExtensionMessageHeader *msg, MessagePackage &msgPkg) {
			assert(msgPkg.tag == omp_get_thread_num());

			WeakElementType element = getSpectrum().getIfExistsWeak( *msg->getKmer() );
			if (element.isValid())
				element.value().trackExtensions(msg->getLeft(), msg->getRight());
			return 0;
		}
	};

	typedef MPIAllToAllMessageBuffer< StoreKmerExtensionMessageHeader, StoreKmerExtensionMessageHeaderProcessor > StoreKmerExtensionMessageBuffer;

	void trackExtensions(const ReadSet &store) {
		int numThreads = omp_get_max_threads();
		int rank = world.rank();
		int worldSize = world.size();

		int messageSize = sizeof(StoreKmerExtensionMessageHeader) + KmerSizer::getByteSize();

		long readSetSize = store.getSize();


		LOG_VERBOSE(2, "starting _buildSpectrumMPI");

		StoreKmerExtensionMessageBuffer *msgBuffers;

		LOG_DEBUG(2, "building kmer extensions using " << numThreads << " threads (" << omp_get_max_threads() << ")");

		ReadSetSizeType globalReadSetOffset = store.getGlobalOffset(world.rank());
		assert( world.rank() == 0 ? (globalReadSetOffset == 0) : (globalReadSetOffset > 0) );
		msgBuffers = new StoreKmerExtensionMessageBuffer(world, messageSize, StoreKmerExtensionMessageHeaderProcessor(*this));

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
				KmerWeights kmers = KmerReadUtils::buildWeightedKmers(read, true, true);
				ReadSetSizeType globalReadIdx = readIdx + globalReadSetOffset;
				LOG_DEBUG(3, "_buildKmerSpectrumMPI(): Read " << readIdx << " (" << globalReadIdx << ") " << kmers.size() );

				for (PositionType readPos = 0 ; readPos < kmers.size(); readPos++) {
					int rankDest, threadDest;
					WeightType weight = kmers.valueAt(readPos);
					if ( TrackingData::isDiscard( (weight<0.0) ? 0.0-weight : weight ) )  {
						LOG_DEBUG(4, "discarded kmer " << readIdx << "@" << readPos << " " << weight << " " << kmers[readPos].toFasta());
					} else {

						this->getThreadIds(kmers[readPos], threadDest, numThreads, rankDest, worldSize, true);

						if (rankDest == rank && threadDest == threadId) {
							this->append(pointers, kmers[readPos], weight, globalReadIdx, readPos, isSolid);
						} else {
							msgBuffers->bufferMessage(rankDest, threadDest)->set(globalReadIdx, readPos, weight, kmers[readPos]);
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

	//TODO
	static StoreKmerExtensionMessageHeader *buildKmerExtensions(const Read &read, bool leastComplement = false, bool leastComplementForNegativeWeight = false) {
		SequenceLengthType readLength = read.getLength();
		bool needMalloc = readLength > MAX_STACK_SIZE;
		bool _bools[ needMalloc ? 0 : readLength ];
		bool *bools = _bools;
		if (needMalloc) {
			bools = new bool[readLength];
		}
		KmerWeights kmers(read.getTwoBitSequence(), readLength, leastComplement, bools);
		std::string quals = read.getQuals();
		size_t markupIdx = 0;

		BaseLocationVectorType markups = read.getMarkups();
		double weight = 0.0;
		double change = 0.0;
		bool isRef = false;
		if (quals.length() > 0 && quals[0] == Read::REF_QUAL)
			isRef = true;

		SequenceLengthType size = (SequenceLengthType) kmers.size();
		for (SequenceLengthType i = 0; i < size; i++) {
			if (isRef) {
				weight = 1.0;
			} else if (i % 1024 == 0 || weight == 0.0) {
				weight = 1.0;
				for (SequenceLengthType j = 0; j
						< KmerSizer::getSequenceLength(); j++)
					weight
							*= Read::qualityToProbability[(unsigned char) quals[i
									+ j]];
			} else {
				change = Read::qualityToProbability[(unsigned char) quals[i
						+ KmerSizer::getSequenceLength() - 1]]
						/ Read::qualityToProbability[(unsigned char) quals[i
								- 1]];
				weight *= change;
			}
			while (markupIdx < markups.size() && markups[markupIdx].second < i)
				markupIdx++;
			if (markupIdx < markups.size() && markups[markupIdx].second < i
					+ KmerSizer::getSequenceLength()) {
				weight = 0.0;
			}
			kmers.valueAt(i) = leastComplementForNegativeWeight && bools[i] ? weight : (0.0-weight);

		}
		if (Log::isDebug(5)) {
			ostream &debug = Log::Debug() << "KmerWeights: idx valueAt toFasta" << std::endl;;
		    for(Kmer::IndexType i = 0 ; i < kmers.size(); i++) {
			  debug << i << " " << kmers.valueAt(i) << " " << kmers[i].toFasta() << std::endl;
		    }
		}
		if (needMalloc) {
			delete [] bools;
		}
		return kmers;
	}


};

#endif /* MERACULOUS_H_ */
