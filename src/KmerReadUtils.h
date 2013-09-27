//
// Kmernator/src/KmerReadUtils.h
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

#ifndef _KMER_READ_UTILS_H
#define _KMER_READ_UTILS_H

#include "config.h"
#include "Sequence.h"
#include "KmerTrackingData.h"
#include "Kmer.h"

#include <boost/shared_ptr.hpp>
#include <boost/unordered_set.hpp>

class KmerReadUtils {
public:
	typedef std::vector< int > PositionVector;
	typedef boost::shared_ptr< PositionVector > PositionVectorPtr;
	typedef KmerArrayPair< PositionVectorPtr > KmerReferenceMap;
	typedef KmerReferenceMap::IndexType IndexType;
	typedef KmerWeightedExtensions::CompareArrayIdx CompareArrayIdx;
private:
	KmerWeightedExtensions kmers;
	KmerReferenceMap kmerMap;
public:
	KmerReadUtils() {
		LOG_DEBUG_OPTIONAL(2, true, "KmerReadUtils()" << &kmers);
		if (! Read::isQualityToProbabilityInitialized() )
			Read::setMinQualityScore();
	}
	~KmerReadUtils() {}

	// reports low with indels & subs
	SequenceLengthType getReferenceMapKmerOverlap(const Read &tgt) {
		SequenceLengthType overlap = 0;
		kmers.build(tgt.getTwoBitSequence(), tgt.getLength(), true);
		for(int i = 0; i < (int) kmers.size(); i++) {
			KmerReferenceMap::Iterator it;
			it = kmerMap.find(kmers.get(i));
			if (it != kmerMap.end()) {
				overlap += it->value()->size();
			}
		}
		return overlap;
	}

	// accounts for every base that is covered by at least one kmer
	SequenceLengthType getReferenceMapOverlap(const Read &tgt) {
		boost::unordered_set<SequenceLengthType> coveredPositions;
		kmers.build(tgt.getTwoBitSequence(), tgt.getLength(), true);
		for(int i = 0; i < (int) kmers.size(); i++) {
			KmerReferenceMap::Iterator it;
			LOG_DEBUG(2, "getReferenceMapOverlap(): finding kmer: " << kmers.get(i).toFasta());
			it = kmerMap.find(kmers.get(i));
			if (it != kmerMap.end()) {
				LOG_DEBUG(2, "getReferenceMapOverlap(): matched kmer: " << it->key().toFasta() << " " << it->value()->front());
				for(PositionVector::iterator it2 = it->value()->begin(); it2 != it->value()->end(); it2++) {
					int startPos, endPos, signedPosition = *it2;
					if (signedPosition < 0) {
						startPos = 0 - signedPosition;
					} else {
						startPos = signedPosition;
					}
					endPos = startPos + KmerSizer::getSequenceLength();
					
					assert(startPos < endPos);
					assert(startPos >= 0);
					assert(endPos < (int) tgt.getLength());
					for(int j = startPos; j < endPos; j++)
						coveredPositions.insert(j);
				}
			}
		}
		LOG_DEBUG(2, "getReferenceMapOverlap(): " << coveredPositions.size());
		return coveredPositions.size();
	}
	SequenceLengthType buildReferenceMap(const Read &read) {
		kmerMap.clear(false);
		SequenceLengthType readLength = read.getLength();
		SequenceLengthType numKmers = readLength - KmerSizer::getSequenceLength() + 1;
		STACK_ALLOC(bool, bools, readLength);

		kmers.build(read.getTwoBitSequence(), readLength, true, bools);
		SequenceLengthType size = kmers.size();
		assert(size == numKmers);

		std::vector< IndexType > sortedIndexes;
		sortedIndexes.reserve(size);
		for(IndexType i = 0; i < size; i++)
			sortedIndexes.push_back(i);
		std::sort(sortedIndexes.begin(), sortedIndexes.end(), CompareArrayIdx(kmers));

		kmerMap.reserve(size);
		SequenceLengthType dupIdx = 0;
		for(SequenceLengthType i = 0; i < sortedIndexes.size(); i++) {
			IndexType sorted = sortedIndexes[i];
			int signedPosition = bools[sorted] ? (int) sorted : 0 - (int) sorted;
			LOG_DEBUG(2, "buildReferenceMap(): mapPos: " << i-dupIdx << ", sorted: " << sorted << ", signedPosition: " << signedPosition << " " << kmers[sorted].toFasta());
		
			if (i > 0 && kmers.get(sorted).compare(kmers.get( sortedIndexes[i-dupIdx-1] )) == 0) {
				dupIdx++;
				assert(kmerMap.size() == i - dupIdx + 1);
				kmerMap.valueAt(i-dupIdx)->push_back( signedPosition );
				continue;
			}
			PositionVectorPtr pvp(new PositionVector());
			pvp->reserve(2);
			pvp->push_back( signedPosition );
			kmerMap.append(kmers.get(sorted), pvp);
			assert(kmerMap.size() == i - dupIdx + 1);
		}
		STACK_DEALLOC(bools);

		kmerMap.setLastSorted();
		return kmers.size();
	}

	KmerWeightedExtensions &buildEmptyWeightedKmers(const Read &read, bool leastComplement = false, bool* bools = NULL) {
		kmers.build(read.getTwoBitSequence(), read.getLength(), leastComplement, bools);
		if (Log::isDebug(5) && bools == NULL) {
			ostream &debug = Log::Debug() << "KmerWeights: idx toFasta" << std::endl;;
			for(Kmer::IndexType i = 0 ; i < kmers.size(); i++) {
				debug << i << " " << kmers[i].toFasta() << std::endl;
			}
		}

		return kmers;
	}

	KmerWeightedExtensions &buildWeightedKmers(const Read &read, bool leastComplement = false, bool leastComplementForNegativeWeight = false) {
		if (read.isDiscarded()) {
			kmers.resize(0);
			return kmers;
		}

		SequenceLengthType readLength = read.getLength();
		STACK_ALLOC(bool, bools, readLength);
		std::string fasta = read.getFastaNoMarkup();
		int kmerLen = KmerSizer::getSequenceLength();

		kmers.build(read.getTwoBitSequence(), readLength, leastComplement, bools);
		std::string quals = read.getQuals();
		size_t markupIdx = 0;

		BaseLocationVectorType markups = read.getMarkups();
		double weight = 0.0;
		double change = 0.0;
		bool isRef = false;
		if (quals.length() > 0 && quals[0] == Read::REF_QUAL)
			isRef = true;

		SequenceLengthType size = (SequenceLengthType) kmers.size();
		assert(size == 0 || size < readLength || (size == readLength && KmerSizer::getSequenceLength() == 1));
		Extension left = Extension('X', ExtensionTracking::getMinQuality()), right;
		for (SequenceLengthType i = 0; i < size; i++) {
			if (isRef) {
				weight = 1.0;
			} else if (i % 1024 == 0 || weight == 0.0) {
				weight = 1.0;
				for (SequenceLengthType j = 0; j
				< KmerSizer::getSequenceLength(); j++)
					weight
					*= Read::qualityToProbability[(unsigned char) quals[i + j]];
			} else {
				change = Read::qualityToProbability[(unsigned char) quals[i + KmerSizer::getSequenceLength() - 1]] / Read::qualityToProbability[(unsigned char) quals[i - 1]];
				weight *= change;
			}
			while (markupIdx < markups.size() && markups[markupIdx].second < i)
				markupIdx++;
			if (markupIdx < markups.size() && markups[markupIdx].second < i
					+ KmerSizer::getSequenceLength()) {
				weight = 0.0;
			}

			// set the weight
			kmers.valueAt(i).setWeight( leastComplementForNegativeWeight && bools[i] ? weight : (0.0-weight) );

			SequenceLengthType rightBase = i + kmerLen;
			if (rightBase < fasta.length())
				right = Extension(fasta[rightBase], ((unsigned char) quals[rightBase]) - Read::FASTQ_START_CHAR);
			else
				right = Extension('X', ExtensionTracking::getMinQuality());

			// set the extensions
			if (bools[i])
				kmers.valueAt(i).setExtensions(left, right);
			else // kmer is reverse complement, so reverse extensions
				kmers.valueAt(i).setExtensions(right.getReverseComplement(), left.getReverseComplement());

			left = Extension(fasta[i], ((unsigned char) quals[i]) - Read::FASTQ_START_CHAR);
		}
		if (Log::isDebug(5)) {
			ostream &debug = Log::Debug() << "KmerWeights: idx valueAt toFasta" << std::endl;;
			for(Kmer::IndexType i = 0 ; i < kmers.size(); i++) {
				debug << i << " " << kmers.valueAt(i).getWeight() << " " << kmers[i].toFasta() << std::endl;
			}
		}

		STACK_DEALLOC(bools);

		return kmers;
	}
};

#endif

