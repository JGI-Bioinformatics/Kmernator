//
// Kmernator/src/ReadSelector.h
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

#ifndef _READ_SELECTOR_H
#define _READ_SELECTOR_H

#include <iostream>
#include <cstdlib>
#include <cstring>

#include <vector>
#include <algorithm>

#include <boost/unordered_set.hpp>
#include <boost/lexical_cast.hpp>

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "Utils.h"
#include "Log.h"

template<typename M>
class ReadSelector {
public:
	typedef Sequence::SequenceLengthType SequenceLengthType;
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	typedef ReadSet::Pair Pair;
	typedef float ScoreType;
	typedef OfstreamMap OFM;
	class ReadTrimType {
	public:
		SequenceLengthType trimOffset;
		SequenceLengthType trimLength;
		ScoreType score;
		std::string label;
		bool isAvailable;
		ReadTrimType() :
			trimOffset(0), trimLength(0), score(0.0), label(), isAvailable(true) {
		}
	};
	class PairScore {
	public:
		Pair pair;
		ScoreType score;
		PairScore() :
			pair(), score(0.0) {
		}
		PairScore(const ReadSet::Pair &_pair, ScoreType _score) :
			pair(_pair), score(_score) {
		}
		PairScore(const PairScore &copy) :
			pair(copy.pair), score(copy.score) {
		}
		bool operator<(const PairScore &cmp) const {
			if (score < cmp.score) {
				return true;
			} else if (score == cmp.score && pair.read1 < cmp.pair.read1) {
				return true;
			} else {
				return false;
			}
		}
	};
	typedef M DataType;
    typedef typename DataType::ReadPositionWeightVector ReadPositionWeightVector;
	typedef KmerMap<DataType> KMType;
	typedef KmerMap<unsigned char> KmerCountMap;
	typedef typename KMType::ConstIterator KMIterator;
	typedef typename KMType::ElementType ElementType;
	typedef KmerMap<unsigned short> KMCacheType;
	typedef std::vector< ReadTrimType > ReadTrimVector;
	typedef typename ReadSet::ReadIdxVector ReadIdxVector;
	typedef ReadIdxVector PicksVector;
	typedef std::vector< PairScore > PairScoreVector;
	typedef boost::unordered_set< std::string > DuplicateSet;
	typedef ReadSet::PairedIndexType PairedIndexType;
	typedef KmerArray<double> KA;

protected:
	const ReadSet &_reads;
	const KMType &_map;
	ReadTrimVector _trims;
	PairedIndexType _picks;
	KmerCountMap _counts;
	bool _needCounts;
	DuplicateSet _duplicateSet;
	bool _needDuplicateCheck;
	ReadSetSizeType _lastSortedPick;
	double _bimodalSigmas;

public:
	ReadSelector(const ReadSet &reads, const KMType &map):
	_reads(reads),
	_map(map),
	_trims(),
	_picks(),
	_counts(),
	_needCounts(false),
	_duplicateSet(),
	_needDuplicateCheck(false),
	_lastSortedPick(0)
	{
		_bimodalSigmas = Options::getBimodalSigmas();
	}

	void clear() {
		resetPicks();
		_trims.clear();
		_counts.clear();
		_needCounts = false;
		_needDuplicateCheck = false;
	}

	void resetPicks() {
		_picks.clear();
		_lastSortedPick = 0;
		_duplicateSet.clear();
	}

	OFM getOFM(std::string outputfile, std::string suffix = FormatOutput::getDefaultSuffix()) {
		return OFM(outputfile, suffix);
	}
	class ScoreCompare : public std::binary_function<ReadSetSizeType,ReadSetSizeType,bool>
	{
	private:
		const ReadSelector &_readSelector;
	public:
		ScoreCompare(const ReadSelector *readSelector) : _readSelector(*readSelector) {}
		inline bool operator()(const ReadSetSizeType &x, const ReadSetSizeType &y) {
			return _readSelector._trims[x].score < _readSelector._trims[y].score;
		}
	};
	class PairScoreCompare : public std::binary_function<ReadSetSizeType,ReadSetSizeType,bool>
	{
	private:
		const ReadSelector &_readSelector;
	public:
		PairScoreCompare(const ReadSelector *readSelector) : _readSelector(*readSelector) {}
		inline bool operator()(const PairScore &x, const PairScore &y) {
			return x.score < y.score;
		}
	};

	void setNeedCounts() {
		if (_needCounts) {
		  return;
		}
		_counts = KmerCountMap(_map.getNumBuckets());
		// TODO verify counts are initialized to 0
		_needCounts = true;
	}
	void setNeedDuplicateCheck() {
		if (_needDuplicateCheck)
			return;
		_needDuplicateCheck = true;
	}

	inline KA getKmersForRead(const Read &read) {
		KA kmers(read.getTwoBitSequence(), read.getLength(), true);
		return kmers;
	}
	inline KA getKmersForRead(ReadSetSizeType readIdx) {
		const Read &read = _reads.getRead(readIdx);
		return getKmersForRead(read);
	}
	inline KA getKmersForTrimmedRead(ReadSetSizeType readIdx) {
		if (_trims[readIdx].trimLength < KmerSizer::getSequenceLength())
			return KA();
		KA kmers(_reads.getRead(readIdx).getTwoBitSequence(), _trims[readIdx].trimOffset + _trims[readIdx].trimLength, true);
		if (_trims[readIdx].trimOffset > 0) {
			return kmers.copyRange(_trims[readIdx].trimOffset, _trims[readIdx].trimLength - KmerSizer::getSequenceLength() + 1);
		} else {
			return kmers;
		}
	}

protected:
	bool _testDup(ReadSetSizeType readIdx) {
		bool isGood = true;
		if (!_reads.isValidRead(readIdx))
			return false;
		std::string str = _reads.getRead(readIdx).getFasta( _trims[readIdx].trimOffset, _trims[readIdx].trimLength );
		DuplicateSet::iterator test;
		test = _duplicateSet.find(str);

		if (test != _duplicateSet.end()) {
			_trims[readIdx].isAvailable = false;
			isGood = false;
		}
		return isGood;
	}
	// Make sure only one thread can call this at a time
	bool _addDup(ReadSetSizeType readIdx) {

		bool isGood = _testDup(readIdx);
		if (! isGood )
			return false;

		std::string str = _reads.getRead(readIdx).getFasta( _trims[readIdx].trimOffset, _trims[readIdx].trimLength );

		_duplicateSet.insert(str);

		return true;
	}
	// Make sure only one thread can call this at a time
	void _storeCounts(ReadSetSizeType readIdx) {
		if (!_reads.isValidRead(readIdx))
			return;
		KA kmers = getKmersForTrimmedRead(readIdx);
		for(Kmer::IndexType j = 0; j < kmers.size(); j++) {
			_counts[ kmers[j] ]++;
		}
	}

public:
	bool isNew(ReadSetSizeType readIdx1, ReadSetSizeType readIdx2 = ReadSet::MAX_READ_IDX) {
		// check readIdx1
		if (readIdx1 != ReadSet::MAX_READ_IDX ) {
			if (! (_reads.isValidRead(readIdx1) && _trims[readIdx1].isAvailable)) {
				readIdx1 = ReadSet::MAX_READ_IDX;
			}
		}
		// check readIdx2
		if (readIdx2 != ReadSet::MAX_READ_IDX ) {
			if (! (_reads.isValidRead(readIdx2) && _trims[readIdx2].isAvailable)) {
				readIdx2 = ReadSet::MAX_READ_IDX;
			}
		}
		// order the two reads, if one is invalid
		if (readIdx1 == ReadSet::MAX_READ_IDX && readIdx2 != ReadSet::MAX_READ_IDX) {
			readIdx1 = readIdx2;
			readIdx2 = ReadSet::MAX_READ_IDX;
		}

		if (readIdx1 != ReadSet::MAX_READ_IDX) {

			// check for duplicates
			if (_needDuplicateCheck) {
				bool isGood = _testDup(readIdx1) || _testDup(readIdx2);
				if (!isGood) {
					return false;
				}
			}

			return true;

		} else {
			return false;
		}
	}
	bool isNew(const ReadSet::Pair &pair) {
		return isNew(pair.read1,pair.read2);
	}

	bool pickIfNew(ReadSetSizeType readIdx1, ReadSetSizeType readIdx2 = ReadSet::MAX_READ_IDX) {
		if (isNew(readIdx1, readIdx2)) {

			if (_needDuplicateCheck) {

				if (_addDup(readIdx1) && (readIdx2 == ReadSet::MAX_READ_IDX || _addDup(readIdx2))) {

				} else
					return false;
			}

			// pick the read
			_picks.push_back(Pair(readIdx1,readIdx2));
			_trims[readIdx1].isAvailable = false;
			if (readIdx2 != ReadSet::MAX_READ_IDX) {
				_trims[readIdx2].isAvailable = false;
			}

			// account for the picked read
			if (_needCounts) {
				_storeCounts(readIdx1);
				_storeCounts(readIdx2);
			}
			return true;

		} else {
			return false;
		}
	}
	bool pickIfNew(const ReadSet::Pair &pair) {
		return pickIfNew(pair.read1,pair.read2);
	}

	bool isPassingRead(ReadSetSizeType readIdx) {
		return _reads.isValidRead(readIdx);
	}
	bool isPassingRead(ReadSetSizeType readIdx, ScoreType minimumScore, SequenceLengthType minimumLength) {
		if (! isPassingRead(readIdx) )
		return false;
		ReadTrimType &trim = _trims[readIdx];
		if (minimumLength > _reads.getMaxSequenceLength()) {
			// use actual length of this sequence
			minimumLength = _reads.getRead(readIdx).getLength();
		}
		return trim.isAvailable && trim.score >= minimumScore && trim.trimLength >= minimumLength;
	}
	bool isPassingPair(const ReadSet::Pair &pair, ScoreType minimumScore, SequenceLengthType minimumLength, bool bothPass) {
		if (bothPass)
		return isPassingRead(pair.read1, minimumScore, minimumLength) && isPassingRead(pair.read2, minimumScore, minimumLength);
		else
		return isPassingRead(pair.read1, minimumScore, minimumLength) || isPassingRead(pair.read2, minimumScore, minimumLength);
	}
	bool isPairAvailable(const ReadSet::Pair &pair, bool bothPass) {
		if (bothPass)
		return isPassingRead(pair.read1) && _trims[pair.read1].isAvailable && isPassingRead(pair.read2) && _trims[pair.read2].isAvailable;
		else
		return (isPassingRead(pair.read1) && _trims[pair.read1].isAvailable) || (isPassingRead(pair.read2) && _trims[pair.read2].isAvailable);
	}

	int pickAllPassingReads(ScoreType minimumScore = 0.0, SequenceLengthType minimumLength = KmerSizer::getSequenceLength()) {
		int picked = 0;
		for(ReadSetSizeType i = 0; i < _reads.getSize(); i++) {
			isPassingRead(i, minimumScore, minimumLength) && pickIfNew(i) && picked++;
		}
		optimizePickOrder();
		return picked;
	}

	ReadSetSizeType pickAllPassingPairs(ScoreType minimumScore = 0.0, SequenceLengthType minimumLength = KmerSizer::getSequenceLength(), bool bothPass = false) {
		ReadSetSizeType picked = 0;
		for(ReadSetSizeType i = 0; i < _reads.getPairSize(); i++) {
			const ReadSet::Pair &pair = _reads.getPair(i);
			if (isPassingPair(pair, minimumScore, minimumLength, bothPass)) {
				pickIfNew(pair.read1) && picked++;
				pickIfNew(pair.read2) && picked++;
			}
		}
		optimizePickOrder();
		return picked;
	}

	bool rescoreByBestCoveringSubset(ReadSetSizeType readIdx, unsigned char maxPickedKmerDepth, ReadTrimType &trim) {
		KA kmers = getKmersForTrimmedRead(readIdx);
		ScoreType score = 0.0;
		for(SequenceLengthType j = 0; j < kmers.size(); j++) {
			ScoreType contribution = getValue(kmers[j]);
			if (contribution > 0) {
				KmerCountMap::ValueType pickedCount = 0;
				if (!_counts.getValueIfExists(kmers[j], pickedCount)) {
					#pragma omp critical (ReadSelector_counts)
					{
						pickedCount =  _counts.getOrSetElement(kmers[j], pickedCount).value();
					}
				}

				if ( pickedCount >= maxPickedKmerDepth ) {
					trim.score = -1.0;
					return false;
				} else {
					score += contribution * (maxPickedKmerDepth - pickedCount);
				}
			} else {
				throw "Reads should have already been trimmed to exclude this kmer: " + kmers[j].toFasta();
			}
		}

		bool hasNotChanged = (score > 0) && ( (score * 1.0001) >= trim.score);
		trim.score = score / (double) trim.trimLength;
		return hasNotChanged;
	}

	bool rescoreByBestCoveringSubset(const ReadSet::Pair &pair, unsigned char maxPickedKmerDepth, ScoreType &score) {
		score = 0.0;
		bool hasNotChanged = true;
		double len = 0;
		if (_reads.isValidRead(pair.read1)) {
			ReadTrimType &trim = _trims[pair.read1];
			hasNotChanged &= rescoreByBestCoveringSubset(pair.read1, maxPickedKmerDepth, trim);
			score = trim.score * trim.trimLength;
			len = trim.trimLength;
		}
		if (_reads.isValidRead(pair.read2)) {
			ReadTrimType &trim = _trims[pair.read2];
			hasNotChanged &= rescoreByBestCoveringSubset(pair.read2, maxPickedKmerDepth, trim);
			if (score > 0) {
				if (trim.score > 0) {
					score += trim.score * trim.trimLength;
					len += trim.trimLength;
				} else {
					score = trim.score;
				}
			}
		}
		if (len > 0.0)
			score /= len;
		return hasNotChanged;
	}

	void _initPickBestCoveringSubset()
	{
		setNeedCounts();
		setNeedDuplicateCheck();
	}

	ReadSetSizeType pickBestCoveringSubsetPairs(unsigned char maxPickedKmerDepth, ScoreType minimumScore = 0.0, SequenceLengthType minimumLength = Options::getMinReadLength(), bool bothPass = false) {
		_initPickBestCoveringSubset();
		ReadSetSizeType picked = 0;

		int numThreads = omp_get_max_threads();

		PairScore bestPairs[numThreads];
		Pair fauxPair(ReadSet::MAX_READ_IDX, ReadSet::MAX_READ_IDX);
		PairScore fauxPairScore( fauxPair, -2.0);


		PairScoreVector heapedPairs[numThreads];
		long heapSize = 0;
		int allAreDone = 0;

		long pairsSize = _reads.getPairSize();
		#pragma omp parallel
		{
			int threadId = omp_get_thread_num();

			bestPairs[threadId].score = -3.0;
			heapedPairs[threadId].resize(0);
			heapedPairs[threadId].reserve(pairsSize / numThreads / 10);

			for(long pairIdx = threadId; pairIdx < pairsSize; pairIdx += numThreads) {
				const ReadSet::Pair &pair = _reads.getPair(pairIdx);
				ScoreType score;
				if (isPairAvailable(pair, bothPass)) {
					rescoreByBestCoveringSubset(pair, maxPickedKmerDepth, score);
					if (score > minimumScore && isPassingPair(pair, minimumScore, minimumLength, bothPass)) {
						heapedPairs[threadId].push_back( PairScore( pair, score ) );
					}
				}
			}

			#pragma omp atomic
			heapSize += heapedPairs[threadId].size();

			#pragma omp barrier
			LOG_VERBOSE_OPTIONAL(1, threadId == 0, "building heap out of " << heapSize << " pairs" );

			std::make_heap(heapedPairs[threadId].begin(), heapedPairs[threadId].end(), PairScoreCompare(this));

			LOG_VERBOSE_OPTIONAL(1, threadId == 0, "picking pairs at depth: " << (int) maxPickedKmerDepth);

			// pick pairs
			long iterations = 0;
			bool threadIsDone = false;
			while (allAreDone < numThreads) {
				if (threadId == 0) {
					if (++iterations % 1000 == 0) {
						LOG_VERBOSE_OPTIONAL(1, true, "Processing heap size " << heapSize << " picked " << picked);
					} else {
						LOG_DEBUG(3, "heap size " << heapSize << " picked " << picked);
					}
				}

				PairScore pairScore;
				bool isEmpty = false;
				if (!heapedPairs[threadId].empty()) {
					pairScore = heapedPairs[threadId].front();
					std::pop_heap(heapedPairs[threadId].begin(), heapedPairs[threadId].end(), PairScoreCompare(this));
					heapedPairs[threadId].pop_back();

					#pragma omp atomic
					heapSize--;

				} else {
					pairScore = fauxPairScore;
					isEmpty = true;
				}


				if ( isEmpty || rescoreByBestCoveringSubset(pairScore.pair, maxPickedKmerDepth, pairScore.score) ) {
					if ( isEmpty || isNew(pairScore.pair)) {

						// spin until all threads are ready
						bestPairs[threadId] = pairScore;
						PairScore bestPair = fauxPairScore;
						LOG_DEBUG(4, "chose " << pairScore.pair.read1 << ", " << pairScore.pair.read2 << ": " << pairScore.score);

						#pragma omp barrier

						for(int i = 0 ; i < numThreads; i++) {
							if (bestPair < bestPairs[i] ) {
								bestPair = bestPairs[i];
							}
						}

						if (bestPair.pair == pairScore.pair && bestPair.pair != fauxPair && !isEmpty) {

							pickIfNew(pairScore.pair);
							#pragma omp atomic
							picked++;
							LOG_DEBUG(4, "Selected pair: " << pairScore.pair.read1 << " " << pairScore.pair.read2);

						}

						if (isEmpty && !threadIsDone ) {
							#pragma omp atomic
							allAreDone++;
							threadIsDone = true;
						}

						#pragma omp barrier

					} else {
						if (pairScore.score > minimumScore && isPassingPair(pairScore.pair, minimumScore, minimumLength, bothPass)) {
							LOG_DEBUG(4, "replacing Pair(" << pairScore.pair.read1 << ", " << pairScore.pair.read2 << "): "
									<< pairScore.score << " " << (isPassingRead(pairScore.pair.read1) ? _trims[pairScore.pair.read1].score : -2.0) << " "
									<< (isPassingRead(pairScore.pair.read2) ? _trims[pairScore.pair.read2].score : -2.0));
							#pragma omp atomic
							heapSize++;

							heapedPairs[threadId].push_back(pairScore);
							std::push_heap(heapedPairs[threadId].begin(), heapedPairs[threadId].end(), PairScoreCompare(this));
						}
					}
				}

			}
			LOG_DEBUG(3, "Finished picking: " << allAreDone);
		}
		LOG_VERBOSE(1, "Picked " << picked);
		optimizePickOrder();
		return picked;
	}

	ReadSetSizeType pickBestCoveringSubsetReads(unsigned char maxPickedKmerDepth, ScoreType minimumScore = 0.0, SequenceLengthType minimumLength = Options::getMinReadLength()) {
		_initPickBestCoveringSubset();
		ReadSetSizeType picked = 0;

		LOG_VERBOSE(1, "initializing reads into a heap");
		// initialize heap of reads
		PicksVector heapedReads;
		for(ReadSetSizeType readIdx = 0; readIdx < _reads.getSize(); readIdx++) {
			ReadTrimType &trim = _trims[readIdx];
			if (trim.isAvailable) {
				rescoreByBestCoveringSubset(readIdx, maxPickedKmerDepth, trim);
				if (isPassingRead(readIdx, minimumScore, minimumLength)) {
					heapedReads.push_back(readIdx);
				}
			}
		}
		LOG_VERBOSE(1, "building heap out of " << heapedReads.size() << " reads");
		std::make_heap( heapedReads.begin(), heapedReads.end(), ScoreCompare(this));

		LOG_VERBOSE(1, "picking reads: ");

		// pick reads
		while (heapedReads.begin() != heapedReads.end()) {
			ReadSetSizeType readIdx = heapedReads.front();
			std::pop_heap(heapedReads.begin(), heapedReads.end(), ScoreCompare(this));
			heapedReads.pop_back();

			ReadTrimType &trim = _trims[readIdx];
			if ( rescoreByBestCoveringSubset(readIdx, maxPickedKmerDepth, trim) ) {
				pickIfNew(readIdx) && picked++;
			} else {
				if (trim.score > minimumScore) {
					heapedReads.push_back(readIdx);
					std::push_heap(heapedReads.begin(), heapedReads.end(), ScoreCompare(this));
				}
			}
		}
		LOG_VERBOSE(1,  picked );
		optimizePickOrder();
		return picked;
	}

	inline ScoreType getValue( const Kmer &kmer ) {
		const ElementType elem = _map.getElementIfExists(kmer);
		if (elem.isValid()) {
			return elem.value().getCount();
		} else {
			return ScoreType(0);
		}
	}

	void trimReadByMarkupLength(const Read &read, ReadTrimType &trim, SequenceLengthType markupLength) {
		  if ( markupLength == 0 ) {
			  trim.trimLength = read.getLength();
		  } else {
			  // TODO find longest stretch...
			  // trim at first N or X markup
			  trim.trimOffset = 0;
			  trim.trimLength = markupLength - 1;
		  }
		  trim.score = trim.trimLength;
	}

	template<typename U>
	void trimReadByMinimumKmerScore(double minimumKmerScore, ReadTrimType &trim, U buffBegin, U buffEnd) {

		ReadTrimType test, best;
		std::stringstream ss;
		U it = buffBegin;
		while (it != buffEnd) {
			ScoreType score = *(it++);
			if (Log::isDebug(2))
				ss << (score > 1 ? (int) log(score) : (int) 0) << " ";
			if (score >= minimumKmerScore) {
				test.trimLength++;
				test.score += score;
			} else {
				if (test.score > best.score) {
					best = test;
				}
				test.score = 0;
				test.trimOffset += test.trimLength + 1;
				test.trimLength = 0;
			}
		}
		if (test.score > best.score) {
			best = test;
		}
		if (Log::isDebug(2))
			trim.label += ss.str();

		if (best.trimLength >= 3 && _bimodalSigmas >= 0.0) {
			Statistics::MeanStdCount f, s;
			U begin = buffBegin + best.trimOffset;
			U end   = begin + best.trimLength;
			U p = Statistics::findBimodalPartition(_bimodalSigmas, f, s, begin, end);
			if (p != end) {
				std::string label("Bimodal@" + boost::lexical_cast<std::string>( (p - begin)+KmerSizer::getSequenceLength() )
						+ ":" + boost::lexical_cast<std::string>( (int) f.mean )
						+ "/" + boost::lexical_cast<std::string>( (int) s.mean ));

				if (f.mean > s.mean) {
					// remove the second partition from the original trim estimate
					for(it = p; it != end; it++) {
						best.score -= (ScoreType) *it;
						best.trimLength--;
					}
					if (!trim.label.empty())
						trim.label += " ";
					trim.label += label;
				} else {
					// second partition is greater than first, remove the first
					for(it = begin; it != p ; it++) {
						best.score -= (ScoreType) *it;
						best.trimLength--;
						best.trimOffset++;
					}
					if (!trim.label.empty())
						trim.label += " ";
					trim.label += "Inv" + label;
				}
			}
		}
		trim.score = best.score;
		trim.trimOffset = best.trimOffset;
		trim.trimLength = best.trimLength;

	};
	void setTrimHeaders(ReadTrimType &trim, bool useKmers) {
		double reportScore;
		if (trim.trimLength > 0) {
			// calculate average score (before adding kmer length)
			reportScore= trim.score /= (ScoreType) trim.trimLength;
			if (useKmers) {
				trim.trimLength += KmerSizer::getSequenceLength() - 1;
			}
		} else {
			// keep available so that pairs will be selected together
			trim.trimOffset = 0;
			trim.score = -1.0;
			reportScore = 0.0;
		}
		std::stringstream ss;
		if (!trim.label.empty())
			ss << " ";
		ss << "Trim:" << trim.trimOffset << "+" << trim.trimLength
				<< " Score:" << std::fixed << std::setprecision(2) << reportScore;
		trim.label += ss.str();
	}

	void _setNumKmers( SequenceLengthType markupLength, SequenceLengthType &numKmers) {
		  if ( markupLength != 0 ) {
			  // find first N or X markup and that is the maximum trim point
			  SequenceLengthType maxTrimPoint = markupLength;
			  if (maxTrimPoint > KmerSizer::getSequenceLength()) {
				  numKmers = maxTrimPoint - KmerSizer::getSequenceLength();
			  } else {
				  numKmers = 0;
			  }
		  }
	}
	void scoreReadByKmers(const Read &read, SequenceLengthType markupLength, ReadTrimType &trim, double minimumKmerScore) {
		KA kmers = getKmersForRead(read);
		SequenceLengthType numKmers = kmers.size();

		_setNumKmers(markupLength, numKmers);
		assert(numKmers <= kmers.size());

		SequenceLengthType j = 0;
		for(; j < numKmers; j++) {
			ScoreType score = getValue(kmers[j]);
			if(score >= minimumKmerScore)
				kmers.valueAt(j) = score;
			else
				kmers.valueAt(j) = 0.0;
		}

		trimReadByMinimumKmerScore(minimumKmerScore, trim, kmers.beginValue(), kmers.beginValue() + numKmers);
	}

	virtual void scoreAndTrimReads(ScoreType minimumKmerScore) {
		_trims.resize(_reads.getSize());
		bool useKmers = Options::getKmerSize() != 0;

		long readsSize = _reads.getSize();
		#pragma omp parallel for schedule(dynamic)
		for(long i = 0; i < readsSize; i++) {
			ReadTrimType &trim = _trims[i];
			const Read &read = _reads.getRead(i);
			if (read.isDiscarded()) {
				continue;
			}
			Sequence::BaseLocationVectorType markups = read.getMarkups();
			SequenceLengthType markupLength = TwoBitSequence::firstMarkupNorX(markups);

			if (useKmers) {
				scoreReadByKmers(read, markupLength, trim, minimumKmerScore);
			} else { // !useKmers
				trimReadByMarkupLength(read, trim, markupLength);
			}
			setTrimHeaders(trim, useKmers);
		}
	}

	ReadSetSizeType pickAllCovering() {
		optimizePickOrder();
		throw;
	}

	ReadSetSizeType pickCoverageNormalizedSubset() {
		optimizePickOrder();
		throw;
	}

	// sorts the latest batch of picks to work best with mmaped ReadSets
	void optimizePickOrder(ReadSetSizeType offset = ReadSet::MAX_READ_IDX) {
	    if (offset == ReadSet::MAX_READ_IDX) {
	        offset = _lastSortedPick;
	    }
		if (offset >= (ReadSetSizeType) _picks.size())
			return;

		std::sort(_picks.begin() + offset, _picks.end());
		_lastSortedPick = _picks.size();
	}

	long _intendedWriteSize(ReadSetSizeType readIdx, const ReadTrimType &trim, int format = Options::getFormatOutput()) const {
		return _reads.getRead(readIdx).getIntendedWriteSize(trim.trimLength, trim.label, format);
	}
	std::ostream &_writePickRead(std::ostream &os, ReadSetSizeType readIdx, int format = Options::getFormatOutput()) const {
		if (readIdx == ReadSet::MAX_READ_IDX)
			return os;
		const ReadTrimType &trim = _trims[readIdx];
		return _writePickRead(os, readIdx, trim, format);
	}
	std::ostream &_writePickRead(std::ostream &os, ReadSetSizeType readIdx, const ReadTrimType &trim, int format = Options::getFormatOutput()) const {
		return _reads.write(os, readIdx, trim.trimOffset, trim.trimLength, trim.label, format);
	}
	std::ostream &writePick(std::ostream &os, ReadSetSizeType pickIdx, int format = Options::getFormatOutput()) const {
		const Pair &pair = _picks[pickIdx];
		_writePickRead(os, pair.read1, format);
		_writePickRead(os, pair.read2, format);
		return os;
	}

	void writePicks(OFM &ofstreamMap, ReadSetSizeType offset = 0, bool byInputFile = (Options::getSeparateOutputs() == 1), int format = Options::getFormatOutput() ) const {
		_writePicks(ofstreamMap, offset, _picks.size() - offset, byInputFile, format);
	}
	void _writePicks(OFM &ofstreamMap, ReadSetSizeType offset, ReadSetSizeType length, bool byInputFile, int format) const {
		for(ReadSetSizeType pickIdx = offset; pickIdx < length + offset; pickIdx++) {
			const Pair &pair = _picks[pickIdx];
			writePick(ofstreamMap, pair.read1, byInputFile, format);
			writePick(ofstreamMap, pair.read2, byInputFile, format);
		}
	}
	void writePick(OFM &ofstreamMap, ReadSetSizeType readIdx, bool byInputFile = (Options::getSeparateOutputs() == 1), int format = Options::getFormatOutput()) const {
		if (readIdx == ReadSet::MAX_READ_IDX)
			return;
		const ReadTrimType &trim = _trims[ readIdx ];
		std::string key;
		if (byInputFile) {
			key += "-" + _reads.getReadFileNamePrefix(readIdx);
		}
		LOG_DEBUG(4, "Writing " << readIdx << " to *" << key);
		_writePickRead(ofstreamMap.getOfstream(key), readIdx, trim, format);
	}
};

#endif


// $Log: ReadSelector.h,v $
// Revision 1.23  2010-08-18 17:50:40  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.22.8.1  2010-08-17 18:07:05  regan
// minor refactor
//
// Revision 1.22  2010-05-24 21:48:46  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.21.2.1  2010-05-20 03:42:46  regan
// optimized
//
// Revision 1.21  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.20.2.4  2010-05-12 20:46:50  regan
// bugfix in names of output files
//
// Revision 1.20.2.3  2010-05-12 18:25:10  regan
// help destructor ordering
//
// Revision 1.20.2.2  2010-05-10 17:57:41  regan
// fixing types
//
// Revision 1.20.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.20  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.19  2010-05-06 22:26:18  regan
// changed debug level for messages
//
// Revision 1.18  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
//
