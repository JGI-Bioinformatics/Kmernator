// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSelector.h,v 1.16 2010-05-01 21:57:53 regan Exp $
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

template<typename M>
class ReadSelector {
public:
	typedef Sequence::SequenceLengthType SequenceLengthType;
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	typedef ReadSet::Pair Pair;
	typedef float ScoreType;
	class ReadTrimType {
	public:
		SequenceLengthType trimLength;
		ScoreType score;
		std::string label;
		bool isAvailable;
		ReadTrimType() :
			trimLength(0), score(0.0), label(), isAvailable(true) {
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

private:
	const ReadSet &_reads;
	const KMType &_map;
	ReadTrimVector _trims;
	PairedIndexType _picks;
	KmerCountMap _counts;
	bool needCounts;
	DuplicateSet _duplicateSet;
	bool needDuplicateCheck;
	ReadSetSizeType _lastSortedPick;

public:
	ReadSelector(const ReadSet &reads, const KMType &map, ScoreType minimumKmerScore = 0.0):
	_reads(reads),
	_map(map),
	_trims(),
	_picks(),
	_counts(),
	needCounts(false),
	_duplicateSet(),
	needDuplicateCheck(false),
	_lastSortedPick(0)
	{
		scoreAndTrimReads(minimumKmerScore);
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
		if (needCounts) {
		  return;
		}
		_counts = KmerCountMap(_map.getNumBuckets());
		// TODO verify counts are initialized to 0
		needCounts = true;
	}
	void setNeedDuplicateCheck() {
		if (needDuplicateCheck)
		return;
		needDuplicateCheck = true;
	}

	inline KmerArray<char> getKmersForRead(ReadSetSizeType readIdx) {
		const Read &read = _reads.getRead(readIdx);
		KmerArray<char> kmers(read.getTwoBitSequence(), read.getLength(), true);
		return kmers;
	}
	inline KmerArray<char> getKmersForTrimmedRead(ReadSetSizeType readIdx) {
		KmerArray<char> kmers(_reads.getRead(readIdx).getTwoBitSequence(), _trims[readIdx].trimLength, true);
		return kmers;
	}

protected:
	bool _testDup(ReadSetSizeType readIdx) {
		bool isGood = true;
		if (!_reads.isValidRead(readIdx))
			return false;
		std::string str = _reads.getRead(readIdx).getFasta( _trims[readIdx].trimLength );
		DuplicateSet::iterator test = _duplicateSet.find(str);
		if (test == _duplicateSet.end()) {
			_duplicateSet.insert(str);
		} else {
			_trims[readIdx].isAvailable = false;
			isGood = false;
		}
		return isGood;
	}
	void _storeCounts(ReadSetSizeType readIdx) {
		if (!_reads.isValidRead(readIdx))
			return;
		KmerArray<char> kmers = getKmersForTrimmedRead(readIdx);
		for(unsigned long j = 0; j < kmers.size(); j++) {
			_counts[ kmers[j] ]++;
		}
	}

public:
	bool pickIfNew(ReadSetSizeType readIdx1, ReadSetSizeType readIdx2 = ReadSet::MAX_READ_IDX) {
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
			if (needDuplicateCheck) {
				bool isGood = _testDup(readIdx1) || _testDup(readIdx2);
				if (!isGood) {
					return false;
				}
			}

			// pick the read
			_picks.push_back(Pair(readIdx1,readIdx2));
			_trims[readIdx1].isAvailable = false;
			if (readIdx2 != ReadSet::MAX_READ_IDX) {
				_trims[readIdx2].isAvailable = false;
			}

			// account for the picked read
			if (needCounts) {
				_storeCounts(readIdx1);
				_storeCounts(readIdx2);
			}
			return true;
		} else
		return false;
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
		for(long i = 0; i < (long) _reads.getSize(); i++) {
			isPassingRead(i, minimumScore, minimumLength) && pickIfNew(i) && picked++;
		}
		optimizePickOrder();
		return picked;
	}

	ReadSetSizeType pickAllPassingPairs(ScoreType minimumScore = 0.0, SequenceLengthType minimumLength = KmerSizer::getSequenceLength(), bool bothPass = false) {
		ReadSetSizeType picked = 0;
		for(long i = 0; i < (long) _reads.getPairSize(); i++) {
			const ReadSet::Pair &pair = _reads.getPair(i);
			if (isPassingPair(pair, minimumScore, minimumLength, bothPass)) {
				pickIfNew(pair.read1) && picked++;
				pickIfNew(pair.read2) && picked++;
			}
		}
		optimizePickOrder();
		return picked;
	}

	bool rescoreByBestCoveringSubset(ReadSetSizeType readIdx, unsigned char maxPickedKmerDepth) {
		ReadTrimType &trim = _trims[readIdx];
		KmerArray<char> kmers = getKmersForTrimmedRead(readIdx);
		ScoreType score = 0.0;
		for(unsigned long j = 0; j < kmers.size(); j++) {
			const ElementType elem = _map.getElementIfExists(kmers[j]);
			if (elem.isValid()) {
				ScoreType contribution = elem.value().getCount();
				int pickedCount = _counts[kmers[j]];
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

		bool hasNotChanged = (score >= trim.score);
		//std::cerr << readIdx << " " << score << " " << trim.score << " " << hasNotChanged << " " << kmers.size() << " " << trim.trimLength << std::endl;
		trim.score = score;
		return hasNotChanged;
	}

	bool rescoreByBestCoveringSubset(const ReadSet::Pair &pair, unsigned char maxPickedKmerDepth, ScoreType &score) {
		score = 0.0;
		bool hasNotChanged = true;
		if (_reads.isValidRead(pair.read1)) {
			hasNotChanged &= rescoreByBestCoveringSubset(pair.read1, maxPickedKmerDepth);
			score += _trims[pair.read1].score;
		}
		if (_reads.isValidRead(pair.read2)) {
			hasNotChanged &= rescoreByBestCoveringSubset(pair.read2, maxPickedKmerDepth);
			score += _trims[pair.read2].score;
		}
		return hasNotChanged;
	}

	void _initPickBestCoveringSubset()
	{
		setNeedCounts();
		setNeedDuplicateCheck();
	}

	ReadSetSizeType pickBestCoveringSubsetPairs(unsigned char maxPickedKmerDepth, ScoreType minimumScore = 0.0, SequenceLengthType minimumLength = KmerSizer::getSequenceLength(), bool bothPass = false) {
		_initPickBestCoveringSubset();
		ReadSetSizeType picked = 0;

		PairScoreVector heapedPairs;
		for(ReadSetSizeType pairIdx = 0; pairIdx < _reads.getPairSize(); pairIdx++) {
			const ReadSet::Pair &pair = _reads.getPair(pairIdx);
			ScoreType score;
			if (isPairAvailable(pair, bothPass)) {
				rescoreByBestCoveringSubset(pair, maxPickedKmerDepth, score);
				if (isPassingPair(pair, minimumScore, minimumLength, bothPass)) {
					heapedPairs.push_back( PairScore( pair, score ) );
				}
			}
		}

		std::cerr << "building heap out of " << heapedPairs.size() << " pairs" << std::endl;
		std::make_heap(heapedPairs.begin(), heapedPairs.end(), PairScoreCompare(this));

		std::cerr << "picking pairs: ";
		// pick pairs
		while (heapedPairs.begin() != heapedPairs.end()) {
			PairScore &pairScore = heapedPairs.front();
			std::pop_heap(heapedPairs.begin(), heapedPairs.end(), PairScoreCompare(this));
			heapedPairs.pop_back();

			if ( rescoreByBestCoveringSubset(pairScore.pair, maxPickedKmerDepth, pairScore.score) ) {
				pickIfNew(pairScore.pair) && picked++;
			} else {
				if (pairScore.score > 0.0) {
					heapedPairs.push_back(pairScore);
					std::push_heap(heapedPairs.begin(), heapedPairs.end(), PairScoreCompare(this));
				}
			}
		}
		std::cerr << picked << std::endl;
		optimizePickOrder();
		return picked;
	}

	ReadSetSizeType pickBestCoveringSubsetReads(unsigned char maxPickedKmerDepth, ScoreType minimumScore = 0.0, SequenceLengthType minimumLength = KmerSizer::getSequenceLength()) {
		_initPickBestCoveringSubset();
		ReadSetSizeType picked = 0;

		std::cerr << "initializing reads into a heap" << std::endl;
		// initialize heap of reads
		PicksVector heapedReads;
		for(ReadSetSizeType readIdx = 0; readIdx < _reads.getSize(); readIdx++) {
			ReadTrimType &trim = _trims[readIdx];
			if (trim.isAvailable) {
				rescoreByBestCoveringSubset(readIdx, maxPickedKmerDepth);
				if (isPassingRead(readIdx, minimumScore, minimumLength)) {
					heapedReads.push_back(readIdx);
				}
			}
		}
		std::cerr << "building heap out of " << heapedReads.size() << " reads" << std::endl;
		std::make_heap( heapedReads.begin(), heapedReads.end(), ScoreCompare(this));

		std::cerr << "picking reads: ";
		// pick reads
		while (heapedReads.begin() != heapedReads.end()) {
			ReadSetSizeType readIdx = heapedReads.front();
			std::pop_heap(heapedReads.begin(), heapedReads.end(), ScoreCompare(this));
			heapedReads.pop_back();

			if ( rescoreByBestCoveringSubset(readIdx, maxPickedKmerDepth) ) {
				pickIfNew(readIdx) && picked++;
			} else {
				if (_trims[readIdx].score > 0.0) {
					heapedReads.push_back(readIdx);
					std::push_heap(heapedReads.begin(), heapedReads.end(), ScoreCompare(this));
				}
			}
		}
		std::cerr << picked << std::endl;
		optimizePickOrder();
		return picked;
	}

	void scoreAndTrimReads(ScoreType minimumKmerScore) {
		_trims.resize(_reads.getSize());
		bool useKmers = Options::getKmerSize() != 0;
#ifdef _USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
		for(long i = 0; i < (long) _reads.getSize(); i++) {
			ReadTrimType &trim = _trims[i];
			const Read read = _reads.getRead(i);
			Sequence::BaseLocationVectorType markups = read.getMarkups();
			if (useKmers) {
			  KmerArray<char> kmers = getKmersForRead(i);
			  SequenceLengthType numKmers = kmers.size();

			  SequenceLengthType markupLength = TwoBitSequence::firstMarkupNorX(markups);
			  if ( markupLength != 0 ) {
			    // find first N or X markup and that is the maximum trim point
				SequenceLengthType maxTrimPoint = markupLength;
				if (maxTrimPoint > KmerSizer::getSequenceLength()) {
					numKmers = maxTrimPoint - KmerSizer::getSequenceLength();
				} else {
					numKmers = 0;
				}
			  }

			  for(unsigned long j = 0; j < numKmers; j++) {
				const ElementType elem = _map.getElementIfExists(kmers[j]);
				if (elem.isValid()) {
					ScoreType score = elem.value().getCount();
					if (score >= minimumKmerScore) {
						trim.trimLength++;
						trim.score += score;
					} else
					break;
				} else
				break;
			  }
			} else { // !useKmers
			  SequenceLengthType markupLength = TwoBitSequence::firstMarkupNorX(markups);
			  if ( markupLength == 0 ) {
				  trim.trimLength = read.getLength();
			  } else {
				  // trim at first N or X markup
				  trim.trimLength = markupLength - 1;
			  }
			  trim.score = trim.trimLength;

			}
			if (trim.trimLength > 0) {
				// calculate average score (before adding kmer length)
				trim.score /= (ScoreType) trim.trimLength;
				if (useKmers) {
				  trim.trimLength += KmerSizer::getSequenceLength() - 1;
				}
				trim.label += " Trim:" + boost::lexical_cast<std::string>( trim.trimLength ) + " Score:" + boost::lexical_cast<std::string>( trim.score );
			} else {
				trim.score = -1.0;
				trim.label += " Trim:" + boost::lexical_cast<std::string>( trim.trimLength ) + " Score: 0";
				// keep available so that pairs will be selected together
			}


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
		if (offset >= _picks.size())
			return;

		std::sort(_picks.begin() + offset, _picks.end());
		_lastSortedPick = _picks.size();
	}

	std::ostream &_writePickRead(std::ostream &os, ReadSetSizeType readIdx, int format = 0) const {
		if (readIdx == ReadSet::MAX_READ_IDX)
			return os;
		const ReadTrimType &trim = _trims[readIdx];
		return _writePickRead(os, readIdx, trim, format);
	}
	std::ostream &_writePickRead(std::ostream &os, ReadSetSizeType readIdx, const ReadTrimType &trim, int format = 0) const {
		return _reads.write(os, readIdx, trim.trimLength, trim.label, format);
	}
	std::ostream &writePick(std::ostream &os, ReadSetSizeType pickIdx, int format = 0) const {
		Pair &pair = _picks[pickIdx];
		_writePickRead(os, pair.read1, format);
		_writePickRead(os, pair.read2, format);
		return os;
	}
	std::ostream &writePicks(std::ostream &os, ReadSetSizeType offset = 0, int format = 0) const {
		return writePicks(os, offset, _picks.size() - offset, format);
	}
	std::ostream &writePicks(std::ostream &os, ReadSetSizeType offset, ReadSetSizeType length, int format = 0) const {
		for(ReadSetSizeType i = offset; i < length + offset; i++) {
			writePick(os, i, format);
		}
		return os;
	}

	void writePicks(OfstreamMap &ofstreamMap, ReadSetSizeType offset = 0, bool byInputFile = true, int format = 0 ) const {
		writePicks(ofstreamMap, offset, _picks.size() - offset, byInputFile, format);
	}
	void writePicks(OfstreamMap &ofstreamMap, ReadSetSizeType offset, ReadSetSizeType length, bool byInputFile = true, int format = 0) const {
		for(ReadSetSizeType pickIdx = offset; pickIdx < length + offset; pickIdx++) {
			const Pair &pair = _picks[pickIdx];
			writePick(ofstreamMap, pair.read1, byInputFile, format);
			writePick(ofstreamMap, pair.read2, byInputFile, format);
		}
	}
	void writePick(OfstreamMap &ofstreamMap, ReadSetSizeType readIdx, bool byInputFile = true, int format = 0) const {
		if (readIdx == ReadSet::MAX_READ_IDX)
			return;
		const ReadTrimType &trim = _trims[ readIdx ];
		std::string key;
		if (byInputFile) {
			key += _reads.getReadFileNamePrefix(readIdx);
		}
		_writePickRead(ofstreamMap.getOfstream(key), readIdx, trim, format);
	}
};

#endif

