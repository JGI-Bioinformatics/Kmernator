//
// Kmernator/src/ReadSelector.h
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

class _ReadSelectorOptions : public OptionsBaseInterface {
public:
	_ReadSelectorOptions() : maxKmerDepth(-1), partitionByDepth(-1), bothPairs(1), remainderTrim(-1), minReadLength(0.45),
	    bimodalSigmas(-1.0), kmerScoringType("MAX"), normalizationMethod("RANDOM"), useLogscaleAboveMax(false), separateOutputs(true)  {
	}
	virtual ~_ReadSelectorOptions() {}
	void _resetDefaults() {
		// Other *::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("Read Selector Options"), expt("Read Selector (Experimental) Options");
		opts.add_options()
				// output read selection

				("separate-outputs", po::value<bool>()->default_value(separateOutputs), "If set, each input (plus consensus) will generate a new outputfile.  If set false, all input files will be merged into one output file.")

				("max-kmer-output-depth", po::value<int>()->default_value(maxKmerDepth), "i.e. targeted read normalization depth.  The maximum number of times a kmer will be output among the selected reads (mutually exclusive with partition-by-depth).  This is not a criteria on the kmer spectrum, just a way to reduce the redundancy of the output.  To get targeted depth, multiply maxKmerOutputDepth by: (readLength - kmerSize - 1) / readLength")

				("use-logscale-above-max", po::value<bool>()->default_value(useLogscaleAboveMax), "if --max-kmer-output-depth is set, then reads above this threshold will be reduced by the log2 kmer abundance")

				("normalization-method", po::value<std::string>()->default_value(normalizationMethod), "If --max-kmer-output-depth is selected, what algorithm to use (RANDOM, OPTIMAL) (optimal is *very* slow and is not implemented in MPI version")

				("partition-by-depth", po::value<int>()->default_value(partitionByDepth), "partition filtered reads by powers-of-two coverage depth (mutually exclusive with max-kmer-depth)")

				("min-passing-in-pair", po::value<int>()->default_value(bothPairs), "1 or 2 reads in a pair must pass filters")

				("min-read-length", po::value<float>()->default_value(minReadLength), "minimum (trimmed) read length of selected reads.  0.0: no minimum - 1.0: full read length, between 0.0 and 1.0, then the fraction of the original read length, >1.0 then the absolute read length")

				("remainder-trim", po::value<float>()->default_value(remainderTrim), "if set >=0, a final round letting single reads and lesser trimmed reads will be selected.  Same rules as --min-read-length")

				("kmer-scoring-type", po::value<std::string>()->default_value(kmerScoringType), "How to evaluate a read's score over all its kmers (SUM, MEDIAN, AVG, MIN, MAX)")

				;

		desc.add(opts);

		expt.add_options()
				("bimodal-sigmas", po::value<float>()->default_value(bimodalSigmas), "(experimental) Detect bimodal kmer-signatures across reads and trim at transition point if the two means are separated by bimodal-sigmas * stdDev (2.0 to 3.0 suggested).  disabled if < 0.0")
				;

		desc.add(expt);
		// Other *::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		setOpt("max-kmer-output-depth", maxKmerDepth);
		setOpt("use-logscale-above-max", useLogscaleAboveMax);
		setOpt("partition-by-depth", partitionByDepth);
		setOpt("min-passing-in-pair", bothPairs);
		setOpt("separate-outputs", separateOutputs);

		// set read length
		setOpt("min-read-length", minReadLength);
		setOpt("remainder-trim", remainderTrim);
		setOpt("kmer-scoring-type", kmerScoringType);
		if (kmerScoringType != "SUM" && kmerScoringType != "MEDIAN" && kmerScoringType != "AVG" && kmerScoringType != "MIN" && kmerScoringType != "MAX") {
			setOptionsErrorMsg("Invalid --kmer-scoring-type: " + kmerScoringType);
			ret = false;
		}
		setOpt("normalization-method", normalizationMethod);
		if (normalizationMethod != "RANDOM" && normalizationMethod != "OPTIMAL") {
			setOptionsErrorMsg("Invalid --normalization-method: " + normalizationMethod);
			ret = false;
		}

		// verify mutually exclusive options are not set
		if ( (getMaxKmerDepth() > 0 && getPartitionByDepth() >  0) )
		{
			ret = false;
			setOptionsErrorMsg("You can not specify both max-kmer-depth and partition-by-depth");
		}
		setOpt("bimodal-sigmas", getBimodalSigmas());

		if (maxKmerDepth > 0 && KmerBaseOptions::getOptions().getKmerSize() == 0) {
			ret = false;
			setOptionsErrorMsg("If you select max-kmer-output-depth, then you must select a valid kmer-size, not: " + boost::lexical_cast<std::string>(KmerBaseOptions::getOptions().getKmerSize()));
		}
		if (partitionByDepth > 0 && KmerBaseOptions::getOptions().getKmerSize() == 0) {
			ret = false;
			setOptionsErrorMsg("If you select partitoin-by-depth, then you must select a valid kmer-size, not: " + boost::lexical_cast<std::string>(KmerBaseOptions::getOptions().getKmerSize()));
		}
		// Other ret &= *::_parseOptions(vm);
		return ret;
	}

	int &getMaxKmerDepth() {
		return maxKmerDepth;
	}
	int &getPartitionByDepth() {
		return partitionByDepth;
	}
	int &getMinPassingInPair() {
		return bothPairs;
	}
	bool getBothPairs() {
		return bothPairs == 2;
	}
	float &getMinReadLength()
	{
		return minReadLength;
	}
	float &getRemainderTrim() {
		return remainderTrim;
	}
	float &getBimodalSigmas()
	{
		return bimodalSigmas;
	}
	std::string &getKmerScoringType() {
		return kmerScoringType;
	}
	std::string &getNormalizationMethod() {
		return normalizationMethod;
	}
	bool &getUseLogscaleAboveMax() {
		return useLogscaleAboveMax;
	}
	bool &getSeparateOutputs()
	{
		return separateOutputs;
	}



private:
	int maxKmerDepth, partitionByDepth, bothPairs;
	float remainderTrim, minReadLength;
	float bimodalSigmas;
	std::string kmerScoringType;
	std::string normalizationMethod;
	bool useLogscaleAboveMax;
	bool separateOutputs;
};
typedef OptionsBaseTemplate< _ReadSelectorOptions > ReadSelectorOptions;

class ReadSelectorUtil {
public:
	typedef Sequence::SequenceLengthType SequenceLengthType;
	static inline bool passesLength(float length, const Read &read, float minimumLength) {
		assert(minimumLength >= 0.0);
		if (length <= 1.0)
			return false;
		if (minimumLength <= 1.0) {
			return read.getLength() * minimumLength <= length;
		} else {
			return minimumLength <= length;
		}
	}
	static inline bool passesLength(float length, SequenceLengthType readLength, float minimumLength) {
		assert(minimumLength >= 0.0);
		if (length <= 1.0)
			return false;
		if (minimumLength <= 1.0) {
			return readLength * minimumLength <= length;
		} else {
			return minimumLength <= length;
		}
	}

};

template<typename Map>
class ReadSelector {
public:
	typedef Sequence::SequenceLengthType SequenceLengthType;
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	typedef ReadSet::Pair Pair;
	typedef float ScoreType;
	typedef OfstreamMap OFM;
	enum KmerScoringType {
		KS_SUM,
		KS_MEDIAN,
		KS_MIN,
		KS_MAX,
		KS_AVG,
		_KS_MAX_SCORING
	} ;
	static const char* getKmerScoringTypeLabel(enum KmerScoringType scoreType) {
		static const char* const _labels[] = {
			"Score",
			"MedianScore",
			"MinScore",
			"MaxScore",
			"AvgScore"
		};
		return _labels[scoreType];
	}

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
		std::string toString() const {
			std::stringstream ss;
			ss << "{ trimOffset: " << trimOffset << ", trimLength: " << trimLength << ", score: " << score << ", label: " << label << ", isAvailable: " << isAvailable << "}";
			return ss.str();
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
	typedef Map KMType;
	typedef typename KMType::ValueType DataType;
	typedef typename DataType::ReadPositionWeightVector ReadPositionWeightVector;
	typedef KmerMapGoogleSparse<unsigned char> KmerCountMap;
	typedef typename KMType::ConstIterator KMIterator;
	typedef typename KMType::ElementType ElementType;
	typedef std::vector< ReadTrimType > ReadTrimVector;
	typedef typename ReadSet::ReadIdxVector ReadIdxVector;
	typedef ReadIdxVector PicksVector;
	typedef std::vector< PairScore > PairScoreVector;
	typedef boost::unordered_set< std::string > DuplicateSet;
	typedef ReadSet::PairedIndexType PairedIndexType;
	typedef KmerArrayPair<double> KA;

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
	enum KmerScoringType _defaultScoringType;


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
		_lastSortedPick(0),
		_defaultScoringType(_KS_MAX_SCORING)
	{
		_bimodalSigmas = ReadSelectorOptions::getOptions().getBimodalSigmas();
		// Let the kernel know how these pages will be used
		if (Options::getOptions().getMmapInput())
			ReadSet::madviseMmapsNormal();
		setDefaultScoringType(ReadSelectorOptions::getOptions().getKmerScoringType());
	}
	void setDefaultScoringType(std::string scoringType) {
		if (scoringType == "SUM")
			_defaultScoringType = KS_SUM;
		if (scoringType == "MEDIAN")
			_defaultScoringType = KS_MEDIAN;
		if (scoringType == "AVG")
			_defaultScoringType = KS_AVG;
		if (scoringType == "MIN")
			_defaultScoringType = KS_MIN;
		if (scoringType == "MAX")
			_defaultScoringType = KS_MAX;
		if (_defaultScoringType == _KS_MAX_SCORING)
			LOG_THROW("Invalid scoring type: " << scoringType);
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

	OFM getOFM(std::string outputFile, std::string suffix = FormatOutput::getDefaultSuffix()) {
		LOG_DEBUG(1, "ReadSelector::getOFM(" << outputFile << ", " << suffix << ")");
		return OFM(outputFile, suffix);
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

	inline void getKmersForRead(const Read &read, KA &kmers) {
		kmers.build(read.getTwoBitSequence(), read.getLength(), true);
	}
	inline void getKmersForRead(ReadSetSizeType readIdx, KA &kmers) {
		const Read &read = _reads.getRead(readIdx);
		getKmersForRead(read, kmers);
	}
	inline void getKmersForTrimmedRead(ReadSetSizeType readIdx, KA &kmers) {
		if (_trims[readIdx].trimLength < KmerSizer::getSequenceLength()) {
			kmers.resize(0);
			return;
		}
		kmers.build(_reads.getRead(readIdx).getTwoBitSequence(), _trims[readIdx].trimOffset + _trims[readIdx].trimLength, true);
		if (_trims[readIdx].trimOffset > 0) {
			kmers.copyRange(_trims[readIdx].trimOffset, _trims[readIdx].trimLength - KmerSizer::getSequenceLength() + 1);
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
		KA kmers;
		getKmersForTrimmedRead(readIdx, kmers);
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
	bool isPairedRead(const ReadSet::Pair &pair) {
		return (_reads.isValidRead(pair.read1) & _reads.isValidRead(pair.read2));
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
			LOG_DEBUG(3, "pickIfNew(): Picked " << (readIdx1 != ReadSet::MAX_READ_IDX ? _reads.getRead(readIdx1).getName() : _reads.getRead(readIdx2).getName()));
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
	bool isPassingRead(ReadSetSizeType readIdx, ScoreType minimumScore, float minimumLength) {
		if (! isPassingRead(readIdx) )
			return false;
		ReadTrimType &trim = _trims[readIdx];
		bool passed = (trim.isAvailable && (trim.score >= minimumScore) && ReadSelectorUtil::passesLength(trim.trimLength, _reads.getRead(readIdx), minimumLength));
		LOG_DEBUG(3, "isPassingRead(" << readIdx << " (" << _reads.getRead(readIdx).getName() << "), " << minimumScore << ", " << minimumLength << "): " << trim.toString() << " " << passed);
		return passed;
	}
	bool isPassingPair(const ReadSet::Pair &pair, ScoreType minimumScore, float minimumLength, bool bothPass) {
		bool passed = true;
		bool r1 = isPassingRead(pair.read1, minimumScore, minimumLength);
		bool r2 = isPassingRead(pair.read2, minimumScore, minimumLength);
		if (isPairedRead(pair) && bothPass)
			passed = r1 & r2;
		else
			passed = r1 | r2;
		LOG_DEBUG(3, "isPassingPair(" << pair.read1 << " / " << pair.read2 << ", " << minimumScore << ", " << minimumLength << ", " << bothPass << "): " << passed);
		return passed;
	}
	bool isPairAvailable(const ReadSet::Pair &pair, bool bothPass) {
		if (isPairedRead(pair) && bothPass)
			return (_trims[pair.read1].isAvailable && _trims[pair.read2].isAvailable && isPassingRead(pair.read1) && isPassingRead(pair.read2));
		else
			return ((_trims[pair.read1].isAvailable && isPassingRead(pair.read1)) || (_trims[pair.read2].isAvailable && isPassingRead(pair.read2)));
	}

	int pickAllPassingReads(ScoreType minimumScore = 0.0, float minimumLength = ReadSelectorOptions::getOptions().getMinReadLength()) {
		int picked = 0;
		for(ReadSetSizeType i = 0; i < _reads.getSize(); i++) {
			isPassingRead(i, minimumScore, minimumLength) && pickIfNew(i) && picked++;
		}
		optimizePickOrder();
		return picked;
	}

	ReadSetSizeType pickAllPassingPairs(ScoreType minimumScore = 0.0, float minimumLength = ReadSelectorOptions::getOptions().getMinReadLength(), bool bothPass = false) {
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

	bool rescoreByBestCoveringSubset(ReadSetSizeType readIdx, unsigned char maxPickedKmerDepth, ReadTrimType &trim, KA &kmers) {
		getKmersForTrimmedRead(readIdx, kmers);
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
				LOG_THROW("rescoreByBestCoveringSubset(): Reads should have already been trimmed to exclude this kmer: " << kmers[j].toFasta());
			}
		}

		bool hasNotChanged = (score > 0) && ( (score * 1.0001) >= trim.score);
		trim.score = score;
		return hasNotChanged;
	}

	bool rescoreByBestCoveringSubset(const ReadSet::Pair &pair, unsigned char maxPickedKmerDepth, ScoreType &score, KA &kmers) {
		score = 0.0;
		bool hasNotChanged = true;
		double len = 0;
		if (_reads.isValidRead(pair.read1)) {
			ReadTrimType &trim = _trims[pair.read1];
			hasNotChanged &= rescoreByBestCoveringSubset(pair.read1, maxPickedKmerDepth, trim, kmers);
			score = trim.score;
			len = trim.trimLength;
		}
		if (_reads.isValidRead(pair.read2)) {
			ReadTrimType &trim = _trims[pair.read2];
			hasNotChanged &= rescoreByBestCoveringSubset(pair.read2, maxPickedKmerDepth, trim, kmers);
			if (score > 0) {
				if (trim.score > 0) {
					score += trim.score;
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

	inline bool chooseRead(long score, long targetDepth, bool useLogscale) {
		if (score <= targetDepth) {
			return true;
		} else {
			long choice = IntRand::rand() % score;
			if (useLogscale) {
				return choice <= targetDepth * log((float) score / (float) targetDepth);
			} else {
				return choice <= targetDepth;
			}
		}
	}
	ReadSetSizeType pickCoverageNormalizedSubset(long targetDepth, ScoreType minimumScore = 0.0, float minimumLength = ReadSelectorOptions::getOptions().getMinReadLength(), bool byPair = false, bool bothPass = false ) {
		LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), "pickCoverageNormalizedSubset(" << targetDepth << ", " << minimumLength << ", " << byPair << ", " << bothPass << ")");
		ReadSetSizeType picked = 0;
		int numThreads = omp_get_max_threads();
		long pairsSize = _reads.getPairSize();
		PairedIndexType myPicks[numThreads];
		const bool useLogscale = ReadSelectorOptions::getOptions().getUseLogscaleAboveMax();

		#pragma omp parallel for schedule(guided)
		for(long pairIdx = 0; pairIdx < pairsSize; pairIdx++) {
			const ReadSet::Pair &pair = _reads.getPair(pairIdx);
			long score1, score2;
			score1 = (long) (isPassingRead(pair.read1, minimumScore, minimumLength) ? _trims[pair.read1].score : -1);
			score2 = (long) (isPassingRead(pair.read2, minimumScore, minimumLength) ? _trims[pair.read2].score : -1);
			LOG_DEBUG(3, "Scores: " << score1 << " , " << score2 << " for " << _reads.getRead(pair.read1).getName());
			if (byPair) {
				if (!isPassingPair(pair, minimumScore, minimumLength, bothPass)) {
					LOG_DEBUG(4, "Both do not pass");
					continue;
				}
			}

			if (byPair) {
				long pairedscore;
				if (bothPass)  {
					if (score1 <= 0 || score2 <= 0) {
						LOG_DEBUG(4, "At least one has score <= 0: " << score1 << ", " << score2);
						continue;
					}
				}
				if (score1 <= 0 && score2 <= 0) {
					LOG_DEBUG(4, "Both do not have >0 scores: " << score1 << ", " << score2)
								continue;
				}
				pairedscore = std::max(score1, score2);
				if (chooseRead(pairedscore, targetDepth, useLogscale)) {
					LOG_DEBUG(3, "pCNS(): Picked Pair " << pairedscore << " " << _reads.getRead(pair.read1).getName());
					myPicks[omp_get_thread_num()].push_back(pair);
				} else {
					LOG_DEBUG(3, "pCNS(): Did not pick pair: " << pairedscore << " , " << score1 << " " << score2);
				}

			} else {
				Pair twoReads;
				LOG_DEBUG(4, "single reads");
				if (score1 > 0.0 && (chooseRead(score1, targetDepth, useLogscale))) {
					LOG_DEBUG(3, "pCNS(): Picked Single " << score1 << " " << _reads.getRead(pair.read1).getName());
					twoReads.read1 = pair.read1;
				} else {
					LOG_DEBUG(3, "pCNS(): Did not pick r1:" << score1);
				}
				if (score2 > 0.0 && (chooseRead(score2, targetDepth, useLogscale))) {
					LOG_DEBUG(3, "pCNS(): Picked Single " << score2 << " " << _reads.getRead(pair.read2).getName());
					twoReads.read2 = pair.read2;
				} else {
					LOG_DEBUG(3, "pCNS(): Did not pick r2: " << score2);
				}

				if (twoReads.hasAValidRead())
					myPicks[omp_get_thread_num()].push_back(twoReads);
			}
		}
		LOG_DEBUG_OPTIONAL(1, true, "Consolidating thread-picked reads");
		long newPicks = 0;
		for(int i = 0; i < numThreads; i++) {
			newPicks += myPicks[i].size();
		}
		_picks.reserve(_picks.size() + newPicks);
		for(int i = 0; i < numThreads; i++) {
			for(PairedIndexType::iterator it = myPicks[i].begin(); it != myPicks[i].end(); it++)
				if (pickIfNew(*it))
					picked++;
			myPicks[i].clear();
		}
		optimizePickOrder();
		return picked;
	}

	ReadSetSizeType pickBestCoveringSubsetPairs(unsigned char maxPickedKmerDepth, ScoreType minimumScore = 0.0,
			float minimumLength = ReadSelectorOptions::getOptions().getMinReadLength(), bool bothPass = false) {
		_initPickBestCoveringSubset();

		ReadSetSizeType picked = 0;

		int numThreads = omp_get_max_threads();

		PairScore bestPairs[numThreads];
		Pair fauxPair(ReadSet::MAX_READ_IDX, ReadSet::MAX_READ_IDX);
		PairScore fauxPairScore( fauxPair, -2.0);

		PairScoreVector heapedPairs[numThreads];
		long heapSize = 0;
		int allAreDone = 0;

		std::vector< KA > _kmers(omp_get_max_threads(), KA());
		long pairsSize = _reads.getPairSize();
#pragma omp parallel
		{
			int threadId = omp_get_thread_num();
			KA &kmers = _kmers[threadId];
			bestPairs[threadId].score = -3.0;
			heapedPairs[threadId].resize(0);
			heapedPairs[threadId].reserve(pairsSize / numThreads / 10);

			for(long pairIdx = threadId; pairIdx < pairsSize; pairIdx += numThreads) {
				const ReadSet::Pair &pair = _reads.getPair(pairIdx);
				ScoreType score;
				if (isPairAvailable(pair, bothPass)) {
					rescoreByBestCoveringSubset(pair, maxPickedKmerDepth, score, kmers);
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


				if ( isEmpty || rescoreByBestCoveringSubset(pairScore.pair, maxPickedKmerDepth, pairScore.score, kmers) ) {
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

	ReadSetSizeType pickBestCoveringSubsetReads(unsigned char maxPickedKmerDepth, ScoreType minimumScore = 0.0, float minimumLength = ReadSelectorOptions::getOptions().getMinReadLength()) {
		_initPickBestCoveringSubset();
		ReadSetSizeType picked = 0;

		LOG_VERBOSE(1, "initializing reads into a heap");
		// initialize heap of reads
		PicksVector heapedReads;
		KA kmers;
		for(ReadSetSizeType readIdx = 0; readIdx < _reads.getSize(); readIdx++) {
			ReadTrimType &trim = _trims[readIdx];
			if (trim.isAvailable) {
				rescoreByBestCoveringSubset(readIdx, maxPickedKmerDepth, trim, kmers);
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
			if ( rescoreByBestCoveringSubset(readIdx, maxPickedKmerDepth, trim, kmers) ) {
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
				ss <<  (char) ('@' + (score > 1 ? (int) log(score) : (int) 0));
			if (score >= minimumKmerScore) {
				test.trimLength++;
				test.score += 1;
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
		if (Log::isDebug(2)) {
			if (!trim.label.empty())
				trim.label += Read::LABEL_SEP;
			trim.label += "LS:" + ss.str();
		}

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
					assert(end - p >= 0);
					// remove the second partition from the original trim estimate
					best.trimLength -= end - p;
					if (!trim.label.empty())
						trim.label += Read::LABEL_SEP;
					trim.label += label;
				} else {
					assert(p - begin >= 0);
					// second partition is greater than first, remove the first
					best.trimLength -= p - begin;
					best.trimOffset += p - begin;
					if (!trim.label.empty())
						trim.label += Read::LABEL_SEP;
					trim.label += "Inv" + label;
				}
			}
		}
		trim.trimOffset = best.trimOffset;
		trim.trimLength = best.trimLength;

	};
	void setTrimHeaders(ReadTrimType &trim, bool useKmers) {
		if (trim.trimLength > 0) {
			if (useKmers) {
				trim.trimLength += KmerSizer::getSequenceLength() - 1;
			}
		} else {
			// keep available so that pairs will be selected together
			trim.trimOffset = 0;
			trim.score = -1.0;
		}
		std::stringstream ss;
		if (!trim.label.empty())
			ss << Read::LABEL_SEP;
		ss << "Trim:" << trim.trimOffset << "+" << trim.trimLength << Read::Read::LABEL_SEP << getKmerScoringTypeLabel(_defaultScoringType) << ":" << (int) (trim.score+0.5);
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
	void scoreReadByKmers(const Read &read, SequenceLengthType markupLength, ReadTrimType &trim, double minimumKmerScore, KA &kmers) {
		getKmersForRead(read, kmers);
		setKmerValues(kmers, minimumKmerScore);
		LOG_DEBUG_OPTIONAL(5, true, "Trim and Score: " << read.getName());

		trimReadByKmers(kmers.beginValue(), kmers.endValue(), markupLength, trim, minimumKmerScore);

		if (trim.trimLength > 0)
			scoreReadByScoringType(kmers.beginValue() + trim.trimOffset, kmers.beginValue() + trim.trimOffset + trim.trimLength, trim, _defaultScoringType);
		else
			trim.score = -1;
	}

	void setKmerValues(KA &kmers, double minimumKmerScore) {
		setKmerValues(kmers.begin(), kmers.end(), minimumKmerScore);
	}
	template<typename IT>
	void setKmerValues(IT begin, IT end, double minimumKmerScore) {
		for(IT it = begin; it != end; it++) {
			ScoreType score = getValue(it->key());
			if (score >= minimumKmerScore)
				it->value() = score;
			else
				it->value() = 0.0;
		}
	}

	template<typename IT>
	void trimReadByKmers(IT begin, IT end, SequenceLengthType markupLength, ReadTrimType &trim, double minimumKmerScore) {

		SequenceLengthType numKmers = end-begin;

		_setNumKmers(markupLength, numKmers);

		trimReadByMinimumKmerScore(minimumKmerScore, trim, begin, begin + numKmers);
	};

	void scoreReadByScoringType(KA &kmers, ReadTrimType &trim, enum KmerScoringType scoringType) {
		scoreReadByScoringType(kmers.beginValue(), kmers.endValue(), trim, scoringType);
	}

	template<typename IT>
	void scoreReadByScoringType(IT begin, IT end,  ReadTrimType &trim, enum KmerScoringType scoringType) {
		if (begin == end) {
			trim.score = -1;
			return;
		}
		switch(scoringType) {
		case KS_SUM:
			scoreReadBySumKmer(begin, end, trim); break;
		case KS_MEDIAN:
			scoreReadByMedianKmer(begin, end, trim); break;
		case KS_AVG:
			scoreReadBySumKmer(begin, end, trim, true); break;
		case KS_MIN:
		case KS_MAX:
			scoreReadByLimit(begin, end, trim, scoringType);
			break;
		default:
			LOG_THROW("Invalid scoring type!");
		}
		if (Log::isDebug(5)) {
			std::stringstream ss;
			ss << "Scores: " << " FinalScore: " << trim.score << ". ";
			for(IT it = begin; it != end ; it++)
				ss << *it << ", ";
			std::string s = ss.str();
			LOG_DEBUG_OPTIONAL(5, true, s);
		}
	};

	template<typename IT>
	void scoreReadByMedianKmer(IT begin, IT end, ReadTrimType &trim) {

		// find the median score
		std::vector< ScoreType > scores;
		scores.reserve(end-begin);
		for(IT it = begin; it != end; it++)
			scores.push_back(*it);
		std::sort(scores.begin(), scores.end());
		trim.score = scores[scores.size()/2];
		if (Log::isDebug(5)) {
			std::stringstream ss;
			ss << "Median Calculation: ";
			int c = 0;
			for(IT it = begin; it != end ; it++)
				ss << c++ << ":" << *it << ", ";
			c = 0;
			ss << std::endl;
			for(std::vector<ScoreType>::iterator it2 = scores.begin(); it2 != scores.end(); it2++)
				ss << c++ << ":" << *it2 << ", ";
			ss << std::endl;
			ss << "Final: " << scores.size() << " " << scores.size()/2 << " " << trim.score;
			std::string s = ss.str();
			LOG_DEBUG_OPTIONAL(5, true, s);
		}

	};

	template<typename IT>
	void scoreReadBySumKmer(IT begin, IT end, ReadTrimType &trim, bool byAvg = false) {

		double sum = 0.0;
		int count = 0;
		for(IT it = begin; it != end; it++) {
			sum += *it;
			count++;
		}

		if (byAvg)
			trim.score = sum / (count > 0 ? count : 1);
	};

	template<typename IT>
	void scoreReadByLimit(IT begin, IT end, ReadTrimType &trim, enum KmerScoringType scoringType) {

		IT it = begin;
		if (it == end)
			return;
		ScoreType limitScore = *it;
		it++;
		for(; it != end; it++) {
			ScoreType val = *it;
			if (scoringType == KS_MAX)
				limitScore = limitScore > val ? limitScore : val;
			else
				limitScore = limitScore < val ? limitScore : val;
		}
		trim.score = limitScore;
	};

	virtual void scoreAndTrimReads(ScoreType minimumKmerScore) {

		_trims.resize(_reads.getSize());
		bool useKmers = KmerBaseOptions::getOptions().getKmerSize() != 0;

		long readsSize = _reads.getSize();
		std::vector< KA > _kmers(omp_get_max_threads(), KA());

#pragma omp parallel for schedule(guided)
		for(long i = 0; i < readsSize; i++) {
			KA &kmers = _kmers[omp_get_thread_num()];
			ReadTrimType &trim = _trims[i];
			const Read &read = _reads.getRead(i);
			if (read.isDiscarded()) {
				continue;
			}
			Sequence::BaseLocationVectorType markups = read.getMarkups();
			SequenceLengthType markupLength = TwoBitSequence::firstMarkupNorX(markups);

			if (useKmers) {
				scoreReadByKmers(read, markupLength, trim, minimumKmerScore, kmers);
			} else { // !useKmers
				trimReadByMarkupLength(read, trim, markupLength);
			}
			setTrimHeaders(trim, useKmers);
		}
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

	long _intendedWriteSize(ReadSetSizeType readIdx, const ReadTrimType &trim, FormatOutput format = FormatOutput::getDefault()) const {
		return _reads.getRead(readIdx).getIntendedWriteSize(trim.trimLength, trim.label, format);
	}
	std::ostream &_writePickRead(std::ostream &os, ReadSetSizeType readIdx, FormatOutput format = FormatOutput::getDefault()) const {
		if (readIdx == ReadSet::MAX_READ_IDX)
			return os;
		const ReadTrimType &trim = _trims[readIdx];
		return _writePickRead(os, readIdx, trim, format);
	}
	std::ostream &_writePickRead(std::ostream &os, ReadSetSizeType readIdx, const ReadTrimType &trim, FormatOutput format = FormatOutput::getDefault()) const {
		return _reads.write(os, readIdx, trim.trimOffset, trim.trimLength, trim.label, format);
	}
	std::ostream &writePick(std::ostream &os, ReadSetSizeType pickIdx, FormatOutput format = FormatOutput::getDefault()) const {
		const Pair &pair = _picks[pickIdx];
		_writePickRead(os, pair.read1, format);
		_writePickRead(os, pair.read2, format);
		return os;
	}

	void writePicks(OFM &ofstreamMap, ReadSetSizeType offset = 0, bool byInputFile = ReadSelectorOptions::getOptions().getSeparateOutputs(), FormatOutput format = FormatOutput::getDefault()) const {
		_writePicks(ofstreamMap, offset, _picks.size() - offset, byInputFile, format);
	}
	void _writePicks(OFM &ofstreamMap, ReadSetSizeType offset, ReadSetSizeType length, bool byInputFile, FormatOutput format = FormatOutput::getDefault()) const {
		for(ReadSetSizeType pickIdx = offset; pickIdx < length + offset; pickIdx++) {
			const Pair &pair = _picks[pickIdx];
			writePick(ofstreamMap, pair.read1, byInputFile, format);
			writePick(ofstreamMap, pair.read2, byInputFile, format);
		}
	}
	void writePick(OFM &ofstreamMap, ReadSetSizeType readIdx, bool byInputFile = ReadSelectorOptions::getOptions().getSeparateOutputs(), FormatOutput format = FormatOutput::getDefault()) const {
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
