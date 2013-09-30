//
// Kmernator/src/ReadSet.h
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

#ifndef _READ_SET_H
#define _READ_SET_H
#include <cstring>
#include <boost/unordered_map.hpp>
#include <sys/mman.h>
#include <deque>

#include "config.h"
#include "Options.h"
#include "Sequence.h"
#include "Utils.h"
#include "ReadFileReader.h"
#include "Log.h"

class ReadSet {
public:
	typedef Kmernator::ReadSetSizeType ReadSetSizeType;
	typedef std::vector< ReadSetSizeType > ReadIdxVector;
	typedef ReadFileReader::MmapSource MmapSource;
	typedef std::pair<MmapSource,MmapSource> MmapSourcePair;
	typedef std::vector< MmapSourcePair > MmapSourceVector;
	typedef ReadFileReader::SequenceStreamParser SequenceStreamParser;
	typedef ReadFileReader::SequenceStreamParserPtr SequenceStreamParserPtr;
	typedef Read::ReadPtr ReadPtr;
	typedef Sequence::RecordPtr RecordPtr;
	typedef std::vector<Read> ReadVector;
	typedef std::vector<ReadSet> ReadSetVector;

	static const ReadSetSizeType MAX_READ_IDX = MAX_READ_SET_SIZE;

	static MmapSourceVector mmapSources;
	static void madviseMmaps(int advise);
	static void madviseMmapsRandom() {
		madviseMmaps(MADV_RANDOM);
	}
	static void madviseMmapsSequential() {
		madviseMmaps(MADV_SEQUENTIAL);
	}
	static void madviseMmapsNormal() {
		madviseMmaps(MADV_NORMAL);
	}
	static void madviseMmapsDontNeed() {
		madviseMmaps(MADV_DONTNEED);
	}

	class Pair {
	public:
		ReadSetSizeType read1;
		ReadSetSizeType read2;

		Pair() :
			read1(MAX_READ_IDX), read2(MAX_READ_IDX) {
		}
		Pair(ReadSetSizeType _read1) :
			read1(_read1), read2(MAX_READ_IDX) {
		}
		Pair(ReadSetSizeType _read1, ReadSetSizeType _read2) :
			read1(_read1), read2(_read2) {
		}
		Pair(const Pair &copy) :
			read1(copy.read1), read2(copy.read2) {
		}
		inline bool operator==(const Pair &other) const {
			return (read1 == other.read1) && (read2 == other.read2);
		}
		inline bool operator!=(const Pair & other) const {
			return !(*this == other);
		}
		inline bool operator<(const Pair &other) const {
			return lesser() < other.lesser();
		}
		inline ReadSetSizeType lesser() const {
			return (read1 < read2 ? read1 : read2);
		}
		inline bool isSingle() const {
			return (read1 != MAX_READ_IDX && read2 == MAX_READ_IDX) || (read1 == MAX_READ_IDX && read2 != MAX_READ_IDX);
		}
		inline bool isPaired() const {
			return (read1 != MAX_READ_IDX && read2 != MAX_READ_IDX);
		}
		inline bool hasAValidRead() const {
			return (read1 != MAX_READ_IDX || read2 != MAX_READ_IDX);
		}
	};

	typedef std::vector<Pair> PairedIndexType;

private:
	ReadVector _reads;

private:
	inline const ReadSet &constThis() const {
		return *this;
	}

protected:
	ReadSetSizeType _baseCount;
	SequenceLengthType _maxSequenceLength;
	ReadSetSizeType _globalSize;
	int _myGlobalRank;
	PartitioningData<ReadSetSizeType> _filePartitions;
	ReadIdxVector _globalOffsets;
	PairedIndexType _pairs;
	std::string previousReadName, previousReadComment; // for fast pairing
	uint8_t inputReadQualityBase; // the scaling of the new reads into this ReadSet (for auto-scaling to Read::FASTQ_START_CHAR)
	bool isReadQualityBaseValidated;

private:
	void addRead(const Read &read);
	void addRead(const Read &read, SequenceLengthType readLength, int rank = -1);
	bool _isSequentialPair(const Read &read);
	inline bool _setMaxSequenceLength(SequenceLengthType len) {
		if (len > _maxSequenceLength) {
#pragma omp critical
			if (len > _maxSequenceLength) {
			_maxSequenceLength = len;
			}
			return true;
		}
		return false;
	}
	void validateFastqStart(const Read &read) {
		if (getSize() < 20000 && read.hasQuals()) {
			if (! read.validateFastqStart() ) {
				if (Read::FASTQ_START_CHAR == Kmernator::FASTQ_START_CHAR_STD) {
					LOG_DEBUG(1, "Detected base 64 quality in: " << read.getName());
					if (getSize() > 10000) {
						Log::Warn() << "expected STD (33) fastq but detected ILLUMINA (64) only very far into the file, please make sure standard fastq and illumina fastq are not mixed" << endl;
					}
					__setFastqStart(Kmernator::FASTQ_START_CHAR_ILLUMINA);
				} else if (Read::FASTQ_START_CHAR == Kmernator::FASTQ_START_CHAR_ILLUMINA) {
					LOG_DEBUG(1, "Detected base 33 quality in: " <<  read.getName());
					if (getSize() > 10000) {
						Log::Warn() << "expected ILLLUMINA (64) fastq but detected STD (33) only very far into the file, please make sure standard fastq and illumina fastq are not mixed" << endl;
					}
					__setFastqStart(Kmernator::FASTQ_START_CHAR_STD);
				} else {
					throw;
				}
			}
		} else if (getSize() >= 20000)
			isReadQualityBaseValidated = true;
	}
	void __setFastqStart(int startChar) {
		if (startChar != inputReadQualityBase) {
#pragma omp critical (_setFastqStart)
		{
			if (startChar != inputReadQualityBase) {
				if (isReadQualityBaseValidated && getSize() > 0) {
					LOG_WARN(1, "Re-scaling fastq quality again!! (from " << (int) inputReadQualityBase << " to " << (int) startChar << ")");
				}
				LOG_VERBOSE(1, "Setting input fastq base quality score to: " << startChar);
				// re-scale any existing reads, if necessary
				rescaleQuality(inputReadQualityBase - startChar);
				inputReadQualityBase = startChar;
				isReadQualityBaseValidated = (getSize() > 0);
			}
		}
		}
	}

	void rescaleQuality(int delta) {
		LOG_DEBUG_OPTIONAL(1, !_reads.empty(), "Rescaling quality for " << _reads.size() << " reads by " << delta);
		for(ReadSetSizeType i = 0 ; i < _reads.size(); i++)
			_reads[i].rescaleQuality(delta);
	}

	void incrementFile(ReadFileReader &reader);
	void incrementFile(SequenceStreamParserPtr parser);

	MmapSource mmapFile(string filePath);

	static void addMmaps(MmapSourcePair mmaps);
	static bool &haveClearedCache() {
		static bool _ = false;
		return _;
	}
	static int32_t &getDefaultInputQualityBase() {
		static int32_t _def = GeneralOptions::getOptions().getFastqBaseQuality();
		return _def;
	}

public:
	ReadSet() :
		_baseCount(0), _maxSequenceLength(0), _globalSize(0), _myGlobalRank(0), inputReadQualityBase(getDefaultInputQualityBase()), isReadQualityBaseValidated(false) {
		if (! Read::isQualityToProbabilityInitialized())
			Read::setMinQualityScore();
		setFastqStart();
		isReadQualityBaseValidated = false;
		// this is needed to fix over subscription of threads where the Sequence::threadCacheSequences is undersized
		if (!omp_in_parallel() && !haveClearedCache()) {
			Sequence::clearCaches();
			haveClearedCache() = true;
		}
	}
	ReadSet(const ReadSet &copy)  :
		_baseCount(0), _maxSequenceLength(0), _globalSize(0), _myGlobalRank(0), inputReadQualityBase(getDefaultInputQualityBase()), isReadQualityBaseValidated(copy.isReadQualityBaseValidated) {
		*this = copy;
	}
	~ReadSet() {
	}

	static ReadSet shred(const Read &read, SequenceLengthType length, SequenceLengthType step) {
		ReadSet shreds;
		SequenceLengthType len = read.getLength();
		assert(len > length);
		assert(step < length);
		for(int i = 0; i < (int) (len - length + step - 1); i += step) {
			Read shreddedRead(read.getName() + "-" + boost::lexical_cast<string>(i) + "-" + boost::lexical_cast<string>(i+length),
					read.getFasta(i, length), read.getQuals(i, length), "shredded");
			shreds.addRead(shreddedRead);
		}
		return shreds;
	}

	void swap(ReadSet &other) {
		_reads.swap(other._reads);
		_filePartitions.swap(other._filePartitions);
		std::swap(_baseCount, other._baseCount);
		std::swap(_maxSequenceLength, other._maxSequenceLength);
		std::swap(_myGlobalRank, other._myGlobalRank);
		std::swap(_globalSize, other._globalSize);
		_globalOffsets.swap(other._globalOffsets);
		_pairs.swap(other._pairs);
		previousReadName.swap(other.previousReadName);
		std::swap(inputReadQualityBase,other.inputReadQualityBase);
		std::swap(isReadQualityBaseValidated,other.isReadQualityBaseValidated);
	}

	void clear() {
		_reads.clear();
		_filePartitions.clear();
		_baseCount = 0;
		_maxSequenceLength = 0;
		_myGlobalRank = 0;
		_globalSize = 0;
		_globalOffsets.clear();
		_pairs.clear();
		previousReadName.clear();
	}

	ReadSet &operator=(const ReadSet &copy) {
		_reads.assign(copy._reads.begin(), copy._reads.end());
		_filePartitions.clear();
		_baseCount = copy._baseCount;
		_maxSequenceLength = copy._maxSequenceLength;
		_globalOffsets.assign(copy._globalOffsets.begin(), copy._globalOffsets.end());
		_globalSize = copy._globalSize;
		_pairs.assign(copy._pairs.begin(), copy._pairs.end());
		previousReadName = copy.previousReadName;
		inputReadQualityBase = copy.inputReadQualityBase;
		isReadQualityBaseValidated = copy.isReadQualityBaseValidated;
		return *this;
	}
	void setFastqStart(char fastqStartChar = GeneralOptions::getOptions().getOutputFastqBaseQuality()) {
		__setFastqStart(fastqStartChar);
	}
	long getStoreSize() const {
		// just store numReads, baseCount, maxSeqLength, readSizeCounts & readData
		if (getSize() == 0)
			return sizeof(ReadSetSizeType);
		long size = sizeof(ReadSetSizeType)*2 + sizeof(SequenceLengthType) * (1+getSize());
		for(ReadSetSizeType i = 0; i < getSize(); i++)
			size += getRead(i).getStoreSize();
		return size;
	}
	long store(void *_dst) const {
		ReadSetSizeType *readSize = (ReadSetSizeType*) _dst;
		if (getSize() == 0) {
			*readSize = 0;
			return sizeof(ReadSetSizeType);
		}
		*(readSize++) = getSize();
		*(readSize++) = _baseCount;
		SequenceLengthType *dstSizes = (SequenceLengthType*) readSize;
		*(dstSizes++) = _maxSequenceLength;
		char *readData = (char*) (dstSizes + getSize());
		for(SequenceLengthType i = 0; i < getSize(); i++) {
			SequenceLengthType readSize = getRead(i).store(readData);
			*(dstSizes++) = readSize;
			readData += readSize;
		}
		return readData - ((char*) _dst);
	}
	void *restore(void *_src) {
		clear();

		ReadSetSizeType size, *rssp = (ReadSetSizeType*) _src;
		size = *(rssp++);
		if (size == 0) {
			return rssp;
		}
		_baseCount = *(rssp++);
		_reads.reserve(size);

		SequenceLengthType *dstSizes = (SequenceLengthType*) rssp;
		_maxSequenceLength = *(dstSizes++);
		char *readData = (char*) (dstSizes + size);
		for(ReadSetSizeType i = 0; i < size; i++) {
			SequenceLengthType readSize = *(dstSizes++);
			Read read;
			char *newReadData = (char*) read.restore(readData, readSize);
			_reads.push_back(read);
			assert(newReadData = readData + readSize);
			readData = newReadData;
		}
		return readData;
	}

	inline SequenceLengthType getMaxSequenceLength() const {
		if (_maxSequenceLength == 0 && _baseCount != 0)
			return getAvgSequenceLength() + 1;
		else
			return _maxSequenceLength;
	}
	inline SequenceLengthType getAvgSequenceLength() const {
		if (getSize() > 0)
			return _baseCount / getSize();
		else
			return 0;
	}
	void circularize(long extraLength);

	void appendAllFiles(OptionsBaseInterface::FileListType &files, int rank = 0, int size = 1);
	SequenceStreamParserPtr appendAnyFile(std::string filePath, std::string filePath2 = "", int rank = 0, int size = 1);
	SequenceStreamParserPtr appendAnyFileMmap(string fastaFilePath, string qualFilePath = "", int rank = 0, int size = 1);
	SequenceStreamParserPtr appendFastaFile(std::string &fastaFile, std::string &qualFile, int rank = 0, int size = 1);
	SequenceStreamParserPtr appendFastaFile(std::string &fastaFile, int rank = 0, int size = 1);
	SequenceStreamParserPtr appendFastaData(std::string &fastaData, int rank = 0, int size = 1);

	void append(const ReadSet &reads);
	void append(const Read &read);

	inline ReadSetSizeType getSize() const {
		return _reads.size();
	}
	inline ReadSetSizeType getBaseCount() const {
		return _baseCount;
	}

	inline ReadSetSizeType getPairSize() const {
		return _pairs.size();
	}

	// returns a copy of the reads from start
	ReadSet subset(ReadSetSizeType startIdx, ReadSetSizeType length) const {
		assert(startIdx + length <= getSize());
		ReadSet subset;
		for(ReadSetSizeType i = startIdx; i < length + startIdx && i < getSize(); i++) {
			subset.append(getRead(i));
		}
		return subset;
	}

	// keeps the first reads and returns the truncated ones
	ReadSet truncate(ReadSetSizeType length) {
		LOG_DEBUG_OPTIONAL(1, true, "ReadSet::truncate(" << length << ") sized:" << getSize());
		assert(length < getSize());
		ReadSet keep = subset(0, length);
		ReadSet rest = subset(length, getSize() - length);
		assert(keep.getSize() + rest.getSize() == getSize());
		*this = keep;
		return rest;
	}

	void setGlobalOffsets(int myRank, ReadIdxVector &globalSizes) {
		assert(myRank < (int) globalSizes.size());
		_myGlobalRank = myRank;
		_globalSize = 0;
		_globalOffsets.clear();
		_globalOffsets.reserve(globalSizes.size());
		for(int i = 0; i < (int) globalSizes.size(); i++) {
			_globalOffsets.push_back(_globalSize);
			_globalSize += globalSizes[i];
		}
	}

	inline bool isGlobal() const {
		return !_globalOffsets.empty();
	}
	inline ReadSetSizeType getGlobalOffset() const {
		return getGlobalOffset(_myGlobalRank);
	}
	inline ReadSetSizeType getGlobalOffset(int rank) const {
		if (isGlobal())
			return _globalOffsets[rank];
		else
			return 0;
	}

	inline ReadSetSizeType getGlobalReadIdx(ReadSetSizeType localReadIdx) const {
		return getGlobalReadIdx(_myGlobalRank, localReadIdx);
	}
	inline ReadSetSizeType getGlobalReadIdx(int rank, ReadSetSizeType localReadIdx) const {
		if (isGlobal())
			return _globalOffsets[rank] + localReadIdx;
		else
			return localReadIdx;
	}

	inline ReadSetSizeType getLocalReadIdx(ReadSetSizeType globalReadIdx) const {
		return getLocalReadIdx(_myGlobalRank, globalReadIdx);
	}
	inline ReadSetSizeType getLocalReadIdx(int rank, ReadSetSizeType globalReadIdx) const {
		if (!isGlobal()) {
			assert(rank == 0);
			return globalReadIdx;
		}
		assert(isLocalRead(rank, globalReadIdx));
		return globalReadIdx - _globalOffsets[rank];
	}
	//returns the rank and rankReadidx (localIdx) for a given globalReadIdx
	void getRankReadForGlobalReadIdx(ReadSetSizeType globalReadIdx, int &rank, ReadSetSizeType &rankReadIdx) const {
		int size = _globalOffsets.size();
		if (size == 0) {
			rank = 0;
			rankReadIdx = globalReadIdx;
		} else {
			for(rank = size - 1 ; rank > 0 ; rank--) {
				if (_globalOffsets[rank] <= globalReadIdx)
					break;
			}
			rankReadIdx = globalReadIdx - _globalOffsets[rank];
			LOG_DEBUG(5, "ReadSet::getRankReadForGlobalReadIdx(" << globalReadIdx << ", " << rank << ", " << rankReadIdx << "): " << _globalOffsets[rank]);
		}
	}

	bool isLocalRead(ReadSetSizeType globalReadIdx) const {
		return isLocalRead(_myGlobalRank, globalReadIdx);
	}
	bool isLocalRead(int rank, ReadSetSizeType globalReadIdx) const {
		if (!isGlobal())
			return true;
		if (globalReadIdx >= _globalSize) {
			LOG_THROW("isLocalRead(" << rank <<", " << globalReadIdx << ") exceeds globalSize: " << _globalSize);
		} else if (rank + 1 < (int) _globalOffsets.size()) {
			if (globalReadIdx >=_globalOffsets[rank+1]) {
				return false;
			}
		}

		if (_globalOffsets[rank] <= globalReadIdx) {
			return true;
		} else {
			return false;
		}
	}

	inline ReadSetSizeType getGlobalSize() const {
		return isGlobal() ? _globalSize : getSize();
	}

	inline bool isValidRead(ReadSetSizeType index) const {
		return index < getSize();
	}

	ReadPtr parseMmapedRead(ReadSetSizeType index) const;
	inline const Read &getRead(ReadSetSizeType index) const {
		assert(isValidRead(index));
		const Read &read = _reads[index];
		return read;
	}
	inline Read &getRead(ReadSetSizeType index) {
		return const_cast<Read&>(constThis().getRead(index));
	}

	inline int getReadFileNum(ReadSetSizeType index) const {
		return _filePartitions.getPartitionIdx(index)+1;
	}
	string _getReadFileNamePrefix(unsigned int filenum) const;
	string getReadFileNamePrefix(ReadSetSizeType index) const;
	// returns the first file to match either read from the pair
	string getReadFileNamePrefix(const Pair &pair) const;

	// by default no pairs are identified
	ReadSetSizeType identifyPairs();
	inline bool hasPairs() const {
		bool hasPairs = (getPairSize() != 0 && getPairSize() < getSize());
		return hasPairs;
	}

	static bool isPairedRead(const std::string readName, const std::string readComment = "");
	static bool isPair(const Read &readA, const Read &readB);
	static bool isPair(const std::string readNameA, const std::string readComment, const Read &readB);
	static bool isPair(const std::string readNameA, const std::string readNameB, const std::string commentA, const std::string commentB);

	// may return either as MAX_READ_IDX
	inline Pair &getPair(ReadSetSizeType pairIndex) {
		return _pairs[pairIndex];
	}
	inline const Pair &getPair(ReadSetSizeType pairIndex) const {
		return _pairs[pairIndex];
	}
	ReadSetSizeType getLocalPairIdx(ReadSetSizeType localIdx) const {
		assert(hasPairs());
		if (!getRead(localIdx).isPaired())
			return MAX_READ_IDX;
		// pairs are at the front of the readset, so jump to the pairIdx that it would be if it existed...
		ReadSetSizeType pairIdx = localIdx / 2;
		if (pairIdx < getPairSize()) {
			const Pair &pair = getPair(pairIdx);
			if (localIdx == pair.read1)
				return pair.read2;
			if (localIdx == pair.read2)
				return pair.read1;
		}
		LOG_WARN(1, "Could not find mate pair for " << localIdx << " pairidx: " << pairIdx << " read: " << getRead(localIdx).getName());
		return MAX_READ_IDX;
	}
	bool isSecondReadOnly() const {
		bool isSecondReadOnly = true;
		for (ReadSetSizeType i = 0 ; isSecondReadOnly & (i < _pairs.size()); i++)
			isSecondReadOnly &= (_pairs[i].read1 == MAX_READ_IDX) & (_pairs[i].read2 != MAX_READ_IDX);
		return isSecondReadOnly;
	}

	const ReadIdxVector getReadIdxVector() const;

	ProbabilityBases getProbabilityBases(unsigned char minQual = Options::getOptions().getMinQuality()) const;
	ReadSetSizeType getCentroidRead() const;
	ReadSetSizeType getCentroidRead(const ProbabilityBases &probs) const;
	Read getConsensusRead(unsigned char minQual = Options::getOptions().getMinQuality()) const;
	static Read getConsensusRead(const ProbabilityBases &probs, std::string name, std::string comment = "");

	inline std::ostream &write(std::ostream &os, ReadSetSizeType readIdx,
			SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH, std::string label = "", FormatOutput format = FormatOutput::getDefault()) const {
		return getRead(readIdx).write(os, trimOffset, trimLength, label, format);
	}
	inline std::ostream &write(OfstreamMap &om, ReadSetSizeType readIdx,
			SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH, std::string label = "", FormatOutput format = FormatOutput::getDefault()) const {
		return write(om.getOfstream( getReadFileNamePrefix(readIdx) ), readIdx, trimOffset, trimLength, label, format);
	}

	inline std::ostream &write(std::ostream &os, const Pair &pair,
			SequenceLengthType trimOffset1 = 0, SequenceLengthType trimLength1 = MAX_SEQUENCE_LENGTH, std::string label1 = "",
			SequenceLengthType trimOffset2 = 0, SequenceLengthType trimLength2 = MAX_SEQUENCE_LENGTH, std::string label2 = "",
			FormatOutput format = FormatOutput::getDefault(), bool forcePair = false) const {
		if (isValidRead(pair.read1))
			write(os, pair.read1, trimOffset1, trimLength1, label1, format);
		else if (forcePair)
			fakePair(getRead(pair.read2)).write(os);

		if (isValidRead(pair.read2))
			write(os, pair.read2, trimOffset2, trimLength2, label2, format);
		else if (forcePair)
			fakePair(getRead(pair.read1)).write(os);
		return os;
	}
	inline std::ostream &write(OfstreamMap &om, const Pair &pair,
			SequenceLengthType trimOffset1 = 0, SequenceLengthType trimLength1 = MAX_SEQUENCE_LENGTH, std::string label1 = "",
			SequenceLengthType trimOffset2 = 0, SequenceLengthType trimLength2 = MAX_SEQUENCE_LENGTH, std::string label2 = "",
			FormatOutput format = FormatOutput::getDefault(), bool forcePair = false) const {
		std::ostream &os = om.getOfstream(getReadFileNamePrefix(pair));

		return write(os, pair, trimOffset1, trimLength1, label1, trimOffset2, trimLength2, label2, format, forcePair);
	}
	inline std::ostream &writeAll(std::ostream &os, FormatOutput format = FormatOutput::getDefault(), bool trimmed = true) const {
		LOG_DEBUG(2, "ReadSet::writeAll()");
		for(ReadSetSizeType i = 0; i < getSize(); i++) {
			if (trimmed)
				write(os, i, 0, getRead(i).getFirstMarkupLength(), "", format);
			else
				write(os, i, 0, MAX_SEQUENCE_LENGTH, "", format);
		}
		return os;
	}
	static Read fakePair(const Read &unPaired);

	ReadSet randomlySample(ReadSetSizeType maxReads) const;

	SequenceStreamParserPtr appendFasta(ReadFileReader &reader, int rank = 0, int size = 1);

	std::string toString() const {
		std::stringstream ss;
		ss << "ReadSet{" <<getSize() << ":";
		for(ReadSetSizeType i = 0 ; i < getSize(); i++)
			ss << getRead(i).getName() << ", ";
		ss << "}";
		return ss.str();
	}

protected:
	SequenceStreamParserPtr appendFasta(std::string fastaFilePath, std::string qualFilePath = "", int rank = 0, int size = 1);
	SequenceStreamParserPtr appendFasta(MmapSource &mmap, int rank = 0, int size = 1);

	SequenceStreamParserPtr appendFastq(MmapSource &mmap);
	//void appendFastq(ReadFileReader &reader);
	SequenceStreamParserPtr appendFastqBlockedOMP(MmapSource &mmap);
	SequenceStreamParserPtr appendFastqBatchedOMP(std::string fastaFilePath,
			std::string qualFilePath = "");

};

class ReadIndexScore {
public:
	ReadSet::ReadSetSizeType readIndex;
	float score;
};

class KmerReadSetStats {
public:
	float kmerScore;
	std::vector<ReadIndexScore> linkedReads;
};


class ReadSetStream {
public:
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	ReadSetStream(const ReadSet &readSet) : _rs(&readSet), _rfr(NULL), _readIdx(0), _rank(0), _size(1) { 
		init();
	}
	ReadSetStream(ReadFileReader &reader) : _rs(NULL), _rfr(&reader), _readIdx(0), _rank(0), _size(1) {
		init();
	}
	ReadSetStream(std::string filename, int rank = 0, int size = 1) : _rs(NULL), _rfr(NULL), _readIdx(0), _rank(rank), _size(size) {
		init();
		_files.push_back(filename);
		setNextFile();
	}
	ReadSetStream(std::vector<std::string> &files, int rank = 0, int size = 1) : _files(files.begin(), files.end()), _rs(NULL), _rfr(NULL), _readIdx(0), _rank(rank), _size(size) {
		init();
		setNextFile();
	}
	void init() {
		if (! Read::isQualityToProbabilityInitialized() )
			Read::setMinQualityScore();
		_inputReadQualityBase = GeneralOptions::getOptions().getOutputFastqBaseQuality();
	}

	bool isReadSet() {
		return _rs == NULL ? false : true;
	}
	bool isReadFileReader() {
		return _rfr == NULL ? false : true;
	}
	bool hasNext() {
		if (isReadSet()) {
			_hasNext = _readIdx < _rs->getSize();
			if (_hasNext)
				_nextRead = _rs->getRead(_readIdx++);
		} else if (isReadFileReader()) {
			_hasNext = _rfr->nextRead(_name, _bases, _quals, _comment);
			if (_hasNext) {
				if (_inputReadQualityBase != Read::FASTQ_START_CHAR)
					Read::rescaleQuality(_quals, Read::FASTQ_START_CHAR - _inputReadQualityBase);
				_nextRead = Read(_name, _bases, _quals, _comment);
				if (!_nextRead.validateFastqStart()) {
					if (_inputReadQualityBase == Kmernator::FASTQ_START_CHAR_STD) {
						_inputReadQualityBase = Kmernator::FASTQ_START_CHAR_ILLUMINA;
						Read::rescaleQuality(_quals, Read::FASTQ_START_CHAR - Kmernator::FASTQ_START_CHAR_ILLUMINA);
					} else {
						_inputReadQualityBase = Kmernator::FASTQ_START_CHAR_STD;
						Read::rescaleQuality(_quals, Read::FASTQ_START_CHAR - Kmernator::FASTQ_START_CHAR_STD);
					}
					if (_readIdx > 0)
						LOG_WARN(1, "Detected different FASTQ quality scaling than expected.  Some data was already processed.  You should re-run with --fastq-base-quality " << _inputReadQualityBase);
					_nextRead = Read(_name, _bases, _quals, _comment);
				}
				_readIdx++;
			} else if (setNextFile())
				return hasNext();
		} else 
			_hasNext = false;
		return _hasNext;
	}
	const Read &getRead() {
		LOG_DEBUG(3, "ReadSetStream::getRead(): " << _nextRead.getName() << " " << getReadCount());
		return _nextRead;
	}
	ReadSetSizeType getReadCount() {
		return _readIdx;
	}
protected:
	bool setNextFile() {
		if (_files.empty()) {
			_rfr = NULL;
			return false;
		}
		LOG_DEBUG(3, "ReadSetStream::setNextFile(): " << _files.front());
		_rfrptr.reset( new ReadFileReader( _files.front(), true ) );
		_rfr = _rfrptr.get();
		_files.pop_front();
		if (_size != 1)
			_rfr->seekToPartition(_rank, _size);
		return true;
	}
private:
	std::deque<std::string> _files;
	boost::shared_ptr<ReadFileReader> _rfrptr;
	const ReadSet *_rs;
	ReadFileReader *_rfr;
	Read _nextRead;
	ReadSetSizeType _readIdx;
	string _name, _bases, _quals, _comment;
	bool _hasNext;
	int _rank, _size;
	uint8_t _inputReadQualityBase;
};

#endif
