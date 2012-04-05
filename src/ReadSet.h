//
// Kmernator/src/ReadSet.h
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

#ifndef _READ_SET_H
#define _READ_SET_H
#include <cstring>
#include <boost/unordered_map.hpp>
#include <sys/mman.h>

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
	std::string previousReadName; // for fast pairing

private:
	void addRead(const Read &read);
	void addRead(const Read &read, SequenceLengthType readLength, int rank = -1);
	void _trackSequentialPair(const Read &read);
	inline bool _setMaxSequenceLength(SequenceLengthType len) {
		if (len > _maxSequenceLength) {
#pragma omp critical
			_maxSequenceLength = len;
			return true;
		}
		return false;
	}
	inline void _setFastqStart(const Read &read) {
		if (read.hasQuals() && Read::FASTQ_START_CHAR != Kmernator::FASTQ_START_CHAR_STD && omp_get_thread_num() == 0) {
			std::string quals = read.getQuals();
			std::string::iterator it = std::min_element(quals.begin(), quals.end());
			if (it != quals.end() && *it < Read::FASTQ_START_CHAR) {
				if (getSize() > 10000) {
					Log::Warn() << "detected standard fastq only very far into the file, please make sure standard fastq and illumina fastq are not mixed" << endl;
				}
				Read::setMinQualityScore(Options::getOptions().getMinQuality(), Kmernator::FASTQ_START_CHAR_STD);
			}
		}
	}

	void incrementFile(ReadFileReader &reader);
	void incrementFile(SequenceStreamParserPtr parser);

	MmapSource mmapFile(string filePath);

	static void addMmaps(MmapSourcePair mmaps);

public:
	ReadSet() :
		_baseCount(0), _maxSequenceLength(0), _globalSize(0), _myGlobalRank(0) {
	}
	ReadSet(const ReadSet &copy)  :
		_baseCount(0), _maxSequenceLength(0), _globalSize(0), _myGlobalRank(0) {
		*this = copy;
	}
	~ReadSet() {
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
		return *this;
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
		assert(length < getSize());
		ReadSet keep = subset(0, length);
		ReadSet rest = subset(length, getSize() - length);
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
		return getPairSize() != 0 && getPairSize() < getSize();
	}

	static bool isPairedRead(const std::string &readName);
	static bool isPair(const Read &readA, const Read &readB);
	static bool isPair(const std::string &readNameA, const Read &readB);
	static bool isPair(const std::string &readNameA, const std::string &readNameB);

	// may return either as MAX_READ_IDX
	inline Pair &getPair(ReadSetSizeType pairIndex) {
		return _pairs[pairIndex];
	}
	inline const Pair &getPair(ReadSetSizeType pairIndex) const {
		return _pairs[pairIndex];
	}
	ReadSetSizeType getGlobalPairIdx(ReadSetSizeType globalIdx) const {
		// FIXME hack!!!
		assert(hasPairs());
		return (globalIdx & 0x1) == 0 ? globalIdx + 1 : globalIdx - 1;
	}

	const ReadIdxVector getReadIdxVector() const;

	ProbabilityBases getProbabilityBases(unsigned char minQual = Options::getOptions().getMinQuality()) const;
	ReadSetSizeType getCentroidRead() const;
	ReadSetSizeType getCentroidRead(const ProbabilityBases &probs) const;
	Read getConsensusRead(unsigned char minQual = Options::getOptions().getMinQuality()) const;
	static Read getConsensusRead(const ProbabilityBases &probs, std::string name);

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

#endif
