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
	PartitioningData<ReadSetSizeType> _filePartitions;
	unsigned long _baseCount;
	SequenceLengthType _maxSequenceLength;
	ReadSetSizeType _globalOffset;
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
				Read::setMinQualityScore(Options::getMinQuality(), Kmernator::FASTQ_START_CHAR_STD);
			}
		}
	}

	void incrementFile(ReadFileReader &reader);
	void incrementFile(SequenceStreamParserPtr parser);

	MmapSource mmapFile(string filePath);

	static void addMmaps(MmapSourcePair mmaps);

public:
	ReadSet() :
		_baseCount(0), _maxSequenceLength(0), _globalOffset(0) {
	}
	~ReadSet() {
	}

	inline SequenceLengthType getMaxSequenceLength() const {
		return _maxSequenceLength;
	}
	void circularize(long extraLength);

	void appendAllFiles(Options::FileListType &files, int rank = 0, int size = 1);
	SequenceStreamParserPtr appendAnyFile(std::string filePath, std::string filePath2 = "", int rank = 0, int size = 1);
	SequenceStreamParserPtr appendAnyFileMmap(string fastaFilePath, string qualFilePath = "", int rank = 0, int size = 1);
	SequenceStreamParserPtr appendFastaFile(std::string &is, int rank = 0, int size = 1);

	void append(const ReadSet &reads);
	void append(const Read &read);

	inline ReadSetSizeType getSize() const {
		return _reads.size();
	}
	inline unsigned long getBaseCount() const {
		return _baseCount;
	}

	inline ReadSetSizeType getPairSize() const {
		return _pairs.size();
	}

	inline void setGlobalOffset(ReadSetSizeType globalOffset) {
		_globalOffset = globalOffset;
	}

	inline ReadSetSizeType getGlobalOffset() const {
		return _globalOffset;
	}

	inline bool isValidRead(ReadSetSizeType index) const {
		return index < getSize();
	}

	ReadPtr parseMmapedRead(ReadSetSizeType index) const;
	inline const Read &getRead(ReadSetSizeType index) const {
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
	inline bool hasPairs() {
		return getPairSize() != 0 && getPairSize() < getSize();
	}
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

	const ReadIdxVector getReadIdxVector() const;

	ProbabilityBases getProbabilityBases(unsigned char minQual = Options::getMinQuality()) const;
	ReadSetSizeType getCentroidRead() const;
	ReadSetSizeType getCentroidRead(const ProbabilityBases &probs) const;
	Read getConsensusRead(unsigned char minQual = Options::getMinQuality()) const;
	static Read getConsensusRead(const ProbabilityBases &probs, std::string name);

	inline std::ostream &write(std::ostream &os, ReadSetSizeType readIdx,
			SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH, std::string label = "", FormatOutput format = FormatOutput::getDefault()) const {
		return getRead(readIdx).write(os, trimOffset, trimLength, label, format);
	}
	inline std::ostream &write(OfstreamMap &om, ReadSetSizeType readIdx,
			SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH, std::string label = "", FormatOutput format = FormatOutput::getDefault()) const {
		return write(om.getOfstream( getReadFileNamePrefix(readIdx) ), readIdx, trimOffset, trimLength, label, format);
	}
	inline std::ostream &write(OfstreamMap &om, const Pair &pair,
			SequenceLengthType trimOffset1 = 0, SequenceLengthType trimLength1 = MAX_SEQUENCE_LENGTH, std::string label1 = "",
			SequenceLengthType trimOffset2 = 0, SequenceLengthType trimLength2 = MAX_SEQUENCE_LENGTH, std::string label2 = "",
			FormatOutput format = FormatOutput::getDefault(), bool forcePair = false) const {
		std::ostream &os = om.getOfstream(getReadFileNamePrefix(pair));
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
	static Read fakePair(const Read &unPaired);


protected:
	SequenceStreamParserPtr appendFasta(std::string fastaFilePath, std::string qualFilePath = "", int rank = 0, int size = 1);
	SequenceStreamParserPtr appendFasta(MmapSource &mmap, int rank = 0, int size = 1);
	SequenceStreamParserPtr appendFasta(ReadFileReader &reader, int rank = 0, int size = 1);

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

//
// $Log: ReadSet.h,v $
// Revision 1.31  2010-08-18 17:50:39  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.30.8.1  2010-07-20 20:02:56  regan
// autodetect fastq quality range
//
// Revision 1.30  2010-05-24 21:48:46  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.29.2.3  2010-05-19 21:53:20  regan
// bugfixes
//
// Revision 1.29.2.2  2010-05-19 21:36:54  regan
// refactored duplicate fragment filter code
// added duplicate fragment on single ended reads
//
// Revision 1.29.2.1  2010-05-19 00:20:46  regan
// refactored fomat output options
// added options to fastq2fasta
//
// Revision 1.29  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.28.2.4  2010-05-12 22:45:00  regan
// added readset circularize method
//
// Revision 1.28.2.3  2010-05-12 18:25:48  regan
// minor refactor
//
// Revision 1.28.2.2  2010-05-10 21:24:29  regan
// minor refactor moved code into cpp
//
// Revision 1.28.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.28  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.27  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
// Revision 1.26.2.2  2010-05-06 18:47:50  regan
// fixed
//
// Revision 1.26.2.1  2010-05-06 18:45:36  regan
// broke it...
//
// Revision 1.26  2010-05-06 16:43:56  regan
// merged changes from ConsensusTesting-20100505
//
// Revision 1.25.2.1  2010-05-05 23:46:23  regan
// checkpoint... seems to compile
//
// Revision 1.25  2010-05-05 06:28:35  regan
// merged changes from FixPairOutput-20100504
//
// Revision 1.24.4.2  2010-05-05 05:57:53  regan
// fixed pairing
// fixed name to exclude labels and comments after whitespace
// applied some performance optimizations from other branch
// created FixPair application
//
// Revision 1.24.4.1  2010-05-04 21:33:42  regan
// checkpoint
//
// Revision 1.24.2.2  2010-05-04 19:49:51  regan
// minor rework on include headers
//
// Revision 1.24.2.1  2010-05-02 05:40:09  regan
// added methods and cache variable for fast sequential read pair identification
//
// Revision 1.24  2010-05-01 21:57:53  regan
// merged head with serial threaded build partitioning
//
// Revision 1.23.2.11  2010-05-01 21:29:00  regan
// fixed naming of output files
//
// Revision 1.23.2.10  2010-05-01 05:56:38  regan
// bugfix
//
// Revision 1.23.2.9  2010-04-30 23:53:14  regan
// attempt to fix a bug.  clearing Sequence caches when it makes sense
//
// Revision 1.23.2.8  2010-04-30 21:53:52  regan
// reuse memory efficiently for cache lookups
//
// Revision 1.23.2.7  2010-04-29 04:26:32  regan
// bugfix in output filenames and content
//
// Revision 1.23.2.6  2010-04-28 22:28:11  regan
// refactored writing routines
//
// Revision 1.23.2.5  2010-04-28 16:57:00  regan
// bugfix in output filenames
//
// Revision 1.23.2.4  2010-04-27 23:17:50  regan
// fixed naming of output files
//
// Revision 1.23.2.3  2010-04-27 22:53:35  regan
// added madvise calls
//
// Revision 1.23.2.2  2010-04-26 23:34:41  regan
// bugfix
//
// Revision 1.23.2.1  2010-04-26 22:53:55  regan
// fixed warnings
//
// Revision 1.23  2010-04-21 00:33:20  regan
// merged with branch to detect duplicated fragment pairs with edit distance
//
// Revision 1.22.2.1  2010-04-20 23:56:43  regan
// added madvise handles for mmap optimization
//
// Revision 1.22  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.21.2.10  2010-04-15 17:29:02  regan
// checkpoint, working with some optimizations
//
// Revision 1.21.2.9  2010-04-14 20:53:49  regan
// checkpoint and passes unit tests!
//
// Revision 1.21.2.8  2010-04-14 17:51:43  regan
// checkpoint
//
// Revision 1.21.2.7  2010-04-14 05:35:37  regan
// checkpoint. compiles but segfaults
//
// Revision 1.21.2.6  2010-04-14 03:51:20  regan
// checkpoint. compiles but segfaults
//
// Revision 1.21.2.5  2010-04-12 22:37:47  regan
// checkpoint
//
// Revision 1.21.2.4  2010-04-12 20:59:45  regan
// mmap checkpoint
//
// Revision 1.21.2.3  2010-04-07 22:33:08  regan
// checkpoint mmaping input files
//
// Revision 1.21.2.2  2010-04-05 05:42:53  regan
// checkpoint mmaping input files
//
// Revision 1.21.2.1  2010-04-05 03:32:10  regan
// moved read file reader
//
// Revision 1.21  2010-03-15 18:07:01  regan
// minor refactor and added consensus read
//
// Revision 1.20  2010-03-14 16:56:38  regan
// added centroid methods
// minor refactor
//
// Revision 1.19  2010-03-03 17:10:26  regan
// added ability to recognize which reads came from which files
//
// Revision 1.18  2010-02-26 13:01:16  regan
// reformatted
//
// Revision 1.17  2010-02-22 14:40:45  regan
// major milestone
//
// Revision 1.16  2010-01-16 01:07:29  regan
// added method
//
// Revision 1.15  2010-01-14 00:50:07  regan
// fixes
//
// Revision 1.14  2010-01-13 23:47:44  regan
// made const class modifications
// fixed identify pairs
//
// Revision 1.13  2010-01-13 00:26:49  regan
// fixed some parallelism
// started pair indentification
//
// Revision 1.12  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.11  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.10  2009-12-23 07:16:52  regan
// fixed reading of fasta files
// parallelized reading of multiple files
//
// Revision 1.9  2009-12-22 18:31:41  regan
// parallelized reading fastq if openmp is enabled
//
// Revision 1.8  2009-11-07 00:28:41  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
// Revision 1.7  2009-11-04 19:32:03  cfurman
// now reads in fasta (with optional qual) files
//
// Revision 1.6  2009-10-31 00:16:35  regan
// minor changes and optimizations
//
// Revision 1.5  2009-10-26 23:02:49  regan
// checkpoint
//
// Revision 1.4  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
// Revision 1.3  2009-10-21 00:00:58  cfurman
// working on kmers....
//
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//
