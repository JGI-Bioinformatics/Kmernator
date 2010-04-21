// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.h,v 1.23 2010-04-21 00:33:20 regan Exp $
//

#ifndef _READ_SET_H
#define _READ_SET_H
#include <string>
#include <boost/unordered_map.hpp>
#include <sys/mman.h>

#include "config.h"
#include "Options.h"
#include "Sequence.h"
#include "Utils.h"
#include "ReadFileReader.h"

class ReadSet {
public:
	typedef unsigned int ReadSetSizeType;
	typedef std::vector< ReadSetSizeType > ReadIdxVector;
	typedef ReadFileReader::MmapSource MmapSource;
	typedef std::pair<MmapSource,MmapSource> MmapSourcePair;
	typedef std::vector< MmapSourcePair > MmapSourceVector;
	typedef ReadFileReader::SequenceStreamParser SequenceStreamParser;
	typedef ReadFileReader::SequenceStreamParserPtr SequenceStreamParserPtr;
	typedef Read::ReadPtr ReadPtr;
	typedef Sequence::RecordPtr RecordPtr;

	static MmapSourceVector mmapSources;
    static void madviseMmaps(int advise) {
		for(MmapSourceVector::iterator it = mmapSources.begin(); it != mmapSources.end(); it++) {
			if (it->first.is_open())
				madvise(const_cast<char*>(it->first.data()), it->first.size(), advise);
			if (it->second.is_open())
				madvise(const_cast<char*>(it->second.data()), it->second.size(), advise);
		}
    }
	static void madviseMmapsRandom() {
		madviseMmaps(MADV_RANDOM);
	}
	static void madviseMmapsSequential() {
		madviseMmaps(MADV_SEQUENTIAL);
	}

	static const ReadSetSizeType MAX_READ_IDX = (unsigned int) -1;
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
		bool operator==(const Pair &other) const {
			return (read1 == other.read1) && (read2 == other.read2);
		}
		bool operator<(const Pair &other) const {
			return lesser() < other.lesser();
		}
		inline ReadSetSizeType lesser() const {
			return (read1 < read2 ? read1 : read2);
		}
	};

	typedef std::vector<Pair> PairedIndexType;

private:
	std::vector<Read> _reads;

private:
	inline const ReadSet &constThis() const {
		return *this;
	}

protected:
	PartitioningData<ReadSetSizeType> _filePartitions;
	unsigned long _baseCount;
	SequenceLengthType _maxSequenceLength;
	PairedIndexType _pairs;

private:
	void addRead(Read &read);
	void addRead(Read &read, SequenceLengthType readLength);
	inline bool setMaxSequenceLength(SequenceLengthType len) {
		if (len > _maxSequenceLength) {
			_maxSequenceLength = len;
			return true;
		}
		return false;
	}

	void incrementFile(ReadFileReader &reader) {
		incrementFile(reader.getParser());
	}
	void incrementFile(SequenceStreamParserPtr parser) {
		_filePartitions.addPartition( _reads.size() );
		addMmaps( MmapSourcePair( parser->getMmap(), parser->getQualMmap() ));
	}

	MmapSource mmapFile(string filePath);

	static void addMmaps(MmapSourcePair mmaps) {
#ifdef _USE_OPENMP
#pragma omp critical
#endif
		mmapSources.push_back(mmaps);
		madviseMmapsSequential();
	}

public:
	ReadSet() :
		_baseCount(0) {
	}
	~ReadSet() {
	}

	inline SequenceLengthType getMaxSequenceLength() const {
		return _maxSequenceLength;
	}

	void appendAllFiles(Options::FileListType &files);
	SequenceStreamParserPtr appendAnyFile(std::string filePath, std::string filePath2 = "");
	SequenceStreamParserPtr appendFastaFile(std::string &is);

	void append(const ReadSet &reads);
	void append(const Read &read) {
		_reads.push_back(read);
	}

	inline ReadSetSizeType getSize() const {
		return _reads.size();
	}
	unsigned long getBaseCount() const {
		return _baseCount;
	}

	inline ReadSetSizeType getPairSize() const {
		return _pairs.size();
	}

	inline bool isValidRead(ReadSetSizeType index) const {
		return index < getSize();
	}

	ReadPtr parseMmapedRead(ReadSetSizeType index) const {
		const Read &read = _reads[index];
		ReadPtr readPtr = read.readMmaped();
//		const SequenceStreamParserPtr parser = getParser(index);
//		string name, bases, qual;
//		parser->readRecord(read.getRecord(), name, bases, qual);
		return readPtr;
	}
	inline const Read &getRead(ReadSetSizeType index) const {
		const Read &read = _reads[index];
//		if (read.isMmaped()) {
//			// assume each thread modifies only one read at a time!!!
//			Read::setThreadCache((Sequence&) read, parseMmapedRead(index));
//		}
        return read;
	}
	inline Read &getRead(ReadSetSizeType index) {
		return const_cast<Read&>(constThis().getRead(index));
	}

	inline int getReadFileNum(ReadSetSizeType index) const {
		return _filePartitions.getPartitionIdx(index)+1;
	}
//	const SequenceStreamParserPtr getParser(ReadSetSizeType index) const {
//		return _fileParsers[getReadFileNum(index)-1];
//	}

	// by default no pairs are identified
	ReadSetSizeType identifyPairs();
	bool hasPairs() {
		return getPairSize() != 0 && getPairSize() < getSize();
	}

	// may return either as MAX_READ_IDX
	inline Pair &getPair(ReadSetSizeType pairIndex) {
		return _pairs[pairIndex];
	}
	inline const Pair &getPair(ReadSetSizeType pairIndex) const {
		return _pairs[pairIndex];
	}

	const ReadIdxVector getReadIdxVector() const {
		ReadIdxVector readIdxs;
		readIdxs.reserve(_reads.size());
		for(ReadSetSizeType i = 0 ; i < _reads.size(); i++) {
			readIdxs.push_back(i);
		}
		return readIdxs;
	}

	ProbabilityBases getProbabilityBases() const;
	ReadSetSizeType getCentroidRead() const;
	ReadSetSizeType getCentroidRead(const ProbabilityBases &probs) const;
	Read getConsensusRead() const;
	static Read getConsensusRead(const ProbabilityBases &probs, std::string name);

protected:
	SequenceStreamParserPtr appendFasta(std::string fastaFilePath, std::string qualFilePath = "");
	SequenceStreamParserPtr appendFasta(MmapSource &mmap);
	SequenceStreamParserPtr appendFasta(ReadFileReader &reader);

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
