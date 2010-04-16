// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.cpp,v 1.32 2010-04-16 22:44:18 regan Exp $
//

#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>

#include <locale>
#include <algorithm>

#include <boost/lexical_cast.hpp>

#include "ReadSet.h"

ReadSet::MmapSourceVector ReadSet::mmapSources;

void ReadSet::addRead(Read &read) {
	SequenceLengthType readLength = read.getLength();
	addRead(read, readLength);
}
void ReadSet::addRead(Read &read, SequenceLengthType readLength) {
	_reads.push_back(read);
	_baseCount += readLength;
	setMaxSequenceLength(readLength);
}
ReadSet::MmapSource ReadSet::mmapFile(string filePath) {
	MmapSource mmap(filePath, ReadFileReader::getFileSize(filePath));
	addMmaps( MmapSourcePair(mmap, MmapSource()) );
	return mmap;
}

ReadSet::SequenceStreamParserPtr ReadSet::appendAnyFile(string fastaFilePath, string qualFilePath) {
    MmapSource mmap = mmapFile(fastaFilePath);

    ReadFileReader reader;
    if (!qualFilePath.empty()) {
    	MmapSource mmap2 = mmapFile(qualFilePath);
    	reader.setReader(mmap, mmap2);
    } else {
	    reader.setReader(mmap);
    }

	switch (reader.getType()) {
	case 0:
		appendFastq(mmap);
		break;
	case 1:
		appendFasta(reader);
	}
	incrementFile(reader);
	return reader.getParser();
}

void ReadSet::appendAllFiles(Options::FileListType &files) {

#ifdef _USE_OPENMP
	ReadSet myReads[ files.size() ];
	SequenceStreamParserPtr parsers[ files.size() ];
	int numThreads = omp_get_max_threads() / MAX_FILE_PARALLELISM;
	if ( numThreads > 1 ) {
		omp_set_nested(1);
	} else {
		numThreads = omp_get_max_threads();
		omp_set_nested(0);
	}

#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
#endif
	for (long i = 0; i < (long) files.size(); i++) {
#ifdef _USE_OPENMP
#pragma omp critical
#endif
		{
			std::cerr << "reading " << files[i] << std::endl;
		}

#ifdef _USE_OPENMP
		// append int this thread's ReadSet buffer (note: line continues)
		parsers[i] = myReads[ i ].appendAnyFile(files[i]);
#else
		SequenceStreamParserPtr parser = appendAnyFile(files[i]);
		incrementFile(parser);
#endif


#ifdef _USE_OPENMP
#pragma omp critical
#endif
		{
			std::cerr << "finished reading " << files[i] << std::endl;
		}
	}
#ifdef _USE_OPENMP
	omp_set_nested(OMP_NESTED_DEFAULT);
	{	std::cerr << "concatenating ReadSet buffers" << std::endl;}
	for(int i = 0; i< (long) files.size(); i++) {
	    append(myReads[i]);
	    incrementFile(parsers[i]);
	}
#endif

}

void ReadSet::append(const ReadSet &reads) {
	unsigned long oldSize = _reads.size();
	unsigned long newSize = oldSize + reads._reads.size();
	_reads.resize(newSize);
#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (long i = 0; i < (long) reads._reads.size(); i++)
		_reads[oldSize + i] = reads._reads[i];

	unsigned long oldPairSize = _pairs.size();
	unsigned long newPairSize = oldPairSize + reads._pairs.size();
	_pairs.resize(newPairSize);

#ifdef _USE_OPENMP
#pragma omp parallel for
#endif
	for (long i = 0; i < (long) reads._pairs.size(); i++) {
		const Pair &tmp = reads._pairs[i];
		_pairs[oldPairSize + i]
				= Pair((tmp.read1 == MAX_READ_IDX ? MAX_READ_IDX : tmp.read1
						+ oldSize), (tmp.read2 == MAX_READ_IDX ? MAX_READ_IDX
						: tmp.read2 + oldSize));
	}
	_baseCount += reads._baseCount;
	setMaxSequenceLength(reads.getMaxSequenceLength());

}

ReadSet::SequenceStreamParserPtr ReadSet::appendFasta(string fastaFilePath, string qualFilePath) {
	ReadFileReader reader(fastaFilePath, qualFilePath);
	appendFasta(reader);
	incrementFile(reader);
	return reader.getParser();
}
ReadSet::SequenceStreamParserPtr ReadSet::appendFasta(ReadSet::MmapSource &mmap) {
	ReadFileReader reader(mmap);
	appendFasta(reader);
	incrementFile(reader);
	return reader.getParser();
}

ReadSet::SequenceStreamParserPtr ReadSet::appendFasta(ReadFileReader &reader) {
	string name, bases, quals;
	if (reader.isMmaped() && Options::getMmapInput() != 0) {
	    RecordPtr recordPtr = reader.getStreamRecordPtr();
	    RecordPtr qualPtr = reader.getStreamQualRecordPtr();
	    RecordPtr nextRecordPtr = recordPtr;
	    std::string name, bases, quals;
	    bool isMultiline;
	    while (reader.nextRead(nextRecordPtr, name, bases, quals, isMultiline)) {
            if (isMultiline) {
            	// store the read in memory
            	Read read(name,bases,quals);
            	addRead(read, bases.length());
            } else {
            	Read read(recordPtr, qualPtr);
            	addRead(read, bases.length());
            }
            recordPtr = nextRecordPtr;
	    	qualPtr = reader.getStreamQualRecordPtr();
	    }
	} else {
	    while (reader.nextRead(name, bases, quals)) {
	        Read read(name, bases, quals);
	        addRead(read, bases.length());
	    }
	}
	return reader.getParser();
}

ReadSet::SequenceStreamParserPtr ReadSet::appendFastaFile(string &str) {
	ReadFileReader reader(str);
	appendFasta(reader);
	incrementFile(reader);
	return reader.getParser();
}

// TODO FIXME does not like mmap..
#ifdef _USE_OPENMP

ReadSet::SequenceStreamParserPtr ReadSet::appendFastqBlockedOMP(ReadSet::MmapSource &mmap)
{

	unsigned long startIdx = _reads.size();
	unsigned long blockSize = 0;

	// set OMP variables
	omp_set_nested(1);
	int numThreads = omp_get_max_threads();
	if (numThreads > MAX_FILE_PARALLELISM)
	numThreads = MAX_FILE_PARALLELISM;
	unsigned long numReads[ numThreads ];
	unsigned long seekPos[ numThreads ];
	for(int i = 0; i < numThreads; i++)
	numReads[i] = 0;
	SequenceStreamParserPtr singleParser;

#pragma omp parallel num_threads(numThreads)
	{

		if (omp_get_num_threads() != numThreads)
		throw "OMP thread count discrepancy!";

		ReadSet myReads;
		ReadFileReader reader(mmap);

		if (omp_get_thread_num() == 0)
			singleParser = reader.getParser();

#pragma omp single
		{
			blockSize = reader.getBlockSize(numThreads);
			if (blockSize < 100)
			blockSize = 100;
			std::cerr << "Reading " << mmap << " with " << numThreads << " threads" << std::endl;
		}
		reader.seekToNextRecord( blockSize * omp_get_thread_num() );
		seekPos[ omp_get_thread_num() ] = reader.getPos();

		string name,bases,quals;
		bool hasNext = omp_get_thread_num() < (long) numThreads;
		if (hasNext)
		hasNext = (reader.getPos() < blockSize * (omp_get_thread_num() +1));

		if (!hasNext)
		seekPos[ omp_get_thread_num() ] = 0;

#pragma omp barrier
		if (hasNext && omp_get_thread_num() != 0 && seekPos[omp_get_thread_num()] == seekPos[omp_get_thread_num()-1])
		hasNext = false;

		//#pragma omp critical
		//{ std::cerr << omp_get_thread_num() << " " << blockSize << " seeked to " << reader.getPos() << " will read: " << hasNext << std::endl; }

		if (reader.isMmaped() && Options::getMmapInput() != 0) {
		    RecordPtr recordPtr = reader.getStreamRecordPtr();
		    RecordPtr qualPtr = reader.getStreamQualRecordPtr();
		    RecordPtr nextRecordPtr = recordPtr;
		    std::string name, bases, quals;
		    bool isMultiline;
		    while (hasNext && reader.nextRead(nextRecordPtr, name, bases, quals, isMultiline)) {
	            if (isMultiline) {
	            	// store the read in memory, as mmap does not currently work
	            	Read read(name,bases,quals);
	            	myReads.addRead(read, bases.length());
	            } else {
	            	Read read(recordPtr, qualPtr);
                    myReads.addRead(read, bases.length());
	            }

		    	recordPtr = nextRecordPtr;
		    	qualPtr = reader.getStreamQualRecordPtr();
		    	hasNext = (reader.getPos() < blockSize * (omp_get_thread_num() +1));
		    }
		} else {
		    while (hasNext && reader.nextRead(name, bases, quals)) {
		        Read read(name, bases, quals);
		        myReads.addRead(read, bases.length());
		        hasNext = (reader.getPos() < blockSize * (omp_get_thread_num() +1));
		    }
		}

		numReads[ omp_get_thread_num() ] = myReads.getSize();

		//#pragma omp critical
		//{ std::cerr << omp_get_thread_num() << " " << blockSize << " finished at " << reader.getPos() << " with " << myReads.getSize() << " " << numReads[ omp_get_thread_num() ]<< std::endl; }

#pragma omp critical
		{
			// set global counters
			_baseCount += myReads._baseCount;
			setMaxSequenceLength(myReads.getMaxSequenceLength());
		}

#pragma omp barrier

#pragma omp single
		{
			unsigned long newReads = 0;
			for(int i=0; i < numThreads; i++) {
				unsigned long tmp = newReads;
				newReads += numReads[i];
				numReads[i] = tmp;
			}

			_reads.resize(startIdx + newReads);
		}

		for(unsigned long j = 0; j < myReads.getSize(); j++)
		_reads[startIdx + numReads[ omp_get_thread_num() ] + j] = myReads.getRead(j);
	}

    incrementFile(singleParser);
	// reset omp variables
	omp_set_nested(OMP_NESTED_DEFAULT);
	return singleParser;
}

ReadSet::SequenceStreamParserPtr ReadSet::appendFastq(ReadSet::MmapSource &mmap)
{
	return appendFastqBlockedOMP(mmap);
}

#else // _USE_OPENMP
ReadSet::SequenceStreamParserPtr ReadSet::appendFastq(ReadSet::MmapSource &mmap) {
	return appendFasta(mmap);
}
#endif // _USE_OPENMP

std::string _commonName(const std::string &readName) {
	return readName.substr(0, readName.length() - 1);
}
int _readNum(const std::string &readName) {
	int retVal = 0;
	int len = readName.length();
	char c = readName[len - 1];
	switch (c) {
	case '1':
		if (readName[len - 2] != '/')
			break;
	case 'A':
	case 'F':
		retVal = 1;
		break;
	case '2':
		if (readName[len - 2] != '/')
			break;
	case 'B':
	case 'R':
		retVal = 2;
		break;
	}
	return retVal;
}

ReadSet::ReadSetSizeType ReadSet::identifyPairs() {
	ReadSetSizeType size = getSize();

	// find the last Pair to be identified
	ReadSetSizeType readIdx = 0;
	if (_pairs.size() > 0) {
		Pair &lastPair = _pairs[_pairs.size() - 1];
		if (isValidRead(lastPair.read1))
			readIdx = lastPair.read1;
		if (isValidRead(lastPair.read2) && readIdx < lastPair.read2)
			readIdx = lastPair.read2;
		readIdx++; // start with the next read
	}
	if (size <= readIdx)
		return _pairs.size();

	boost::unordered_map<std::string, ReadSetSizeType> unmatchedNames;
	boost::unordered_map<std::string, ReadSetSizeType>::iterator unmatchedIt;

	// build the unmatchedNames map for existing identified pairs
	for (ReadSetSizeType i = 0; i < _pairs.size(); i++) {
		Pair &pair = _pairs[i];
		if ((pair.read1 == MAX_READ_IDX || pair.read2 == MAX_READ_IDX)
				&& pair.read1 != pair.read2) {
			ReadSetSizeType idx = pair.lesser();
			unmatchedNames[_commonName(_reads[idx].getName())] = idx;
		}
	}

	std::string lastName, common;
	int readNum = 0;
	bool isPairable = true;

	while (readIdx < size) {

		string name = getRead(readIdx).getName();
		readNum = _readNum(name);
		common = _commonName(name);

		unmatchedIt = unmatchedNames.find(common);
		if (unmatchedIt != unmatchedNames.end()) {
			Pair &test = _pairs[unmatchedIt->second];
			if (readNum == 2) {
				if (test.read2 != MAX_READ_IDX) {
					isPairable = false;
					std::cerr
							<< "Detected a conflicting read2. Aborting pair identification: "
							<< name << std::endl;
					break;
				}
				_pairs[unmatchedIt->second].read2 = readIdx;
			} else {
				if (test.read1 != MAX_READ_IDX) {
					isPairable = false;
					std::cerr
							<< "Detected a conflicting read1. Aborting pair identification: "
							<< name << std::endl;
					break;
				}
				_pairs[unmatchedIt->second].read1 = readIdx;
			}
			unmatchedNames.erase(unmatchedIt);
		} else {
			Pair pair;
			if (readNum == 2)
				pair.read2 = readIdx;
			else
				pair.read1 = readIdx;
			if (readNum > 0)
				unmatchedNames[common] = _pairs.size();
			_pairs.push_back(pair);
		}

		readIdx++;
		if (readIdx % 10000000 == 0)
			std::cerr << "Processed " << readIdx << " reads for pairing"
					<< std::endl;
	}

	if (!isPairable) {
		// Pair identification was aborted, re-assigning all 'pairs' to be single reads
		_pairs.clear();
		_pairs.resize(size);
		readIdx = 0;
		while (readIdx < size) {
			_pairs[readIdx] = Pair(readIdx);
			readIdx++;
		}
	}

	return _pairs.size();
}

ProbabilityBases ReadSet::getProbabilityBases() const {
	ProbabilityBases probs(0);
	for(ReadSetSizeType readIdx = 0 ; readIdx < getSize(); readIdx++) {
	    const Read &read = getRead(readIdx);
		probs += read.getProbabilityBases();
    }
	probs *= 1.0 / (double) getSize();
	return probs;
}
ReadSet::ReadSetSizeType ReadSet::getCentroidRead() const {
	ProbabilityBases probs = getProbabilityBases();
	return getCentroidRead(probs);
}

ReadSet::ReadSetSizeType ReadSet::getCentroidRead(const ProbabilityBases &probs) const {
	double bestScore = 0.0;
	ReadSetSizeType bestRead = MAX_READ_IDX;
	for(ReadSetSizeType readIdx = 0 ; readIdx < getSize(); readIdx++) {
		double score = getRead(readIdx).scoreProbabilityBases(probs);
	    if (bestRead == MAX_READ_IDX || bestScore < score ) {
	    	bestRead = readIdx;
	    	bestScore = score;
	    }
	}
	return bestRead;
}

Read ReadSet::getConsensusRead() const {
	ProbabilityBases probs = getProbabilityBases();
    return getConsensusRead(probs, string("C") + boost::lexical_cast<std::string>(getSize()) + string("-") + getRead(0).getName());
}
Read ReadSet::getConsensusRead(const ProbabilityBases &probs, std::string name) {
    stringstream fasta;
    stringstream qual;
    for(size_t i = 0 ; i < probs.size(); i++) {
    	BaseQual base = probs[i].getBaseQual();
    	fasta << base.base;
    	qual << base.qual;
    }
    return Read(name, fasta.str(), qual.str());
}

//
// $Log: ReadSet.cpp,v $
// Revision 1.32  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.31.2.13.2.2  2010-04-16 21:38:39  regan
// addressed part of hack where multi-line records are dangerous to read mmaped
//
// Revision 1.31.2.13.2.1  2010-04-16 17:42:34  regan
// fixed parallelism
//
// Revision 1.31.2.13  2010-04-15 21:31:50  regan
// bugfix in markups and duplicate fragment filter
//
// Revision 1.31.2.12  2010-04-15 20:48:04  regan
// honor mmap-input option in parallel read
//
// Revision 1.31.2.11  2010-04-15 20:42:35  regan
// bugfix in parallel fasta/fastq read
//
// Revision 1.31.2.10  2010-04-15 17:59:52  regan
// made mmap optional
//
// Revision 1.31.2.9  2010-04-14 22:36:06  regan
// round of bugfixes
//
// Revision 1.31.2.8  2010-04-14 20:53:49  regan
// checkpoint and passes unit tests!
//
// Revision 1.31.2.7  2010-04-14 17:51:43  regan
// checkpoint
//
// Revision 1.31.2.6  2010-04-14 05:35:37  regan
// checkpoint. compiles but segfaults
//
// Revision 1.31.2.5  2010-04-12 22:37:47  regan
// checkpoint
//
// Revision 1.31.2.4  2010-04-12 20:59:45  regan
// mmap checkpoint
//
// Revision 1.31.2.3  2010-04-07 22:33:08  regan
// checkpoint mmaping input files
//
// Revision 1.31.2.2  2010-04-05 05:42:53  regan
// checkpoint mmaping input files
//
// Revision 1.31.2.1  2010-04-05 03:32:10  regan
// moved read file reader
//
// Revision 1.31  2010-03-15 18:06:50  regan
// minor refactor and added consensus read
//
// Revision 1.30  2010-03-15 07:44:30  regan
// better logging to track down a non-existent bug (file was corrupted)
//
// Revision 1.29  2010-03-14 17:16:39  regan
// bugfix in reading fastq in parallel
//
// Revision 1.28  2010-03-14 16:55:55  regan
// added centroid methods
//
// Revision 1.27  2010-03-10 13:17:53  regan
// fixed quality ignoring
//
// Revision 1.26  2010-03-03 17:10:26  regan
// added ability to recognize which reads came from which files
//
// Revision 1.25  2010-02-26 12:57:26  regan
// added progress to pair id routine
//
// Revision 1.24  2010-02-26 12:53:05  regan
// reformatted
//
// Revision 1.23  2010-02-22 15:18:37  regan
// bugfix
//
// Revision 1.22  2010-02-22 14:39:30  regan
// optimized to not open too many filesystem threads
// fixed pairing to abort identification when there are duplicate names
//
// Revision 1.21  2010-01-14 00:50:07  regan
// fixes
//
// Revision 1.20  2010-01-13 23:47:44  regan
// made const class modifications
// fixed identify pairs
//
// Revision 1.19  2010-01-13 00:26:49  regan
// fixed some parallelism
// started pair indentification
//
// Revision 1.18  2010-01-08 06:23:07  regan
// fixed to support nested parallel reading of files
//
// Revision 1.17  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.16  2010-01-05 06:44:39  regan
// fixed warnings
//
// Revision 1.15  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.14  2009-12-23 07:16:52  regan
// fixed reading of fasta files
// parallelized reading of multiple files
//
// Revision 1.13  2009-12-22 18:31:41  regan
// parallelized reading fastq if openmp is enabled
//
// Revision 1.12  2009-12-21 22:04:38  regan
// minor optimization
//
// Revision 1.11  2009-12-21 07:54:18  regan
// minor parallelization of reading files step
//
// Revision 1.10  2009-11-28 01:00:07  regan
// fixed bugs and warnings
//
// Revision 1.9  2009-11-21 15:58:29  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.8  2009-11-07 00:28:41  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
// Revision 1.7  2009-11-04 20:14:46  cfurman
// added conversion to uppercase
//
// Revision 1.6  2009-11-04 19:32:03  cfurman
// now reads in fasta (with optional qual) files
//
// Revision 1.5  2009-10-31 00:16:35  regan
// minor changes and optimizations
//
// Revision 1.4  2009-10-23 07:06:59  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
// Revision 1.3  2009-10-21 00:00:58  cfurman
// working on kmers....
//
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//

