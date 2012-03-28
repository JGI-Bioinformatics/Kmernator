//
// Kmernator/src/ReadSet.cpp
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

#include <exception>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>

#include <locale>
#include <algorithm>

#include <boost/lexical_cast.hpp>

#include "ReadSet.h"
#include "Log.h"

ReadSet::MmapSourceVector ReadSet::mmapSources;
void ReadSet::madviseMmaps(int advise) {

#pragma omp critical (readset_madvise_map)
	{
		for(MmapSourceVector::iterator it = mmapSources.begin(); it != mmapSources.end(); it++) {
			if (it->first.is_open())
				madvise(const_cast<char*>(it->first.data()), it->first.size(), advise);
			if (it->second.is_open())
				madvise(const_cast<char*>(it->second.data()), it->second.size(), advise);
		}
	}
}
void ReadSet::addMmaps(MmapSourcePair mmaps) {

#pragma omp critical (readset_madvise_map)
	{
		mmapSources.push_back(mmaps);
	}
	madviseMmapsSequential();
}

void ReadSet::incrementFile(ReadFileReader &reader) {
	incrementFile(reader.getParser());
}
void ReadSet::incrementFile(SequenceStreamParserPtr parser) {
	_filePartitions.addPartition( _reads.size() );
	addMmaps( MmapSourcePair( parser->getMmap(), parser->getQualMmap() ));
}

void ReadSet::_trackSequentialPair(const Read &read) {
	std::string readName = read.getName();
	if (!isPairedRead(readName))
		return;
	if (_reads.size() > 1 && !previousReadName.empty()) {
		if (isPair(previousReadName, read)) {
			// mark the previous read, but not this one
			// signaling, in conjunction with NOT populating this entry in _pairs
			// that these two are sequential and paired...
			// this gets around the parallel nature of building ReadSets in multiple threads
			// as the pairs are identified after all ReadSets are consolidated
			Read &lastRead = _reads[_reads.size() - 2];
			lastRead.markPaired();
			previousReadName.clear();
		} else {
			previousReadName = read.getName();
		}
	} else {
		previousReadName = read.getName();
	}

}

void ReadSet::circularize(long extraLength) {
	for(ReadSetSizeType idx = 0 ; idx < getSize() ; idx++) {
		Read &read = getRead(idx);
		std::string name, fasta, quals;
		name = read.getName();
		fasta = read.getFasta();
		quals = read.getQuals(0, MAX_SEQUENCE_LENGTH, false, true);
		read = Read(name, fasta + fasta.substr(0, extraLength), quals + quals.substr(0,extraLength));
	}
}

void ReadSet::addRead(const Read &read) {
	SequenceLengthType readLength = read.getLength();
	addRead(read, readLength);
}
void ReadSet::addRead(const Read &read, SequenceLengthType readLength, int rank) {
	_reads.push_back(read);
	_baseCount += readLength;
	_setMaxSequenceLength(readLength);
	_trackSequentialPair(read);
	_setFastqStart(read);
	long countReads = getSize();
	LOG_VERBOSE_OPTIONAL(2, ((countReads & 131071) == 0 && (rank & 15) == 0), "Just read " << countReads << " reads");
}
ReadSet::MmapSource ReadSet::mmapFile(string filePath) {
	MmapSource mmap(filePath, FileUtils::getFileSize(filePath));
	addMmaps( MmapSourcePair(mmap, MmapSource()) );
	return mmap;
}

ReadSet::SequenceStreamParserPtr ReadSet::appendAnyFile(string filePath, string filePath2, int rank, int size) {
	if (size == 1 && Options::getOptions().getMmapInput() > 0)
		return appendAnyFileMmap(filePath, filePath2);
	LOG_DEBUG(2, "appendAnyFile(" << filePath << ", " << filePath2 << ", " << rank << ", " << size << ")");
	ReadFileReader reader(filePath, filePath2);
	appendFasta(reader, rank, size);
	incrementFile(reader);
	return reader.getParser();
}

ReadSet::SequenceStreamParserPtr ReadSet::appendAnyFileMmap(string fastaFilePath, string qualFilePath, int rank, int size) {
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
		if (size > 1)
			appendFasta(reader, rank, size);
		else
			appendFastq(mmap);
		break;
	case 1:
		appendFasta(reader, rank, size);
	}
	incrementFile(reader);
	return reader.getParser();
}

void ReadSet::appendAllFiles(OptionsBaseInterface::FileListType &files, int rank, int size) {

#ifdef _USE_OPENMP
	int fileCount = files.size();
	ReadSet myReads[ fileCount ];
	SequenceStreamParserPtr parsers[ fileCount ];
	int numThreads = omp_get_max_threads() / MAX_FILE_PARALLELISM;
	if ( numThreads >= 1 ) {
		omp_set_nested(1);
	} else {
		numThreads = omp_get_max_threads();
		omp_set_nested(0);
	}
	LOG_DEBUG(2, "reading " << fileCount << " file(s) using " << numThreads << " files at a time");

#endif

	long filesSize = files.size();
#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
	for (long i = 0; i < filesSize; i++) {

		LOG_DEBUG(2, "reading " << files[i] << " using " << omp_get_max_threads() << " threads per file");

		string qualFile;
		if (!Options::getOptions().getIgnoreQual()) {
			// test for an implicit qual file
			qualFile = files[i] + ".qual";
			ifstream _qs;
			_qs.open(qualFile.c_str());
			if (! _qs.good()) {
				qualFile.clear();
			} else {
				LOG_DEBUG(2, "detected qual file: " << qualFile);
			}
		}
#ifdef _USE_OPENMP
		// append int this thread's ReadSet buffer (note: line continues)
		parsers[i] = myReads[ i ].appendAnyFile(files[i], qualFile, rank, size);
#else
		SequenceStreamParserPtr parser = appendAnyFile(files[i], qualFile, rank, size);
		incrementFile(parser);
#endif

		LOG_DEBUG(2, "finished reading " << files[i]);

	}
#ifdef _USE_OPENMP
	omp_set_nested(OMP_NESTED_DEFAULT);
	LOG_DEBUG(2,"concatenating ReadSet buffers");
	for(int i = 0; i< (long) files.size(); i++) {
		append(myReads[i]);
		incrementFile(parsers[i]);
	}
#endif

	// Let the kernel know how these pages will be used
	if (Options::getOptions().getMmapInput() == 0)
		madviseMmapsDontNeed();
	else
		madviseMmapsNormal();
}

void ReadSet::append(const Read &read) {
	addRead(read);
}

void ReadSet::append(const ReadSet &reads) {
	unsigned long oldSize = _reads.size();
	unsigned long newSize = oldSize + reads._reads.size();
	_reads.resize(newSize);

	long readSize = reads._reads.size();
#pragma omp parallel for
	for (long i = 0; i < readSize; i++)
		_reads[oldSize + i] = reads._reads[i];

	unsigned long oldPairSize = _pairs.size();
	unsigned long newPairSize = oldPairSize + reads._pairs.size();
	_pairs.resize(newPairSize);

	long pairsSize = reads._pairs.size();
#pragma omp parallel for
	for (long i = 0; i < pairsSize; i++) {
		const Pair &tmp = reads._pairs[i];
		_pairs[oldPairSize + i]
		       = Pair((tmp.read1 == MAX_READ_IDX ? MAX_READ_IDX : tmp.read1
		    		   + oldSize), (tmp.read2 == MAX_READ_IDX ? MAX_READ_IDX
		    				   : tmp.read2 + oldSize));
	}
	_baseCount += reads._baseCount;
	_setMaxSequenceLength(reads.getMaxSequenceLength());

}

ReadSet::SequenceStreamParserPtr ReadSet::appendFasta(string fastaFilePath, string qualFilePath, int rank, int size) {
	ReadFileReader reader(fastaFilePath, qualFilePath);
	appendFasta(reader, rank, size);
	incrementFile(reader);
	return reader.getParser();
}
ReadSet::SequenceStreamParserPtr ReadSet::appendFasta(ReadSet::MmapSource &mmap, int rank, int size) {
	ReadFileReader reader(mmap);
	appendFasta(reader, rank, size);
	incrementFile(reader);
	return reader.getParser();
}

ReadSet::SequenceStreamParserPtr ReadSet::appendFasta(ReadFileReader &reader, int rank, int size) {
	string name, bases, quals;
	LOG_DEBUG(2, "appendFasta(reader, " << rank << ", " << size << ")");
	reader.seekToPartition(rank,size);
	unsigned long firstPos = reader.getPos();
	if (reader.isMmaped() && Options::getOptions().getMmapInput() != 0) {
		RecordPtr recordPtr = reader.getStreamRecordPtr();
		RecordPtr qualPtr = reader.getStreamQualRecordPtr();
		RecordPtr nextRecordPtr = recordPtr;
		std::string name, bases, quals;
		bool isMultiline;
		LOG_DEBUG(3, "Reading mmap file");
		while (reader.nextRead(nextRecordPtr, name, bases, quals, isMultiline)) {

			if (isMultiline) {
				// store the read in memory
				Read read(name,bases,quals);
				addRead(read, bases.length(), rank);
			} else {
				Read read(recordPtr, qualPtr);
				addRead(read, bases.length(), rank);
			}
			recordPtr = nextRecordPtr;
			qualPtr = reader.getStreamQualRecordPtr();

		}
	} else {
		LOG_DEBUG(3, "Reading file stream");
		while (reader.nextRead(name, bases, quals)) {
			Read read(name, bases, quals);
			addRead(read, bases.length(), rank);
		}
	}
	unsigned long lastPos = reader.getPos();
	LOG_DEBUG(2, "Finished reading " << (lastPos - firstPos)/1024 << " KB, " << getSize() << " reads");
	return reader.getParser();
}

ReadSet::SequenceStreamParserPtr ReadSet::appendFastaFile(string &fastaFile, string &qualFile, int rank, int size) {
	LOG_DEBUG(2, "ReadSet::appendFastaFile(" << fastaFile << ", " << qualFile << ", " << rank << ", " << size << ")");
	ReadFileReader reader(fastaFile, qualFile);
	appendFasta(reader, rank, size);
	incrementFile(reader);
	return reader.getParser();
}

ReadSet::SequenceStreamParserPtr ReadSet::appendFastaFile(string &fastaFile, int rank, int size) {
	LOG_DEBUG(2, "ReadSet::appendFastaFile(" << fastaFile << ", " << rank << ", " << size << ")");
	ReadFileReader reader(fastaFile, false);
	appendFasta(reader, rank, size);
	incrementFile(reader);
	return reader.getParser();
}

ReadSet::SequenceStreamParserPtr ReadSet::appendFastaData(string &fastaData, int rank, int size) {
	LOG_DEBUG(2, "ReadSet::appendFastaData(" << fastaData.size() << ", " << rank << ", " << size << ")");
	ReadFileReader reader(fastaData);
	appendFasta(reader, rank, size);
	incrementFile(reader);
	return reader.getParser();
}

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

	LOG_DEBUG(2, "appendFastqBlockedOMP(mmap): " << mmap << " with " << numThreads << " threads");

#pragma omp parallel num_threads(numThreads)
	{

		if (omp_get_num_threads() != numThreads)
			throw "OMP thread count discrepancy!";
		int threadId = omp_get_thread_num();

		ReadSet myReads;
		ReadFileReader reader(mmap);

		if (threadId == 0)
			singleParser = reader.getParser();

		unsigned long lastPos = MAX_UI64;
#pragma omp single
		{
			lastPos = reader.getFileSize();
			blockSize = lastPos / numThreads;

			if (blockSize < 100)
				blockSize = 100;
			LOG_DEBUG(2, "Reading " << mmap << " with " << numThreads << " threads" );
		}


		string name,bases,quals;
		bool hasNext = true;

		unsigned long minimumPos = blockSize * threadId;
		if (minimumPos >= lastPos) {
			seekPos[ threadId ] = lastPos;
		} else {
			reader.seekToNextRecord( minimumPos );
			seekPos[ threadId ] = reader.getPos();
		}
		hasNext = seekPos[ threadId ] < lastPos;

#pragma omp barrier

		if (hasNext && threadId + 1 < numThreads )
			lastPos = seekPos[ threadId + 1 ];

		if (hasNext)
			hasNext = (seekPos[ threadId ] < lastPos);

#pragma omp barrier

		if (!hasNext)
			seekPos[ threadId ] = 0;

#pragma omp barrier

		if (hasNext && threadId != 0 && seekPos[ threadId ] == seekPos[ threadId-1 ])
			hasNext = false;

		if (!hasNext) {
			seekPos[ threadId ] = 0;
			lastPos = 0;
		}

		if (reader.isMmaped() && Options::getOptions().getMmapInput() != 0) {
			RecordPtr recordPtr = reader.getStreamRecordPtr();
			RecordPtr qualPtr = reader.getStreamQualRecordPtr();
			RecordPtr nextRecordPtr = recordPtr;
			std::string name, bases, quals;
			bool isMultiline;
			while (hasNext && reader.nextRead(nextRecordPtr, name, bases, quals, isMultiline)) {
				if (isMultiline) {
					// store the read in memory, as mmap does not currently work
					Read read(name,bases,quals);
					myReads.addRead(read, bases.length(), threadId);
				} else {
					Read read(recordPtr, qualPtr);
					myReads.addRead(read, bases.length(), threadId);
				}

				recordPtr = nextRecordPtr;
				qualPtr = reader.getStreamQualRecordPtr();
				hasNext = (reader.getPos() < lastPos);
			}
		} else {
			while (hasNext && reader.nextRead(name, bases, quals)) {
				Read read(name, bases, quals);
				myReads.addRead(read, bases.length(), threadId);
				hasNext = (reader.getPos() < lastPos);
			}
		}

		numReads[ omp_get_thread_num() ] = myReads.getSize();

#pragma omp critical (readsetGlobals)
		{
			// set global counters
			_baseCount += myReads._baseCount;
			_setMaxSequenceLength(myReads.getMaxSequenceLength());
		}

#pragma omp barrier

#pragma omp single
		{
			ReadSetSizeType newReads = 0;
			for(int i=0; i < numThreads; i++) {
				ReadSetSizeType tmp = newReads;
				newReads += numReads[i];
				numReads[i] = tmp;
			}

			_reads.resize(startIdx + newReads);
		}

		for(ReadSetSizeType j = 0; j < myReads.getSize(); j++)
			_reads[startIdx + numReads[ omp_get_thread_num() ] + j] = myReads.getRead(j);
	}

	incrementFile(singleParser);
	// reset omp variables
	omp_set_nested(OMP_NESTED_DEFAULT);
	return singleParser;
}

ReadSet::SequenceStreamParserPtr ReadSet::appendFastq(ReadSet::MmapSource &mmap)
{
#ifdef _USE_OPENMP
	return appendFastqBlockedOMP(mmap);
#else
	return appendFasta(mmap);
#endif
}

ReadSet::ReadPtr ReadSet::parseMmapedRead(ReadSetSizeType index) const {
	const Read &read = _reads[index];
	ReadPtr readPtr = read.readMmaped();
	return readPtr;
}

string ReadSet::_getReadFileNamePrefix(unsigned int filenum) const {
	if (filenum > Options::getOptions().getInputFiles().size()) {
		return std::string("consensus-") + boost::lexical_cast<std::string>(filenum);
	} else {
		return Options::getOptions().getInputFileSubstring(filenum-1);
	}
}
string ReadSet::getReadFileNamePrefix(ReadSetSizeType index) const {
	unsigned int filenum = getReadFileNum(index);
	return _getReadFileNamePrefix(filenum);
}
// returns the first file to match either read from the pair
string ReadSet::getReadFileNamePrefix(const Pair &pair) const {
	unsigned int filenum1 = -1;
	unsigned int filenum2 = -1;
	if (isValidRead(pair.read1))
		filenum1 = getReadFileNum(pair.read1);
	if (isValidRead(pair.read2))
		filenum2 = getReadFileNum(pair.read2);
	return _getReadFileNamePrefix( filenum1 < filenum2 ? filenum1 : filenum2);
}

const ReadSet::ReadIdxVector ReadSet::getReadIdxVector() const {
	ReadIdxVector readIdxs;
	readIdxs.reserve(_reads.size());
	for(size_t i = 0 ; i < _reads.size(); i++) {
		readIdxs.push_back(i);
	}
	return readIdxs;
}

bool ReadSet::isPairedRead(const std::string &readName) {
	return SequenceRecordParser::isPairedRead(readName);
}
bool ReadSet::isPair(const std::string &readNameA, const std::string &readNameB) {
	return SequenceRecordParser::isPair(readNameA, readNameB);
}
bool ReadSet::isPair(const std::string &readNameA, const Read &readB) {
	std::string readNameB = readB.getName();
	return isPair(readNameA, readNameB);
}
bool ReadSet::isPair(const Read &readA, const Read &readB) {
	std::string readNameA = readA.getName();
	return isPair(readNameA, readB);
}

Read ReadSet::fakePair(const Read &unPaired) {
	std::string name = unPaired.getName();
	int readNum = SequenceRecordParser::readNum(name);
	std::string newName = SequenceRecordParser::commonName(name);
	if (readNum == 1) {
		newName += '2';
	} else if (readNum == 2) {
		newName += '1';
	} else if (readNum == 0)
		throw std::invalid_argument( (std::string("Can not fake pair reads that were not paired end to start with: ") + name).c_str() );
	return Read(newName, "N", "A", true);
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

	long newPairs = 0;
	// first scan the reads adding in sequential pairs which
	// were flagged during addRead()
	long sequentialPairs = 0;
	for(ReadSetSizeType spIdx = readIdx; spIdx < size - 1; spIdx++) {
		Read &read1 = _reads[spIdx];
		Read &read2 = _reads[spIdx+1];
		if (read1.isPaired() && ! read2.isPaired()) {
			// special signal that these two are sequential pairs
			_pairs.push_back( Pair(spIdx, spIdx+1) );
			read2.markPaired();
			spIdx++;
			LOG_DEBUG(4, "Paired sequential reads: " << read1.getName() << " " << read2.getName());
			sequentialPairs++;
		}
	}

	LOG_DEBUG(2, "Paired sequential reads (fast): " << sequentialPairs);

	newPairs += sequentialPairs;

	boost::unordered_map<std::string, ReadSetSizeType> unmatchedNames;
	boost::unordered_map<std::string, ReadSetSizeType>::iterator unmatchedIt;

	// build the unmatchedNames map for existing identified pairs
	for (size_t i = 0; i < _pairs.size(); i++) {
		Pair &pair = _pairs[i];
		if ((pair.read1 == MAX_READ_IDX || pair.read2 == MAX_READ_IDX)
				&& pair.read1 != pair.read2) {
			ReadSetSizeType idx = pair.lesser();
			std::string readName = _reads[idx].getName();
			if (isPairedRead(readName)) {
				unmatchedNames[SequenceRecordParser::commonName( readName )] = idx;
			}
		}
	}

	std::string lastName, common;
	int readNum = 0;
	bool isPairable = true;

	long countNewPaired = 0;
	while (readIdx < size) {

		const Read &read = getRead(readIdx);
		if (read.isPaired()) {
			readIdx++;
			continue;
		}

		string name = read.getName();
		readNum = SequenceRecordParser::readNum(name);
		common =  SequenceRecordParser::commonName(name);

		unmatchedIt = unmatchedNames.find(common);
		if (unmatchedIt != unmatchedNames.end()) {
			Pair &test = _pairs[unmatchedIt->second];
			if (readNum == 2) {
				if (test.read2 != MAX_READ_IDX) {
					if (isPairable) {
						Log::Warn() << "Detected a conflicting read2. Skipping pair identification: "
								<< name << std::endl;
					}
					isPairable = false;
					unmatchedNames.erase(unmatchedIt);
					_pairs.push_back(Pair(MAX_READ_IDX, readIdx));
					readIdx++;
					continue;
				}
				test.read2 = readIdx;
			} else {
				if (test.read1 != MAX_READ_IDX) {
					if (isPairable) {
						Log::Warn()
						<< "Detected a conflicting read1. Skipping pair identification: "
						<< name << std::endl;
					}
					isPairable = false;
					unmatchedNames.erase(unmatchedIt);
					_pairs.push_back(Pair(readIdx, MAX_READ_IDX));
					readIdx++;
					continue;
				}
				test.read1 = readIdx;
			}
			unmatchedNames.erase(unmatchedIt);
			getRead(test.read1).markPaired();
			getRead(test.read2).markPaired();
			newPairs++;

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
		if (++countNewPaired % 10000000 == 0)
			LOG_VERBOSE(3, "Processed " << countNewPaired << " pairs for pairing");
	}

	LOG_VERBOSE(2, "Identified new pairs: " << newPairs << " (" << sequentialPairs << " sequential)");

	return _pairs.size();
}

ProbabilityBases ReadSet::getProbabilityBases(unsigned char minQual) const {
	ProbabilityBases probs(0);
	for(ReadSetSizeType readIdx = 0 ; readIdx < getSize(); readIdx++) {
		const Read &read = getRead(readIdx);
		probs += read.getProbabilityBases(minQual);
	}
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
	assert(bestRead != MAX_READ_IDX);
	return bestRead;
}

Read ReadSet::getConsensusRead(unsigned char minQual) const {
	ProbabilityBases probs = getProbabilityBases(minQual);
	std::string consensusName = string("C") + boost::lexical_cast<std::string>(getSize()) + string("-") + getRead(0).getName();
	Read consensus = getConsensusRead(probs, consensusName);
	if (Log::isDebug(2)) {

		{
			std::stringstream ss;
			ss << "Consensus details for: " << consensusName << std::endl;
			for(ReadVector::const_iterator it = _reads.begin(); it != _reads.end(); it++) {
				ss << it->toString() << std::endl;
			}
			ss << consensus.toString() << std::endl;
			ss << probs.toString() << std::endl;

			std::string s = ss.str();
			LOG_DEBUG(2, s);
		}
	}
	return consensus;
}
Read ReadSet::getConsensusRead(const ProbabilityBases &probs, std::string name) {
	std::string fasta(probs.size(), ' ');
	std::string qual(probs.size(), ' ');
	for(size_t i = 0 ; i < probs.size(); i++) {
		BaseQual base = probs[i].getBaseQual();
		fasta[i] = base.base;
		qual[i] = base.qual;
	}
	return Read(name, fasta, qual);
}

ReadSet ReadSet::randomlySample(ReadSet::ReadSetSizeType maxReads) const {
	ReadSetSizeType size = getSize();
	if (size <= maxReads)
		return *this;

	Random<ReadSetSizeType>::Set readIds = Random<ReadSetSizeType>::sample(size, maxReads);

	ReadSet sampledReadSet;
	for(Random<ReadSetSizeType>::SetIterator it = readIds.begin(); it != readIds.end(); it++) {
		sampledReadSet.append( getRead( *it ));
	}
	LOG_DEBUG(4, "ReadSet::randomlySample(): sampled " << size << " to: " << sampledReadSet.getSize());
	assert( sampledReadSet.getSize() == std::min(size, maxReads));
	return sampledReadSet;
}


