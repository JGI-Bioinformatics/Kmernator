//
// Kmernator/src/ReadSet.cpp
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

bool ReadSet::_isSequentialPair(const Read &read) {
	std::string readName = read.getName();
	std::string comment = read.getComment();
	LOG_DEBUG(5, "_isSequentialPair(" << readName << ")");
	if (!isPairedRead(readName, comment)) {
		previousReadName.clear();
		previousReadComment.clear();
		return false;
	}
	if (!previousReadName.empty()) {
		if (isPair(previousReadName, previousReadComment, read)) {
			previousReadName.clear();
			previousReadComment.clear();
			return true;
		} else {
			previousReadName = readName;
			previousReadComment = comment;
			return false;
		}
	} else {
		previousReadName = readName;
		previousReadComment= comment;
		return false;
	}
}

void ReadSet::circularize(long extraLength) {
	for(ReadSetSizeType idx = 0 ; idx < getSize() ; idx++) {
		Read &read = getRead(idx);
		std::string name, fasta, quals, comment;
		name = read.getName();
		comment = read.getComment();
		fasta = read.getFasta();
		quals = read.getQuals(0, MAX_SEQUENCE_LENGTH, false, true);
		read = Read(name, fasta + fasta.substr(0, extraLength), quals + quals.substr(0,extraLength), comment);
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
	_setFastqStart(read);
	long countReads = getSize();
	LOG_DEBUG(5, "addRead(" << read.getName() << ", " << readLength << ", " << rank << ") len: " << read.getLength() << " size: " << countReads);
	LOG_VERBOSE_OPTIONAL(2, ((countReads & 131071) == 0 && (rank & 15) == 0), "Just read " << countReads << " reads");
}
ReadSet::MmapSource ReadSet::mmapFile(string filePath) {
	MmapSource mmap(filePath, FileUtils::getFileSize(filePath));
	addMmaps( MmapSourcePair(mmap, MmapSource()) );
	return mmap;
}

ReadSet::SequenceStreamParserPtr ReadSet::appendAnyFile(string filePath, string filePath2, int rank, int size) {
	if (size == 1 && Options::getOptions().getMmapInput())
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

	// TODO fix this assertion
	// assert(files == Options::getOptions().getInputFiles() || files == Options::getOptions().getReferenceFiles());
	int fileCount = files.size();
	ReadSet myReads[ fileCount ];
	SequenceStreamParserPtr parsers[ fileCount ];
	int numThreads = std::min(omp_get_max_threads(), MAX_FILE_PARALLELISM);
	LOG_DEBUG(2, "reading " << fileCount << " file(s) using " << numThreads << " files at a time");

	OptionsBaseInterface::FileListType reorder;
	long filesSize = files.size();
#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
	for (long i = 0; i < filesSize; i++) {

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
		// append int this thread's ReadSet buffer (note: line continues)
		parsers[i] = myReads[ i ].appendAnyFile(files[i], qualFile, rank, size);
		LOG_DEBUG_OPTIONAL(2, true, "finished reading " << files[i]);
		myReads[i].identifyPairs();
		if (myReads[i].hasPairs() && myReads[i].getPairSize() != myReads[i].getSize() / 2)
			LOG_WARN(1, "Paired file: " << files[i] << " has incomplete number of pairs in my slice: " << myReads[i].getPairSize() << " vs reads: " <<  myReads[i].getSize());

	}
	LOG_DEBUG(2,"concatenating ReadSet buffers");
	// First add any files with paired reads
	for(int i = 0; i< (long) files.size(); i++) {
		if (myReads[i].hasPairs()) {
			append(myReads[i]);
			myReads[i].clear();
			incrementFile(parsers[i]);
			reorder.push_back( files[i] );
		}
	}
	// Then add any unpaired files
	for(int i = 0; i< (long) files.size(); i++) {
		if (myReads[i].getSize() > 0) {
			append(myReads[i]);
			myReads[i].clear();
			incrementFile(parsers[i]);
			reorder.push_back( files[i] );
		}
	}
	LOG_DEBUG_OPTIONAL(1, Logger::isMaster(), "Reordered input files to scan paired files first: " << Log::toString<OptionsBaseInterface::FileListType>(reorder) );
	files = reorder;

	// Let the kernel know how these pages will be used
	if (!Options::getOptions().getMmapInput())
		madviseMmapsDontNeed();
	else
		madviseMmapsNormal();

	if (Log::isDebug(4)) {
		std::stringstream ss;
		for(ReadSetSizeType i = 0; i < getPairSize(); i++) {
			const Pair &pair = getPair(i);
			ss << "Pair: " << i << " (" << pair.read1 << ", " << pair.read2 << "): " << (isValidRead(pair.read1)?getRead(pair.read1).getName(): "x1") << " " << (isValidRead(pair.read2)?getRead(pair.read2).getName(): "x2") << "\n";
		}
		std::string s = ss.str();
		LOG_DEBUG(4, "AppendAllFiles(): IdentifedPairs:" << getPairSize() << "\n" << s);
	}
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
	unsigned long newPairSize = oldPairSize;
	if (reads.isSecondReadOnly())
		identifyPairs();
	else {
		newPairSize += reads._pairs.size();
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
	string name, bases, quals, comment;
	LOG_DEBUG(2, "appendFasta(reader, " << rank << ", " << size << ")");
	reader.seekToPartition(rank,size);
	unsigned long firstPos = reader.getPos();
	if (reader.isMmaped() && Options::getOptions().getMmapInput()) {
		RecordPtr recordPtr = reader.getStreamRecordPtr();
		RecordPtr qualPtr = reader.getStreamQualRecordPtr();
		RecordPtr nextRecordPtr = recordPtr;
		bool isMultiline;
		LOG_DEBUG(3, "Reading mmap file");
		while (reader.nextRead(nextRecordPtr, name, bases, quals, comment, isMultiline)) {

			if (isMultiline) {
				// store the read in memory
				Read read(name,bases,quals,comment);
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
		while (reader.nextRead(name, bases, quals, comment)) {
			Read read(name, bases, quals, comment);
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

ReadSet::SequenceStreamParserPtr ReadSet::appendFastq(ReadSet::MmapSource &mmap)
{
	return appendFasta(mmap);
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

bool ReadSet::isPairedRead(const std::string readName, const std::string comment) {
	return SequenceRecordParser::isPairedRead(readName, comment);
}
bool ReadSet::isPair(const std::string readNameA, const std::string readNameB, const std::string commentA, const std::string commentB) {
	return SequenceRecordParser::isPair(readNameA, readNameB, commentA, commentB);
}
bool ReadSet::isPair(const std::string readNameA, const std::string commentA, const Read &readB) {
	return isPair(readNameA, readB.getName(), commentA, readB.getComment());
}
bool ReadSet::isPair(const Read &readA, const Read &readB) {
	return isPair(readA.getName(), readB.getName(), readA.getComment(), readB.getComment());
}

Read ReadSet::fakePair(const Read &unPaired) {
	std::string name = unPaired.getName();
	std::string comment = unPaired.getComment();
	int readNum = SequenceRecordParser::readNum(name, comment);
	std::string newName = SequenceRecordParser::commonName(name);
	if (SequenceRecordParser::isCommentCasava18(comment)) {
		if (readNum == 1)
			comment[0] = '2';
		else if (readNum == 2)
			comment[0] = '1';
		else
			LOG_THROW( "ReadSet::fakePair(): Can not fake pair reads that were not paired end to start with: " << name << " " << comment);
	} else {
		if (newName[newName.length()-1] != '/')
			newName += '/';

		if (readNum == 1) {
			newName += '2';
		} else if (readNum == 2) {
			newName += '1';
		} else if (readNum == 0)
			LOG_THROW( "ReadSet::fakePair(): Can not fake pair reads that were not paired end to start with: " << name);
	}
	return Read(newName, "N", "A", comment + ":fakePair", true);
}

ReadSet::ReadSetSizeType ReadSet::identifyPairs() {
	ReadSetSizeType size = getSize();
	LOG_DEBUG_OPTIONAL(2, true, "ReadSet::identifyPairs(): " << size << " reads, " << _pairs.size() << " pairs starting");

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
	for(ReadSetSizeType spIdx = readIdx; spIdx < size; spIdx++) {
		Read &read = _reads[spIdx];
		if (_isSequentialPair(read)) {
			assert(spIdx > 0);
			read.markPaired();
			Read &readMate = _reads[spIdx-1];
			readMate.markPaired();
			_pairs.push_back( Pair(spIdx-1, spIdx) );
			sequentialPairs++;
			LOG_DEBUG(4, "Paired sequential reads: " << spIdx-1 << ", " << spIdx << ": " << read.getName() << " " << readMate.getName());
		}
	}

	LOG_DEBUG_OPTIONAL(2, true, "Paired sequential reads (fast): " << sequentialPairs);

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
		string comment = read.getComment();
		readNum = SequenceRecordParser::readNum(name, comment);
		common =  SequenceRecordParser::commonName(name);

		unmatchedIt = unmatchedNames.find(common);
		if (unmatchedIt != unmatchedNames.end()) {
			Pair &test = _pairs[unmatchedIt->second];
			if (readNum == 2) {
				if (test.read2 != MAX_READ_IDX) {
					if (isPairable) {
						LOG_WARN(1, "Detected a conflicting read2. Skipping pair identification: " << name << " common: " << common << " readNum: " << readNum); // << "\n" << this->toString());
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
						LOG_WARN(1, "Detected a conflicting read1. Skipping pair identification: " << name << " common: " << common << " readNum: " << readNum);
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

	LOG_VERBOSE_OPTIONAL(2, true, "Identified new pairs: " << newPairs << " (" << sequentialPairs << " sequential)");

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
Read ReadSet::getConsensusRead(const ProbabilityBases &probs, std::string name, std::string comment) {
	std::string fasta(probs.size(), ' ');
	std::string qual(probs.size(), ' ');
	for(size_t i = 0 ; i < probs.size(); i++) {
		BaseQual base = probs[i].getBaseQual();
		fasta[i] = base.base;
		qual[i] = base.qual;
	}
	return Read(name, fasta, qual, comment);
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


