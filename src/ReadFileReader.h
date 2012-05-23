//
// Kmernator/src/ReadFileReader.h
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


// added support for new Illumina 1.8 fasta/q file header & comment format
// @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
// 1: member of a pair
// Y: read fails filter
// control bits (0 is no control, even is something else)
// ATCACG barcode
// virtually translate to
// @EAS139:136:FC706VJ:2:2104:15343:197393/1 Y:18:ATCACG
// (and set the read discard flag if filtered 'Y')

#ifndef _READ_FILE_READER_H
#define _READ_FILE_READER_H
#include <cstring>
#include <cstdlib>

#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include "config.h"
#include "Log.h"
#include "Utils.h"

using namespace std;

class ReadFileReader {
public:
	typedef Kmernator::RecordPtr RecordPtr;
	typedef Kmernator::MmapSource MmapSource;
	typedef Kmernator::MmapIStream MmapIStream;
	typedef Kmernator::FilteredIStream FilteredIStream;

	class SequenceStreamParser;
	typedef boost::shared_ptr< SequenceStreamParser > SequenceStreamParserPtr;

private:
	SequenceStreamParserPtr _parser;
	string _path;
	ifstream _ifs;
	ifstream _qs;
	istringstream _iss;
	Kmernator::FilteredIStream _fis;
	int _streamType;

public:
	ReadFileReader(): _parser() {}

	ReadFileReader(string fastaFilePath, bool autoFindQual) :
		_parser(), _path(fastaFilePath), _streamType(0) {

		openFastaFile(fastaFilePath);

		if (autoFindQual) {
			std::string qualFilePath;
			openQualFile(qualFilePath);
		} else {
			_qs.close();
		}
		setParser(_ifs, _qs);

	}
	ReadFileReader(string fastaFilePath, string qualFilePath, bool autoFindQual = GeneralOptions::getOptions().getIgnoreQual()) :
		_parser(), _path(fastaFilePath), _streamType(0) {

		openFastaFile(fastaFilePath);

		if (!qualFilePath.empty() || autoFindQual) {
			openQualFile(qualFilePath);
		} else {
			qualFilePath.clear();
			_qs.close();
		}

		setParser(_ifs, _qs);
	}

	ReadFileReader(string &fasta) :
		_iss(fasta), _streamType(1) {
		LOG_DEBUG(2, "ReadFileReader(" << fasta.length() << "): Constructor on a string-based fasta file");
		_parser = SequenceStreamParserPtr(new FastaStreamParser(_iss));
	}

	ReadFileReader(MmapSource &mmap) : _streamType(2) {
		setReader(mmap);
	}

	ReadFileReader(MmapSource &mmap1, MmapSource &mmap2) : _streamType(2) {
		setReader(mmap1,mmap2);
	}

	~ReadFileReader() {
		switch (_streamType) {
		case(0) : _ifs.close();	_qs.close(); break;
		case(1) : break;
		case(2) : break;
		}
	}

	void openFastaFile(string fastaFilePath) {

		_ifs.open(_path.c_str());

		if (_ifs.fail())
			LOG_THROW("ReadFileReader::openFastaFile(): Could not open : " << _path);

	}
	void openQualFile(string qualFilePath) {
		if (!qualFilePath.empty()) {
			_qs.open(qualFilePath.c_str());
			if (_qs.fail())
				LOG_THROW("ReadFileReader::openQualFile(): Could not open : " << qualFilePath);
		} else {
			_qs.open((_path + ".qual").c_str());
			LOG_DEBUG(3, "ReadFileReader() opened " << _path << ".qual " << _qs.good());
		}
		LOG_DEBUG(2, "ReadFileReader(" << _path << ", " << qualFilePath << ")");

	}

	std::string getFilePath() {
		return _path;
	}

	void setReader(MmapSource &mmap) {
		assert(mmap.is_open());
		LOG_DEBUG(3, "setReader(mmap):" << (void*)mmap.data());
		setParser(mmap, *mmap.data());
	}
	void setReader(MmapSource &mmap1, MmapSource &mmap2) {
		assert(mmap1.is_open());
		assert(mmap2.is_open());
		LOG_DEBUG(3, "setReader(mmap, mmap):" << (void*)mmap1.data() << " " << (void*)mmap2.data());
		setParser(mmap1,mmap2);
	}
	void setParser(istream &fs1) {
		if (fs1.fail() || !fs1.good())
			LOG_THROW("ReadFileRader::setParser(): istream fail() or !good()");
		setParser(fs1, fs1.peek());
	}
	template<typename U> void setParser(U &data, char marker) {
		LOG_DEBUG(3, "setParser(U)");
		if (marker == '@')
			_parser = SequenceStreamParserPtr(new FastqStreamParser(data));
		else if (marker == '>')
			_parser = SequenceStreamParserPtr(new FastaStreamParser(data));
		else
			LOG_THROW("ReadFileReader::setParser(): Unknown file format: marker == " << marker);
	}
	SequenceStreamParserPtr getParser() {
		return _parser;
	}
	void setParser(istream &fs1, istream &fs2) {
		if (fs2.fail() || fs2.eof() || !fs2.good() ) {
			setParser(fs1);
		} else {
			LOG_DEBUG(3, "setParser(istream,istream) FastaQualStreamParser");
			_parser = SequenceStreamParserPtr(new FastaQualStreamParser(fs1, fs2));
		}
	}
	void setParser(MmapSource &mmap1, MmapSource &mmap2) {
		_parser = SequenceStreamParserPtr(new FastaQualStreamParser(mmap1, mmap2));
		LOG_DEBUG(3, "setParser(mmap,mmap) FastaQualStreamParser");
	}

	bool isMmaped() {
		return _parser->isMmaped();
	}
	RecordPtr getStreamRecordPtr() const {
		return _parser->getStreamRecordPtr();
	}
	RecordPtr getStreamQualRecordPtr() const {
		return _parser->getStreamQualRecordPtr();
	}
	MmapSource getMmap() const {
		return _parser->getMmap();
	}
	MmapSource getQualMmap() const {
		return _parser->getQualMmap();
	}
	const string &getLastName() const {
		return _parser->getName();
	}
	const string &getLastBases() const {
		return _parser->getBases();
	}
	const string &getLastQuals() const {
		return _parser->getQuals();
	}
	bool getLastMultiline() const {
		return _parser->isMultiline();
	}
	bool nextRead(RecordPtr &recordStart) {
		LOG_DEBUG(5, "nextRead(recordStart = " << (void*) recordStart << ")");
		try {
			assert(recordStart == _parser->getStreamRecordPtr());
			RecordPtr endPtr = _parser->readRecord();
			while (_parser->isDiscardFiltered() && endPtr != NULL && (!_parser->getName().empty()) && (!_parser->isGood()) ) {
				LOG_DEBUG_OPTIONAL(2, true, "nextRead(): Skipping failed-filter read: " << _parser->getName());
				endPtr = _parser->readRecord();
			}
			if (endPtr == NULL || _parser->getName().empty()) {
				_parser->resetBuffers();
				recordStart = NULL;
				return false;
			} else {
				recordStart = endPtr;
				return true;
			}
		} catch (runtime_error &e) {
			stringstream error;
			error << e.what() << " in file '" << _path << "' at line "
					<< _parser->lineNumber() << " position:" << _parser->tellg() << " " \
					<< _parser->getName() << _parser->getBases() << _parser->getQuals() \
					<< " '" << _parser->getLineBuffer() << _parser->nextLine() << "'";
			LOG_THROW("ReadFileReader::nextRead(RecordPtr): " << error.str());
		}
	}
	bool nextRead(RecordPtr &recordStart, string &name, string &bases, string &quals, bool &isMultiline) {
		bool passed = nextRead(recordStart);
		if (passed) {
			name = getLastName();
			bases = getLastBases();
			quals = getLastQuals();
			isMultiline = getLastMultiline();
		}
		return passed;
	}
	bool nextRead(string &name, string &bases, string &quals) {
		try {
			_parser->readRecord();
			name = _parser->getName();
			while (_parser->isDiscardFiltered() && (!name.empty()) && (!_parser->isGood())) {
				LOG_DEBUG_OPTIONAL(2, true, "nextRead(,,): Skipping failed-filter read: " << _parser->getName());
				_parser->readRecord();
				name = _parser->getName();
			}
			if (name.empty())
				return false;

			bases = _parser->getBases();
			std::transform(bases.begin(), bases.end(), bases.begin(), ::toupper);

			quals = _parser->getQuals();

			if (quals.length() != bases.length())
				if (!(quals.length() == 1 && quals[0] == Read::REF_QUAL))
					LOG_THROW("ReadFileReader::nextRead(name,base,qual): Number of bases and quals not equal: " << bases << " " << quals);
			LOG_DEBUG(5, "nextRead(name, bases, quals) " << name);
			return true;
		}

		catch (runtime_error &e) {
			stringstream error;
			error << e.what() << " in file '" << _path << "' at line "
					<< _parser->lineNumber() << " position:" << _parser->tellg() << " " \
					<< _parser->getName() << _parser->getBases() << _parser->getQuals() \
					<< " '" << _parser->getLineBuffer() << _parser->nextLine() << "'";
			LOG_THROW("ReadFileReader::nextRead(name,base,qual): " << error.str());
		}
	}

	unsigned long getFileSize() {
		if (_parser->isMmaped()) {
			return _parser->getMmapFileSize();
		} else {
			long size = 0;
			switch(_streamType) {
			case(0) : size = FileUtils::getFileSize(_ifs); break;
			case(1) : size = FileUtils::getFileSize(_iss); break;
			case(2) : size = 0; break;
			}
			return size;
		}
	}

	unsigned long getBlockSize(unsigned int numThreads) {
		unsigned long fileSize = getFileSize();
		unsigned long blockSize = fileSize / numThreads;
		LOG_DEBUG(3, "ReadFileReader(" << _path << ")::getBlockSize(" << numThreads <<"): " << blockSize << " out of " << fileSize);
		return blockSize;
	}

	void setPos(unsigned long pos) {
		_parser->setPos(pos);
	}
	inline unsigned long getPos() {
		return _parser->getPos();
	}
	inline unsigned long getLastPos() {
		return _parser->getLastPos();
	}
	void setLastPos(unsigned long pos) {
		_parser->setLastPos(pos);
	}

	bool seekToNextRecord(unsigned long minimumPos) {
		return seekToNextRecord(minimumPos, true);
	}
	bool seekToNextRecord(unsigned long minimumPos, bool byPair) {
		return _parser->seekToNextRecord(minimumPos, byPair);
	}
	int getType() const {
		return _parser->getType();
	}

	void seekToPartition(int rank, int size) {
		unsigned long lastPos = getFileSize();
		unsigned long firstPos = 0;
		LOG_DEBUG(2, "ReadFileReader(" << _path <<")::seekToPartition(" << rank << ", " << size << ")");
		setLastPos(lastPos);
		if (size > 1) {
			unsigned long blockSize = getBlockSize(size);
			if (rank + 1 != size ) {
				seekToNextRecord( blockSize * (rank+1) );
				lastPos = getPos();
			}
			if (seekToNextRecord( blockSize * rank )) {
				firstPos = getPos();
			} else {
				firstPos = lastPos;
			}
		}
		setLastPos(lastPos);
		LOG_DEBUG(2, "ReadFileReader(" << _path <<")::seekToPartition(" << rank << ", " << size << ") " << firstPos << ", reading until " << lastPos);
	}

	bool eof() const {
		return _parser->endOfStream();
	}

public:

	class SequenceStreamParser {
	private:
		istream * _stream;
		unsigned long _line;
		mutable unsigned long _pos;
		unsigned long _lastPos;
		bool _discardFiltered;

	protected:
		char _marker;
		ReadFileReader::MmapSource _mmap;
		ReadFileReader::RecordPtr _lastPtr;
		bool _freeStream;
		mutable std::vector<string> _nameBuffer;
		mutable std::vector<string> _basesBuffer;
		mutable std::vector<string> _qualsBuffer;
		mutable std::vector<string> _lineBuffer;
		mutable std::vector<bool>   _isMultiline;
		mutable std::vector<bool>   _isGood;

	public:
		// returns a termination pointer (the start of next record or lastByte+1 at end)
		virtual RecordPtr readRecord(RecordPtr recordPtr) const = 0;
		virtual RecordPtr readRecord() = 0;
		RecordPtr readRecord(RecordPtr recordPtr, string &name, string &bases, string &quals) const {
			RecordPtr end = readRecord(recordPtr);
			if (end != NULL) {
				name = getName();
				bases = getBases();
				quals = getQuals();
			} else {
				LOG_DEBUG(5, "SequenceStreamParser()::readRecord(" << recordPtr <<",,,) found the end");
				resetBuffers();
				return NULL;
			}
			return end;
		}

		static inline string &nextLine(string &buffer, RecordPtr &recordPtr) {
			return SequenceRecordParser::nextLine(buffer, recordPtr);
		}
		inline string &nextLine(string &buffer) {
			_line++;
			buffer.clear();
			if (isPastPartition())
				return buffer;
			getline(*_stream, buffer);
			_pos += buffer.length() + 1;
			LOG_DEBUG(6, "SequenceStreamParser::nextLine(" << buffer << ") position at: " << _pos);
			return buffer;
		}
		inline string &nextLine() {
			int threadNum = omp_get_thread_num();
			nextLine(_lineBuffer[threadNum]);
			return _lineBuffer[threadNum];
		}

		SequenceStreamParser(istream &stream, char marker) :
			_stream(&stream), _line(0), _pos(0), _lastPos(-1), _discardFiltered(true), _marker(marker), _mmap(), _lastPtr(NULL), _freeStream(false) {
			if (!_stream->good() || _stream->fail())
				LOG_THROW("SequenceStreamParser(" << _stream << ", " << marker << ") has an invalid stream!");
			LOG_DEBUG(4, "SequenceStreamParser(istream, " << marker << ")");
			setBuffers();
		}
		SequenceStreamParser(ReadFileReader::MmapSource &mmap, char marker) :
			_stream(NULL), _line(0), _pos(0), _lastPos(-1), _discardFiltered(true), _marker(marker), _mmap( mmap ), _lastPtr(NULL), _freeStream(false) {
			if (!_mmap.is_open() || !_mmap.size()>0 || _mmap.data() == NULL)
				LOG_THROW("SequenceStreamParser(" << _mmap << ", " << marker << ") has an invalid memory map!");

			_stream = new MmapIStream(_mmap);
			_lastPtr = _mmap.data() + _mmap.size();
			_freeStream = true;
			LOG_DEBUG(4, "SequenceStreamParser(MmapSource, " << marker << ")");
			setBuffers();
		}
		SequenceStreamParser(ReadFileReader::FilteredIStream &stream, char marker) :
			_line(0), _pos(0), _lastPos(-1), _discardFiltered(true), _marker(marker), _mmap(), _lastPtr(NULL), _freeStream(false) {
			_stream = (istream *) &stream;
			if (!_stream->good() || _stream->fail())
				LOG_THROW("SequenceStreamParser(FilteredIStream " << _stream << ", " << marker << ") has an invalid stream!");

			LOG_DEBUG(4, "SequenceStreamParser(FilteredIStream, " << marker << ")");
			setBuffers();
		}
		void setBuffers() {
			_nameBuffer.assign(OMP_MAX_THREADS_DEFAULT, string());
			_basesBuffer.assign(OMP_MAX_THREADS_DEFAULT, string());
			_qualsBuffer.assign(OMP_MAX_THREADS_DEFAULT, string());
			_lineBuffer.assign(OMP_MAX_THREADS_DEFAULT, string());
			_isMultiline.assign(OMP_MAX_THREADS_DEFAULT, false);
			_isGood.assign(OMP_MAX_THREADS_DEFAULT, false);
		}

		void resetBuffers() const {
			int threadNum = omp_get_thread_num();
			_nameBuffer[threadNum].clear();
			_basesBuffer[threadNum].clear();
			_qualsBuffer[threadNum].clear();
			_lineBuffer[threadNum].clear();
			_isMultiline[threadNum] = false;
			_isGood[threadNum] = false;
		}

		virtual ~SequenceStreamParser() {
			// TODO fix this -- open a new mmap that can be closed!
			// Do not free/close istream as that closes the mmap too!!!
			//if (_freeStream)
			//	delete _stream;
		}

		void setLastPos(unsigned long lastPos) {
			_lastPos = lastPos;
		}
		unsigned long lineNumber() {
			return _line;
		}
		unsigned long tellg() {
			return _pos = _stream->tellg();
		}

		void reset() {
			_stream->clear(_stream->rdstate() & ( ~ios_base::eofbit ) );
			_stream->clear(_stream->rdstate() & ( ~ios_base::failbit ) );
		}
		void seekg(unsigned long pos) {
			_pos = pos;
			reset(); // reset a potential eof...
			_stream->seekg(_pos);
		}
		bool endOfStream() {
			return _stream->eof();
		}
		int peek() {
			return _stream->peek();
		}

		string &getLineBuffer() const {
			int threadNum = omp_get_thread_num();
			return _lineBuffer[threadNum];
		}
		string &getName() const {
			int threadNum = omp_get_thread_num();
			return _nameBuffer[threadNum];
		}
		string &getBases() const {
			int threadNum = omp_get_thread_num();
			return _basesBuffer[threadNum];
		}
		string &getQuals() const {
			int threadNum = omp_get_thread_num();
			return _qualsBuffer[threadNum];
		}
		bool isMultiline() const {
			int threadNum = omp_get_thread_num();
			return _isMultiline[threadNum];
		}
		bool isGood() const {
			int threadNum = omp_get_thread_num();
			return _isGood[threadNum];
		}
		void setIsGood(bool isGood) const {
			int threadNum = omp_get_thread_num();
			_isGood[threadNum] = isGood;
		}
		inline bool &isDiscardFiltered() {
			return _discardFiltered;
		}


		virtual string &readName() {
			std::string &name = getName();
			nextLine( name );

			int count = 0;
			while (name.length() == 0 || name[0] != _marker) // skip empty lines at end of stream
			{
				LOG_DEBUG(5, "SequenceStreamParser::readName() found nothing on the last line: " << name << " looking for '" << _marker << "' at " << getPos());
				if (isPastPartition() || endOfStream()) {
					name.clear();
					return name;
				}
				nextLine( name );
				if (++count > 100000)
					break;
			}
			LOG_DEBUG(5, "SequenceStreamParser::readName() found " << name << " at " << getPos());

			if (name.length() == 0 || name[0] != _marker)
				LOG_THROW("Missing name marker '" << _marker << "'" << " in '" << name << "' at " << getPos() << " after " << count << " attempts to find the marker");

			// remove marker and any extra comments or fields
			bool isGood = SequenceRecordParser::trimName( name );
			setIsGood(isGood);

			return name;
		}


		virtual int getType() const = 0;

		unsigned long getPos() const {
			return _pos;
		}
		void setPos(unsigned long newPos) const {
			_pos = newPos;
		}
		unsigned long getLastPos() const {
			return _lastPos;
		}
		bool isPastPartition() const {
			return _pos >= _lastPos;
		}
		bool isMmaped() const {
			return _mmap.is_open();
		}
		unsigned long getMmapFileSize() const {
			return _mmap.size();
		}
		RecordPtr getStreamRecordPtr() const {
			if (_mmap.is_open()) {
				return _mmap.data() + getPos();
			} else {
				return NULL;
			}
		}
		MmapSource getMmap() const {
			return _mmap;
		}
		virtual MmapSource getQualMmap() const = 0;
		RecordPtr getLastRecordPtr() const {
			return _lastPtr;
		}
		virtual RecordPtr getStreamQualRecordPtr() const = 0;
		virtual RecordPtr getLastQualRecordPtr() const = 0;

		virtual bool seekToNextRecord(unsigned long minimumPos, bool byPair) {
			int threadNum = omp_get_thread_num();
			LOG_DEBUG(2, "seekToNextRecord(" << minimumPos << ", " << byPair << ")");
			// get to the first line after the pos
			if (minimumPos > 0) {
				seekg(minimumPos - 1);
				if (endOfStream())
					return false;
				if (peek() == '\n') {
					seekg(minimumPos);
				} else
					nextLine();
			} else {
				seekg(0);
				assert(!endOfStream());
				return true;
			}
			LOG_DEBUG(2, "seeked to " << tellg() );

			if (endOfStream()) {
				LOG_DEBUG(3, "At endofstream already1");
				return false;
			}

			while ( (!isPastPartition()) && (!endOfStream()) && _stream->peek() != _marker) {
				std::string line = nextLine();
				LOG_DEBUG(4, "reading line (" << tellg() << "): " << line);
			}
			if (isPastPartition()) {
				LOG_WARN(1, "seekToNextRecord(" << minimumPos << "): is past partition: " << tellg());
				return false;
			}

			if (_marker == '@' && (!endOfStream())) {
				// since '@' is a valid quality character in a FASTQ (sometimes)
				// verify that the next line is not also starting with '@', as that would be the true start of the record
				unsigned long tmpPos = tellg();
				string current = _lineBuffer[threadNum];
				string tmp = nextLine();
				if (endOfStream() || _stream->peek() != _marker) {
					seekg(tmpPos);
					_lineBuffer[threadNum] = current;

					LOG_DEBUG(3, "Correctly picked record break point to read fastq in parallel: "
							<< current << "\n" << tmp << " " << tellg());


				} else {
					LOG_DEBUG(3, "Needed to skip an extra line to read fastq in parallel: "
							<< current << "\n" << tmp << " " << tellg());
				}
			}

			if (endOfStream()) {
				LOG_DEBUG(3, "At endofstream already2");
				return false;
			}

			if (byPair) {
				// now verify that we have not split a sequential pair
				unsigned long here1 = tellg();
				string name1;
				readRecord();
				if ( (name1 = getName()).empty() || endOfStream()) {
					LOG_DEBUG(3, "Found endofstream 0 records in leaving at eof to preserve pair");
					return false;
				} else {
					LOG_DEBUG(3, "Checking pairing against " << name1);
				}

				unsigned long here2 = tellg();
				string name2;
				readRecord();
				if ( (name2 = getName()).empty() || endOfStream()) {
					LOG_DEBUG(3, "Found endofstream 1 record in leaving at eof to preserve pair " << name1);
					return false;
				}

				string name3;
				readRecord();
				if ( (name3 = getName()).empty() || endOfStream()) {
					seekg(here1);
					LOG_DEBUG(3, "Found endofstream two records in, rewinding to initial boundary " << here1 << " : " << name1 << " & " << name2);
					return true;
				}

				if (SequenceRecordParser::isPair(name1,name2)) {
					seekg(here1);
					LOG_DEBUG(3, "Found natural pair at boundary, rewinding to " << here1);
				} else if (SequenceRecordParser::isPair(name2, name3)) {
					seekg(here2);
					LOG_DEBUG(3, "Found split pair at boundary, incrementing one record to " << here2);
				} else {
					LOG_DEBUG(3, "Found no pairs at boundary, rewinding to " << here1);
					seekg(here1);
				}
				return true;
			} else {
				return true;
			}
		}

	};

	class FastqStreamParser: public SequenceStreamParser {

	public:
		FastqStreamParser(istream &s) :
			SequenceStreamParser(s, '@') {
		}
		FastqStreamParser(ReadFileReader::MmapSource &mmap) :
			SequenceStreamParser(mmap, '@') {

		}
		RecordPtr readRecord() {

			int threadNum = omp_get_thread_num();
			if (readName().empty()) {
				resetBuffers();
				LOG_DEBUG(5, "FastqStreamParser::readRecord() found the end");
				return NULL;
			} else {
				readBases();
				readQuals();
				LOG_DEBUG(5, "FastqStreamParser::readRecord() found " << _nameBuffer[threadNum]);
			}

			return getStreamRecordPtr();
		}
		RecordPtr readRecord(RecordPtr _recordPtr) const {
			int threadNum = omp_get_thread_num();
			RecordPtr recordPtr = _recordPtr;
			if (isPastPartition()) {
				resetBuffers();
				return NULL;
			}
			if (*(recordPtr++) != _marker) {
				// skip the first character
				LOG_THROW("Detected invalid line marker: at " << getPos() << ": " << *recordPtr);
			}
			std::string &name = _nameBuffer[threadNum];
			nextLine(name, recordPtr);
			bool isGood = SequenceRecordParser::trimName( name );  // name
			setIsGood(isGood);

			nextLine(_basesBuffer[threadNum], recordPtr); // fasta
			nextLine(_lineBuffer[threadNum], recordPtr);  // qual name
			nextLine(_qualsBuffer[threadNum], recordPtr); // quals
			setPos( getPos() + (recordPtr - _recordPtr));
			return recordPtr;
		}
		string &readBases() {
			int threadNum = omp_get_thread_num();
			nextLine(_basesBuffer[threadNum]);

			if (_basesBuffer[threadNum].empty() || _basesBuffer[threadNum].length()
					== _basesBuffer[threadNum].max_size())
				LOG_THROW("FastqStreamParser::readBases(): Missing or too many bases");

			string &qualName = nextLine();
			if (qualName.empty() || qualName[0] != '+')
				LOG_THROW("FastqStreamParser::readBases(): Missing '+' in fastq, got: " <<  _basesBuffer[threadNum] << " " << qualName);

			_isMultiline[threadNum] = false;
			return _basesBuffer[threadNum];
		}

		string &readQuals() {
			int threadNum = omp_get_thread_num();
			nextLine(_qualsBuffer[threadNum]);
			if (Options::getOptions().getIgnoreQual()) {
				_qualsBuffer[threadNum].assign(_qualsBuffer[threadNum].length(), Read::REF_QUAL);
			}
			return _qualsBuffer[threadNum];
		}
		int getType() const {
			return 0;
		}
		RecordPtr getStreamQualRecordPtr() const { return NULL; }
		RecordPtr getLastQualRecordPtr() const { return NULL; }
		MmapSource getQualMmap() const { return MmapSource(); }
	};

	class FastaStreamParser: public SequenceStreamParser {
		static const char fasta_marker = '>';

	public:

		FastaStreamParser(istream &s) :
			SequenceStreamParser(s, fasta_marker) {
		}
		FastaStreamParser(ReadFileReader::MmapSource &mmap) :
			SequenceStreamParser(mmap, fasta_marker) {
		}

		RecordPtr readRecord() {
			int threadNum = omp_get_thread_num();
			if (readName().empty()) {
				resetBuffers();
				LOG_DEBUG(5, "FastaStreamParser::readRecord() found the end at " << tellg());
				return NULL;
			} else {
				readBases();
				LOG_DEBUG(5, "FastaStreamParser::readRecord() found " << _nameBuffer[threadNum]);
			}
			return getStreamRecordPtr();
		}
		RecordPtr readRecord(RecordPtr _recordPtr) const {
			int threadNum = omp_get_thread_num();
			if (isPastPartition()) {
				resetBuffers();
				return NULL;
			}
			RecordPtr recordPtr = _recordPtr;
			if (*recordPtr != _marker) {
				LOG_THROW("FastaStreamParser::readRecord(): Could not FastaStreamParser::readRecord()");
			}
			std::string &name = _nameBuffer[threadNum];
			nextLine(name, recordPtr);  // name
			bool isGood = SequenceRecordParser::trimName(name);
			setIsGood(isGood);

			_basesBuffer[threadNum].clear();
			_isMultiline[threadNum] = false;
			long count = 0;
			while (recordPtr < _lastPtr && *recordPtr != fasta_marker) {
				_basesBuffer[threadNum] += nextLine(_lineBuffer[threadNum], recordPtr);
				if (++count > 1)
					_isMultiline[threadNum] = true;
			}
			_qualsBuffer[threadNum].assign(_basesBuffer[threadNum].length(), Read::REF_QUAL);
			setPos(getPos() + (recordPtr - _recordPtr));

			return recordPtr;
		}

		string &readBases() {
			int threadNum = omp_get_thread_num();
			string &bases = getBasesOrQuals();
			_qualsBuffer[threadNum].assign(bases.length(), Read::REF_QUAL);
			return bases;
		}


		virtual string &getBasesOrQuals() {
			int threadNum = omp_get_thread_num();
			_basesBuffer[threadNum].clear();
			_isMultiline[threadNum] = false;
			long count = 0;
			while (!endOfStream()) {
				nextLine();

				if (_lineBuffer[threadNum].empty())
					break;

				_basesBuffer[threadNum] += _lineBuffer[threadNum];
				if (++count > 1)
					_isMultiline[threadNum] = true;
				if (_basesBuffer[threadNum].size() == _basesBuffer[threadNum].max_size())
					LOG_THROW("FastaStreamParser::getBasesOrQuals(): Sequence/Qual too large to read");

				if (peek() == fasta_marker)
					break;
			}
			return _basesBuffer[threadNum];
		}
		int getType() const {
			return 1;
		}
		RecordPtr getStreamQualRecordPtr() const { return NULL; }
		RecordPtr getLastQualRecordPtr() const { return NULL; }
		MmapSource getQualMmap() const { return MmapSource(); }
	};

	class FastaQualStreamParser: public FastaStreamParser {
		FastaStreamParser _qualParser; // odd, but it works
		boost::unordered_map<RecordPtr, RecordPtr> translate;

	public:
		FastaQualStreamParser(istream &fastaStream, istream &qualStream) :
			FastaStreamParser(fastaStream), _qualParser(qualStream) {
		}
		FastaQualStreamParser(ReadFileReader::MmapSource &mmap, ReadFileReader::MmapSource &qualMmap) :
			FastaStreamParser(mmap), _qualParser(qualMmap) {
		}

		RecordPtr readRecord() {
			int threadNum = omp_get_thread_num();
			RecordPtr fastaPtr = getStreamRecordPtr();
			RecordPtr qualPtr = _qualParser.getStreamRecordPtr();
			translate[fastaPtr] = qualPtr;
			if (readName().empty()) {
				resetBuffers();
				LOG_DEBUG(5, "FastaQualStreamParser::readRecord() found the end");
				return NULL;
			} else {
				readBases();
				readQuals();
				_isMultiline[threadNum] = _isMultiline[threadNum] || _qualParser.isMultiline();
				LOG_DEBUG(5, "FastaQualStreamParser::readRecord() found " << _nameBuffer[threadNum]);
			}
			return getStreamRecordPtr();
		}
		RecordPtr readRecord(RecordPtr _recordPtr) const {
			int threadNum = omp_get_thread_num();
			if (isPastPartition()) {
				resetBuffers();
				return NULL;
			}
			RecordPtr recordPtr = _recordPtr;
			RecordPtr qualPtr = translate.find(recordPtr)->second;
			FastaStreamParser::readRecord(recordPtr);
			_qualParser.readRecord( qualPtr );
			_qualsBuffer[threadNum] = _qualParser.getQuals();
			_isMultiline[threadNum] = _isMultiline[threadNum] || _qualParser.isMultiline();
			setPos(getPos() + (recordPtr - _recordPtr));
			return recordPtr;
		}

		string &readName() {
			string &qualName = _qualParser.readName();
			string &name = FastaStreamParser::readName();
			if (name != qualName)
				LOG_THROW("FastaQualStreamParser::readName(): fasta and quals have different names");
			return name;
		}

		string &readBases() {
			return getBasesOrQuals();
		}

		string &readQuals() {
			int threadNum = omp_get_thread_num();
			// odd, but it works
			string &qualInts = _qualParser.getBasesOrQuals();
			_qualsBuffer[threadNum] = SequenceRecordParser::convertQualIntsToChars(qualInts, Read::FASTQ_START_CHAR);
			return _qualsBuffer[threadNum];
		}
		int getType() const {
			return 1;
		}
		RecordPtr getStreamQualRecordPtr() const { return _qualParser.getStreamRecordPtr(); }
		RecordPtr getLastQualRecordPtr() const { return _qualParser.getLastRecordPtr(); }
		MmapSource getQualMmap() const { return _qualParser.getMmap(); }
	};


};

#endif

