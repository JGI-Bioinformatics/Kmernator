//
// Kmernator/src/Utils.h
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

#ifndef _UTILS_H
#define _UTILS_H

#include <ostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <set>
#include <cstring>
#include <cmath>
#include <memory>
#include <csignal>
#include <sys/wait.h>

#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>
#include <dirent.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <cstdlib>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>


#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

#define foreach BOOST_FOREACH

#include "config.h"
#include "Options.h"
#include "Log.h"
#include "lookup3.h"


class FormatOutput
{
public:
	enum FormatType {FASTQ, FASTA, FASTQ_UNMASKED, FASTA_UNMASKED};

	FormatOutput(const FormatOutput &copy) : _type(copy._type) {}
	FormatOutput &operator=(const FormatOutput copy) {
		_type = copy._type;
		return *this;
	}
	static const FormatOutput Fastq() { static FormatOutput singleton = FormatOutput(FASTQ); return singleton; }
	static const FormatOutput Fasta() { static FormatOutput singleton = FormatOutput(FASTA); return singleton; }
	static const FormatOutput FastqUnmasked() { static FormatOutput singleton = FormatOutput(FASTQ_UNMASKED); return singleton; }
	static const FormatOutput FastaUnmasked() { static FormatOutput singleton = FormatOutput(FASTA_UNMASKED); return singleton; }

	static inline FormatOutput getDefault() {
		return FormatOutput(Options::getOptions().getFormatOutput());
	}
	static std::string getDefaultSuffix() {
		return getDefault().getSuffix();
	}

	inline FormatType getType() const {
		return _type;
	}
	bool inline operator==(const FormatOutput &other) const {
		return _type == other._type;
	}
	inline std::string getSuffix() const {
		return getSuffix(*this);
	}
	static inline std::string getSuffix(const FormatOutput &format) {
		return getSuffix(format._type);
	}
	static inline std::string getSuffix(int type) {
		std::string ret;
		switch(type) {
		case 0:
		case 2: ret = std::string(".fastq"); break;
		case 1:
		case 3: ret = std::string(".fasta"); break;
		default : ret = std::string(".txt");
		}
		return ret;
	}

	std::string toString(std::string &name, std::string &fasta, std::string &qual) {
		std::stringstream ss;
		switch(_type) {
		case 0:
		case 2: ss << "@" << name << "\n" << fasta << "\n+\n" << qual << "\n"; break;
		case 1:
		case 3: ss << ">" << name << "\n" << fasta << "\n"; break;
		}
		return ss.str();
	}
	std::string toString(std::string &name, std::string &fasta, std::string &qual, std::string &comment) {
		std::stringstream ss;
		std::string nameAndComment = name;
		if (!comment.empty())
			nameAndComment += " " + comment;
		switch(_type) {
		case 0:
		case 2: ss << "@" << nameAndComment << "\n" << fasta << "\n+\n" << qual << "\n"; break;
		case 1:
		case 3: ss << ">" << nameAndComment << "\n" << fasta << "\n"; break;
		}
		return ss.str();
	}

private:
	FormatType _type;
	FormatOutput(int type) {
		switch(type) {
		case(0) : _type = FASTQ; break;
		case(1) : _type = FASTA; break;
		case(2) : _type = FASTQ_UNMASKED; break;
		case(3) : _type = FASTA_UNMASKED; break;
		default: LOG_THROW("Invalid format type: " << type);
		}
	}

};

class UniqueName {
	static int getUnique() {
		static int id = 0;
#pragma omp atomic
		id++;
		return id;
	}
	static int getGlobalUnique() {
		static int id = 0;
#pragma omp master
		{
			id++;
		}
#pragma omp barrier
		return id;
	}

public:
	typedef OptionsBaseInterface::StringListType StringListType;
	static std::string generateUniqueGlobalName(std::string filename = "") {
		filename += "-" + boost::lexical_cast<std::string>( getGlobalUnique() );
		return filename;
	}
	static std::string generateUniqueName(std::string filename = "") {
		filename += "-" + boost::lexical_cast<std::string>( getpid() );
		filename += "-" + boost::lexical_cast<std::string>( getUnique() );
		filename += "-" + boost::lexical_cast<std::string>( omp_get_thread_num() );
		filename += OptionsBaseInterface::getHostname();
		return filename;
	}
	static std::string generateHashName(std::string filename) {
		uint64_t hash = 0xDEADBEEF;
		uint32_t *pc, *pb;
		pc = (uint32_t*) &hash;
		pb = pc+1;
		Lookup3::hashlittle2(filename.c_str(), filename.size(), pc, pb);
		std::stringstream ss;
		ss << std::hex;
		ss << hash;
		return ss.str();
	}
	static std::string generateHashName(StringListType fileNames) {
		std::string hash;
		for(StringListType::iterator it = fileNames.begin(); it != fileNames.end(); it++)
			hash += *it + "\n";
		return generateHashName(hash);
	}

};

class OfstreamMap {
public:
	typedef std::set< std::string > KeySet;
	class OStreamPtr {
	private:
		boost::shared_ptr< std::ofstream > of;
		boost::shared_ptr< std::stringstream> ss;
		std::string filePath;
	public:
		OStreamPtr() {
			assert(empty());
		}
		OStreamPtr(std::string _filePath, bool append) : filePath(_filePath) {
			std::ios_base::openmode mode = std::ios_base::out;
			if (append)
				mode |= std::ios_base::app;
			else
				mode |= std::ios_base::trunc;

			LOG_DEBUG_OPTIONAL(2, true, "OfstreamMap::OStreamPtr(): Writing to " << filePath);
			of.reset(new std::ofstream(filePath.c_str(), mode));
			assert(isFileStream());
			if (of->fail() || !of->is_open() || !of->good())
				LOG_THROW("Could not open " << filePath << " for writing! " << of->rdstate());
		}
		OStreamPtr(std::string _filePath) : filePath(_filePath) {
			LOG_DEBUG_OPTIONAL(2, true, "OfstreamMap::OStreamPtr(): In-memory writing to " << filePath);
			ss.reset(new std::stringstream());
			assert(isStringStream());
		}
		void close() {
			if (isFileStream()) {
				LOG_DEBUG_OPTIONAL(2, true, "OfstreamMap::OStreamPtr::close(): Closing " << getFilePath());
				of->flush();
				of->close();
				if (of->fail())
					LOG_WARN(1, "Could not properly close: " << getFilePath());
			}
			if (isStringStream()) {
				LOG_DEBUG_OPTIONAL(2, true, "OfstreamMap::OStreamPtr::close(): Writing out in-memory : " << getFilePath());
				OStreamPtr osp(getFilePath(), false);
				assert(osp.isFileStream());
				*osp << *ss;
				osp.close();
			}
			reset();
		}
		void reset() {
			of.reset();
			ss.reset();
			filePath.clear();
			assert(empty());
		}
		~OStreamPtr() {
			reset(); // do not close(), as this violates behavior when in a container
		}

		bool isFileStream() const {
			return of.get() != NULL && ss.get() == NULL;
		}
		bool isStringStream() const {
			return of.get() == NULL && ss.get() != NULL;
		}
		bool empty() const {
			return !(isFileStream() | isStringStream());
		}
		std::ostream &operator*() {
			if (isFileStream())
				return *of;
			else if (isStringStream())
				return *ss;
			else
				LOG_THROW("Invalid state for OfstreamMap - not file, not memory!");
		}
		std::string getFilePath() const {
			return filePath;
		}
		std::string getFinalString() {
			assert(isStringStream());
			int64_t bytes = ss->tellp();
			LOG_DEBUG(3, "OfstreamMap::OStreamPtr::getFinalString(): Writing out " << bytes << " bytes in-memory for virtual file: " << getFilePath());
			std::string s = ss->str();
			reset();
			return s;
		}
	};
	typedef boost::unordered_map< std::string, OStreamPtr > Map;
	typedef Map::iterator Iterator;
	typedef boost::shared_ptr< Map > MapPtr;
protected:
	MapPtr _map;
	std::string _outputFilePathPrefix;
	std::string _suffix;
	bool _append;
	bool _isStdout;
	bool _buildInMemory;
#ifdef _USE_MPI
	mpi::communicator *_world;
#endif
	virtual void close() {
		LOG_DEBUG_OPTIONAL(2, true, "Calling OfstreamMap::close()");
		for(Iterator it = _map->begin() ; it != _map->end(); it++) {
			OStreamPtr &_osp = it->second;
			_osp.close();
		}
	}

public:
	static bool &getDefaultAppend() {
		static bool _defaultAppend = false;
		return _defaultAppend;
	}

	OfstreamMap(std::string outputFilePathPrefix = Options::getOptions().getOutputFile(), std::string suffix = FormatOutput::getDefaultSuffix())
	: _map(new Map()), _outputFilePathPrefix(outputFilePathPrefix), _suffix(suffix), _append(false), _isStdout(false), _buildInMemory(false) {
		_append = getDefaultAppend();
		if (Options::getOptions().getOutputFile() == std::string("-")) {
			_isStdout = true;
			LOG_VERBOSE_OPTIONAL(1, true, "Writing output(s) to stdout");
		}
#ifdef _USE_MPI
		_world = NULL;
#endif
	}
	virtual ~OfstreamMap() {
		LOG_DEBUG_OPTIONAL(2, true, "~OfstreamMap():");
		this->clear();
	}
	std::string strip(std::string filename, std::string suffix) {
		if (!suffix.empty()) {
			filename.substr(0, filename.find(suffix));
		}
		return filename;
	}
	KeySet getKeySet() {
		KeySet keys;
		for(Iterator it = _map->begin() ; it != _map->end(); it++) {
			std::string key = it->first;
			keys.insert(key);
		}
		return keys;
	}
	void setBuildInMemory(bool buildInMemory = true) {
		_buildInMemory = buildInMemory;
	}
	bool isBuildInMemory() const {
		return _buildInMemory;
	}
	bool &getAppend() {
		return _append;
	}
	const bool &getAppend() const {
		return _append;
	}
	virtual void clear() {
		LOG_DEBUG_OPTIONAL(2, true, "Calling OfstreamMap::clear() " << _map->size() << " " << _outputFilePathPrefix);
		this->close();
		_clear();
	}
	void _clear() {
		_map->clear();
	}
	virtual std::string getRank() const {
		return std::string();
	}
	std::string getOutputPrefix() const {
		return _outputFilePathPrefix;
	}
	std::string getSuffix() const {
		return _suffix;
	}
	std::string getFilename(std::string key) const {
		return key + getSuffix() + getRank();
	}
	std::string getFilePath(std::string key) const {
		return getOutputPrefix() + getFilename(key);
	}
	virtual std::string getRealFilePath(std::string key) const {
		return getFilePath(key);
	}
	std::ostream &getOfstream(std::string key) {
		if (_isStdout)
			return std::cout;
		// lockless lookup
		MapPtr thisMap = _map;
		Iterator it = thisMap->find(key);
		if (it == thisMap->end()) {

			// lock if not found and map needs to be updated
#ifdef _USE_OPENMP
#pragma omp critical (ofStreamMap)
#endif
			{
				// re-copy the map first
				thisMap = _map;

				// recheck map
				it = thisMap->find(key);
				if (it == thisMap->end()) {
					OStreamPtr osp;
					if (_buildInMemory)
						osp = OStreamPtr(getFilePath(key));
					else
						osp = OStreamPtr(getFilePath(key), getAppend());

					MapPtr copy = MapPtr(new Map(*thisMap));
					it = copy->insert( copy->end(), Map::value_type(key, osp) );
					_map = thisMap = copy;
				}
			}
		}
		return *(it->second);
	}

};

template<typename S>
class PartitioningData {
public:
	typedef S DataType;
	typedef std::vector< DataType > Partitions;

private:
	Partitions _partitions;

public:
	PartitioningData() : _partitions() {}

	void swap(PartitioningData &other) {
		_partitions.swap(other._partitions);
	}
	void clear() {
		_partitions.clear();
	}
	inline bool hasPartitions() const {
		return ! _partitions.empty();
	}
	inline int getPartitionIdx(DataType score) const {
		// TODO binary search?
		for(unsigned int i = 0; i < _partitions.size(); i++)
			if (score < _partitions[i])
				return i;
		return _partitions.size();
	}
	inline Partitions getPartitions() const {
		return _partitions;
	}

	inline int addPartition(DataType partition) {
		_partitions.push_back(partition);
		std::sort(_partitions.begin(), _partitions.end());
		return _partitions.size();
	}


};

template<typename Raw, typename Store>
class BucketedData {
public:
	typedef Store StoreType;
	typedef Raw RawType;

private:
	RawType _minValue, _maxValue;
	StoreType steps;

public:
	// TODO inclusive permutations: (), [), (], []
	//      presently it truncates fraction, so [)
	// TODO implement log scale
	BucketedData(RawType minValue, RawType maxValue) :
		_minValue(minValue), _maxValue(maxValue) {
		// assumes unsigned Store...
		steps = ((StoreType) -1);
	}
	~BucketedData() {
	}

	StoreType getStore(RawType value) {
		if (value < _minValue || value > _maxValue)
			LOG_THROW( "InvalidArgument: getStore() out of range" );
		return ((RawType) steps) * (value - _minValue)
				/ (_maxValue - _minValue);
	}

	RawType getValue(StoreType store) {
		return (_maxValue - _minValue) * ((RawType) store) / ((RawType) steps);
	}
};


typedef BucketedData<double, Kmernator::UI8> DoubleToOneByte;
typedef BucketedData<double, Kmernator::UI16> DoubleToTwoByte;

class SequenceRecordParser
{
public:
	static inline std::string &nextLine(std::string &buffer, Kmernator::RecordPtr &recordPtr) {
		Kmernator::RecordPtr nextPtr = strchr(recordPtr, '\n');
		long len = nextPtr - recordPtr;
		if (len > 0) {
			buffer.assign(recordPtr, len);
		} else {
			buffer.clear();
		}
		recordPtr = ++nextPtr;
		return buffer;
	}

	// returns true if the read is to be kept, based on the name and comment
	static bool trimName(std::string &nameLine, std::string &comment) {

		if (nameLine.length() == 0)
			return false;

		bool isGood = true;

		char &marker = nameLine[0];
		if (marker != '>' && marker != '@') {
			LOG_THROW( "Can not parse name without a marker: " <<  nameLine );
		} else {
			// remove marker
			nameLine.erase(0,1);
		}

		// trim at first whitespace
		size_t pos = nameLine.find_first_of(" \t\r\n");
		if (pos != std::string::npos) {
			if (nameLine.length() >= pos + 2) {
				comment = nameLine.substr(pos+1);
				// detect Illumina 1.8 naming scheme in comment '[12]:[YN]:.*'
				if (isCommentCasava18(comment) && (pos <= 2 || nameLine[pos-2] != '/')) {
					if (!GlobalOptions::isCommentStored()) {
						// rewrite "name [12]:Y:..." to: "name/[12]" to detect pairs in old format, as comment is not preserved this run...
						nameLine[pos] = '/';
						pos += 2;
					}
					if (nameLine[pos+3] == 'Y')
						isGood = false;
				}
			} else
				comment.clear();
			nameLine.erase(pos);
		} else {
			comment.clear();
		}
		return isGood;
	}

	// returns true if the read is to be kept, based on the name and comment
	static bool parse(Kmernator::RecordPtr record, Kmernator::RecordPtr lastRecord,
			std::string &name, std::string &bases, std::string &quals, std::string &comment,
			Kmernator::RecordPtr qualRecord = NULL, Kmernator::RecordPtr lastQualRecord = NULL,
			char fastqStartChar = Kmernator::FASTQ_START_CHAR_ILLUMINA) {
		std::string buf;
		bool isGood = true;
		if (*record == '@') {
			// FASTQ
			nextLine(name,  record);
			isGood = trimName(name, comment);
			nextLine(bases, record);
			nextLine(buf,   record);
			nextLine(quals, record);
			if (buf[0] != '+') {
				LOG_THROW( "Invalid FASTQ record!" );
			}
		} else if (*record == '>') {
			// FASTA
			nextLine(name, record);
			isGood = trimName(name, comment);

			// TODO FIXME HACK!
			if (lastRecord == NULL) // only read one line if no last record is given
				lastRecord = record+1;

			while (record < lastRecord && *record != '>') {
				bases += nextLine(buf, record);
			}
			if (qualRecord != NULL) {
				nextLine(buf, qualRecord);
				trimName(buf, comment);
				if (buf != name) {
					LOG_THROW( "fasta and qual do not match names!" );
				}

				// TODO FIXME HACK!
				if (lastQualRecord == NULL) // only read one line if no last record is given
					lastQualRecord = qualRecord+1;
				while (qualRecord < lastQualRecord && *qualRecord != '>') {
					quals += nextLine(buf, qualRecord);
				}
				quals = convertQualIntsToChars(quals, fastqStartChar);
			} else {
				quals.assign(bases.length(), Kmernator::REF_QUAL);
			}
		} else {
			LOG_THROW( "Do not know how to parse this file!" );
		}
		return isGood;
	}

	static std::string convertQualIntsToChars(const std::string &qualInts, char fastqStartChar) {
		std::istringstream ss(qualInts);
		std::ostringstream oss;
		int maxQual = Kmernator::REF_QUAL - fastqStartChar - 1;
		while (!ss.eof()) {
			int qVal;
			ss >> qVal;
			if (ss.fail())
				break;
			if (qVal > maxQual)
				qVal = maxQual;
			qVal += fastqStartChar;
			oss << (char) qVal;
		}
		return oss.str();
	}

	static std::string commonName(const std::string readName) {
		if (readName.length() <= 2)
			return readName;
		else if (readName[readName.length()-2] == '/')
			return readName.substr(0, readName.length() - 1);
		else
			return readName;
	}

	static bool isCommentCasava18(const std::string comment) {
		if (comment.empty() || comment.length() < 6)
			return false;
		else if ( (':' == comment[1] && ':' == comment[3] && ':' == comment[5]) && (comment[0] == '1' || comment[0] == '2') && (comment[2] == 'Y' || comment[2] == 'N'))
			return true;
		else
			return false;
	}

	// returns 0 for unpaired reads, 1 or 2 if the read is paired.
	// supports the following three patterns: xxx/1 xxx/2, xxx/A xxx/B, or xxx/F xxx/R
	static int readNum(const std::string readName, const std::string comment) {
		if (isCommentCasava18(comment)) {
			return comment[0] == '2' ? 2 : 1;
		}
		int retVal = 0;
		int len = readName.length();
		if (len < 2)
			return retVal;
		if (readName[len - 2] != '/')
			return retVal;
		char c = readName[len - 1];
		switch (c) {
		case '1':
		case 'A':
		case 'F':
			retVal = 1;
			break;
		case '2':
		case 'B':
		case 'R':
			retVal = 2;
			break;
		}
		return retVal;
	}

	static bool isPairedRead(const std::string readName, const std::string comment) {
		return readNum(readName, comment) != 0;
	}

	static bool isPair(const std::string readNameA, const std::string readNameB, const std::string commentA, const std::string commentB) {
		std::string commonA = commonName(readNameA);
		std::string commonB = commonName(readNameB);
		if (commonA.compare(commonB) == 0) {
			int readNumA = readNum(readNameA, commentA);
			int readNumB = readNum(readNameB, commentB);
			if (readNumA != 0 && readNumB != 0 && readNumA != readNumB) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
};

// does not work across pipe!
class StdErrInterceptor2 {
public:
	StdErrInterceptor2(bool _intercept = false)
	: _stderrBuffer(), _olderr(_intercept ? std::cerr.flush().rdbuf(_stderrBuffer.rdbuf()) : NULL)
	{
	}
	virtual ~StdErrInterceptor2() {
		reset();
	}
	inline bool isIntercepted() const {
		return _olderr != NULL;
	}
	void reset() {
		if (isIntercepted()) {
			std::cerr.rdbuf( _olderr );
		}
		_olderr = NULL;
	}
	std::string getStdErr () {
		return _stderrBuffer.str();
	}
private:
	std::stringstream _stderrBuffer;
	std::streambuf *_olderr;
};

// works across a pipe!
class StdErrInterceptor {
public:
	StdErrInterceptor(bool _intercept = false)
	: _errFileName( _intercept ? getNewFile() : "")
	{
	}
	virtual ~StdErrInterceptor() {
		reset();
	}
	inline bool isIntercepted() const {
		return !_errFileName.empty();
	}
	void setErrorFile(std::string errFileName) {
		_errFileName = errFileName;
	}
	void reset() {
		if (isIntercepted())
			unlink(_errFileName.c_str());
		_errFileName.clear();
	}
	std::string getErrIntercept() const {
		if (isIntercepted())
			return " 2> " + _errFileName + " ";
		else
			return "";
	}
	std::string getStdErr () {
		std::stringstream ss;
		std::fstream strace(_errFileName.c_str(), std::ios_base::in);
		std::string buffer;
		while (strace.good()) {
			getline(strace, buffer);
			ss << buffer << std::endl;
		}
		reset();
		return ss.str();
	}
	static std::string getNewFile() {
		return Options::getOptions().getTmpDir() + UniqueName::generateUniqueName("/.tmp-stderr-");
	}
private:
	std::string _errFileName;
};


class OPipestream : public StdErrInterceptor, public boost::iostreams::stream< boost::iostreams::file_descriptor_sink >
{
public:
	typedef boost::iostreams::stream<  boost::iostreams::file_descriptor_sink > base ;
	explicit OPipestream() : _pipe(NULL), _exitStatus(0) {}
	explicit OPipestream( const std::string command, bool interceptStderr = false)
	: StdErrInterceptor(interceptStderr),
	  base( fileno( _pipe = popen( (command+getErrIntercept()).c_str(), "w" ) ), boost::iostreams::never_close_handle), _cmd(command+getErrIntercept()), _exitStatus(0) {
		assert(_pipe != NULL);
		assert(fileno(_pipe) >= 0);
		assert(is_open());
	}
	void close() {
		if (_pipe == NULL)
			return;
		try {
			this->flush();
			base::close();
			_exitStatus = pclose(_pipe);
			_pipe = NULL;
			if (_exitStatus != 0) {
				LOG_DEBUG(3, "OPipestream::close() '" << _cmd << "' closed with an error: " << _exitStatus << "." << getStdErr());
			}
		} catch(...) {
			// ignoring this pipe closure error.
			LOG_WARN(1, "OPipestream::close(): Potentially failed to close pipe properly." << getStdErr() );
		}
	}
	int getExitStatus() { return _exitStatus; }
	virtual ~OPipestream() {
		close();
	}
private :
	FILE* _pipe ;
	std::string _cmd;
	int _exitStatus;
};

class IPipestream : public StdErrInterceptor, public boost::iostreams::stream< boost::iostreams::file_descriptor_source >
{
public:
	typedef boost::iostreams::stream<  boost::iostreams::file_descriptor_source > base ;
	explicit IPipestream() : StdErrInterceptor(), _pipe(NULL), _exitStatus(0) {}
	explicit IPipestream( const std::string command, bool interceptStderr = false )
	: StdErrInterceptor( interceptStderr ),
	  base( fileno( _pipe = popen( (command+getErrIntercept()).c_str(), "r" ) ), boost::iostreams::never_close_handle ),
	  _cmd(command+getErrIntercept()) , _exitStatus(0) {

		assert(_pipe != NULL);
		assert(fileno(_pipe) >= 0);
		assert(is_open());

	}
	void close() {
		if (_pipe == NULL)
			return;
		this->set_auto_close(true);
		_exitStatus = 0;
		try {
			_exitStatus = pclose(_pipe);
			_pipe = NULL;
			if (_exitStatus != 0) {
				LOG_WARN(1, "IPipestream::close() '" << _cmd << "' closed with an error: " << _exitStatus << "." << getStdErr());
			}
		} catch(...) {
			// ignoring this pipe closure error.
			LOG_WARN(1, "IPipestream::close(): Potentially failed to close pipe properly." << getStdErr());

		}
	}
	int getExitStatus() { return _exitStatus; }
	virtual ~IPipestream() {
		close();
	}
private :
	FILE* _pipe ;
	std::string _cmd;
	int _exitStatus;
};


class FileUtils
{
public:

	static std::string getDirname(const std::string filePath) {
		char buf[filePath.size() + 1];
		memcpy(buf, filePath.c_str(), filePath.size());
		std::string dirPath(dirname(buf));
		return dirPath;
	}
	static void syncDir(const std::string filePath) {
		std::string dirPath = getDirname(filePath);
		DIR *d = opendir(dirPath.c_str());
		if (d == NULL) {
			LOG_ERROR(1, "Could not opendir " << dirPath << "! " << strerror(errno));
		} else {
			if (fsync(dirfd(d)) != 0)
				LOG_ERROR(1, "Could not dirsync " << dirPath << "! " << strerror(errno));
			closedir(d);
		}
	}
	static void syncFile(const std::string filePath) {
		int fd = open(filePath.c_str(), std::ios_base::in | std::ios_base::out);
		if (fd < 0) {
			LOG_ERROR(1, "Could not open " << filePath << "! " << strerror(errno));
		} else {
			if (fsync(fd) != 0)
				LOG_ERROR(1, "Could not fsync " << filePath << "! " << strerror(errno));
		}
		if (close(fd) != 0)
			LOG_ERROR(1, "Could not close " << filePath << "! " << strerror(errno));
	}
	static std::auto_ptr<struct stat> statFile(const std::string filePath, bool dofsync = false, bool dodirsync = false) {
		if (dodirsync) {
			syncDir(filePath);
		}
		if (dofsync) {
			syncFile(filePath);
		}
		std::auto_ptr<struct stat> fileStat(new struct stat);
		if (stat(filePath.c_str(), fileStat.get()) != 0) {
			LOG_ERROR(1, "Could not stat " << filePath << "! " << strerror(errno));
		} else {
			LOG_DEBUG_OPTIONAL(2, true, "stat of " << filePath << ": dev " << fileStat->st_dev << ", ino " << fileStat->st_ino << ", mode " << fileStat->st_mode << ", nlink " << fileStat->st_nlink << ", size " << fileStat->st_size)
		}
		return fileStat;
	}
	static unsigned long getFileSize(const std::string filePath) {
		std::ifstream ifs(filePath.c_str());
		if (ifs.good())
			return getFileSize(ifs);
		else
			return 0;
	}
	static unsigned long getFileSize(std::istream &is) {
		assert( !is.eof() );
		assert( is.good() );
		assert( !is.fail() );
		std::ifstream::streampos current = is.tellg();
		is.seekg(0, std::ios_base::end);
		unsigned long size = is.tellg();
		is.seekg(current);
		return size;
	}
	static bool fileExists(std::string filePath) {
		std::ifstream ifs(filePath.c_str());
		return fileExists(ifs);
	}
	static bool fileExists(std::ifstream &ifs) {
		if (ifs.fail())
			return false;
		else
			return true;
	}
	static std::string dumpFile(std::string filePath) {
		std::stringstream ss;
		std::fstream f(filePath.c_str(), std::ios_base::in);
		std::string buf;
		while (f.good()) {
			f >> buf;
			ss << buf;
		}
		return ss.str();
	}
	static std::string getBasePath(std::string path) {
		std::string basePath;
		if (path[0] == '/') {
			basePath = path.substr(0, path.find_last_of('/'));
		} else {
			std::string PATH = getenv("PATH");
			std::string name = path;
			size_t pos = path.find_last_of('/');
			if (pos != std::string::npos)
				name = path.substr( pos + 1);
			boost::char_separator<char> sep(":");
			boost::tokenizer< boost::char_separator<char> > tok(PATH, sep);
			for(boost::tokenizer< boost::char_separator<char> >::iterator it = tok.begin(); it != tok.end() ; it++) {
				if (fileExists( *it + "/" + name )) {
					basePath = *it;
					break;
				}
			}
		}
		return basePath;
	}
};

class Statistics
{
public:
	class MeanStdCount
	{
	public:
		double mean;
		double stdDev;
		long count;
		MeanStdCount() : mean(0.0), stdDev(0.0), count(0) {}

		template<typename forwardIterator>
		MeanStdCount(forwardIterator begin, forwardIterator end, bool isPoisson = true) : mean(0.0), stdDev(0.0), count(0) {
			for(forwardIterator it = begin; it != end ; it++) {
				count++;
				mean += *it;
			}
			if (count > 1) {
				mean /= (double) count;
				for(forwardIterator it = begin; it != end; it++) {
					double diff = ((double)*it) - mean;
					stdDev += diff*diff;
				}
				stdDev /= (double) (count-1);
			}
			if (isPoisson && mean > 0.0) {
				double poissonStdDev = sqrt(mean);
				if (stdDev < poissonStdDev) {
					stdDev = poissonStdDev;
				}
			}
		}
	};

	template<typename forwardIterator>
	static forwardIterator findBimodalPartition(double numSigmas, MeanStdCount &firstMSC, MeanStdCount &secondMSC, forwardIterator begin, forwardIterator end, bool isPoisson = true) {
		if (begin == end)
			return end;
		forwardIterator bestPart = end;
		double bestDiff = 0.0;
		forwardIterator part = begin;
		while (++part != end) {
			MeanStdCount first(begin, part, isPoisson);
			MeanStdCount second(part, end, isPoisson);
			if (first.count == 1 && second.count == 1)
				continue; // can not evaluate unless one sample is >1
			double diff = abs(first.mean - second.mean);
			// use the larger of the two standard deviations
			double stdDev = std::max(first.stdDev, second.stdDev);
			if (diff > numSigmas * stdDev) {	
				if (diff > bestDiff) {
					bestDiff = diff;
					bestPart = part;
					firstMSC = first;
					secondMSC = second;
				}
			}
		}
		return bestPart;
	};
	template<typename forwardIterator>
	static forwardIterator findBimodalPartition(double numSigmas, forwardIterator begin, forwardIterator end, bool isPoisson = true) {
		MeanStdCount f, s;
		return findBimodalPartition(numSigmas, f, s, begin, end, isPoisson);
	};
};

class LongRand {
public:
	// THREAD SAFE
	static unsigned long rand() {
		return rand( getSeed() );
	}
	static unsigned long rand(unsigned int &seed) {
		return
				(((unsigned long)(rand_r(&seed) & 0xFFFF)) << 48) |
				(((unsigned long)(rand_r(&seed) & 0xFFFF)) << 32) |
				(((unsigned long)(rand_r(&seed) & 0xFFFF)) << 16) |
				(((unsigned long)(rand_r(&seed) & 0xFFFF)) );
	}
	// NOT THREAD SAFE!
	static unsigned long rand2() {
		assert(omp_get_num_threads() == 1 || !omp_in_parallel());
		return (((unsigned long)(std::rand() & 0xFFFF)) << 48) |
				(((unsigned long)(std::rand() & 0xFFFF)) << 32) |
				(((unsigned long)(std::rand() & 0xFFFF)) << 16) |
				(((unsigned long)(std::rand() & 0xFFFF)) );

	}
	static void srand(unsigned int seed) {
		assert(!omp_in_parallel());
		for(int i = 0 ; i < omp_get_max_threads(); i++) {
			getSeed(i) = (seed+i) ^ i;
		}
	}
private:
	typedef std::vector< unsigned int > SeedVector;
	typedef boost::shared_ptr< SeedVector > SeedPtr;
	static unsigned int &getSeed(int threadId) {
		assert(threadId < omp_get_max_threads());
		static SeedPtr threadSeedsPtr;

		if (threadSeedsPtr.get() == NULL || (int) threadSeedsPtr->size() <= threadId) {
#pragma omp critical
			{
				unsigned int baseSeed = 0;
				if (threadSeedsPtr.get() == NULL || (int) threadSeedsPtr->size() <= threadId) {
					SeedPtr newSeeds;
					if (threadSeedsPtr.get() != NULL)
						newSeeds.reset( new SeedVector(*threadSeedsPtr) );
					else
						newSeeds.reset( new SeedVector );
					int numSeeds = std::max(threadId+1, omp_get_max_threads());
					newSeeds->reserve(numSeeds);

					while ((int) newSeeds->size() < numSeeds) {
						baseSeed += (time(NULL) + 2*time(NULL)) ^ newSeeds->size();
						newSeeds->push_back( rand_r(&baseSeed) );
					}
					threadSeedsPtr = newSeeds;
				}
			}
		}
		return (*threadSeedsPtr)[threadId];
	}
	static unsigned int &getSeed() {
		return getSeed(omp_get_thread_num());
	}

};

class Timer {
public:
	typedef std::pair< std::string, double > LabeledTime;
	typedef std::vector< LabeledTime > LabeledTimes;
	Timer() : _times() {
	}
	void resetTimes() {
		_times.clear();
	}
	void resetTimes(std::string label, double time) {
		resetTimes();
		recordTime(label,time);
	}
	void recordTime(std::string label, double time) {
		_times.push_back( LabeledTime(label, time));
	}
	std::string getTimes(std::string label) {
		std::stringstream ss;
		ss << label;
		double last = 0.0;
		double total = 0.0;
		if (!_times.empty())
			last = _times.begin()->second;
		ss << std::fixed << std::setprecision(2);
		for(LabeledTimes::iterator it = _times.begin(); it != _times.end(); it++) {
			double diff = it->second - last;
			total += diff;
			ss << " " << it->first << ": " << diff;
			last = it->second;
		}
		ss << " Total: " << total;
		return ss.str();
	}
private:
	LabeledTimes _times;
};

template<typename T>
class Random {
public:
	typedef typename boost::unordered_set<T> Set;
	typedef typename Set::iterator SetIterator;
	static Set sample(T existingSize, T sampleSize) {
		Set ids;
		if (existingSize > (3 * sampleSize / 2)) {
			// select ids to include
			while (ids.size() < sampleSize) {
				T includeCandidate = LongRand::rand() % existingSize;
				ids.insert(includeCandidate);
			}
		} else if (existingSize > sampleSize){
			// select ids to exclude
			T excludeCount = existingSize - sampleSize;
			Set exclude;
			while (exclude.size() < excludeCount) {
				T excludeCandidate = LongRand::rand() % existingSize;
				exclude.insert(excludeCandidate);
			}
			for(T i = 0; i < existingSize; i++) {
				if (exclude.find(i) == exclude.end())
					ids.insert(i);
			}
		} else {
			// really no need to sample!
			for(T i = 0; i < existingSize; i++)
				ids.insert(i);
		}
		assert( ids.size() == std::min(existingSize, sampleSize));
		return ids;
	}
};


class Cleanup {
public:
	typedef std::set< std::string > FileSet;
	typedef std::set< pid_t > PidSet;

	static Cleanup &getInstance() {
		static Cleanup cleanup;
		return cleanup;
	}
	static void cleanup(int param = 0) {
		if (param != 0)
			LOG_WARN(1, "Cleanup::cleanup(): Caught signal: " << param);
		LOG_DEBUG_OPTIONAL(4, true, "Entered Cleanup::cleanup()");
		getInstance()._clean(param);
	}
	~Cleanup() {
		// only allow parent process that first called Cleanup() to cleanup (if fork() was ever called)
		_unsetSignals();
		if (myPid == getpid()) {
			_clean(0);
		}
	}
	static std::string makeTempDir(std::string dir = "", std::string prefix = "") {
		if (dir.empty())
			dir = Options::getOptions().getTmpDir();
		if (prefix.empty())
			prefix = ".tmp-";
		std::string tempDir = dir + "/" + UniqueName::generateUniqueName(prefix);
		if (mkdir(tempDir.c_str(), 0700) != 0)
			LOG_THROW("Could not mkdir: " << tempDir);
		getInstance().tempDirs.insert(tempDir);
		return tempDir;
	}
	static void removeTempDir(std::string tempDir) {
		if (getInstance().tempDirs.find(tempDir) != getInstance().tempDirs.end()){
			getInstance().tempDirs.erase(tempDir);
			std::string cmd = std::string("rm -r " + tempDir);
			if (!getInstance().keepTempDir.empty()) {
				if (rename(tempDir.c_str(), getInstance().keepTempDir.c_str()) != 0) {
					cmd = std::string("mv " + tempDir + " " + getInstance().keepTempDir);
				}
			}
			system(cmd.c_str()); // belt & suspenders
		}
	}

	static void addTemp(std::string tempFile) {
		getInstance().tempFiles.insert(tempFile);
	}
	static void removeTemp(std::string tempFile) {
		if (getInstance().tempFiles.find(tempFile) != getInstance().tempFiles.end()) {
			getInstance().tempFiles.erase(tempFile);
			if (!getInstance().keepTempDir.empty()) {
				if (rename(tempFile.c_str(), getInstance().keepTempDir.c_str()) != 0) {
					LOG_DEBUG_OPTIONAL(1, true, "Could not move " << tempFile << " to " << getInstance().keepTempDir);
				}
			}
			unlink(tempFile.c_str()); // belt & suspenders
		}
	}

	static void trackChild(pid_t pid) {
		getInstance().children.insert(pid);
		LOG_DEBUG_OPTIONAL(2, true, "Tracking child pid: " << pid);
	}
	static void releaseChild(pid_t pid) {
		PidSet &children = getInstance().children;
		PidSet::iterator it = children.find(pid);
		if (it != children.end()) {
			LOG_DEBUG_OPTIONAL(2, true, "Released child pid: " << *it);
			children.erase(*it);
		}
	}
	void setKeepTempDir(std::string _keepTempDir) {
		keepTempDir = _keepTempDir;
		if (!_keepTempDir.empty()) {
			mkdir(_keepTempDir.c_str(), 0700);
			LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), "Copying / preserving all temporary files to: " << _keepTempDir);
		}
	}

protected:
	typedef std::vector< struct sigaction > SigHandlers;
	static SigHandlers _initHandlers() {
		int numSignals = 16;
		SigHandlers sh;
		sh.resize(numSignals);
		for(int sig = 0; sig < numSignals ; sig++) {
			struct sigaction &oact = sh[sig];
			bzero(&oact, sizeof(oact));
			sigaction(sig, NULL, &oact);
			if ( oact.sa_handler != SIG_ERR && oact.sa_handler != SIG_IGN && oact.sa_handler != NULL ) {
			      LOG_DEBUG_OPTIONAL(2, true, "set fall through signal handler for " << sig << " to " << oact.sa_handler );
			}
		}
		return sh;
	}
	static SigHandlers &getSigHandlers() {
		static SigHandlers _sh = _initHandlers();
		return _sh;
	}
	static struct sigaction &getSigHandle(int sig) {
		return getSigHandlers()[sig];
	}
	static void _setSignals() {
		_setSignal(SIGHUP, Cleanup::cleanup);
		_setSignal(SIGINT, Cleanup::cleanup);
		_setSignal(SIGQUIT, Cleanup::cleanup);
		_setSignal(SIGABRT, Cleanup::cleanup);
		_setSignal(SIGUSR1, Cleanup::cleanup);
		_setSignal(SIGSEGV, Cleanup::cleanup);
		_setSignal(SIGUSR2, Cleanup::cleanup);
		_setSignal(SIGPIPE, Cleanup::cleanup);
		_setSignal(SIGTERM, Cleanup::cleanup);
	}
	static void _unsetSignals() {
		LOG_DEBUG_OPTIONAL(2, true, "Unsetting cleanup signal handlers");
		_setSignal(SIGHUP, getSigHandle(SIGHUP).sa_handler, false);
		_setSignal(SIGINT, getSigHandle(SIGINT).sa_handler, false);
		_setSignal(SIGQUIT, getSigHandle(SIGQUIT).sa_handler, false);
		_setSignal(SIGABRT, getSigHandle(SIGABRT).sa_handler, false);
		_setSignal(SIGUSR1, getSigHandle(SIGUSR1).sa_handler, false);
		_setSignal(SIGSEGV, getSigHandle(SIGSEGV).sa_handler, false);
		_setSignal(SIGUSR2, getSigHandle(SIGUSR2).sa_handler, false);
		_setSignal(SIGPIPE, getSigHandle(SIGPIPE).sa_handler, false);
		_setSignal(SIGTERM, getSigHandle(SIGTERM).sa_handler, false);
	}
	static void _setSignal(int sig, sighandler_t handler, bool keepOld = true) {
		//void (*prev_fn)(int);
		//sighandler_t prev_fn;
		//prev_fn = signal(sig, handler);
		//if (prev_fn == SIG_IGN) signal(sig, SIG_IGN); // no change if it was ignoring
		struct sigaction act;
		bzero(&act, sizeof(act));
		act.sa_handler = handler;
		act.sa_flags = SA_NOCLDSTOP | SA_RESTART;
		sigaction(sig, &act, keepOld ? &getSigHandle(sig) : NULL);
	}
	void _clean(int param) {
		// save this signals next action, and reset the rest.
		struct sigaction &oact = getSigHandle(param);
		if (param != 0 && param != SIGTERM) {
			LOG_VERBOSE_OPTIONAL(1, true, "Sending SIGTERM to master pid: " << getpid());
			kill(getpid(), SIGTERM);
		}
		sleep(2);
		_unsetSignals();
		FileSet _copy;

		std::stringstream ss;
		if (param != 0) {
			// kill children with same signal
			if (!children.empty())
				ss << "\tTERMinating: ";
			for(PidSet::iterator it = children.begin(); it != children.end(); it++) {
				kill(*it, SIGTERM);
				ss << *it << " , ";
			}
		}
		children.clear();

		if (!tempFiles.empty())
			ss << "\tUnlinked: ";
		_copy = tempFiles;
		for(FileSet::iterator it = _copy.begin(); it != _copy.end(); it++) {
			ss << *it << " , ";
			removeTemp(*it);
		}
		tempFiles.clear();

		if (!tempDirs.empty())
			ss << "rm -rf : ";
		_copy = tempDirs;
		for(FileSet::iterator it = _copy.begin(); it != _copy.end(); it++) {
			ss << *it << " , ";
			removeTempDir(*it);
		}
		tempDirs.clear();

		std::string mesg = ss.str();
		if (param != 0) {
			LOG_ERROR(1, "Cleaned up after signal: " << param << "... " << mesg);
			if ( oact.sa_handler != SIG_ERR && oact.sa_handler != SIG_IGN && oact.sa_handler != NULL ) {
				LOG_DEBUG_OPTIONAL(2, true, "Now running default signal handler" );
				(oact.sa_handler)( param );
			}
		} else {
			LOG_DEBUG_OPTIONAL(2, true, "Cleanup: " << mesg);
		}

		int waitCount = 0;
		bool childrenLive = true;
		while (waitCount++ < 3 && childrenLive) {
			childrenLive = false;
			int status = 0;
			for(PidSet::iterator it = children.begin(); it != children.end(); it++) {
				int wpid = waitpid(*it, &status, WNOHANG);
				if (wpid == 0 || !WIFEXITED(status)) {
					childrenLive = true;
				}
			}
			sleep(1);
		}
		if (childrenLive) {
			int status = 0;
			for(PidSet::iterator it = children.begin(); it != children.end(); it++) {
				int wpid = waitpid(*it, &status, WNOHANG);
				if (wpid == 0 || !WIFEXITED(status)) {
					LOG_WARN(1, "KILLing: " << *it);
					kill(*it, SIGKILL);
				}
			}
		}
		if (param != 0) {
			Logger::getAbortFlag() = true;
			if (param > 0)
				LOG_THROW("Terminating because of signal " << param);
		}
	}

	Cleanup() : tempFiles(), tempDirs(), myPid(0), keepTempDir() {
		myPid = getpid();
		_setSignals();
		setKeepTempDir(	GeneralOptions::getOptions().getKeepTempDir() );
	}

private:
	FileSet tempFiles;
	FileSet tempDirs;
	pid_t myPid;
	PidSet children;
	std::string keepTempDir;
};

class Fork {
public:
	static pid_t forkCommand(std::string command) {
		std::string dir = Options::getOptions().getTmpDir();
		std::string temp = dir + UniqueName::generateUniqueName("/.tmpScript-");
		int child = fork();
		if (child < 0) {
			LOG_THROW("Could not fork a child process");
		} else if (child == 0) {
			// child
			std::fstream cmdFile(temp.c_str(), std::ios_base::out);
			if (!(cmdFile.is_open() && cmdFile.good()))
				LOG_THROW("Could not open " << temp << " for writing!");
			cmdFile << "#!/bin/bash" << std::endl << command << std::endl << "exit $?" << std::endl;
			cmdFile.close();
			if (cmdFile.fail())
				LOG_THROW("Could not close " << temp << " for writing!");
			chmod(temp.c_str(), 0700);
			if (setpgrp() != 0)
				LOG_WARN(1, "Child could not set a new process group");
			execv(temp.c_str(), NULL);
			std::cerr << "ERROR, you should never get here!" << std::endl;
			exit(1);
		} else {
			// parent
			child = 0 - child;
			Cleanup::addTemp(temp);
			Cleanup::trackChild(child);
			LOG_DEBUG_OPTIONAL(1, true, "Executing '" << command << "' in pid " << child << " through script: " << temp);
			return child;
		}
	}
	static int wait(pid_t pid) {
		bool exited = false;
		int status = 0;
		int wpid = 0;
		int exitStatus = -1;

		while (!exited) {
			wpid = waitpid(pid, &status, WNOHANG);
			if (wpid != 0 && WIFEXITED(status)) {
				exitStatus = WEXITSTATUS(status);
				exited = true;
				LOG_DEBUG_OPTIONAL(1, true, "waitpid(" << pid << ") returned: " << wpid << ", " << exitStatus);
			} else if (wpid < 0) {
				LOG_DEBUG_OPTIONAL(1, true, "waitpid(" << pid << ") returned: " << wpid);
				exitStatus = 1;
				exited = true;
			} else {
				LOG_DEBUG_OPTIONAL(2, true, "Waiting for pid: " << pid);
				sleep(1);
			}
		}
		Cleanup::releaseChild(pid);
		LOG_DEBUG_OPTIONAL(1, wpid != 0, "waitpid(" << pid << ") returned: " << exitStatus << " for " << wpid);
		return exitStatus;
	}
};

template<typename T>
class Partition {
public:
	typedef typename T::iterator I;
	Partition(T& container, int num, int size) {
		assert(num < size);
		assert(num >= 0);
		I c_begin = container.begin(), c_end = container.end();
		long batch = container.size() / size + 1;
		assert(c_begin + container.size() == c_end);
		begin = c_begin;
		if (batch * num > (int) container.size())
			begin = c_end;
		else
			begin += batch * num;
		end = c_begin;
		if (batch * (num + 1) > (int) container.size())
			end = c_end;
		else
			end += batch * (num + 1);
		LOG_DEBUG_OPTIONAL(2, true, "Partition(" << container.size() << "," << num << "," << size <<"): " << (begin-c_begin) << ", " << (end-c_begin));
	}
	I begin, end;
};


template<typename RandomAccessIterator, typename Comp >
class MergeSortedRanges {
public:
	class SortedRange {
	public:
		SortedRange(RandomAccessIterator _b, RandomAccessIterator _e) : begin(_b), end(_e) {}
		RandomAccessIterator begin, end;
	};
	typedef std::vector< SortedRange > SortedRangeVector;
	class Compare {
	public:
		Compare(const Comp &c) : comp(c) {}
		inline bool operator()(SortedRange &a, SortedRange &b) {
			return comp(*a.begin, *b.begin);
		}
	private:
		Comp comp;
	};
	MergeSortedRanges(Comp _c) : isHeap(false), comp(_c) {}
	~MergeSortedRanges() {}
	void addSortedRange(const SortedRange &sr) {
		if (sr.begin != sr.end) {
			sortedRanges.push_back(sr);
			if (isHeap)
				std::push_heap(sortedRanges.begin(), sortedRanges.end());
		}
	}
	void addSortedRange(RandomAccessIterator begin, RandomAccessIterator end) {
		SortedRange sr(begin, end);
		addSortedRange(sr);
	}

	RandomAccessIterator getNext() {
		if (!isHeap) {
			std::make_heap(sortedRanges.begin(), sortedRanges.end(), comp);
			isHeap = true;
		}
		SortedRange f = sortedRanges.front();
		std::pop_heap(sortedRanges.begin(), sortedRanges.end(), comp);
		sortedRanges.pop_back();
		RandomAccessIterator nextVal = f.begin++;
		addSortedRange(f);
		return nextVal;
	}

	bool hasNext() {
		return !sortedRanges.empty();
	}

private:
	SortedRangeVector sortedRanges;
	bool isHeap;
	Compare comp;

};
#endif

