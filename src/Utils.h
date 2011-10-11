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

#include <sys/types.h>
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
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>

#define foreach BOOST_FOREACH

#include "config.h"
#include "Options.h"
#include "Log.h"
#include "lookup3.h"


class FormatOutput
{
public:
	enum FormatType {FASTQ, FASTA, FASTQ_UNMASKED, FASTA_UNMASKED};

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
private:
	FormatType _type;
	FormatOutput(int type) {
		switch(type) {
		case(0) : _type = FASTQ; break;
		case(1) : _type = FASTA; break;
		case(2) : _type = FASTQ_UNMASKED; break;
		case(3) : _type = FASTA_UNMASKED; break;
		default: throw;
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
public:
	typedef OptionsBaseInterface::StringListType StringListType;
	static std::string generateUniqueName(std::string filename = "") {
		filename += boost::lexical_cast<std::string>( getpid() );
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
			assert(!isFileStream() && !isStringStream());
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
			}
			reset();
		}
		void reset() {
			of.reset();
			ss.reset();
			filePath.clear();
		}
		~OStreamPtr() {
			reset();
		}

		bool isFileStream() const {
			return of.get() != NULL && ss.get() == NULL;
		}
		bool isStringStream() const {
			return of.get() == NULL && ss.get() != NULL;
		}
		std::ostream &operator*() {
			if (isFileStream())
				return *of;
			else if (isStringStream())
				return *ss;
			else
				throw;
		}
		std::string getFilePath() const {
			return filePath;
		}
		std::string getFinalString() {
			assert(isStringStream());
			long long int bytes = ss->tellp();
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
	~OfstreamMap() {
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
			throw std::invalid_argument("getStore() out of range");
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
	static std::string &trimName(std::string &nameLine) {
		if (nameLine.length() == 0)
			return nameLine;

		char &marker = nameLine[0];
		if (marker != '>' && marker != '@') {
			throw std::invalid_argument( (std::string("Can not parse name without a marker: ") + nameLine).c_str() );
		} else {
		    // remove marker
		    nameLine.erase(0,1);
		}

		// trim at first whitespace
		size_t pos = nameLine.find_first_of(" \t\r\n");
		if (pos != std::string::npos) {
			nameLine.erase(pos);
		}
		return nameLine;
	}
	static void parse(Kmernator::RecordPtr record, Kmernator::RecordPtr lastRecord,
			          std::string &name, std::string &bases, std::string &quals,
			          Kmernator::RecordPtr qualRecord = NULL, Kmernator::RecordPtr lastQualRecord = NULL,
			          char fastqStartChar = Kmernator::FASTQ_START_CHAR_ILLUMINA) {
		std::string buf;
		if (*record == '@') {
			// FASTQ
			nextLine(name,  record);
			trimName(name);
			nextLine(bases, record);
			nextLine(buf,   record);
			nextLine(quals, record);
			if (buf[0] != '+') {
				throw "Invalid FASTQ record!";
			}
		} else if (*record == '>') {
			// FASTA
			nextLine(name, record);
			trimName(name);

			// TODO FIXME HACK!
			if (lastRecord == NULL) // only read one line if no last record is given
				lastRecord = record+1;

			while (record < lastRecord && *record != '>') {
				bases += nextLine(buf, record);
			}
			if (qualRecord != NULL) {
				nextLine(buf, qualRecord);
				trimName(buf);
				if (buf != name) {
					throw "fasta and qual do not match names!";
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
			throw "Do not know how to parse this file!";
		}
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

	static std::string commonName(const std::string &readName) {
		return readName.substr(0, readName.length() - 1);
	}
	static int readNum(const std::string &readName) {
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

	static bool isPair(const std::string &readNameA, const std::string readNameB) {
		std::string commonA = commonName(readNameA);
		std::string commonB = commonName(readNameB);
		if (commonA == commonB) {
			int readNumA = readNum(readNameA);
			int readNumB = readNum(readNameB);
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
	      base( fileno( _pipe = popen( (command+getErrIntercept()).c_str(), "w" ) ) ), _cmd(command+getErrIntercept()), _exitStatus(0) {
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
				LOG_WARN(1, "OPipestream::close() '" << _cmd << "' closed with an error: " << _exitStatus << "." << getStdErr());
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
	      base( fileno( _pipe = popen( (command+getErrIntercept()).c_str(), "r" ) ) ),
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
	static bool fileExists(std::string &filePath) {
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

		if (threadSeedsPtr.get() == NULL || (int) threadSeedsPtr->size() < omp_get_max_threads()) {
            #pragma omp critical
			{
				unsigned int baseSeed = 0;
				if (threadSeedsPtr.get() == NULL || (int) threadSeedsPtr->size() < omp_get_max_threads()) {
					SeedPtr newSeeds;
					if (threadSeedsPtr.get() != NULL)
						newSeeds.reset( new SeedVector(*threadSeedsPtr) );
					else
						newSeeds.reset( new SeedVector );
					newSeeds->reserve(omp_get_max_threads());

					while ((int) newSeeds->size() < omp_get_max_threads()) {
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
#endif

