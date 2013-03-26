/*
 * Kmernator/src/Log.h
 *
 *  Created on: Sep 1, 2010
 *      Author: regan
 *
 *****************

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

 *****************
 *
 *
 * Log guildines:
 * Verbose 1 - Major milestones, information
 * Verbose 2 - Progress/batch milestones
 * Verbose 3 - more information
 *
 * Debug 1 - debug major steps
 * Debug 2 - debug progress/batch steps
 * Debug 3 - detail progress/batch
 * Debug 4 - per outer loop diagnostics
 * Debug 5 - per inner loop (per-read / kmer) diagnostics
 */

#ifndef LOG_H_
#define LOG_H_

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <unistd.h>
#include <sys/syscall.h>
#include <csignal>
#include <set>
#include <execinfo.h>

#include "config.h"
#include <boost/lexical_cast.hpp>

class Logger
{
	typedef std::ostream* OstreamPtr;
	OstreamPtr _os;
	std::string _attribute;
	unsigned int _level;
	mutable int _thisLevel;

	static bool &_debugGather() {
		static bool _debugGather = false;
		return _debugGather;
	}
	static void **_getWorld() {
		static void * _world = NULL;
		return &_world;
	}
	static int &_getWorldRank() {
		static int _rank = -1;
		return _rank;
	}
	static std::string &getRank() {
		static std::string _rankHeader;
		return _rankHeader;
	}

	static inline std::string getThread() {
#ifdef _USE_OPENMP
		return " T" + boost::lexical_cast<std::string>( omp_get_thread_num() );
#else
		return " T" + boost::lexical_cast<std::string>( syscall(SYS_gettid) ); // gettid() );
#endif
	}
	static inline std::string getTime() {
		time_t rawtime;
		struct tm *timeinfo;
		time (&rawtime);
		timeinfo = localtime( &rawtime );
		char buf[32];
		strftime(buf, 32, "%x %X", timeinfo);
		return std::string( buf );
	}
	std::ostream &misuseWarning() {
		return std::cerr << "WARNING: Using Logger(" << _attribute << ") when it was explicitly disconnected!" << std::endl;
	}
public:
	Logger(std::ostream &os, std::string attr, unsigned int level) : _os(&os), _attribute(attr), _level(level), _thisLevel(-1) {
		prepBuffer();
	}
	Logger(const Logger &copy) {
		*this = copy;
		prepBuffer();
	}
	~Logger() {
		releaseBuffer();
	}
	static void prepBuffer() {
		if (*getReservedBuffer() == NULL && ! getAbortFlag()) {
#pragma omp critical
			{
				if (*getReservedBuffer() == NULL) {
					*getReservedBuffer() = new char[ 64 * 1024 ];
				}
			}
		}
	}
	static void releaseBuffer() {
		if (*getReservedBuffer() != NULL) {
#pragma omp critical
			{
				if (*getReservedBuffer() != NULL) {
					delete [] *getReservedBuffer();
					*getReservedBuffer() = NULL;
				}
			}

		}
	}
private:
	static char **getReservedBuffer() {
		static char *_reservedBuffer = NULL;
		return &_reservedBuffer;
	}

public:

	static bool &getAbortFlag() {
		static bool _ = false;
		return _;
	}
	static void setWorld(void *_w, bool debugGather = false) {
		*_getWorld() = _w;
		_debugGather() = debugGather;
		if (_w == NULL)
			return;
#ifdef _USE_MPI
		_getWorldRank() = (*((mpi::communicator**)_getWorld()))->rank();
		getRank() = " R" + boost::lexical_cast<std::string>(_getWorldRank());
#endif
	}
#ifdef _USE_MPI
	std::string gatherMessages(std::string msg) {
		if (*_getWorld() != NULL ) {
			mpi::communicator &w = *(*((mpi::communicator**) _getWorld()));
			std::string out[w.size()];
			if (!msg.empty()) {
				msg = getStamp("M") + msg + "\n";
				if (_debugGather())
					*_os << "--DEBUG-GATHER--" << msg << std::endl;
			}

			try {
				mpi::gather(w, msg, out, 0);
			} catch (...) {
				std::string errMsg("ERROR: Failed to gather all messages: " + msg);
				*_os << errMsg << std::endl;
				throw(errMsg);
			}
			std::stringstream ss;
			for(int i = 0 ; i < w.size(); i++)
				ss << out[i];
			msg = ss.str();
			if (!msg.empty() && _getWorldRank() == 0)
				return "\tMPI Gathered Log Entries:\n" + msg;
			else
				return std::string();
		} else {
			return msg;
		}
	}
#endif
	Logger &operator=(const Logger &copy) {
		if (this == &copy)
			return *this;
		_os = copy._os;
		_attribute = copy._attribute;
		_level = copy._level;
		return *this;
	}
	static inline bool isMaster() {
#ifdef _USE_OPENMP
		if ( omp_get_thread_num() != 0 )
			return false;
#endif
		if (*_getWorld() == NULL)
			return true;
		else
#ifdef  _USE_MPI
			return _getWorldRank() == 0;
#else
		    return true;
#endif
	}
	inline bool isActive(unsigned int level = 1) const {
		bool isActive = _os != NULL && _level >= level;
		if (isActive && _thisLevel < (int) level)
			_thisLevel = level;
		return isActive;
	};
	inline void setOstream(std::ostream &os) {
		_os = &os;
	}
	inline unsigned int &setLevel(unsigned int level) {
		_level = level;
		return _level;
	}
	inline unsigned int &setLevel() {
		return _level;
	}
	inline unsigned int getLevel() const {
		return _level;
	}
	inline void setThisLevel(int level) {
		_thisLevel = level;
	}
	inline std::string getThisLevel() const {
		std::string val;
		if (_thisLevel >= 0) {
			val = boost::lexical_cast<std::string>(_thisLevel);
			_thisLevel = -1;
		}
		return val;
	}
	inline void unsetOstream() {
		_os = NULL;
	}
	inline static std::string getCommonStamp() {
		return getTime() + " " + getRank() + getThread() + ": ";
	}
	inline std::string getStamp(std::string attribLabel = "") const {
		return getTime() + " " + _attribute + getThisLevel() + attribLabel + getRank() + getThread() + ": ";
	}
	template<typename T>
	inline std::ostream &operator<<(T log) {
		if (isActive()) {
			return *_os << toString(log);
		} else
			return misuseWarning();
	}
	template<typename T>
	inline std::string toString(T log) const {
		std::stringstream ss;
		ss << getStamp() << log;
		return ss.str();
	}
	// cast to ostream
	inline operator std::ostream&() {
		if (isActive())
			return *_os;
		else
			return misuseWarning();
	}
	inline OstreamPtr getOstreamPtr() {
		return &( (std::ostream&) *this);
	}

	class ScopedLevel {
	public:
		ScopedTempValue<unsigned int> _s;
		ScopedLevel(Logger &log, unsigned int newLevel) : _s(log.setLevel(), newLevel) {}
	};
};

class Log
{
public:
	static std::string &getErrorMessages() {
		static std::string errorMessages;
		return errorMessages;
	}
	template<typename T> static std::string toString(const T &container) {
		std::stringstream ss;
		for(typename T::const_iterator it = container.begin() ; it != container.end(); it++)
			ss << *it << ", ";
		return ss.str();
	};

private:

	static void setErrorMessage(std::string &msg) {
		getErrorMessages() += msg + "\n";
	}

	static inline Logger &getDebugOstream() {
		static Logger _debugOstream =  Logger( std::cerr, "DEBUG", 0 );
		return _debugOstream;
	}
	static inline Logger &getVerboseOstream() {
		static Logger _verboseOstream = Logger( std::cerr, "INFO", 1 );
		return _verboseOstream;
	}
	static inline Logger &getWarningOstream() {
		static Logger _warningOstream = Logger( std::cerr, "WARNING", 1 );
		return _warningOstream;
	}
public:
	static inline Logger &getErrorOstream() {
		static Logger _errorOstream = Logger( std::cerr, "ERROR", 1 );
		return _errorOstream;
	}

private:
	static inline Logger &_log(Logger &log, std::string &msg, bool bypassmpi = true) {
#ifdef _USE_MPI
		if (!bypassmpi)
			msg = log.gatherMessages(msg);
#endif
		if (!msg.empty()) {
			std::string s(msg + std::string("\n"));
#ifdef _USE_OPENMP
#pragma omp critical(Log)
#endif
			{
				log << s;
			}
		}
		return log;
	}
public:
	static inline const bool printOptions() {
		return Logger::isMaster() && ((isVerbose(1) || isDebug(1)));
	}
	static inline Logger &getDebug() {
		return getDebugOstream();
	}
	static inline Logger &getVerbose() {
		return getVerboseOstream();
	}
	static inline Logger &getWarning() {
		return getWarningOstream();
	}
	static inline Logger &getError() {
		return getErrorOstream();
	}

	static inline void setVerboseOstream(std::ostream &os) {
		getVerboseOstream().setOstream(os);
	}
	static inline void setDebugOstream(std::ostream &os) {
		getDebugOstream().setOstream(os);
	}
	static inline bool isVerbose(unsigned int level = 1) {
		return getVerboseOstream().isActive(level);
	}
	static inline bool isDebug(unsigned int level = 1) {
		return getDebugOstream().isActive(level);
	}
	static inline bool isWarn(unsigned int level = 1) {
		return getWarningOstream().isActive(level);
	}
	static inline bool isError(unsigned int level = 1) {
		return getErrorOstream().isActive(level);
	}
	static inline Logger &Verbose() {
		std::string msg;
		return Verbose(msg);
	}
	static inline Logger &Verbose(std::string msg, bool bypassMPI = true) {
		return _log(getVerboseOstream(), msg, bypassMPI);
	}
	static inline Logger &Debug() {
		std::string msg;
		return Debug(msg);
	}
	static inline Logger &Debug(std::string msg, bool bypassMPI = true) {
		return _log(getDebugOstream(), msg, bypassMPI);
	}
	static inline Logger &Warn() {
		std::string msg;
		return Warn(msg);
	}
	static inline Logger &Warn(std::string msg, bool bypassMPI = true) {
		return _log(getWarningOstream(), msg, bypassMPI);
	}
	static inline Logger &Error() {
		std::string msg;
		return Error(msg);
	}
	static inline Logger &Error(std::string msg, bool bypassMPI = true) {
		setErrorMessage(msg);
		return _log(getErrorOstream(), msg, bypassMPI);
	}
	template<typename T>
	static std::string toString(T begin, T end) {
		std::stringstream ss;
		ss << (end-begin) << ": ";
		for(T it = begin; it != end; it++)
			ss << *it << ", ";
		return ss.str();
	}

};

#define LOG_VERBOSE(level, log) if ( Log::isVerbose(level)) { std::stringstream ss ; ss << log; Log::Verbose(ss.str(), true); }
#define LOG_DEBUG(level,   log) if ( Log::isDebug(level)  ) { std::stringstream ss ; ss << log; Log::Debug(ss.str(), true); }
#define LOG_WARN(level,    log) if ( Log::isWarn(level)   ) { std::stringstream ss ; ss << log; Log::Warn(ss.str(), true); }
#define LOG_ERROR(level,   log) if ( Log::isError(level)  ) { std::stringstream ss ; ss << log; Log::Error(ss.str(), true); }

#define LOG_VERBOSE_OPTIONAL(level, test, log) if ( test && Log::isVerbose(level)) { std::stringstream ss ; ss << log; Log::Verbose(ss.str(), true); }
#define LOG_DEBUG_OPTIONAL(level,   test, log) if ( test && Log::isDebug(level)  ) { std::stringstream ss ; ss << log; Log::Debug(ss.str(), true); }

#define LOG_VERBOSE_GATHER(level, log) if ( Log::isVerbose(level)) { std::stringstream ss ; ss << log; Log::Verbose(ss.str(), false); }
#define LOG_DEBUG_GATHER(level,   log) if ( Log::isDebug(level)  ) { std::stringstream ss ; ss << log; Log::Debug(ss.str(), false); }

class LoggedException : public std::exception {
public:
	LoggedException(std::string msg) throw() : _msg(msg) {}
	LoggedException(const std::exception &e) throw() {
		_msg = e.what();
	}
	LoggedException &operator=(const LoggedException &copy) throw() {
		_msg = copy._msg;
		return *this;
	}
	virtual ~LoggedException() throw() {}
	const char* what() const throw() { return _msg.c_str(); }
private:
	std::string _msg;
};

class StackTrace {
public:
	static void printStackTrace() {
		void *array[200];
		size_t size;

		size = backtrace(array, 200);
		backtrace_symbols_fd(array, size, STDERR_FILENO);
	}
	static std::string getStackTrace() {
		std::stringstream ss;
		void *array[200];
		size_t size, i;
		char **strings;

		size = backtrace(array, 200);
		strings = backtrace_symbols(array, size);

		for (i = 0; i < size; i++)
			ss << strings[i] << std::endl;

		free(strings);
		return ss.str();
	}
};

#define LOG_THROW(log) { Logger::getAbortFlag() = true; Logger::releaseBuffer(); std::cerr << Logger::getCommonStamp() << "Exception!: " << log << std::endl; StackTrace::printStackTrace(); { std::stringstream ss ; ss << log; Log::Error(ss.str(), true); } throw LoggedException(Log::getErrorMessages()); }

#endif /* LOG_H_ */
