/*
 * Kmernator/src/Log.h
 *
 *  Created on: Sep 1, 2010
 *      Author: regan
 *
 * Copyright 2010 The Regents of the University of California.
 * All rights reserved.
 *
 * The United States Government has rights in this work pursuant
 * to contracts DE-AC03-76SF00098, W-7405-ENG-36 and/or
 * W-7405-ENG-48 between the United States Department of Energy
 * and the University of California.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that: (1) source distributions retain this entire
 * copyright notice and comment, and (2) distributions including
 * binaries display the following acknowledgement:  "This product
 * includes software developed by the University of California,
 * JGI-PSF and its contributors" in the documentation or other
 * materials provided with the distribution and in all advertising
 * materials mentioning features or use of this software.  Neither the
 * name of the University nor the names of its contributors may be
 * used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE.
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

	inline std::string getThread() const {
#ifdef _USE_OPENMP
		return " T" + boost::lexical_cast<std::string>( omp_get_thread_num() );
#else
		return std::string();
#endif
	}
	inline std::string getTime() const {
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
	Logger(std::ostream &os, std::string attr, unsigned int level) : _os(&os), _attribute(attr), _level(level), _thisLevel(-1) {}
	Logger(const Logger &copy) {
		*this = copy;
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
	static Logger verboseOstream;
	static Logger debugOstream;
	static Logger warningOstream;
	static Logger errorOstream;

	static void setErrorMessage(std::string &msg) {
		getErrorMessages() += msg + "\n";
	}

	static inline Logger &getDebugOstream() {
		return debugOstream;
	}
	static inline Logger &getVerboseOstream() {
		return verboseOstream;
	}
	static inline Logger &getWarningOstream() {
		return warningOstream;
	}
	static inline Logger &getErrorOstream() {
		return errorOstream;
	}
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
		return debugOstream;
	}
	static inline Logger &getVerbose() {
		return verboseOstream;
	}
	static inline Logger &getWarning() {
		return warningOstream;
	}
	static inline Logger &getError() {
		return errorOstream;
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

#define LOG_THROW(log) {  std::stringstream ss ; ss << log; Log::Error(ss.str() , true); throw LoggedException(Log::getErrorMessages()); }
#define LOG_VERBOSE(level, log) if ( Log::isVerbose(level)) { std::stringstream ss ; ss << log; Log::Verbose(ss.str(), level >= 3); }
#define LOG_DEBUG(level,   log) if ( Log::isDebug(level)  ) { std::stringstream ss ; ss << log; Log::Debug(ss.str(), level >= 2); }
#define LOG_WARN(level,    log) if ( Log::isWarn(level)   ) { std::stringstream ss ; ss << log; Log::Warn(ss.str(), true); }
#define LOG_ERROR(level,   log) if ( Log::isError(level)  ) { std::stringstream ss ; ss << log; Log::Error(ss.str(), true); }

#define LOG_VERBOSE_OPTIONAL(level, test, log) if ( test && Log::isVerbose(level)) { std::stringstream ss ; ss << log; Log::Verbose(ss.str(), true); }
#define LOG_DEBUG_OPTIONAL(level,   test, log) if ( test && Log::isDebug(level)  ) { std::stringstream ss ; ss << log; Log::Debug(ss.str(), true); }

#endif /* LOG_H_ */
