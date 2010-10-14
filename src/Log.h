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
 */

#ifndef LOG_H_
#define LOG_H_

#include <cstdlib>
#include <iostream>
#include <ctime>

class Logger
{
	typedef std::ostream* OstreamPtr;
	OstreamPtr _os;
	std::string _attribute;
	unsigned int _level;
	inline std::string getTime() {
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
	Logger(std::ostream &os, std::string attr, unsigned int level) : _os(&os), _attribute(attr), _level(level) {}
	Logger(const Logger &copy) {
		*this = copy;
	}
	Logger &operator=(const Logger &copy) {
		if (this == &copy)
			return *this;
		_os = copy._os;
		_attribute = copy._attribute;
		_level = copy._level;
		return *this;
	}
	inline bool isActive(unsigned int level = 1) {
		return _os != NULL && _level >= level;
	};
	inline void setOstream(std::ostream &os) {
		_os = &os;
	}
	inline void setLevel(unsigned int level) {
		_level = level;
	}
	inline unsigned int &getLevel() {
		return _level;
	}
	inline void unsetOstream() {
		_os = NULL;
	}
	template<typename T>
	inline std::ostream &operator<<(T log) {
		if (isActive())
			return *_os << getTime() << " " << _attribute << ": " << log;
		else
			return misuseWarning();
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
};

class Log
{
	static Logger verboseOstream;
	static inline Logger &getVerboseOstream() {
		return verboseOstream;
	}
	static Logger debugOstream;
	static inline Logger &getDebugOstream() {
		return debugOstream;
	}
	static Logger warningOstream;
	static inline Logger &getWarningOstream() {
		return warningOstream;
	}
	static Logger errorOstream;
	static inline Logger &getErrorOstream() {
		return errorOstream;
	}
public:
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
	static inline Logger &Verbose(std::string msg) {
		Logger &log = getVerboseOstream();
		if (!msg.empty()) {
			#pragma omp critical(Log)
			log << msg << std::endl;
		}
		return log;
	}
	static inline Logger &Debug() {
		std::string msg;
		return Debug(msg);
	}
	static inline Logger &Debug(std::string msg) {
		Logger &log =  getDebugOstream();
		if (!msg.empty()) {
			#pragma omp critical(Log)
			log << msg << std::endl;
		}
		return log;
	}
	static inline Logger &Warn() {
		std::string msg;
		return Warn(msg);
	}
	static inline Logger &Warn(std::string msg) {
		Logger &log =  getWarningOstream();
		if (!msg.empty()) {
			#pragma omp critical (Log)
			log << msg << std::endl;
		}
		return log;
	}
	static inline Logger &Error() {
		std::string msg;
		return Error(msg);
	}
	static inline Logger &Error(std::string msg) {
		Logger &log =  getErrorOstream();
		if (!msg.empty()) {
			#pragma omp critical (Log)
			log << msg << std::endl;
		}
		return log;
	}
};

#ifdef NDEBUG
// Higher possible verbosity and debug levels with NDEBUG set, optimized out otherwise
#define LOG_VERBOSE(level, log) if ( Log::isVerbose(level)) { Log::Verbose() << log << std::endl; }
#define LOG_DEBUG(level,   log) if ( Log::isDebug(level)  ) { Log::Debug()   << log << std::endl; }
#else
#define LOG_VERBOSE(level, log) if ( level <=2 && Log::isVerbose(level)) { Log::Verbose() << log << std::endl; }
#define LOG_DEBUG(level,   log) if ( level <=2 && Log::isDebug(level)  ) { Log::Debug()   << log << std::endl; }
#endif
#define LOG_WARN(level,    log) if ( Log::isWarn(level)   ) { Log::Warn()    << log << std::endl; }
#define LOG_ERROR(level,   log) if ( Log::isError(level)  ) { Log::Error()   << log << std::endl; }

#define LOG_VERBOSE_MT(level, log) if ( Log::isVerbose(level) ) { std::stringstream ss ; ss << log << std::endl; Log::Verbose(ss.str()); }
#define LOG_DEBUG_MT(level,   log) if ( Log::isDebug(level)   ) { std::stringstream ss ; ss << log << std::endl; Log::Debug(ss.str()); }

#endif /* LOG_H_ */
