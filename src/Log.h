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
	int _thisLevel;
	static void *world;
	inline std::string getRank() {
#ifdef _USE_MPI
		if (world != NULL)
			return " R" + boost::lexical_cast<std::string>( ((mpi::communicator*)world)->rank() );
		else
#endif
			return std::string();
	}
	inline std::string getThread() {
#ifdef _USE_OPENMP
		return " T" + boost::lexical_cast<std::string>( omp_get_thread_num() );
#else
		return std::string();
#endif
	}
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
	Logger(std::ostream &os, std::string attr, unsigned int level) : _os(&os), _attribute(attr), _level(level), _thisLevel(-1) {}
	Logger(const Logger &copy) {
		*this = copy;
	}
	static void setWorld(void *_w) {
		world = _w;
	}
#ifdef _USE_MPI
	std::string gatherMessages(std::string msg) {
		if (world != NULL) {
			//*_os << getStamp() << "Entering gatherMessages: " << msg << std::endl;
			mpi::communicator &w = *((mpi::communicator*) world);
			std::string out[w.size()];
			if (!msg.empty())
				msg = getStamp("M") + msg + "\n";
			try {
				mpi::gather(w, msg, out, 0);
			} catch (...) {
				*_os << "ERROR: Failed to gather all messages: " << msg;
				throw;
			}
			std::stringstream ss;
			for(int i = 0 ; i < w.size(); i++)
				ss << out[i];
			msg = ss.str();
			if (!msg.empty() && w.rank() == 0)
				return "MPI Gathered Log Entries:\n" + msg;
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
		if (world == NULL)
			return true;
		else
#ifdef  _USE_MPI
			return ((mpi::communicator*)world)->rank() == 0;
#else
			return true;
#endif
	}
	inline bool isActive(unsigned int level = 1) {
		bool isActive = _os != NULL && _level >= level;
		if (isActive && _thisLevel < (int) level)
			_thisLevel = level;
		return isActive;
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
	inline void setThisLevel(int level) {
		_thisLevel = level;
	}
	inline std::string getThisLevel() {
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
	inline std::string getStamp(std::string attribLabel = "") {
		return getTime() + " " + _attribute + getThisLevel() + attribLabel + getRank() + getThread() + ": ";
	}
	template<typename T>
	inline std::ostream &operator<<(T log) {
		if (isActive())
			return *_os << getStamp() << log;
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
	static inline Logger &_log(Logger &log, std::string &msg, bool bypassmpi = true) {
#ifdef _USE_MPI
		if (!bypassmpi)
			msg = log.gatherMessages(msg);
#endif
		if (!msg.empty()) {
			#pragma omp critical(Log)
			{
				log << msg << std::endl;
			}
		}
		return log;
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
		return _log(getErrorOstream(), msg, bypassMPI);
	}
};

#define LOG_VERBOSE(level, log) if ( Log::isVerbose(level)) { std::stringstream ss ; ss << log; Log::Verbose(ss.str(), level >= 3); }
#define LOG_DEBUG(level,   log) if ( Log::isDebug(level)  ) { std::stringstream ss ; ss << log; Log::Debug(ss.str(), level >= 2); }
#define LOG_WARN(level,    log) if ( Log::isWarn(level)   ) { std::stringstream ss ; ss << log; Log::Warn(ss.str(), true); }
#define LOG_ERROR(level,   log) if ( Log::isError(level)  ) { std::stringstream ss ; ss << log; Log::Error(ss.str(), true); }

#define LOG_VERBOSE_OPTIONAL(level, test, log) if ( test && Log::isVerbose(level)) { std::stringstream ss ; ss << log; Log::Verbose(ss.str(), true); }
#define LOG_DEBUG_OPTIONAL(level,   test, log) if ( test && Log::isDebug(level)  ) { std::stringstream ss ; ss << log; Log::Debug(ss.str(), true); }

#endif /* LOG_H_ */
