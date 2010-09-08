/*
 * Log.h
 *
 *  Created on: Sep 1, 2010
 *      Author: regan
 */

#ifndef LOG_H_
#define LOG_H_

#include <cstdlib>
#include <ostream>
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
	static inline Logger &getVerboseOstream() {
		static Logger verboseOstream = Logger( std::cerr, "INFO", 1 );
		return verboseOstream;
	}
	static inline Logger &getDebugOstream() {
		static Logger debugOstream = Logger( std::cerr, "DEBUG", 0 );
		return debugOstream;
	}
	static inline Logger &getWarningOstream() {
		static Logger warningOstream = Logger( std::cerr, "WARNING", 1 );
		return warningOstream;
	}
	static inline Logger &getErrorOstream() {
		static Logger errorOstream = Logger( std::cerr, "ERROR", 1 );
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

#define LOG_VERBOSE(level, log) if ( Log::isVerbose(level)) { Log::Verbose() << log << std::endl; }
#define LOG_DEBUG(level,   log) if ( Log::isDebug(level)  ) { Log::Debug()   << log << std::endl; }
#define LOG_WARN(level,    log) if ( Log::isWarn(level)   ) { Log::Warn()    << log << std::endl; }
#define LOG_ERROR(level,   log) if ( Log::isError(level)  ) { Log::Error()   << log << std::endl; }

#define LOG_VERBOSE_MT(level, log) if ( Log::isVerbose(level) ) { std::stringstream ss ; ss << log << std::endl; Log::Verbose(ss.str()); }
#define LOG_DEBUG_MT(level,   log) if ( Log::isDebug(level)   ) { std::stringstream ss ; ss << log << std::endl; Log::Debug(ss.str()); }

#endif /* LOG_H_ */
