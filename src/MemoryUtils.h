// $Header: /repository/PI_annex/robsandbox/KoMer/src/MemoryUtils.h,v 1.12 2010-05-06 22:55:05 regan Exp $
//

#ifndef _MEMORY_UTILS_H
#define _MEMORY_UTILS_H

#include <cstdlib>
#include <vector>
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdexcept>
#include <iomanip>
#include <unistd.h>
#include <fstream>

#include <execinfo.h>

#include <boost/unordered_map.hpp>
#include <boost/pool/pool.hpp>
#include <boost/shared_ptr.hpp>


#include "config.h"
#include "Sequence.h"

class StackTrace {
public:
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

class PoolManager {
public:
	virtual void *malloc(unsigned long size) = 0;
	virtual void free(void* memory, unsigned long size) = 0;
};

class ClassicMemory: public PoolManager {
private:
	static ClassicMemory singleton;
public:
	inline void *malloc(unsigned long size) {
		return std::malloc(size);
	}
	inline void free(void* memory, unsigned long size) {
		std::free(memory);
	}

	inline static PoolManager &get() {
		return singleton;
	}
};

class BoostPoolManager: public PoolManager {
private:
	static BoostPoolManager singleton;
public:
	typedef boost::pool<> Pool;
	typedef boost::shared_ptr<Pool> PoolPtr;
	typedef std::vector<PoolPtr> SizePools;

	BoostPoolManager() :
		mallocs(0), frees(0) {
	}
	~BoostPoolManager() {
		purgePools();
	}

	void *malloc(unsigned long size) {
		mallocs++;
		return getPool(size).ordered_malloc();
	}
	void free(void *memory, unsigned long size) {
		getPool(size).free(memory);
		frees++;
		if (mallocs > 10000 && frees > mallocs / 4) {
			releasePools();
			mallocs -= frees;
			frees = 0;
		}
	}
	inline static PoolManager &get() {
		return ClassicMemory::get();
	}
	//inline static PoolManager &get() { return singleton; }

private:
	SizePools pools;
	unsigned long mallocs, frees;

	Pool &getPool(unsigned long poolByteSize) {
		if (pools.size() <= poolByteSize)
			pools.resize(poolByteSize + 1);
		if (pools[poolByteSize].get() == NULL) {
			PoolPtr pool(new Pool(poolByteSize));
			pools[poolByteSize] = pool;
		}
		return *(pools[poolByteSize]);
	}

public:

	// frees ALL memory that EVERY pool has ever allocated
	// Use at end of program or before calling KmerSizer::set()...
	void purgePools() {
		pools.clear();
	}

	// frees unused memory in existing pools
	// use anytime you want
	void releasePools() {
		for (SizePools::iterator it = pools.begin(); it != pools.end(); it++)
			if (*it != NULL)
				(*it)->release_memory();
	}
};

class MemoryUtils {
public:

	static time_t getLastCall() {
		static time_t last = 0;
		time_t copy = last;
		time(&last);
		return copy;
	}
	static time_t getTime() {
		time_t copy;
		time(&copy);
		return copy - getLastCall();
	}

	static std::string getMemoryUsage() {
		std::stringstream ss;
		rusage usage;
		if (getrusage(RUSAGE_SELF, &usage) != 0)
			throw std::runtime_error("Could not probe memory");

		double t = (double) usage.ru_utime.tv_sec
				+ (double) usage.ru_utime.tv_usec / 1000000.0;
		ss << "time: " << getTime() << " utime: " << std::fixed
				<< std::setprecision(1) << t;

		t = usage.ru_stime.tv_sec + (double) usage.ru_stime.tv_usec / 1000000.0;
		ss << " stime: " << std::fixed << std::setprecision(1) << t;

		if (usage.ru_maxrss != 0) {
			ss << " maxrss: " << usage.ru_maxrss;
			ss << " ixrss: " << usage.ru_ixrss;
			ss << " idrss: " << usage.ru_idrss;
			ss << " isrss: " << usage.ru_isrss;
		} else {
			pid_t pid = getpid();
			char buffer[1024];
			sprintf(buffer, "/proc/%d/statm", pid);
			std::fstream statm(buffer, std::fstream::in);
			statm.getline(buffer, 1024);
			statm.close();
			ss << " statm: " << buffer;
		}

		return ss.str();
	}

};

#endif

//
// $Log: MemoryUtils.h,v $
// Revision 1.12  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.11  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
// Revision 1.10.6.1  2010-05-04 19:49:51  regan
// minor rework on include headers
//
// Revision 1.10.14.1  2010-05-06 18:45:36  regan
// broke it...
//
// Revision 1.10  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.9.2.2  2010-04-15 17:29:02  regan
// checkpoint, working with some optimizations
//
// Revision 1.9.2.1  2010-04-04 15:58:02  regan
// fixed assertion code to obey debug rules
//
// Revision 1.9  2010-03-04 06:37:31  regan
// bugfix
//
// Revision 1.8  2010-02-26 13:01:16  regan
// reformatted
//
// Revision 1.7  2009-12-22 23:12:00  regan
// added wall time to stats
//
// Revision 1.6  2009-12-18 19:04:20  regan
// helpe function to print a stack trace
//
// Revision 1.5  2009-11-11 07:57:23  regan
// built framework for autoPromote (not working) - make_heap is broken
//
// Revision 1.4  2009-11-04 18:23:14  regan
// added a memory usage reporting function
//
// Revision 1.3  2009-11-02 21:19:25  regan
// fixed types and boundary tests
//
// Revision 1.2  2009-11-02 18:27:43  regan
// refactor memory pools (out)
//
// Revision 1.1  2009-10-31 23:44:17  regan
// fixed bug in KmerArray::remove
// refactored memory pool out of KmerArray
//
//
