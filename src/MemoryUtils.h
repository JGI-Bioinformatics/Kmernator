//
// Kmernator/src/MemoryUtils.h
//
// Author: Rob Egan
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
		pid_t pid = getpid();
		char filename[128];
		sprintf(filename, "/proc/%d/statm", pid);
		std::fstream statm(filename, std::fstream::in);
		std::string buffer;
		statm >> buffer;
		statm.close();
		return std::string("statm: " + buffer);
	}
	static std::string getMemoryUsage2() {
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
			ss << " " << getMemoryUsage();
		}
		return ss.str();
	}

	static std::string getMmapUsage() {
		std::stringstream ss;
		std::string buffer;
		pid_t pid = getpid();
		ss << "smaps for " << pid << ": " << std::endl;
		std::fstream meminfo("/proc/meminfo", std::fstream::in);
		while ( !meminfo.eof() ) {
			getline(meminfo, buffer);
			if (buffer.length() == 0)
				break;
			ss << buffer << std::endl;
		}
		meminfo.close();

		char filename[128];
		sprintf(filename, "/proc/%d/smaps", pid);
		std::fstream smaps(filename, std::fstream::in);
		std::string record;
		while ( !smaps.eof() ) {
			getline(smaps, buffer);
			if (buffer.length() == 0)
				break;
			if (buffer.find(" kB") == std::string::npos) {
				if (record.find("/") != std::string::npos && record.find("lib") == std::string::npos)
					ss << record.substr(record.find("/")) << std::endl;
				record = buffer + "\t";
			} else {
				if (buffer.find(" 0 kB") == std::string::npos && buffer.find("PageSize:") == std::string::npos)
					record = record + buffer + "\t";
			}
		}
		smaps.close();
		if (record.find("/") != std::string::npos && record.find("lib") == std::string::npos)
			ss << record.substr(record.find("/")) << std::endl;

		return ss.str();
	}

};

#endif

