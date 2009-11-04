// $Header: /repository/PI_annex/robsandbox/KoMer/src/MemoryUtils.h,v 1.4 2009-11-04 18:23:14 regan Exp $
//

#ifndef _MEMORY_UTILS_H
#define _MEMORY_UTILS_H

#include <boost/unordered_map.hpp>
#include <boost/pool/pool.hpp>

#include <tr1/memory>
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

class PoolManager
{
public:
  virtual void *malloc(unsigned long size) = 0;
  virtual void free(void* memory, unsigned long size) = 0;
};

class ClassicMemory : public PoolManager
{
private:
  static ClassicMemory singleton;
public:
  inline void *malloc(unsigned long size) { return std::malloc(size); }
  inline void free(void* memory, unsigned long size) { std::free(memory); }
  
  inline static PoolManager &get() { return singleton; }
};

class BoostPoolManager : public PoolManager
{
private:
  static BoostPoolManager singleton;
public:
  typedef boost::pool< > Pool;
  typedef std::tr1::shared_ptr<Pool> PoolPtr;
  typedef std::vector< PoolPtr > SizePools;

  BoostPoolManager() : mallocs(0), frees(0) {}
  ~BoostPoolManager() { purgePools(); };
  
  void *malloc(unsigned long size) {
  	mallocs++;
  	return getPool(size).ordered_malloc();
  }
  void free(void *memory, unsigned long size) {
  	getPool(size).free(memory);
  	frees++;
  	if ( mallocs > 10000 && frees > mallocs / 4) {
  	  releasePools();
  	  mallocs-=frees;
  	  frees=0;
  	}
  }
  inline static PoolManager &get() {return ClassicMemory::get();}
  //inline static PoolManager &get() { return singleton; } 
   
private:
  SizePools pools;
  unsigned long mallocs,frees;

  Pool &getPool(unsigned long poolByteSize) {
     if (pools.size() <= poolByteSize)
       pools.resize(poolByteSize+1);
     if ( pools[poolByteSize].get() == NULL ) {
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
  	 for(SizePools::iterator it = pools.begin() ; it != pools.end(); it++)
  	   if (*it != NULL)
         (*it)->release_memory();
  }
};


class MemoryUtils
{
public:

	static std::string getMemoryUsage()
	{
	  std::stringstream ss;
	  rusage usage;
	  if (getrusage(RUSAGE_SELF, &usage) != 0)
	    throw std::runtime_error("Could not probe memory");
	    
	  double t = (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1000000.0;     
	  ss << " utime: " << std::fixed << std::setprecision(1) << t;
	  
	  t = usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec / 1000000.0;       
	  ss << " stime: " << std::fixed << std::setprecision(1) << t;
	  
	  if (usage.ru_maxrss != 0) { 
	    ss << " maxrss: " << usage.ru_maxrss;
	    ss << " ixrss: "  << usage.ru_ixrss;
	    ss << " idrss: "  << usage.ru_idrss;
	    ss << " isrss: "  << usage.ru_isrss;
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