// $Header: /repository/PI_annex/robsandbox/KoMer/src/MemoryUtils.h,v 1.1 2009-10-31 23:44:17 regan Exp $
//

#ifndef _MEMORY_UTILS_H
#define _MEMORY_UTILS_H

#include <boost/unordered_map.hpp>
#include <boost/pool/pool.hpp>

#include <tr1/memory>
#include <cstdlib>
#include <vector>

class PoolManager
{
public:
  virtual void *malloc(unsigned long size) = 0;
  virtual void free(void* memory, unsigned long size) = 0;
};

class ClassicMemory : public PoolManager
{
public:
  void *malloc(unsigned long size) { return std::malloc(size); }
  void free(void* memory, unsigned long size) { std::free(memory); }
  
  inline static PoolManager &get() { static ClassicMemory singleton; return singleton; }
};

class BoostPoolManager : public PoolManager
{
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
  // inline static PoolManager &get() {return ClassicMemory::get();}
  inline static PoolManager &get() { static BoostPoolManager singleton; return singleton; } 
   
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


#endif

//
// $Log: MemoryUtils.h,v $
// Revision 1.1  2009-10-31 23:44:17  regan
// fixed bug in KmerArray::remove
// refactored memory pool out of KmerArray
//
//