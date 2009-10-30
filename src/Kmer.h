// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.33 2009-10-30 00:51:40 regan Exp $
//

#ifndef _KMER_H
#define _KMER_H
#include <tr1/memory>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <stdexcept>

#include <boost/unordered_map.hpp>
#include <boost/pool/pool.hpp>
#include <boost/functional/hash.hpp>

#include "TwoBitSequence.h"

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 1024
#endif

typedef std::tr1::shared_ptr<TwoBitEncoding> KmerSharedPtr;

// TODO remove ExtraBytes from this!!
class KmerSizer
{
private:
  KmerSizer() : _sequenceLength(21), _extraBytes(0) {}

  SequenceLengthType _sequenceLength;
  unsigned long _extraBytes;

  SequenceLengthType _twoBitLength;
  unsigned long _totalSize;
public:

  static KmerSizer &getSingleton() { /* TODO make thread safe */ static KmerSizer singleton; return singleton; }
  static void set(SequenceLengthType sequenceLength, unsigned long extraBytes=0)
  {
    KmerSizer &singleton = getSingleton();
    singleton._sequenceLength = sequenceLength;
    singleton._extraBytes = extraBytes;
    singleton._twoBitLength =  TwoBitSequence::fastaLengthToTwoBitLength(singleton._sequenceLength);
    singleton._totalSize = singleton._twoBitLength + extraBytes;
  }

  static SequenceLengthType getSequenceLength()  {
    return getSingleton()._sequenceLength;
  }
  static unsigned long getExtraBytes()  {
    return getSingleton()._extraBytes;
  }
  static SequenceLengthType getTwoBitLength()   {
    return getSingleton()._twoBitLength;
  }
  static unsigned long getByteSize()  {
    return getSingleton()._totalSize;
  }
};


typedef void *VoidPtr;

class KmerPtr {

private:
	class Kmer
	{
	public:
	   typedef Kmer *RawKmerPtr;
	   
    private:
       static boost::hash<std::string> &getHasher() { static boost::hash<std::string> hasher; return hasher; }
	   Kmer(); // never construct, just use as cast
	
	#ifdef STRICT_MEM_CHECK
	   TwoBitEncoding _someData[MAX_KMER_SIZE]; // need somedata to hold a pointer and a large amount to avoid memory warnings
	public:
	   const void *_data() const { return _someData;}
	   void *_data()  { return _someData;}
	#else
	   // No data for you!!!
	public:
	   const void *_data() const { return this;}
	   void *_data()  { return this;}
	#endif
	
	public:
	   
	   Kmer(const Kmer &copy) {
	   	  *this = copy;
	   }
	   
	   int compare(const Kmer &other) const
	   {
	     //return memcmp(_data(), other._data(), getTwoBitLength());
         return (toFasta().compare(other.toFasta()));
	   }
	
	   Kmer &operator=(const Kmer &other)
	   {
	   	  if (this == &other)
	   	    return *this;
	   	    
	      memcpy(_data(), other._data(), getTwoBitLength());
	      return *this;
	   }
	
	   KmerPtr operator&()
       {
         return KmerPtr(_data());
	   }

       Kmer *get() {
         return this;
       }
       
	   bool operator ==(const Kmer &other) const
	   {
	      return compare(other) == 0;
	   }
	   bool operator !=(const Kmer &other) const
	   {
	   	  return compare(other) != 0;
	   }
	   bool operator <(const Kmer &other) const
	   {
	   	  return compare(other) < 0;
	   }
	   bool operator <=(const Kmer &other) const
	   {
	   	  return compare(other) <= 0;
	   }
	   bool operator >(const Kmer &other) const
	   {
	   	  return compare(other) > 0;
	   }
	   bool operator >=(const Kmer &other) const
	   {
	   	  return compare(other) >= 0;
	   }
	
	   void swap(Kmer &other)
	   {
	      TwoBitEncoding buffer[getByteSize()];
	      Kmer &temp = (Kmer &)buffer;
	      temp = other;
	      other = *this;
	      *this = temp;
	   }
	   void swap(KmerPtr &other)
	   {
	   	  swap(*(other._me));
	   }
	
	   TwoBitEncoding *getTwoBitSequence() const
	   {
	     return (TwoBitEncoding *)_data();
	   }
	   SequenceLengthType getTwoBitLength() const
	   {
	   	 return KmerSizer::getTwoBitLength();
	   }
	   SequenceLengthType getByteSize() const
	   {
	   	 return KmerSizer::getByteSize();
	   }
	   SequenceLengthType getLength() const
	   {
	   	 return KmerSizer::getSequenceLength();
	   }
	
	   void buildReverseComplement(Kmer &output) const
	   {
	     TwoBitSequence::reverseComplement((TwoBitEncoding*)_data(), (TwoBitEncoding*)output._data(), getLength());
	   }
	
	   std::string toFasta() const
	   {
	      return TwoBitSequence::getFasta(getTwoBitSequence(), getLength());
	   }
	   std::string toFastaFull() const
	   {
	   	  return TwoBitSequence::getFasta(getTwoBitSequence(), getTwoBitLength()*4);
	   }
	   long hash() const
	   {
	     return getHasher()(std::string((const char *)getTwoBitSequence(), getTwoBitLength()));
	   }
  };
   
private:
   Kmer *_me;
   

public:
   KmerPtr():
   _me(NULL)
   {  }
   
   KmerPtr( void *in):
   _me((Kmer *)in)
   { }

   void *get() const  { return _me; };
   Kmer &operator*() const  { return *_me; }   
   Kmer *operator->() const { return _me;  } 
    
   int compare(const KmerPtr &other) {
      return _me->compare(*(other._me));
   }
   bool equals(const KmerPtr &other) {
   	  return compare(other) == 0;
   }
  
   KmerPtr &operator=(const void *right)    { _me = (Kmer *)right; return *this; }
   KmerPtr &operator=(const KmerPtr &right) { _me = right._me ;    return *this; }
   
   KmerPtr  operator+ (unsigned long right) const { return KmerPtr((Kmer *)((char *)_me + right * KmerSizer::getByteSize())); }
   KmerPtr  operator- (unsigned long right) const { return *this + (-right); }
   
   KmerPtr &operator+=(unsigned long right)       { *this = *this + right; return *this; }
   KmerPtr &operator-=(unsigned long right)       { *this = *this - right; return *this; }

   KmerPtr &operator++()           { return *this += 1;}
   KmerPtr operator++(int unused)  { KmerPtr saved = *this; ++(*this); return saved; }

   KmerPtr &operator--()           { return *this -= 1;}
   KmerPtr operator--(int unused)  { KmerPtr saved = *this; --(*this); return saved; }

   const Kmer &operator[](unsigned long index) const { return *(*this + index); }
   Kmer       &operator[](unsigned long index)       { return *(*this + index); }
   
   // cast operator
   operator VoidPtr() { return (VoidPtr)_me ; }
   
   std::string toFasta() const { return _me->toFasta(); }
};

class SolidKmerTag {};
class WeakKmerTag {};

template<typename Tag, typename Value>
class KmerValue { public: Value value; };

static KmerPtr NullKmerPtr(NULL);

template<typename Value = SolidKmerTag>
class KmerArray
{

public:
  
  typedef boost::pool< > Pool;
  typedef std::tr1::shared_ptr<Pool> PoolPtr;
  typedef std::vector< PoolPtr > SizePools;
  typedef Value ValueType;

private:
  KmerPtr _begin;
  unsigned long _size;

  static SizePools &getPools() { static SizePools pools; return pools; }

public:

  static Pool &getPool(unsigned long poolByteSize) {
     SizePools &pools = getPools();
     if (pools.size() <= poolByteSize)
       pools.resize(poolByteSize+1);
     if ( pools[poolByteSize].get() == NULL ) {
       PoolPtr pool(new Pool(poolByteSize));
       pools[poolByteSize] = pool;
     }
     return *(pools[poolByteSize]);
  }

  // frees ALL memory that EVERY pool has ever allocated
  // Use at end of program or before calling KmerSizer::set()...
  static void purgePools() {
  	 SizePools &pools = getPools();
  	 pools.clear();
  }
  // frees unused memory in existing pools
  // use anytime you want
  //static void releasePools() {
  // 	 SizePools &pools = getPools();
  //	 for(SizePools::iterator it = pools.begin() ; it != pools.end(); it++)
  //	   (*it)->release_memory();
  //}

public:
 
  KmerArray(unsigned long size = 0):
   _size(0), _begin(NULL)
  {
    resize(size);
  }

  KmerArray(TwoBitEncoding *twoBit, SequenceLengthType length):
  _size(0),
  _begin(NULL)
  {
    SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
    resize(numKmers);
    build(twoBit,length);
  }

  KmerArray(const KmerArray &copy)
  {
    *this = copy;
  }
   
  ~KmerArray()
  {
     reset();
  }

  KmerArray &operator=(const KmerArray &other)
  {
    if (this == &other)
  	  return *this;
    reset();
    resize(other.size());
    if (size() == 0)
      return *this;
    if (_begin.get() == NULL)
       throw new std::runtime_error("Could not allocate memory");
    
    memcpy(_begin.get(),other._begin.get(),_size*getElementByteSize());
    return *this;
  }

  const KmerPtr::Kmer &operator[](unsigned long index) const
  {
    if (index >= _size)
       throw new std::invalid_argument("attempt to access index greater than size"); 
    return get(index);
  }
  KmerPtr::Kmer &operator[](unsigned long index)
  {
    if (index >= _size)
       throw new std::invalid_argument("attempt to access index greater than size"); 
    return get(index);
  }
  
  ValueType *getValueStart() const {
  	if (size() > 0)
  	  return (ValueType*) (_begin + size() ).get();
  	else
  	  return NULL;
  }
  const ValueType &valueAt(unsigned long index) const
  {
    if (index >= _size)
    {
      throw std::invalid_argument("attempt to access index greater than size");
  	}
    return *( getValueStart() + index );
  }
  ValueType &valueAt(unsigned long index)
  {
    if (index >= _size)
    {
      throw std::invalid_argument("attempt to access index greater than size");
  	}
    return *( getValueStart() + index );
  }

  const KmerPtr::Kmer &get(unsigned long index) const
  {
    return *(_begin + index);
  }
  KmerPtr::Kmer &get(unsigned long index)
  {
    return *(_begin + index);
  }

  unsigned int size() const { return _size; }

  unsigned int getElementByteSize() const { 
  	return (KmerSizer::getByteSize() + sizeof(ValueType)); 
  }
  
  void reset()
  {
    void *test = (void *)(_begin.get());
    if (test != NULL) {
      getPool( size() * getElementByteSize() ).free(test); 
      //free(test);
    } 
    _begin = NULL;
    _size = 0;
  }

  void resize(unsigned long size) {
  	resize(size, -1);
  }
  void resize(unsigned long size, unsigned long idx)
  {
    if (size == _size)
      return;
    unsigned long oldSize = _size;
    
    // alloc / realloc memory
    _setMemory(size, idx);
      
    if(_begin.get() == NULL) {
       throw new std::runtime_error("Could not allocate memory");
    }

    if (size > oldSize && idx == -1) {
       // zero fill remainder
       char *start = (char*) (_begin + oldSize).get();
       memset(start, 0, KmerSizer::getByteSize()*(size-oldSize));
       start = (char*) (getValueStart() + oldSize);
       memset(start, 0, sizeof(ValueType) * (size - oldSize));
    }
    
  }

  void _setMemory(unsigned long size, unsigned long idx)
  {
    void *old = _begin.get();
    void *memory = NULL;
    
    if (size != 0 ) {
    	// allocate new memory
        boost::pool<> &newPool = getPool( size * getElementByteSize() );
        memory = newPool.malloc();
    	//memory = malloc( size * getElementByteSize() );

        if(memory == NULL) {
           throw new std::runtime_error("Could not allocate memory");
        }        
    }
    unsigned long oldSize = _size;
    unsigned long lesserSize = std::min(size,oldSize);
    ValueType *oldValue;
    if (oldSize > 0) {
    	oldValue = getValueStart();
    }
    KmerPtr oldBegin(_begin.get());
    _begin = KmerPtr( memory );
    _size = size;
     
    // TODO refactor this block
    if (old != NULL && memory != NULL && oldValue != NULL && lesserSize > 0) {
      // copy the old contents
      if (idx == -1 || idx >= lesserSize) {
      	// copy all records in order (default ; end is trimmed or expanded)
        memcpy(memory, old, lesserSize*KmerSizer::getByteSize());
        memcpy(getValueStart(), oldValue, lesserSize*sizeof(ValueType));
      } else {
      	
      	// shrink or expand the first records will be copied
      	if (idx > 0) {
      	  memcpy(memory, old, (idx)*KmerSizer::getByteSize());
          memcpy(getValueStart(), oldValue, (idx)*sizeof(ValueType));
      	}
      	
      	if (lesserSize == size) {
      	  // shrink, removing record at idx
      	  if (idx < lesserSize-1) {
      	  	void *memory2 = _begin[idx].get();
      	  	void *old2 = oldBegin[idx+1].get();
      	  	unsigned long length = (lesserSize-idx)*KmerSizer::getByteSize();
      	  	memcpy(memory2, old2, length);
      	  	
      	  	length = (lesserSize-idx)*sizeof(ValueType);
            memcpy(getValueStart()+idx, oldValue+idx+1, length);
      	  }
      	} else {
      	  // expand, leaving new (uninitialized) record at idx
      	  if (idx < lesserSize) {
      	  	void *memory2 = _begin[idx+1].get();
      	  	void *old2 = oldBegin[idx].get();
      	  	unsigned long length = (lesserSize-idx)*KmerSizer::getByteSize();
      	  	memcpy(memory2, old2, length);
      	  	
      	  	length = (lesserSize-idx)*sizeof(ValueType);
            memcpy(getValueStart()+idx+1, oldValue+idx, length);
      	  }
      	}
      }
    }
    if (old != NULL) {
      // free old memory
      boost::pool<> &oldPool = getPool( oldSize * getElementByteSize() );
      oldPool.free(old);
      //free(old);
    }
  }
    


  void build(TwoBitEncoding *twoBit, SequenceLengthType length)
  {
    SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
    if (_size != numKmers)
      throw new std::invalid_argument("attempt to build an incorrectly sized KmerArray"); ;

    KmerArray &kmers = *this;
    for(SequenceLengthType i=0; i < numKmers ; i+=4) {
      TwoBitEncoding *ref = twoBit+i/4;
      for (int bitShift=0; bitShift < 4 && i+bitShift < numKmers; bitShift++) {
        TwoBitSequence::shiftLeft(ref, kmers[i+bitShift].get(), KmerSizer::getTwoBitLength(), bitShift, bitShift != 0);
        TwoBitEncoding *lastByte = kmers[i+bitShift].getTwoBitSequence()+KmerSizer::getTwoBitLength()-1;
        switch (KmerSizer::getSequenceLength() % 4) {
          case 1: *lastByte &= 0xc0; break;
          case 2: *lastByte &= 0xf0; break;
          case 3: *lastByte &= 0xfc; break;
        }
      }
    }
  }
  unsigned long find(const KmerPtr &target) const {
  	return find(*target);
  }
  unsigned long find(const KmerPtr::Kmer &target) const {
    for(unsigned long i=0; i<_size; i++)
      if (target.compare(_begin[i]) == 0)
        return i;
    return -1;
  }
  unsigned long findSorted(const KmerPtr &target, bool &targetIsFound) const {
  	return findSorted(*target, targetIsFound);
  }
  unsigned long findSorted(const KmerPtr::Kmer &target, bool &targetIsFound) const {
  	// binary search
  	unsigned long min = 0;
  	unsigned long max = size();
 
  	if (max == 0)
    {
       targetIsFound = false;
       return 0;
    }
  	unsigned long mid;
  	int comp;
  	do {
  		mid = (min+max) / 2;
  		comp = target.compare(_begin[mid]);
  		if (comp > 0)
  		  min = mid+1;
  		else if (comp < 0)
  		  max = mid-1;
  	} while (comp != 0 && max != -1 && min <= max);
  	if (comp == 0)
  	  targetIsFound = true;
  	else
  	  targetIsFound = false;
  	return mid + (comp>0 && size()>mid?1:0);
  }
  void insertAt(unsigned long idx, const KmerPtr &target) {
  	insertAt(idx, *target);
  }
  void insertAt(unsigned long idx, const KmerPtr::Kmer &target) {
  	if (idx > size())
  	  throw new std::invalid_argument("attempt to access index greater than size");
  	resize(size() + 1, idx);
  	_begin[idx] = target;
  }
  unsigned long append(const KmerPtr &target) {
    return append(*target);	
  }
  unsigned long append(const KmerPtr::Kmer &target) {
  	unsigned long idx = size();
   	insertAt(size(), target);
   	return idx;
  }
  unsigned long insertSorted(const KmerPtr &target) {
  	return insertSorted(*target);
  }
  unsigned long insertSorted(const KmerPtr::Kmer &target) {
  	bool isFound;
  	unsigned long idx = findSorted(target, isFound);
  	if (!isFound)
  	  insertAt(idx, target);
  	return idx;
  }
  void remove(unsigned long idx) {
    resize(size()-1,idx);
  }
  void swap(unsigned long idx1, unsigned long idx2) {
  	if (idx1 == idx2)
  	  return;
  	if (idx1 >= size() || idx2 >= size())
  	  throw new std::invalid_argument("attempt to access index greater than size");
  	  
  	get(idx1).swap(get(idx2));
  	if (sizeof(ValueType) > 0) {
  	  ValueType tmp = valueAt(idx1);
  	  valueAt(idx1) = valueAt(idx2);
  	  valueAt(idx2) = tmp; 
  	}
  }

  std::string toString() {
  	std::stringstream ss;
  	ss <<  "{";
  	for(unsigned long idx=0; idx<size(); idx++) {
  		ss << get(idx).toFasta() << ":" << valueAt(idx) << ", ";
  	} 
  	ss << "}";
  	return ss.str();
  }
};



template<typename Value>
class KmerMap
{

public:
   typedef KmerPtr::Kmer KeyType;
   typedef Value ValueType;
   typedef KmerArray<Value> BucketType;
   typedef std::vector< BucketType > BucketsVector;

private:
   BucketsVector _buckets;
   
public:
   KmerMap(unsigned long bucketCount = 1024*1024) {
     _buckets.resize(bucketCount);
   }
   ~KmerMap() 
   {
   	 clear();
   }
   
   void clear() {
     for(int i=0; i< _buckets.size(); i++)
       _buckets[i].reset();
   }
   BucketType &getBucket(long hash) {
   	return _buckets[hash % _buckets.size()];
   }
   const BucketType &getBucket(long hash) const {
   	return _buckets[hash % _buckets.size()];
   }
   BucketType &getBucket(const KmerPtr &key) {
     return getBucket(*key);
   }
   const BucketType &getBucket(const KmerPtr &key) const {
   	 return getBucket(*key);
   }
   BucketType &getBucket(const KeyType &key)  {
     return getBucket(key.hash());
   }
   const BucketType &getBucket(const KeyType &key) const {
     return getBucket(key.hash());
   }

   ValueType &insert(const KmerPtr &key, const ValueType &value, BucketType *bucketPtr = NULL) {
     return insesrt(*key, value, bucketPtr);
   }
   ValueType &insert(const KeyType &key, const ValueType &value, BucketType *bucketPtr = NULL) {
   	  if (bucketPtr == NULL)
   	    bucketPtr = &getBucket(key);
   	    
   	  unsigned long idx = bucketPtr->insertSorted(key);
   	  return bucketPtr->valueAt(idx) = value;
   }
   
   bool remove(const KmerPtr &key, BucketType *bucketPtr = NULL) {
   	  return remove(*key, bucketPtr);
   }
   bool remove(const KeyType &key, BucketType *bucketPtr = NULL) {
   	  if (bucketPtr == NULL)
   	    bucketPtr = &getBucket(key);
   	  bool isFound;
   	  unsigned long idx = bucketPtr->findSorted(key, isFound);
   	  if (isFound && idx != -1)
   	    bucketPtr->remove(idx);
   	  return isFound;
   }
   
   bool exists(const KmerPtr &key, BucketType *bucketPtr = NULL) const {
     return exists(*key, bucketPtr);
   }
   bool exists(const KeyType &key, BucketType *_bucketPtr = NULL) const {
   	 const BucketType *bucketPtr = _bucketPtr;
     if (bucketPtr == NULL)
   	   bucketPtr = &getBucket(key);
   	 bool isFound;
   	 bucketPtr->findSorted(key, isFound);
   	 return isFound;
   }
   
   ValueType &operator[](const KmerPtr &key) {
   	 return operator[](*key);
   }
   ValueType &operator[](const KeyType &key) {
     BucketType &bucket = getBucket(key);
     bool isFound;
     unsigned long idx = bucket.findSorted(key, isFound);
     if (isFound && idx != -1)
       return bucket.valueAt(idx);
     else 
       return insert(key, Value(), &bucket);       
   }
   
   unsigned long size() const {
   	unsigned long size = 0;
   	for(int i = 0; i<_buckets.size() ; i++)
   	  size += _buckets[i].size();
   	return size;
   } 
   
   std::string toString() {
  	std::stringstream ss;
  	ss << this << "[";
  	for(unsigned long idx=0; idx<_buckets.size(); idx++) {
  		ss << "bucket:" << idx << ' ' << _buckets[idx].toString() << ", ";
  	} 
  	ss << "]";
  	return ss.str();
  }

};

#endif


 

//
// $Log: Kmer.h,v $
// Revision 1.33  2009-10-30 00:51:40  regan
// bug fix and working on executable
//
// Revision 1.32  2009-10-30 00:10:32  regan
// cleaned up a bit
//
// Revision 1.31  2009-10-30 00:07:59  regan
// bugfix on KmerArray.build trailing bits
//
// Revision 1.30  2009-10-29 23:30:01  regan
// checkpoint
//
// Revision 1.29  2009-10-29 23:04:49  regan
// works
//
// Revision 1.28  2009-10-29 20:59:23  cfurman
// fixed testing bugs
//
// Revision 1.27  2009-10-29 19:01:33  regan
// checkpoint
//
// Revision 1.26  2009-10-29 17:00:58  regan
// checkpoint (with bugs)
//
// Revision 1.25  2009-10-29 07:03:33  regan
// fixed some bugs , added others
// KmerArray is working, *Sorted methods are untested
//
// Revision 1.24  2009-10-28 18:50:57  regan
// made KmerArray behave properly and not like a KmerPtrArray
//
// Revision 1.23  2009-10-28 18:42:59  regan
// added debug flags, fixed tests, bugs
//
// Revision 1.22  2009-10-28 02:29:55  cfurman
// fixed KmerArray  bugs
//
// Revision 1.21  2009-10-28 00:00:41  regan
// added more bugs
//
// Revision 1.20  2009-10-27 07:16:09  regan
// checkpoint
// defined KmerMap and KmerArray lookup methods
//
// Revision 1.19  2009-10-26 23:04:33  regan
// checkpoint make Kmer private inner class
//
// Revision 1.18  2009-10-26 17:50:54  regan
// templated KmerArray; added boost pool allocation
//
// Revision 1.17  2009-10-24 00:32:46  regan
// added bugs
//
// Revision 1.16  2009-10-24 00:03:49  regan
// checkpoint
//
// Revision 1.15  2009-10-23 23:22:41  regan
// checkpoint
//
// Revision 1.14  2009-10-23 21:54:46  regan
// checkpoint
//
// Revision 1.13  2009-10-23 20:32:50  cfurman
// more kmer changes
//
// Revision 1.12  2009-10-23 17:22:39  regan
// added more tests
//
// Revision 1.11  2009-10-23 07:06:59  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
// Revision 1.10  2009-10-23 01:24:53  cfurman
// ReadSet test created
//
// Revision 1.9  2009-10-22 01:39:43  cfurman
// bug fix in kmer.h
//
// Revision 1.8  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
// Revision 1.7  2009-10-21 18:58:44  regan
// checkpoint
//
// Revision 1.6  2009-10-21 18:44:20  regan
// checkpoint
//
// Revision 1.5  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
//
