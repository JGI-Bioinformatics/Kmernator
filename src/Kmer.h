// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.20 2009-10-27 07:16:09 regan Exp $
//

#ifndef _KMER_H
#define _KMER_H
#include <tr1/memory>
#include <cstring>
#include <cstdlib>
#include <vector>
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
	   //Kmer(); // never construct, just use as cast
	
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
	   
	   int compare(const Kmer &other) const
	   {
	     return memcmp(_data(), other._data(), getTwoBitLength());
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
	   long hash() const
	   {
	   	  return getHasher()(std::string((const char *)getTwoBitSequence(), getLength()));
	   }
  };
   
private:
   Kmer *_me;
   

public:
   KmerPtr():
   _me(NULL)
   {  }
   
   KmerPtr(Kmer &in):
   _me((&in)._me)
   { }

   KmerPtr( void *in):
   _me((Kmer *)in)
   { }

   void *get() const  { return _me; };
   void *operator*() const  { return _me; }   
   Kmer *operator->() const { return _me;  } 
    
   int compare(const KmerPtr &other) {
      return _me->compare(*(other._me));
   }
   bool equals(const KmerPtr &other) {
   	  return compare(other) == 0;
   }
  
   KmerPtr &operator=(const void *right)    { _me = (Kmer *)right; return *this; }
   KmerPtr &operator=(const KmerPtr &right) { _me = right._me ;    return *this; }
   
   //bool operator==(const KmerPtr &right) const { return _me == right._me; }
   //bool operator!=(const KmerPtr &right) const { return _me != right._me; }
   
   KmerPtr  operator+ (unsigned long right) const { return KmerPtr((Kmer *)((char *)_me + right * KmerSizer::getByteSize())); }
   KmerPtr  operator- (unsigned long right) const { return *this + (-right); }
   
   KmerPtr &operator+=(unsigned long right)       { *this = *this + right; return *this; }
   KmerPtr &operator-=(unsigned long right)       { *this = *this - right; return *this; }

   KmerPtr &operator++()           { return *this += 1;}
   KmerPtr operator++(int unused)  { KmerPtr saved = *this; ++(*this); return saved; }

   KmerPtr &operator--()           { return *this -= 1;}
   KmerPtr operator--(int unused)  { KmerPtr saved = *this; --(*this); return saved; }

   const KmerPtr operator[](unsigned long index) const { return (*this + index); }
   KmerPtr       operator[](unsigned long index)       { return (*this + index); }
   
   // cast operator
   operator VoidPtr() { return (VoidPtr)_me ; }
   
   std::string toFasta() const { return _me->toFasta(); }
};

   
/* class KmerInstance : public KmerPtr
{

private:
   TwoBitEncoding _somedata[1024];
public:

   KmerInstance()
   {
   	 if ((void*)this != (void*) _somedata)
   	   throw;
   }
   KmerInstance &operator=(const KmerPtr &other)
   {
   	  if (this == &other)
   	    return *this;
   	    
      memcpy(_data(), other._data(), getTwoBitLength());
      return *this;
   }
   KmerInstance &operator=(const KmerPtr &other);
};
 */

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

  static void purgePools() {
  	 SizePools &pools = getPools();
  	 pools.clear();
  }
  static void releasePools() {
  	 SizePools &pools = getPools();
  	 for(SizePools::iterator it = pools.begin() ; it != pools.end(); it++)
  	   (*it)->release_memory();
  }

public:
 
  KmerArray(unsigned long size = 0):
   _size(0), _begin(NULL)
  {
    resize(size);
    if (_begin.get() == NULL)
       throw;
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
    memcpy(_begin.get(),other._begin.get(),_size*getElementByteSize());
    return *this;
  }

  const KmerPtr operator[](unsigned long index) const
  {
    if (index >= _size)
       throw; 
    return _begin + index;
  }
  KmerPtr operator[](unsigned long index)
  {
    if (index >= _size)
       throw; 
    return _begin + index;
  }
  
  ValueType *valueAt(unsigned long index) const
  {
  	if (sizeof(ValueType) > 0 && _size > 0)
      return ((ValueType*) (_begin + _size ).get()) + index;
    else
      return NULL;
  }

  const KmerPtr get(unsigned long index) const
  {
    return _begin + index;
  }

  unsigned int size() const { return _size; }

  unsigned int getElementByteSize() const { 
  	return (KmerSizer::getByteSize() + sizeof(ValueType)); 
  }
  
  void reset()
  {
    void *test = (void *)_begin.get();
    if (test != NULL) {
      getPool( size() * getElementByteSize() ).free(test); 
    }
    _begin = NULL;
    _size = 0;
  }

  void resize(unsigned long size)
  {
    if (size == _size)
      return;
    unsigned long oldSize = _size;
    
    // alloc / realloc memory
    _setMemory(size);
      
    if(_begin.get() == NULL) {
       throw;
    }

    if (size > oldSize) {
       // zero fill remainder
       char *start = (char*) (_begin + oldSize).get();
       memset(start, 0, KmerSizer::getByteSize()*(size-oldSize));
       if (sizeof(ValueType) > 0) {
         start = (char*) (valueAt(0) + oldSize);
         memset(start, 0, sizeof(ValueType) * (size - oldSize));
       }
    }
  }

  void _setMemory(unsigned long size)
  {
    void *old = _begin.get();

    boost::pool<> &newPool = getPool( size * getElementByteSize() );

    void *memory = newPool.malloc();
    if(memory == NULL) {
       throw;
    }
    unsigned long oldSize = _size;
    ValueType *oldValue = valueAt(0);
    _begin = KmerPtr( memory );
    _size = size;
     
    if (old != NULL && _size > 0) {      
      // copy the old contents
      memcpy(memory, old, (size < oldSize ? size : oldSize)*KmerSizer::getByteSize());
      if (sizeof(ValueType)>0) 
        memcpy(valueAt(0), oldValue, (size < oldSize ? size : oldSize)*sizeof(ValueType));
      boost::pool<> &oldPool = getPool( _size * getElementByteSize() );
      oldPool.free(old);
    }
  }
    
  KmerArray(TwoBitEncoding *twoBit, SequenceLengthType length):
  _size(0),
  _begin(NULL)
  {
    SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
    resize(numKmers);
    build(twoBit,length);
  }

  void build(TwoBitEncoding *twoBit, SequenceLengthType length)
  {
    SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
    if (_size != numKmers)
      throw;

    KmerArray &kmers = *this;
    for(SequenceLengthType i=0; i < numKmers ; i+=4) {
      TwoBitEncoding *ref = twoBit+i/4;
      for (int bitShift=0; bitShift < 4 && i+bitShift < numKmers; bitShift++)
        TwoBitSequence::shiftLeft(ref, kmers[i+bitShift].get(), KmerSizer::getTwoBitLength(), bitShift, bitShift != 0);
    }
  }
  
  unsigned long find(const KmerPtr &target) {
    for(unsigned long i=0; i<_size; i++)
      if (target.compare(_begin[i]) == 0)
        return i;
    return -1;
  }
  unsigned long insert(const KmerPtr &target) {
   	unsigned long idx = size();
   	resize(idx + 1);
  	_begin[idx] = target;
  	return idx; 
  }
  void erase(unsigned long idx) {
    swap(idx, size()-1);
    resize(size()-1);
  }
  void swap(unsigned long idx1, unsigned long idx2) {
  	if (idx1 == idx2)
  	  return;
  	if (idx1 >= size() || idx2 >= size())
  	  throw;
  	  
  	get(idx1)->swap(get(idx2));
  	if (sizeof(ValueType) > 0) {
  	  ValueType tmp = *valueAt(idx1);
  	  *valueAt(idx1) = *valueAt(idx2);
  	  *valueAt(idx2) = tmp; 
  	}
  }

};



template<typename Value>
class KmerMap
{

public:
   typedef KmerPtr KeyType;
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
   	 //for(BucketsVector::iterator iter = _buckets.begin(); iter != _buckets.end(); iter++)
   	 //  iter->reset();
   	 for(int i=0; i< _buckets.size(); i++)
   	   _buckets[i].reset();
   	 _buckets.clear();
   	 BucketType::releaseMemory();
   }
   
   BucketType &getBucket(const KeyType &key) {
     return _buckets[key->hash() % _buckets.size()];
   }

   ValueType &insert(const KeyType &key, const ValueType &value, BucketType *bucketPtr = NULL) {
   	  if (bucketPtr == NULL)
   	    bucketPtr = &getBucket(key);
   	    
   	  unsigned long idx = bucketPtr->insert(key);
   	  ValueType *val = bucketPtr->valueAt(idx);
   	  *val = value;

   }
   
   void erase(const KeyType &key, BucketType *bucketPtr = NULL) {
   	  if (bucketPtr == NULL)
   	    bucketPtr = &getBucket(key);
   	  unsigned long idx = bucketPtr->find(key);
   	  if (idx != -1)
   	    bucketPtr->erase(idx);
   }
   
   ValueType &operator[](const KmerPtr &key) {
     BucketType &bucket = getBucket(key);
     unsigned long idx = bucket.find(key);
     if (idx == -1)
       return insert(key, Value(), &bucket);
     else 
       return *( bucket.valueAt(idx) );
   }
    

};

#endif


 

//
// $Log: Kmer.h,v $
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
