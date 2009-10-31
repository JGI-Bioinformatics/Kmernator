// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.37 2009-10-31 23:44:17 regan Exp $
//

#ifndef _KMER_H
#define _KMER_H
#include <tr1/memory>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <stdexcept>

#include <boost/functional/hash.hpp>

#include "TwoBitSequence.h"
#include "MemoryUtils.h"

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 1024
#endif

typedef std::tr1::shared_ptr<TwoBitEncoding> KmerSharedPtr;

class KmerSizer
{
private:
  KmerSizer() : _sequenceLength(21) {}
  static KmerSizer singleton;

  SequenceLengthType _sequenceLength;

  SequenceLengthType _twoBitLength;
  unsigned long _totalSize;
public:

  static inline KmerSizer &getSingleton() { return singleton; }
  static void set(SequenceLengthType sequenceLength, unsigned long extraBytes=0)
  {
    KmerSizer &singleton = getSingleton();
    singleton._sequenceLength = sequenceLength;
    singleton._twoBitLength =  TwoBitSequence::fastaLengthToTwoBitLength(singleton._sequenceLength);
    singleton._totalSize = singleton._twoBitLength;
  }

  static inline SequenceLengthType getSequenceLength()  {
    return getSingleton()._sequenceLength;
  }
  static inline SequenceLengthType getTwoBitLength()   {
    return getSingleton()._twoBitLength;
  }
  static inline unsigned long getByteSize()  {
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


template<typename Value>
class KmerValue { public: Value value; };

typedef KmerValue<unsigned char> KmerValueByte;
typedef KmerValue<unsigned short> KmerValue2Byte;
typedef KmerValue<unsigned int> KmerValue4Byte;


class SolidKmerTag : public KmerValue<unsigned short> {};
class WeakKmerTag  : public KmerValue<unsigned char> {};

static KmerPtr NullKmerPtr(NULL);

template<typename Value = SolidKmerTag>
class KmerArray
{

public:
  
  typedef Value ValueType;
  class ElementType { 
  private:
  	  KmerPtr _key; 
  	  ValueType *_value;
  public:
      ElementType(): _key(NULL), _value(NULL) {}
  	  ElementType(KmerPtr::Kmer &key, ValueType &value): _key(&key), _value(&value) {}
  	  ElementType(const ElementType &copy) {
  	  	*this = copy;
  	  }
  	  ~ElementType() {}
  	  ElementType &operator=(const ElementType &other) {
  	  	_key = other._key;
  	  	_value = other._value;
  	  	return *this;
  	  }
  	  KmerPtr::Kmer &key() { return *_key; }
  	  ValueType &value()   { return *_value; }
  };

private:
  KmerPtr _begin;
  unsigned long _size;

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
       throw std::runtime_error("Could not allocate memory in KmerArray operator=()");
    
    memcpy(_begin.get(),other._begin.get(),_size*getElementByteSize());
    return *this;
  }

  const KmerPtr::Kmer &operator[](unsigned long index) const
  {
    if (index >= _size)
       throw std::invalid_argument("attempt to access index greater than size in KmerArray operator[] const"); 
    return get(index);
  }
  KmerPtr::Kmer &operator[](unsigned long index)
  {
    if (index >= _size)
       throw std::invalid_argument("attempt to access index greater than size in KmerArray operator[]"); 
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
      throw std::invalid_argument("attempt to access index greater than size in KmerArray valueAt() const");
  	}
    return *( getValueStart() + index );
  }
  ValueType &valueAt(unsigned long index)
  {
    if (index >= _size)
    {
      throw std::invalid_argument("attempt to access index greater than size in KmerArray valueAt()");
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

  const ElementType getElement(unsigned long idx) const
  {
  	return ElementType( get(idx), valueAt(idx) );
  }
  ElementType getElement(unsigned long idx)
  {
  	return ElementType( get(idx), valueAt(idx) );
  }
  
  unsigned int size() const { return _size; }

  unsigned int getElementByteSize() const { 
  	return (KmerSizer::getByteSize() + sizeof(ValueType)); 
  }
  
  void reset(PoolManager &pool = BoostPoolManager::get())
  {
    void *test = (void *)(_begin.get());
    if (test != NULL) {
      pool.free(test, size() * getElementByteSize());
      //getPool( size() * getElementByteSize() ).free(test); 
      //free(test);
    } 
    _begin = NULL;
    _size = 0;
  }

  void resize(unsigned long size) {
  	resize(size, -1);
  }
  void resize(unsigned long size, unsigned long idx, PoolManager &pool = BoostPoolManager::get())
  {
    if (size == _size)
      return;
    unsigned long oldSize = _size;
    
    // alloc / realloc memory
    _setMemory(size, idx, pool);
      
    if(_begin.get() == NULL) {
       throw std::runtime_error("Could not allocate memory in KmrArray resize()");
    }

    if (size > oldSize && idx == -1) {
       // zero fill remainder
       char *start = (char*) (_begin + oldSize).get();
       memset(start, 0, KmerSizer::getByteSize()*(size-oldSize));
       start = (char*) (getValueStart() + oldSize);
       memset(start, 0, sizeof(ValueType) * (size - oldSize));
    }    
  }
    
  void _copyRange(KmerPtr &srcKmer, ValueType *srcValue, unsigned long idx, unsigned long srcIdx, unsigned long count) {
	memcpy((_begin + idx).get(), (srcKmer + srcIdx).get(), count * KmerSizer::getByteSize()); 
	memcpy(getValueStart() + idx, srcValue + srcIdx, count * sizeof(ValueType));
  }

  void _setMemory(unsigned long size, unsigned long idx, PoolManager &pool = BoostPoolManager::get())
  {
    void *oldMemory = _begin.get();
    void *memory    = NULL;
    
    KmerPtr oldBegin(_begin.get());
    
    if (size != 0 ) {
    	// allocate new memory
        //boost::pool<> &newPool = getPool( size * getElementByteSize() );
        //memory = newPool.ordered_malloc();
    	//memory = malloc( size * getElementByteSize() );
    	memory = pool.malloc( size * getElementByteSize() );

        if(memory == NULL) {
           throw std::runtime_error("Could not allocate memory in KmerArray _setMemory()");
        }        
    }
    unsigned long oldSize = _size;
    unsigned long lesserSize = std::min(size,oldSize);
    ValueType *oldValueStart;
    if (oldSize > 0) {
    	oldValueStart = getValueStart();
    }
    _begin = KmerPtr( memory );
    _size = size;
    
    if (oldMemory != NULL && memory != NULL && oldValueStart != NULL && lesserSize > 0) {
      // copy the old contents
      if (idx == -1 || idx >= lesserSize) {
      	// copy all records in order (default ; end is trimmed or expanded)
      	_copyRange(oldBegin,oldValueStart, 0, 0, lesserSize);
      } else {      	
      	// shrink or expand.
      	
      	// the first record(s) leading to idx will be copied
      	if (idx > 0) {
      	  _copyRange(oldBegin, oldValueStart, 0, 0, idx);
      	}
      	
      	if (lesserSize == size) {
      	  // shrink: skipping the old record at idx
      	  if (idx < size) {
      	  	_copyRange(oldBegin, oldValueStart, idx, idx+1, lesserSize-idx);
      	  }
      	} else {
      	  // expand: leaving new (uninitialized) record at idx
      	  if (idx < oldSize) {
      	  	_copyRange(oldBegin, oldValueStart, idx+1, idx, lesserSize-idx);
      	  }
      	}
      }
    }
    if (oldMemory != NULL) {
      // free old memory
      //boost::pool<> &oldPool = getPool( oldSize * getElementByteSize() );
      //oldPool.free(old);
      //free(old);
      pool.free(oldMemory, oldSize * getElementByteSize());
    }
  }
    


  void build(TwoBitEncoding *twoBit, SequenceLengthType length)
  {
    SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
    if (_size != numKmers)
      throw std::invalid_argument("attempt to build an incorrectly sized KmerArray in KmerArray build()"); ;

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
  	  throw std::invalid_argument("attempt to access index greater than size in KmerArray insertAt");
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
  	  throw std::invalid_argument("attempt to access index greater than size in KmerArray swap()");
  	  
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

public:
  class Iterator : public std::iterator<std::forward_iterator_tag, KmerArray>
  {
    private:
      KmerArray *_tgt;
      unsigned long _idx;
      ElementType elementCopy;
    public:
      Iterator(KmerArray *target, unsigned long idx = 0): _tgt(target), _idx(idx) { }
      Iterator(const Iterator &copy) { *this = copy; }
      ~Iterator() {}
      Iterator& operator=(const Iterator& other) {
        _tgt = other._tgt;
        _idx = other._idx;
        return *this;
      }
      bool operator==(const Iterator& other) { return _idx == other._idx && _tgt == other._tgt; }
      bool operator!=(const Iterator& other) { return _idx != other._idx || _tgt != other._tgt; }
      Iterator& operator++() { if ( !isEnd() ) _idx++; return *this; }
      Iterator operator++(int unused) { Iterator tmp(*this) ; ++(*this); return tmp; }
      ElementType &operator*() { return elementCopy = _tgt->getElement(_idx); }
      KmerPtr::Kmer &key() { return _tgt->get(_idx); }
      Value &value() { return _tgt->valueAt(_idx); }
      ElementType *operator->() { elementCopy = _tgt->getElement(_idx); return &elementCopy; }
      bool isEnd() { return _idx == _tgt->size(); }
  };

  Iterator begin() { return Iterator(this, 0); }
  Iterator end()   { return Iterator(this, size()); }

};



template<typename Value>
class KmerMap
{

public:
   typedef KmerPtr::Kmer KeyType;
   typedef Value ValueType;
   typedef KmerArray<Value> BucketType;
   typedef typename BucketType::Iterator BucketTypeIterator;
   typedef typename BucketType::ElementType ElementType;
   
   typedef std::vector< BucketType > BucketsVector;
   typedef typename BucketsVector::iterator BucketsVectorIterator;

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
     //BucketType::releasePools();
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
     return insert(*key, value, bucketPtr);
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

public:
  class Iterator : public std::iterator<std::forward_iterator_tag, KmerMap>
  {
    private:
      KmerMap *_tgt;
      BucketsVectorIterator _pos; 
      BucketTypeIterator    _pos2;
    public:
      Iterator(KmerMap *map, const BucketsVectorIterator &pos, unsigned long idx = 0): _tgt(map), _pos(pos), _pos2( (*pos).begin() ) {
        if ( ! isEnd() ) {
          setPos2(*_pos, idx);
          if ( _pos2.isEnd() )
            ++(*this);
        }
      }
      Iterator(const Iterator &copy): _pos2( (*copy._pos).begin() ){
        *this = copy;
      }
      ~Iterator() {}
      void setPos2(BucketType &bucket, unsigned long idx = 0) {
        _pos2 = bucket.begin();
        for(unsigned long i=0; i<idx; i++)
          _pos2++;
      }
      Iterator& operator=(const Iterator& other) {
      	_tgt = other._tgt;
        _pos = other._pos;
        _pos2 = other._pos2;
        return *this;
      }
      bool operator==(const Iterator& other) { return _pos == other._pos && _pos2 == other._pos2; }
      bool operator!=(const Iterator& other) { return _pos != other._pos || _pos2 != other._pos2; }
      Iterator& operator++() { 
        bool movedBucket = false;
        while ( !isEnd() ) {
          if ( _pos2.isEnd() ) {
            setPos2(*(++_pos));
            movedBucket = true;
          } else if (movedBucket) {
            break;
          } else {
            _pos2++;
            if ( !_pos2.isEnd() )
              break;
          }
        }
        return *this;
      }
      Iterator operator++(int unused) { Iterator tmp(*this) ; ++(*this); return tmp; }
      ElementType &operator*() { return *_pos2; }
      KmerPtr::Kmer &key() {return _pos2.key(); }
      Value &value() { return _pos2.value(); }
      BucketType &bucket() { return *_pos; }
      unsigned long bucketIndex() { return (_pos - _tgt->_buckets.begin()); }
      ElementType *operator->() { return &(*_pos2); }
      bool isEnd() { return _pos == _tgt->_buckets.end(); }
  };

  Iterator begin() { return Iterator( this, _buckets.begin(), 0); }
  Iterator end()   { return Iterator( this, _buckets.end()  , 0); }

};

typedef KmerMap<SolidKmerTag>   KmerSolidMap;
typedef KmerMap<WeakKmerTag>    KmerWeakMap;
typedef KmerMap<unsigned short> KmerCountMap;

#endif


 

//
// $Log: Kmer.h,v $
// Revision 1.37  2009-10-31 23:44:17  regan
// fixed bug in KmerArray::remove
// refactored memory pool out of KmerArray
//
// Revision 1.36  2009-10-31 00:16:35  regan
// minor changes and optimizations
//
// Revision 1.35  2009-10-30 20:56:27  regan
// fixed kmermap iterator
//
// Revision 1.34  2009-10-30 19:27:46  regan
// added iterator goodness, but KmerMap::Iterator still does not work
//
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
