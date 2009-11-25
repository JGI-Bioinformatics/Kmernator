// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.54 2009-11-25 18:39:08 regan Exp $
//

#ifndef _KMER_H
#define _KMER_H
#include <tr1/memory>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>

#include <boost/functional/hash.hpp>

#include "TwoBitSequence.h"
#include "MemoryUtils.h"

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 1024
#endif

typedef std::tr1::shared_ptr<TwoBitEncoding> KmerSharedPtr;
typedef unsigned long IndexType;

const IndexType MAX_INDEX = (IndexType) -1;

class KmerSizer
{
private:
  KmerSizer() : _sequenceLength(21) {}
  static KmerSizer singleton;

  SequenceLengthType _sequenceLength;

  SequenceLengthType _twoBitLength;
  IndexType _totalSize;
public:

  static inline KmerSizer &getSingleton() { return singleton; }
  static void set(SequenceLengthType sequenceLength, IndexType extraBytes=0)
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
  static inline IndexType getByteSize()  {
    return getSingleton()._totalSize;
  }
};


#define TEMP_KMER(name)  TwoBitEncoding _stack_##name[KmerSizer::getByteSize()]; Kmer &name = (Kmer &)(_stack_##name); 


class Kmer
	{
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
          TEMP_KMER(temp);
	      temp = other;
	      other = *this;
	      *this = temp;
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
	   
	   bool buildLeastComplement(Kmer &output) const
	   {
	   	 buildReverseComplement(output);
	   	 if ( *this <= output ) {
	   	   output = *this;
	   	   return true;
	   	 } else {
	   	   return false;
	   	 }
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
	   // check for trivial patterns AAAAA... GGGGG.... etc
	   bool isTrivial() const
	   {
	     const TwoBitEncoding *firstByte = (const TwoBitEncoding *) this;
	     if (*firstByte == 0x00 || *firstByte == 0x55 || *firstByte == 0xaa || *firstByte == 0xff) {
	     	const TwoBitEncoding *nextByte = firstByte;
	     	for(int i=1; i<getTwoBitLength()-1; i++)
	     	  if (*firstByte != *(++nextByte))
	     	    return false;
	     	return true;
	     } else
	        return false;
	   }
  };

template<typename Value>
class KmerValue { 
public: 
  Value value;
  KmerValue() : value() {}
  ~KmerValue() {}
  KmerValue(const KmerValue &copy) {
  	*this = copy;
  }
  KmerValue &operator=(const KmerValue &other) {
  	this->value = other.value;
  }
  inline bool operator==(const KmerValue &other) const {
  	return this->value == other.value;
  }
  inline bool operator<(const KmerValue &other) const {
  	return this->value < other.value;
  }
  inline bool operator>(const KmerValue &other) const {
  	return this->value > other.value;
  }
  // cast operator
  operator Value() { return value; }
};

template<typename Value>
std::ostream &operator<<(std::ostream &stream, KmerValue<Value> &ob)
{
  stream << ob.value;
  return stream;
};


typedef KmerValue<unsigned char> KmerValueByte;
typedef KmerValue<unsigned short> KmerValue2Byte;
typedef KmerValue<unsigned int> KmerValue4Byte;


class TrackingData 
{
public:
  typedef unsigned int   CountType;
  typedef float          WeightType;
  
  static WeightType minimumWeight;
  static CountType minimumDepth;
  static const CountType MAX_COUNT = (CountType) -1;
  static unsigned long discarded;
  static unsigned long singletonCount;
  
  static CountType  maxCount;
  static WeightType maxWeightedCount;
  
  static void resetGlobalCounters() {
  	discarded = 0;
  	singletonCount = 0;
  	maxCount = 0;
  	maxWeightedCount = 0.0;
  }
  
protected:
  CountType           count;
  CountType           directionBias;
  WeightType          weightedCount;

public:
  TrackingData() : count(0), directionBias(0), weightedCount(0.0) {}
  
  void reset() {
  	count = 0;
  	directionBias = 0;
  	weightedCount = 0.0;
  }
  inline bool operator==(const TrackingData &other) const {
  	return count == other.count;// && weightedCount == other.weightedCount;
  }
  inline bool operator<(const TrackingData &other) const {
  	return count < other.count;// || (count == other.count && weightedCount < other.weightedCount);
  }
  inline bool operator>(const TrackingData &other) const {
  	return count > other.count;// || (count == other.count && weightedCount > other.weightedCount);
  }
  bool track(double weight, bool forward) {
    if (weight < minimumWeight) {
      discarded++;
  	  return false;
    }
    if (count < MAX_COUNT) {
      count++;
      if (forward) 
        directionBias++;
      weightedCount += weight;
      
      if (count > maxCount) 
        maxCount = count;
      if (weightedCount > maxWeightedCount)
        maxWeightedCount = weightedCount;
      
      if (count == 1)
        singletonCount++;
      else if (count == 2)
        singletonCount--;
      
      return true;
    } else
      return false;
  }

  inline CountType getCount() { return count; }
  inline CountType getDirectionBias() { return directionBias; }
  inline WeightType getWeightedCount() { return weightedCount; }
  inline double getNormalizedDirectionBias() { return (double) directionBias / (double)count ; }

  
  std::string toString() {
    std::stringstream ss;
    ss << count << ":" << std::fixed << std::setprecision(2) << getNormalizedDirectionBias();
    ss << ':' << std::fixed << std::setprecision(2) << ((double)weightedCount / (double)count);
    return ss.str();
  }

};

class TrackingDataWithLastRead : public TrackingData
{
public:
  typedef std::pair< unsigned int, unsigned short > ReadAndPositionType;
  typedef std::vector< ReadAndPositionType > ReadAndPositionVectorType;
  
protected:
  ReadAndPositionType lastReadAndPosition;

public:
  bool track(double weight, bool forward, unsigned int readIdx, unsigned short readPos) { 
  	bool ret = TrackingData::track(weight,forward);
  	if (ret)
  	  lastReadAndPosition = ReadAndPositionType(readIdx,readPos);
    return ret;
  }
};

class TrackingDataWithAllReads : public TrackingData
{
public:
  typedef std::pair< unsigned int, unsigned short > ReadAndPositionType;
  typedef std::vector< ReadAndPositionType > ReadAndPositionVectorType;
  
protected:
  ReadAndPositionVectorType readsAndPositions;

public:
  bool track(double weight, bool forward, unsigned int readIdx, unsigned short readPos) { 
  	bool ret = TrackingData::track(weight,forward);
  	if (ret)
  	  readsAndPositions.push_back( ReadAndPositionType(readIdx,readPos) );
    return ret;
  }
};

class SolidTrackingData : public TrackingData
{
public:
  bool track(double weight, bool forward) {
  	return TrackingData::track(weight, forward);
  }
  SolidTrackingData &operator=(TrackingData &other) {
  	TrackingData::operator=(other);
  	return *this;
  }
};

//typedef TrackingData SolidTrackingData;

std::ostream &operator<<(std::ostream &stream, TrackingData &ob);


class SolidKmerTag : public KmerValue<SolidTrackingData> {};
class WeakKmerTag  : public KmerValue<TrackingDataWithLastRead> {};



template<typename Value = WeakKmerTag>
class KmerArray
{

public:
  
  typedef Value ValueType;
  
  class ElementType { 
    private:
  	  Kmer *_key; 
  	  ValueType *_value;
    public:
      ElementType(): _key(NULL), _value(NULL) {}
  	  ElementType(Kmer &key, ValueType &value): _key(&key), _value(&value) {}

  	  Kmer &key() { return *_key; }
  	  ValueType &value()   { return *_value; }
  	  inline bool operator==(const ElementType &other) const {
  	  	if (this->_value != NULL && other._value != NULL)
  	  	  return *(this->_value) == *(other._value);
  	  	else
  	  	  return true;
  	  }
  	  inline bool operator<(const ElementType &other) const {
  	  	if (this->_value != NULL && other._value != NULL)
  	  	  return *(this->_value) < *(other._value);
  	  	else
  	  	  return false;
  	  }
  	  inline bool operator>(const ElementType &other) const {
  	  	if (this->_value != NULL && other._value != NULL)
  	  	  return *(this->_value) > *(other._value);
  	  	else
  	  	  return false;
  	  }
  };

private:
  void * _begin; // safer than Kmer *: prevents incorrect ptr arithmetic: _begin+1 , use _add instead
  IndexType _size;

private:
  static void *_add(void *ptr,IndexType i)  {return ((char *)ptr + i * KmerSizer::getByteSize());}
  
  const Kmer &get(IndexType index) const    { return *(Kmer *)_add(_begin,index); }
  Kmer &get(IndexType index)                { return *(Kmer *)_add(_begin,index); }

public:
 
  KmerArray(IndexType size = 0):
   _size(0), _begin(NULL)
  {
    resize(size);
  }

  KmerArray(TwoBitEncoding *twoBit, SequenceLengthType length):
  _size(0),_begin(NULL)
  {
    SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
    resize(numKmers);
    build(twoBit,length);
  }

  KmerArray(const KmerArray &copy):
  _size(0), _begin(NULL)
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
    if (_begin == NULL)
       throw std::runtime_error("Could not allocate memory in KmerArray operator=()");
    
    memcpy(_begin,other._begin,_size*getElementByteSize());
    return *this;
  }

  const Kmer &operator[](IndexType index) const
  {
    if (index >= _size)
       throw std::invalid_argument("attempt to access index greater than size in KmerArray operator[] const"); 
    return get(index);
  }
  Kmer &operator[](IndexType index)
  {
    if (index >= _size)
       throw std::invalid_argument("attempt to access index greater than size in KmerArray operator[]"); 
    return get(index);
  }
  
  ValueType *getValueStart() const {
  	if (size() > 0)
  	  return (ValueType*) _add(_begin,size());
  	else
  	  return NULL;
  }
  const ValueType &valueAt(IndexType index) const
  {
    if (index >= _size)
    {
      throw std::invalid_argument("attempt to access index greater than size in KmerArray valueAt() const");
  	}
    return *( getValueStart() + index );
  }
  ValueType &valueAt(IndexType index)
  {
    if (index >= _size)
    {
      throw std::invalid_argument("attempt to access index greater than size in KmerArray valueAt()");
  	}
    return *( getValueStart() + index );
  }


  const ElementType getElement(IndexType idx) const
  {
  	return ElementType( get(idx), valueAt(idx) );
  }
  ElementType getElement(IndexType idx)
  {
  	return ElementType( get(idx), valueAt(idx) );
  }
  
  IndexType size() const { return _size; }

  IndexType getElementByteSize() const { 
  	return (KmerSizer::getByteSize() + sizeof(ValueType)); 
  }
  
  void reset()//PoolManager &pool = BoostPoolManager::get())
  {
    if (_begin != NULL) {
      //pool.free(test, size() * getElementByteSize());
      std::free(_begin);
    } 
    _begin = NULL;
    _size = 0;
  }

  void resize(IndexType size) {
  	resize(size, MAX_INDEX);
  }
  void resize(IndexType size, IndexType idx)//, PoolManager &pool = BoostPoolManager::get())
  {
    if (size == _size)
      return;
    IndexType oldSize = _size;
    
    // alloc / realloc memory
    _setMemory(size, idx);//, pool);
      
    if(_begin == NULL && size > 0) {
       throw std::runtime_error("Could not allocate memory in KmerArray resize()");
    }

    if (size > oldSize && idx == MAX_INDEX) {
       // zero fill remainder
       void *start = _add(_begin,oldSize); 
       memset(start, 0, KmerSizer::getByteSize()*(size-oldSize));
       start =  (getValueStart() + oldSize);
       memset(start, 0, sizeof(ValueType) * (size - oldSize));
    }    
  }
    
  void _copyRange(void * srcKmer, ValueType *srcValue, IndexType idx, IndexType srcIdx, IndexType count) {
	memcpy(_add(_begin,idx), _add(srcKmer,srcIdx), count * KmerSizer::getByteSize()); 
	memcpy(getValueStart() + idx, srcValue + srcIdx, count * sizeof(ValueType));
  }

  void _setMemory(IndexType size, IndexType idx)//, PoolManager &pool = BoostPoolManager::get())
  {
    void *memory    = NULL;
    
    void  *oldBegin = _begin;
    
    if (size != 0 ) {
    	// allocate new memory
    	//memory = pool.malloc( size * getElementByteSize() );
        memory = std::malloc( size * getElementByteSize() );
        
        if(memory == NULL) {
           throw std::runtime_error("Could not allocate memory in KmerArray _setMemory()");
        }        
    }
    IndexType oldSize = _size;
    IndexType lesserSize = std::min(size,oldSize);
    ValueType *oldValueStart;
    if (oldSize > 0) {
    	oldValueStart = getValueStart();
    }
    _begin =  memory ;
    _size = size;
    
    if (oldBegin != NULL && memory != NULL && oldValueStart != NULL && lesserSize > 0) {
      // copy the old contents
      if (idx == MAX_INDEX || idx >= lesserSize) {
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
    if (oldBegin != NULL) {
      // free old memory
      //pool.free(oldBegin, oldSize * getElementByteSize());
      std::free(oldBegin);
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
  // return a KmerArray that has one entry for each possible single-base substitution
  static KmerArray permuteBases(const Kmer &kmer, bool leastComplement = false ) {
    KmerArray kmers( KmerSizer::getSequenceLength() * 3 );
   
    TEMP_KMER(tmp);
    for(SequenceLengthType byteIdx=0; byteIdx<KmerSizer::getByteSize(); byteIdx++) { 
  	  int max = 12;
  	  if (byteIdx+1 == KmerSizer::getByteSize())
  	    max = 3*(KmerSizer::getSequenceLength()%4);
  	  if (max == 0)
  	    max = 12;
  	  for(SequenceLengthType j=0; j<max; j++) {
  	    SequenceLengthType kmerIdx = byteIdx*12+j;
        tmp = kmer;
        TwoBitEncoding *ptr = (TwoBitEncoding*) &tmp;
        ptr += byteIdx;
        *ptr = TwoBitSequence::permutations[ ((TwoBitEncoding)*ptr)*12 + j ];
        if (leastComplement)
          tmp.buildLeastComplement( kmers[kmerIdx] );
        else
          kmers[kmerIdx] = tmp;
      }
    }
    return kmers;
  }
  
  

  IndexType find(const Kmer &target) const {
    for(IndexType i=0; i<_size; i++)
      if (target.compare(get(i)) == 0)
        return i;
    return MAX_INDEX;
  }

  IndexType findSorted(const Kmer &target, bool &targetIsFound) const {
  	// binary search
  	IndexType min = 0;
  	IndexType max = size();
 
  	if (max == 0)
    {
       targetIsFound = false;
       return 0;
    }
  	IndexType mid;
  	int comp;
  	do {
  		mid = (min+max) / 2;
  		comp = target.compare(get(mid));
  		if (comp > 0)
  		  min = mid+1;
  		else if (comp < 0)
  		  max = mid-1;
  	} while (comp != 0 && max != MAX_INDEX && min <= max);
  	if (comp == 0)
  	  targetIsFound = true;
  	else
  	  targetIsFound = false;
  	return mid + (comp>0 && size()>mid?1:0);
  }

  void insertAt(IndexType idx, const Kmer &target) {
  	if (idx > size())
  	  throw std::invalid_argument("attempt to access index greater than size in KmerArray insertAt");
  	resize(size() + 1, idx);
  	get(idx) = target;
  }
  
  IndexType append(const Kmer &target) {
  	IndexType idx = size();
   	insertAt(size(), target);
   	return idx;
  }

  IndexType insertSorted(const Kmer &target) {
  	bool isFound;
  	IndexType idx = findSorted(target, isFound);
  	if (!isFound)
  	  insertAt(idx, target);
  	return idx;
  }

  void remove(const Kmer &target) {
  	bool isFound;
  	IndexType idx = find(target, isFound);
  	if (isFound)
  	  remove(idx);
  }
  void remove(IndexType idx) {
    resize(size()-1,idx);
  }
  void swap(IndexType idx1, IndexType idx2) {
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
  	IndexType idx=0;
  	for(idx=0; idx<size() && idx < 30; idx++) {
  		ss << get(idx).toFasta() << ":" << valueAt(idx) << ", ";
  	} 
  	if (idx < size())
  	    ss << " ... " << size() - idx << " more ";
  	ss << "}";
  	return ss.str();
  }

public:
  class Iterator : public std::iterator<std::forward_iterator_tag, KmerArray>
  {
    private:
      KmerArray *_tgt;
      IndexType _idx;
      ElementType elementCopy;
    public:
      Iterator(KmerArray *target, IndexType idx = 0): _tgt(target), _idx(idx) { }
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
      ElementType &operator*() { return elementCopy =   _tgt->getElement(_idx); }
      Kmer &key() { return _tgt->get(_idx); }
      Value &value() { return _tgt->valueAt(_idx); }
      ElementType *operator->() {  elementCopy = _tgt->getElement(_idx); return &elementCopy; }
      bool isEnd() { return _idx == _tgt->size(); }
  };

  Iterator begin() { return Iterator(this, 0); }
  Iterator end()   { return Iterator(this, size()); }

};



template<typename Value>
class KmerMap
{

public:
   typedef Kmer KeyType;
   typedef Value ValueType;
   typedef KmerArray<Value> BucketType;
   typedef typename BucketType::Iterator BucketTypeIterator;
   typedef typename BucketType::ElementType ElementType;
   
   typedef std::vector< BucketType > BucketsVector;
   typedef typename BucketsVector::iterator BucketsVectorIterator;

private:
   BucketsVector _buckets;
   
public:
   KmerMap(IndexType bucketCount = 1024) {
     _buckets.resize(bucketCount);
   }
   ~KmerMap() 
   {
     clear();
   }
   
   void reset() {
     for(int i=0; i< _buckets.size(); i++)
       _buckets[i].reset();
   }
   void clear() {
   	 reset();
     _buckets.resize(1); // iterators require at least 1
   }
   BucketType &getBucket(long hash) {
   	return _buckets[hash % _buckets.size()];
   }
   const BucketType &getBucket(long hash) const {
   	return _buckets[hash % _buckets.size()];
   }
   
   BucketType &getBucket(const KeyType &key)  {
     return getBucket(key.hash());
   }
   const BucketType &getBucket(const KeyType &key) const {
     return getBucket(key.hash());
   }

   ValueType &insert(const KeyType &key, const ValueType &value, BucketType *bucketPtr = NULL) {
   	  if (bucketPtr == NULL)
   	    bucketPtr = &getBucket(key);
   	    
   	  IndexType idx = bucketPtr->insertSorted(key);
   	  return bucketPtr->valueAt(idx) = value;
   }

   bool remove(const KeyType &key, BucketType *bucketPtr = NULL) {
   	  if (bucketPtr == NULL)
   	    bucketPtr = &getBucket(key);
   	  bool isFound;
   	  IndexType idx = bucketPtr->findSorted(key, isFound);
   	  if (isFound && idx != MAX_INDEX)
   	    bucketPtr->remove(idx);
   	  return isFound;
   }
   
   bool exists(const KeyType &key, BucketType *_bucketPtr = NULL) const {
   	 const BucketType *bucketPtr = _bucketPtr;
     if (bucketPtr == NULL)
   	   bucketPtr = &getBucket(key);
   	 bool isFound;
   	 bucketPtr->findSorted(key, isFound);
   	 return isFound;
   }
   ValueType &operator[](const KeyType &key) {
     BucketType &bucket = getBucket(key);
     bool isFound;
     IndexType idx = bucket.findSorted(key, isFound);
     if (isFound && idx != MAX_INDEX)
       return bucket.valueAt(idx);
     else 
       return insert(key, Value(), &bucket);       
   }
   
   IndexType size() const {
   	IndexType size = 0;
   	for(int i = 0; i<_buckets.size() ; i++)
   	  size += _buckets[i].size();
   	return size;
   } 
   
   std::string toString() {
  	std::stringstream ss;
  	ss << this << "[";
  	IndexType idx=0;
  	for(; idx<_buckets.size() && idx < 30; idx++) {
  		ss << "bucket:" << idx << ' ' << _buckets[idx].toString() << ", ";
  	} 
  	if (idx < _buckets.size())
  	  ss << " ... " << _buckets.size() - idx << " more ";
  	ss << "]";
  	return ss.str();
  }

public:
  class Iterator : public std::iterator<std::forward_iterator_tag, KmerMap>
  {
    friend class KmerMap;
    private:
     const KmerMap *_target;
     BucketsVectorIterator _iBucket;
     BucketTypeIterator  _iElement;

     void _moveToNextValidElement()  {
      while(_iElement == _iBucket->end() && ++_iBucket != _target->_buckets.end())
        _iElement = _iBucket->begin();
     }
          
      Iterator(KmerMap *target):    
      _target(target),
      _iBucket(target->_buckets.begin()),
      _iElement(_iBucket->begin())
     {
       _moveToNextValidElement();
     }

      Iterator(KmerMap *target,bool endConstructor):
      _target(target),
      _iBucket(target->_buckets.end()), 
      _iElement((_iBucket-1)->end())    // Assumes at least 1 bucket.
     {

     }

     public:
        
      bool operator==(const Iterator& other) 
      { return _iBucket == other._iBucket && _iElement == other._iElement; }
      
      bool operator!=(const Iterator& other) 
      { return !(*this == other); }

      Iterator& operator++() 
      { 
        ++_iElement;
        _moveToNextValidElement();
        return *this;
      }

      Iterator operator++(int unused) 
      { Iterator tmp(*this) ; ++(*this); return tmp; }

      ElementType &operator*()  { return *_iElement;        }
      ElementType *operator->() { return &(*_iElement);     }

      Kmer &key()               { return _iElement.key();   }
      Value &value()            { return _iElement.value(); }

      BucketType &bucket() { return *_iBucket; }
      IndexType bucketIndex() { return (_iBucket - _target->_buckets.begin()); }

  };

  Iterator begin()  { return Iterator(this);     }
  Iterator end()    { return Iterator(this,true);}

};

typedef KmerArray<double> KmerWeights;

typedef KmerMap<SolidKmerTag>   KmerSolidMap;
typedef KmerMap<WeakKmerTag>    KmerWeakMap;
typedef KmerMap<unsigned short> KmerCountMap;



#endif


 

//
// $Log: Kmer.h,v $
// Revision 1.54  2009-11-25 18:39:08  regan
// fixed bug across versions of compiler
//
// Revision 1.53  2009-11-24 13:35:29  cfurman
// removed KmerPtr class.
//
// Revision 1.52  2009-11-22 08:16:41  regan
// some fixes some bugs... optimized vs debug vs deb4/5 give different results
//
// Revision 1.51  2009-11-21 18:46:53  regan
// added bugs
//
// Revision 1.50  2009-11-21 15:58:29  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.49  2009-11-12 17:01:51  regan
// checkpoint
//
// Revision 1.48  2009-11-11 17:23:24  regan
// fixed bugs in heap generation
// solid picking logic needs work
//
// Revision 1.47  2009-11-11 07:57:23  regan
// built framework for  (not working) - make_heap is broken
//
// Revision 1.46  2009-11-09 19:37:17  regan
// enhanced some debugging / analysis output
//
// Revision 1.45  2009-11-07 00:26:13  cfurman
// minor formatting
//
// Revision 1.44  2009-11-06 16:59:11  regan
// added base substitution/permutations table and build function
//
// Revision 1.43  2009-11-06 04:08:23  regan
// minor changes
//
// Revision 1.42  2009-11-04 18:24:25  regan
// reworked tracking data
//
// Revision 1.41  2009-11-03 17:15:40  regan
// minor refactor
//
// Revision 1.40  2009-11-02 21:19:25  regan
// fixed types and boundary tests
//
// Revision 1.39  2009-11-02 20:02:48  cfurman
// KmerMap::iterator refactor
//
// Revision 1.38  2009-11-02 18:48:18  regan
// minor refactor and performance tweaks
//
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
