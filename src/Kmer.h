// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.58 2009-11-28 01:00:07 regan Exp $
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
#include <boost/shared_ptr.hpp>

#include "TwoBitSequence.h"
#include "MemoryUtils.h"

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 1024
#endif

typedef std::tr1::shared_ptr<TwoBitEncoding> KmerSharedPtr;
typedef unsigned long IndexType;

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
	     	for(unsigned int i=1; i<getTwoBitLength()-1; i++)
	     	  if (*firstByte != *(++nextByte))
	     	    return false;
	     	return true;
	     } else
	        return false;
	   }
  };


class TrackingDataSingleton;
class TrackingDataWithAllReads;

class TrackingData 
{
public:
  typedef unsigned short CountType;
  typedef float          WeightType;
  typedef unsigned int   ReadIdType;
  typedef unsigned short PositionType;
    
  class ReadPosition
  {
   public:
  	ReadIdType   readId;
  	PositionType position;
  	
  	ReadPosition(ReadIdType _read = 0, PositionType _pos = 0) : readId(_read), position(_pos) {}
  };
  	
  class ReadPositionWeight : public ReadPosition
  {
  public:
  	WeightType   weight;
  	
  	ReadPositionWeight(ReadIdType _read = 0, PositionType _pos = 0, WeightType _weight = 0.0) : ReadPosition(_read,_pos), weight(_weight) {}
  };
  typedef std::vector< ReadPositionWeight > ReadPositionWeightVector;
  
  
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
  static void setGlobals(CountType count, WeightType weightedCount) {
  	  if (count > maxCount) 
        maxCount = count;
      if (weightedCount > maxWeightedCount)
        maxWeightedCount = weightedCount;
      
      if (count == 1)
        singletonCount++;
      else if (count == 2)
        singletonCount--;
  }
  static inline bool isDiscard(WeightType weight) {
  	if ( weight < minimumWeight ) {
  		discarded++;
  		return true;
  	} else
  	    return false;
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
  	return getCount() == other.getCount();// && weightedCount == other.weightedCount;
  }
  inline bool operator<(const TrackingData &other) const {
  	return getCount()  < other.getCount();// || (count == other.count && weightedCount < other.weightedCount);
  }
  bool track(double weight, bool forward, ReadIdType readIdx = 0, PositionType readPos = 0) {
    if (isDiscard(weight))
  	  return false;
    
    if (count < MAX_COUNT) {
      count++;
      if (forward) 
        directionBias++;
      weightedCount += weight;
      
      setGlobals(getCount(), getWeightedCount());
      
      return true;
    } else
      return false;
  }

  inline unsigned long getCount() const { return count; }
  inline unsigned long getDirectionBias() const { return directionBias; }
  inline double getWeightedCount() const { return weightedCount; }
  inline double getNormalizedDirectionBias() const { return (double) directionBias / (double)count ; }

  ReadPositionWeightVector getEachInstance() const { 
  	//returns count entries from read -1, position 0
  	ReadPositionWeight dummy1((ReadIdType) -1, 0, getWeightedCount()/getCount());
  	ReadPositionWeightVector dummy(getCount(), dummy1); 
    return dummy;
  }
  
  std::string toString() const {
    std::stringstream ss;
    ss << count << ":" << std::fixed << std::setprecision(2) << getNormalizedDirectionBias();
    ss << ':' << std::fixed << std::setprecision(2) << ((double)weightedCount / (double)count);
    return ss.str();
  }

  TrackingData &operator=(const TrackingDataSingleton &other);
  TrackingData &operator=(const TrackingDataWithAllReads &other);
};

class TrackingDataWithLastRead : public TrackingData
{
public:
  
protected:
  ReadPosition readPosition;

public:
  TrackingDataWithLastRead() : TrackingData(), readPosition() {}
  
  bool track(double weight, bool forward, ReadIdType readIdx, PositionType readPos) { 
  	bool ret = TrackingData::track(weight,forward);
  	if (ret)
  	  readPosition = ReadPosition(readIdx,readPos);
    return ret;
  }
  ReadPositionWeightVector getEachInstance() { 
  	//returns count entries from read -1, position 0
  	ReadPositionWeight dummy1(readPosition.readId, readPosition.position, getWeightedCount()/getCount());
  	ReadPositionWeightVector dummy(getCount(), dummy1); 
    return dummy;
  }

};

class TrackingDataSingleton
{
public:
  typedef TrackingData::ReadPositionWeight       ReadPositionWeight;
  typedef TrackingData::ReadPositionWeightVector ReadPositionWeightVector;
  typedef TrackingData::CountType                CountType;
  typedef TrackingData::WeightType               WeightType;
  typedef TrackingData::ReadIdType               ReadIdType;
  typedef TrackingData::PositionType             PositionType;

protected:
  ReadPositionWeight instance; // save space and store direction within sign of weight.
  
public:
  TrackingDataSingleton() : instance(0,0,0.0) {}
  ~TrackingDataSingleton() { reset(); }
  void reset() {
  	// construction does not set singletonCount
  	if (instance.weight != 0.0)
  	  TrackingData::singletonCount--;
  	instance = ReadPositionWeight(0,0,0.0);
  }
  inline bool operator==(const TrackingDataSingleton &other) const {
  	return getCount() == other.getCount();// && weightedCount == other.weightedCount;
  }
  inline bool operator<(const TrackingDataSingleton &other) const {
  	return getCount() < other.getCount();// || (count == other.count && weightedCount < other.weightedCount);
  }
  bool track(double weight, bool forward, ReadIdType readIdx = 0, PositionType readPos = 0) {
    if (TrackingData::isDiscard(weight)) 
  	  return false;
    
    instance = ReadPositionWeight(readIdx, readPos, forward ? weight : -1.0*weight);
    
    TrackingData::setGlobals(1, weight);

    return true;
  }

  inline unsigned long getCount() const { return instance.weight == 0.0 ? 0 : 1; }
  inline unsigned long getDirectionBias() const { return instance.weight > 0 ? 1 : 0; }
  inline double getWeightedCount() const { 
    return instance.weight > 0 ? instance.weight : -1.0 * instance.weight;
  }
  inline double getNormalizedDirectionBias() const { return (double) getDirectionBias() / (double) getCount() ; }

  inline ReadPositionWeightVector getEachInstance() const {
  	return ReadPositionWeightVector(1, instance);
  }
  
  inline ReadIdType   getReadId()  { return instance.readId;   }
  inline PositionType getPosition() { return instance.position; }
  
  std::string toString() const {
    std::stringstream ss;
    ss << getCount() << ":" << std::fixed << std::setprecision(2) << getNormalizedDirectionBias();
    ss << ':' << std::fixed << std::setprecision(2) << ((double)getWeightedCount() / (double)getCount());
    return ss.str();
  }
  	
};

class TrackingDataWithAllReads
{
public:
  typedef TrackingData::ReadPositionWeight       ReadPositionWeight;
  typedef TrackingData::ReadPositionWeightVector ReadPositionWeightVector;
  typedef TrackingData::CountType                CountType;
  typedef TrackingData::WeightType               WeightType;
  typedef TrackingData::ReadIdType               ReadIdType;
  typedef TrackingData::PositionType             PositionType;
  
protected:
  ReadPositionWeightVector instances;
  unsigned int directionBias;

public:
  TrackingDataWithAllReads() : instances(), directionBias(0) { }
  
  void reset() {
  	instances.clear();
  	directionBias = 0;
  }
  inline bool operator==(const TrackingDataWithAllReads &other) const {
  	return getCount() == other.getCount();// && weightedCount == other.weightedCount;
  }
  inline bool operator<(const TrackingDataWithAllReads &other) const {
  	return getCount() < other.getCount();// || (count == other.count && weightedCount < other.weightedCount);
  }
  bool track(double weight, bool forward, ReadIdType readIdx = 0, PositionType readPos = 0) {
    if (TrackingData::isDiscard(weight))
  	  return false;

    if (forward) 
      directionBias++;
    
    ReadPositionWeight rpw(readIdx, readPos, weight);
    
    instances.push_back( rpw );
    
    TrackingData::setGlobals(getCount(), getWeightedCount());

    return true;
  }

  inline unsigned long getCount() const { return instances.size(); }
  inline unsigned long getDirectionBias() const { return directionBias; }
  inline double getWeightedCount() const { 
    double weightedCount = 0.0;
    for(ReadPositionWeightVector::const_iterator it = instances.begin(); it != instances.end(); it++)
      weightedCount += it->weight;
    return weightedCount;
  }
  inline double getNormalizedDirectionBias() const { return (double) getDirectionBias() / (double) getCount() ; }

  inline ReadPositionWeightVector getEachInstance() const {
  	return instances;
  }
  
  std::string toString() const {
    std::stringstream ss;
    ss << getCount() << ":" << std::fixed << std::setprecision(2) << getNormalizedDirectionBias();
    ss << ':' << std::fixed << std::setprecision(2) << ((double)getWeightedCount() / (double)getCount());
    return ss.str();
  }
  TrackingDataWithAllReads &operator=(const TrackingDataSingleton &other) {
  	this->reset();
  	ReadPositionWeightVector rpw = other.getEachInstance();
  	track(other.getWeightedCount(), other.getDirectionBias() == 1, rpw[0].readId, rpw[0].position);
  	return *this;
  }
};


std::ostream &operator<<(std::ostream &stream, TrackingData &ob);
std::ostream &operator<<(std::ostream &stream, TrackingDataSingleton &ob);
std::ostream &operator<<(std::ostream &stream, TrackingDataWithAllReads &ob);




template<typename Value>
class KmerArray
{

public:
  
  typedef Value ValueType;
  static const IndexType MAX_INDEX = (IndexType) -1;  
  
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
   _begin(NULL),_size(0)
  {
    resize(size);
  }

  KmerArray(TwoBitEncoding *twoBit, SequenceLengthType length):
  _begin(NULL),_size(0)
  {
    SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
    resize(numKmers);
    build(twoBit,length);
  }

  KmerArray(const KmerArray &copy):
  _begin(NULL),_size(0)
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

       // Values should already have been constructed
    }    
  }
    
  void _copyRange(void * srcKmer, ValueType *srcValue, IndexType idx, IndexType srcIdx, IndexType count) {
	memcpy(_add(_begin,idx), _add(srcKmer,srcIdx), count * KmerSizer::getByteSize()); 
	//assignment copy Values
	Value *startValue=getValueStart();
	for(IndexType i=0; i<count; i++)
	  *(startValue + idx + i) = *(srcValue + srcIdx + i);
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
    
    // construct all unintialized Values
    for(IndexType i = 0 ; i<size ; i++)
  	  ::new((void*) (getValueStart() + i)) Value();

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
      //   pool.free(oldBegin, oldSize * getElementByteSize());
      // destruct old Values
      for(IndexType i = 0; i < oldSize; i++)
        (oldValueStart + i)->~Value();
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
  	  unsigned int max = 12;
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
      Iterator(KmerArray *target, IndexType idx = 0): _tgt(target), _idx(idx) {
      	setElementCopy();
      }
      Iterator(const Iterator &copy) { *this = copy; }
      ~Iterator() {}
      Iterator& operator=(const Iterator& other) {
        _tgt = other._tgt;
        _idx = other._idx;
        elementCopy = other.elementCopy;
        return *this;
      }
      bool operator==(const Iterator& other) const { return _idx == other._idx && _tgt == other._tgt; }
      bool operator!=(const Iterator& other) const { return _idx != other._idx || _tgt != other._tgt; }
      Iterator& operator++() { if ( !isEnd() ) { ++_idx; setElementCopy(); } return *this; }
      Iterator operator++(int unused) { Iterator tmp(*this) ; ++(*this); return tmp; }
      ElementType &operator*() { return elementCopy; }
      Kmer &key() { return elementCopy.key(); }
      Value &value() { return elementCopy.value(); }
      ElementType *operator->() { return &elementCopy; }
      bool isEnd() const { return _idx >= _tgt->size(); }
      void setElementCopy() { if ( !isEnd() ) elementCopy = _tgt->getElement(_idx); }
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
     for(IndexType i=0; i< _buckets.size(); i++)
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
   	  if (isFound && idx != BucketType::MAX_INDEX)
   	    bucketPtr->remove(idx);
   	  return isFound;
   }
   
   bool _exists(const KeyType &key, IndexType &existingIdx, const BucketType *bucketPtr = NULL) const {
     if (bucketPtr == NULL)
   	   bucketPtr = &getBucket(key);
   	 bool isFound;
   	 existingIdx = bucketPtr->findSorted(key, isFound);
   	 return isFound;
   }
   bool exists(const KeyType &key, BucketType *bucketPtr = NULL) const {
   	 IndexType dummy;
   	 return _exists(key, dummy, bucketPtr);
   }
   
   Value *getIfExists(const KeyType &key, BucketType *bucketPtr = NULL) {
     if (bucketPtr == NULL)
   	   bucketPtr = &getBucket(key);
   	 IndexType existingIdx;
   	 if (_exists(key, existingIdx, bucketPtr))
   	 	return &(bucketPtr->valueAt(existingIdx));
   	 else
   	   return NULL;
   }
   ValueType &operator[](const KeyType &key) {
     BucketType &bucket = getBucket(key);
     IndexType existingIdx;
   	 if (_exists(key, existingIdx, &bucket))
   	   return bucket.valueAt(existingIdx);
     else 
       return insert(key, Value(), &bucket);       
   }
   
   IndexType size() const {
   	IndexType size = 0;
   	for(IndexType i = 0; i<_buckets.size() ; i++)
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

typedef KmerArray<char>          Kmers;
typedef KmerArray<double>        KmerWeights;
typedef KmerArray<unsigned long> KmerCounts;

#endif


 

//
// $Log: Kmer.h,v $
// Revision 1.58  2009-11-28 01:00:07  regan
// fixed bugs and warnings
//
// Revision 1.57  2009-11-27 23:16:58  regan
// refactored and broke it
//
// Revision 1.56  2009-11-27 01:53:40  regan
// refactored and got first pass at error rate by position
//
// Revision 1.55  2009-11-26 09:03:29  regan
// refactored and stuff
//
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
