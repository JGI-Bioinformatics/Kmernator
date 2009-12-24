// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.64 2009-12-24 00:39:22 cfurman Exp $
//

#ifndef _KMER_H
#define _KMER_H
#include <tr1/memory>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <boost/functional/hash.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>

#include "config.h"
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
	public:
	   typedef unsigned long NumberType;
	   
    private:
       inline static boost::hash<NumberType> &getHasher() { static boost::hash<NumberType> hasher; return hasher; }
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
	
       static NumberType toNumber(const Kmer &kmer) {
       	  NumberType val;
       	  switch(KmerSizer::getTwoBitLength()) {
       	  	case 1:  val =   (NumberType) *((boost::uint8_t *)   kmer.getTwoBitSequence()); 
       	  	         break;
       	  	case 2:  val =   (NumberType) *((boost::uint16_t *)  kmer.getTwoBitSequence()); 
       	  	         break;
       	  	case 3:  val = (((NumberType) *((boost::uint16_t *)  kmer.getTwoBitSequence()))<<8)
       	  	              +  (NumberType) *((boost::uint8_t *)  (kmer.getTwoBitSequence()+2)); 
       	  	         break;
       	  	case 4:  val =   (NumberType) *((boost::uint32_t *)  kmer.getTwoBitSequence()); 
       	  	         break;
       	  	case 5:  val = (((NumberType) *((boost::uint32_t *)  kmer.getTwoBitSequence()))<<8)
       	  	              +  (NumberType) *((boost::uint8_t *)  (kmer.getTwoBitSequence()+4));
       	  	         break;
       	  	case 6:  val = (((NumberType) *((boost::uint32_t *)  kmer.getTwoBitSequence()))<<16)
       	  	              +  (NumberType) *((boost::uint16_t *) (kmer.getTwoBitSequence()+4)); 
       	  	         break;
       	  	case 7:  val = (((NumberType) *((boost::uint32_t *)  kmer.getTwoBitSequence()))<<24)
       	  	              +(((NumberType) *((boost::uint16_t *) (kmer.getTwoBitSequence()+4)))<<8)
       	  	              +  (NumberType) *((boost::uint8_t *)  (kmer.getTwoBitSequence()+6)); 
       	  	         break;
       	  	default: val =  (NumberType) *((boost::uint64_t *) kmer.getTwoBitSequence());
       	  }
       	  return val;
       }
	
	public:
	   
	   Kmer(const Kmer &copy) {
	   	  *this = copy;
	   }
	   
	   inline int compare(const Kmer &other) const
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


       inline Kmer *get() {
         return this;
       }
       
	   inline bool operator ==(const Kmer &other) const
	   {
	      return compare(other) == 0;
	   }
	   inline bool operator !=(const Kmer &other) const
	   {
	   	  return compare(other) != 0;
	   }
	   inline bool operator <(const Kmer &other) const
	   {
	   	  return compare(other) < 0;
	   }
	   inline bool operator <=(const Kmer &other) const
	   {
	   	  return compare(other) <= 0;
	   }
	   inline bool operator >(const Kmer &other) const
	   {
	   	  return compare(other) > 0;
	   }
	   inline bool operator >=(const Kmer &other) const
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

	
	   inline TwoBitEncoding *getTwoBitSequence() const
	   {
	     return (TwoBitEncoding *)_data();
	   }
	   inline SequenceLengthType getTwoBitLength() const
	   {
	   	 return KmerSizer::getTwoBitLength();
	   }
	   inline SequenceLengthType getByteSize() const
	   {
	   	 return KmerSizer::getByteSize();
	   }
	   inline SequenceLengthType getLength() const
	   {
	   	 return KmerSizer::getSequenceLength();
	   }
	
	   void buildReverseComplement(Kmer &output) const
	   {
	     TwoBitSequence::reverseComplement((TwoBitEncoding*)_data(), (TwoBitEncoding*)output._data(), getLength());
	   }
	   
	   // returns true if this is the least complement, false otherwise
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
	   inline NumberType toNumber() const
	   {
	   	 return Kmer::toNumber(*this);
	   }
	   inline NumberType hash() const
	   {
	     return getHasher()(toNumber(*this));
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
  	  if (count > maxCount) {
        maxCount = count;
  	  }
      if (weightedCount > maxWeightedCount) {
        maxWeightedCount = weightedCount;
      }
      
      if (count == 1) {
        #ifdef _USE_THREADSAFE_KMER
        #pragma omp atomic
        #endif
        singletonCount++;
      }
      else if (count == 2) {
        #ifdef _USE_THREADSAFE_KMER
        #pragma omp atomic
        #endif
        singletonCount--;
      }
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
      #ifdef _USE_THREADSAFE_KMER
      #pragma omp critical
	  #endif
      {    	
      count++;
      if (forward) 
        directionBias++;
      weightedCount += weight;
      }
      
      setGlobals(getCount(), getWeightedCount());
      
      return true;
    } else
      return false;
   
  }

  inline unsigned long getCount() const { return count; }
  inline unsigned long getDirectionBias() const { return directionBias; }
  inline double getWeightedCount() const { return weightedCount; }
  inline double getAverageWeight() const { return weightedCount/count; }
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

  void reset() {
  	// destructor should not call this
  	if (getCount() == 1)
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
  	#ifdef _USE_THREADSAFE_KMER
    #pragma omp critical
	#endif
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
  WeightType weightedCount; // for performance reasons

public:
  TrackingDataWithAllReads() : instances(), directionBias(0), weightedCount(0.0) { }
  
  void reset() {
  	instances.clear();
  	directionBias = 0;
  	weightedCount = 0.0;
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
  	#ifdef _USE_THREADSAFE_KMER
    #pragma omp critical
	#endif
	{
    if (forward) 
      directionBias++;
    
    ReadPositionWeight rpw(readIdx, readPos, weight);
    
    instances.push_back( rpw );
    weightedCount += weight;
	}
    TrackingData::setGlobals(getCount(), getWeightedCount());

    return true;
  }

  inline unsigned long getCount() const { return instances.size(); }
  inline unsigned long getDirectionBias() const { return directionBias; }
  inline double getWeightedCount() const { return weightedCount; }
  double _getWeightedCount() const { 
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
  	if (other.getCount() > 0) {
  	  directionBias = other.getDirectionBias();
  	  instances = other.getEachInstance();
  	  weightedCount = _getWeightedCount();
  	}
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
  	  KmerArray *_array;
    public:
      ElementType(): _key(NULL), _value(NULL), _array(NULL) { }
  	  ElementType(Kmer &key, ValueType &value, KmerArray &array): _key(&key), _value(&value), _array(&array) { setLock(); }
  	  ElementType(const ElementType &copy): _key(NULL), _value(NULL), _array(NULL) {
  	  	*this = copy;
  	  }
  	  ~ElementType() { reset(); }
  	  ElementType &operator=(const ElementType &copy) {
  	  	reset();
  	  	_key = copy._key;
  	  	_value = copy._value;
  	  	_array = copy._array;
  	  	setLock();
  	  	return *this;
  	  }
  	  void reset() {
  	  	unsetLock();
  	  	_key = NULL;
  	  	_value = NULL;
  	  	_array = NULL;
  	  }
  	  inline bool isValid() const { return _array != NULL; }
  	  void setLock() {
  	  	if (isValid())
  	  	  _array->setSharedLock();
  	  }
  	  void unsetLock() {
  	  	if (isValid())
  	  	  _array->unsetSharedLock();
  	  }

  	  const Kmer &key() const { return *_key; }
  	  Kmer &key() { return *_key; }
  	  const ValueType &value() const { return *_value; }
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
  void                   *_begin; // safer than Kmer *: prevents incorrect ptr arithmetic: _begin+1 , use _add instead
  IndexType               _size;
  IndexType               _capacity;
  

#ifdef _USE_THREADSAFE_KMER

  typedef std::vector< omp_nest_lock_t > LockVector;
  mutable omp_nest_lock_t                _lock;
  mutable LockVector                     _sharedLocks;

private:
 
  // thread locking routines
  void initLock()         const { omp_init_nest_lock(&_lock);     initSharedLocks(); }
  void destroyLock()      const { omp_destroy_nest_lock(&_lock);  destroySharedLocks(); }
  void initSharedLocks()  const { 
  	_sharedLocks.clear();
  	_sharedLocks.resize( omp_get_num_procs() );
  	for(unsigned int i=0; i < _sharedLocks.size(); i++)
  	  omp_init_nest_lock( &(_sharedLocks[i]) );
  }
  void destroySharedLocks() const {
  	for(unsigned int i=0; i < _sharedLocks.size(); i++)
  	  omp_destroy_nest_lock( &(_sharedLocks[i]) );
  	_sharedLocks.clear();
  }
  
  inline void setLock()   const { omp_set_nest_lock(&_lock);         }
  inline bool testLock()  const { return omp_test_nest_lock(&_lock); }
  inline void unsetLock() const { omp_unset_nest_lock(&_lock);       }


public:
  // TODO make a true shared lock without a severe performance penalty

  void setSharedLock() const {
  	if (_sharedLocks.empty())
  	  return;
  	unsigned int myThread = omp_get_thread_num();
  	bool gotShared = false;
  	while (!gotShared) {
  		
  		bool gotExclusive = false;
  		if (testLock()) {
  		    gotShared = omp_test_nest_lock( &( _sharedLocks[myThread] ) );
  			unsetLock();
  			gotExclusive = true;
  		}
  		if (! gotShared ) {
  			
  				//std::stringstream ss;
  	            //ss  << "Waiting for shared lock: " << this << " " << myThread << "/" << _sharedLocks.size() << " " << gotExclusive << std::endl;
  		        //#pragma omp critical
  		        //{ std::cerr << ss.str(); }
  		        usleep(1);
  			  			
  		} 		
    }

  	//std::stringstream ss2;
  	//ss2 << "Got shared lock: " << this << " " << myThread << "/" << _sharedLocks.size() <<std::endl << StackTrace::getStackTrace();
  	//#pragma omp critical
  	//{ std::cerr <<  ss2.str(); }
  	
  }
  
  inline void unsetSharedLock() const { 
  	if (_sharedLocks.empty())
  	  return;
  	unsigned int myThread = omp_get_thread_num();  
  	omp_unset_nest_lock( &( _sharedLocks[ myThread ] ) );
    //std::stringstream ss;
  	//ss  << "Rel shared lock: " << this << " " << myThread <<std::endl;
    
    //#pragma omp critical
    //{ std::cerr << ss.str(); } 
  }

  inline void setExclusiveLock() const {  
  	if (_sharedLocks.empty())
  	  return;
  	//std::stringstream ss;
  	//ss << "Trying exclusive lock: " << this << " " << omp_get_thread_num() << "/" << _sharedLocks.size() <<std::endl << StackTrace::getStackTrace(); 
  	//#pragma omp critical
    //{ std::cerr << ss.str(); }
    bool gotExclusive = false;
  	while(!gotExclusive) {
  	
  	   gotExclusive = true;
  	   if (testLock()) {
  	   	 unsigned int i = 0;
  	     for(; i < _sharedLocks.size() ; i++) {
  	        if (! omp_test_nest_lock( &( _sharedLocks[i] ) )) {
  	        	gotExclusive = false;
  				break;
  			}
  	     }
  		 if (! gotExclusive) {
  		 	// unset acquired locks
  		 	unsetLock();
  		 	for(unsigned int j = 0; j<i ; j++)
  		 	  omp_unset_nest_lock( &( _sharedLocks[j] ) );
  		 	//std::stringstream ss2;
  	        //ss2 << "Waiting for exclusive lock: " << this << " " << omp_get_thread_num() << "/" << _sharedLocks.size() << " failed to get shared " << i << std::endl;
  	        //#pragma omp critical
  	        //{ std::cerr << ss2.str(); }
  	        usleep(1);
  		 }
  	   } else {
  	   	 gotExclusive = false;
  	   }
  	   
  	}
  	//std::stringstream ss3;
  	//ss3 <<"Got exclusive lock: " << this << " " << omp_get_thread_num() << "/" << _sharedLocks.size() <<std::endl;
  	//#pragma omp critical
    //{ std::cerr << ss3.str(); }
  }
  inline void unsetExclusiveLock() const { 
  	if (_sharedLocks.empty())
  	  return;
  	  //throw std::runtime_error("unsetExclusiveLock() called while under readOnlyOptimization mode");
  	
  	unsetLock();
  	// unset acquired shared locks
  	for(unsigned int j = 0 ; j < _sharedLocks.size() ; j++)
        omp_unset_nest_lock( &( _sharedLocks[j] ) );
        
  	//std::stringstream ss;
  	//ss << "Rel exclusive lock: " << this << " " << omp_get_thread_num() << "/" << _sharedLocks.size() <<std::endl;
  	//#pragma omp critical
  	//{ std::cerr << ss.str(); }
  }
  
#else // _USE_THREADSAFE_KMER

private:
  void initLock()            const { }
  void destroyLock()         const { }
  void initSharedLocks()     const { }
  void destroySharedLocks()  const { }
public:
  inline void setSharedLock() const { }
  inline void unsetSharedLock() const { }
  inline void setExclusiveLock()const { }
  inline void unsetExclusiveLock() const { }

#endif //  _USE_THREADSAFE_KMER


private:
  static inline const void *_add(const void *ptr, IndexType i) 
                                                   { return ((char *)ptr + i * KmerSizer::getByteSize()); }
  static inline void *_add(void *ptr,IndexType i)  { return ((char *)ptr + i * KmerSizer::getByteSize()); }
  
  // these are never thread safe!
  inline const Kmer &get(IndexType index) const    { return *((Kmer *) _add(_begin,index)); }
  inline Kmer &get(IndexType index)                { return *((Kmer *) _add(_begin,index)); }


public:
  void setReadOnlyOptimization() const {
  	destroyLock();
  }
  void unsetReadOnlyOptimization() {
  	initLock();
  }
  
  
public:
 
  KmerArray(IndexType size = 0):
   _begin(NULL),_size(0),_capacity(0)
  {
  	initLock();
    resize(size);
  }

  KmerArray(TwoBitEncoding *twoBit, SequenceLengthType length, bool leastComplement = false):
  _begin(NULL),_size(0),_capacity(0)
  {
  	initLock();
    SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
    resize(numKmers, MAX_INDEX, false);
    build(twoBit,length, leastComplement);
  }

  KmerArray(const KmerArray &copy):
  _begin(NULL),_size(0),_capacity(0)
  {
  	initLock();
    *this = copy;
  }
   
  ~KmerArray()
  {
     reset();
     destroyLock();
  }

  KmerArray &operator=(const KmerArray &other)
  {
    if (this == &other)
  	  return *this;
    reset();
  	setExclusiveLock();
    reserve(other.size());
    if (size() == 0) {
      unsetExclusiveLock();
      return *this;
    }
    if (_begin == NULL)
       throw std::runtime_error("Could not allocate memory in KmerArray operator=()");
  
    _copyRange(other._begin,other.getValueStart(),0,0,_size,false);
    
    unsetExclusiveLock();
    return *this;
  }
  
protected:
  // never thread safe!
  const ValueType *getValueStart() const {
  	if (capacity() > 0) 
  	  return (ValueType*) _add(_begin,capacity());
  	else
  	  return NULL;
  }
  
  // never thread safe!
  ValueType *getValueStart() {
  	if (capacity() > 0) 
  	  return (ValueType*) _add(_begin,capacity());
  	else
  	  return NULL;
  }

public:
  // never thread safe!
  const Kmer &operator[](IndexType index) const
  {
    if (index >= size())
       throw std::invalid_argument("attempt to access index greater than size in KmerArray operator[] const"); 
    return get(index);
  }
  
  // never thread safe!
  Kmer &operator[](IndexType index)
  {
    if (index >= size())
       throw std::invalid_argument("attempt to access index greater than size in KmerArray operator[]"); 
    return get(index);
  }
  
  
  // never thread safe!
  const ValueType &valueAt(IndexType index) const
  {
    if (index >= _size)
    {
      throw std::invalid_argument("attempt to access index greater than size in KmerArray valueAt() const");
  	}
    return *( getValueStart() + index );
  }
  
  // never thread safe!
  ValueType &valueAt(IndexType index)
  {
    if (index >= _size)
    {
      throw std::invalid_argument("attempt to access index greater than size in KmerArray valueAt()");
  	}
    return *( getValueStart() + index );
  }

public:
  const ElementType getElement(IndexType idx) const
  {
  	return ElementType( get(idx), valueAt(idx), *this );
  }
  ElementType getElement(IndexType idx)
  {
  	return ElementType( get(idx), valueAt(idx), *this );
  }
  
  static inline IndexType getElementByteSize() { 
  	return (KmerSizer::getByteSize() + sizeof(ValueType)); 
  }

  void reset(bool releaseMemory = true)
  {
  	setExclusiveLock();
  	
    if (_begin != NULL && releaseMemory) {
      // destruct old Values
      for(IndexType i = 0; i < _capacity; i++)
        (getValueStart() + i)->~Value();
      // free memory    	
      std::free(_begin);
      
      _begin = NULL;
      _capacity = 0;
    }
    _size = 0;
    
    unsetExclusiveLock();
  }

  inline IndexType size() const { return _size; }
  inline IndexType capacity() const { return _capacity; }
  inline bool empty() const { return _size == 0; }
  void reserve(IndexType size) { resize(size, MAX_INDEX, false); };
  
  
  void resize(IndexType size) {
  	resize(size, MAX_INDEX, false);
  }
  void resize(IndexType size, IndexType idx, bool reserveExtra = true)
  {
    if (size == _size)
      return;
    IndexType oldSize = _size;
    
    setExclusiveLock();
    // alloc / realloc memory
    _setMemory(size, idx, reserveExtra);
      
    if(_begin == NULL && size > 0) {
       throw std::runtime_error("Could not allocate memory in KmerArray resize()");
    }

    if (size > oldSize && idx == MAX_INDEX) {
       // zero fill remainder
       void *start = _add(_begin,oldSize); 
       memset(start, 0, KmerSizer::getByteSize()*(size-oldSize));

       // Values should already have been constructed
    }
    unsetExclusiveLock();
  }
    
  void _copyRange(const void * srcKmer, const ValueType *srcValue, IndexType idx, IndexType srcIdx, IndexType count, bool isOverlapped = false) {
	
	if (isOverlapped)
	  memmove(_add(_begin,idx), _add(srcKmer,srcIdx), count * KmerSizer::getByteSize()); 
	else
	  memcpy(_add(_begin,idx), _add(srcKmer,srcIdx), count * KmerSizer::getByteSize()); 

	//assignment copy Values
	Value *dstValue=getValueStart();	
	if (isOverlapped) {// && dstValue + idx > srcValue + srcIdx) {
	  Value *tmp = new Value[count];  
	  for(IndexType i=0; i<count; i++)
	    tmp[i] = *(srcValue + srcIdx + i);
	  for(IndexType i=0; i<count; i++)
	    *(dstValue + idx + i) = tmp[i];
	  delete [] tmp;
	} else {
	  // copy forward, without a temp buffer, is safe
	  for(IndexType i=0; i<count; i++)
	    *(dstValue + idx + i) = *(srcValue + srcIdx + i);
	}
  }

  void _setMemory(IndexType size, IndexType idx, bool reserveExtra = true)
  {    
    // preserve old pointers and metrics
    void  *oldBegin = _begin;
    IndexType oldSize = _size;
    IndexType oldCapacity = _capacity;
    IndexType lesserSize = std::min(size,oldSize);
    ValueType *oldValueStart = NULL;
    if (oldCapacity > 0) {
    	oldValueStart = getValueStart();
    }
    
    bool memChanged = false;
    if (size == 0) {
    	// ignore reserveExtra for zero resize
    	if ( _capacity > 0 ) {
    	    _capacity = 0;
    	    _begin = NULL;
    		memChanged = true;
    	}
    } else if ((size > _capacity) || (_capacity > size && !reserveExtra))  {
    	// allocate new memory
    	if (reserveExtra)
    	  _capacity = std::max( (IndexType) (size*1.2), (IndexType) (size+10));
    	else
    	  _capacity = size;
        void *memory = std::malloc( _capacity * getElementByteSize() );
        if(memory == NULL) {
           std::cerr << "Attempt to malloc " << _capacity * getElementByteSize() << std::endl;
           throw std::runtime_error("Could not allocate memory in KmerArray _setMemory()");
        }        
        _begin = memory;
        memChanged = true;
    }   
    
    _size = size;
    
    if (memChanged && _begin != NULL) {
      // construct all values that are newly allocated
      for(IndexType i = 0 ; i < _capacity ; i++)
  	    ::new((void*) (getValueStart() + i)) Value();
    }
    
    if (_begin != NULL && oldBegin != NULL && oldValueStart != NULL && lesserSize > 0) {
      // copy the old contents
      if (idx == MAX_INDEX || idx >= lesserSize) {
      	if (memChanged) {
      	  // copy all records in order (default ; end is trimmed or expanded)
      	  _copyRange(oldBegin,oldValueStart, 0, 0, lesserSize);
        } else {
      	  // noop
      	}
      } else {      	
      	// shrink or expand.
      	
      	if (memChanged && idx > 0) {
      	  // the first record(s) leading to idx will be copied	
      	  _copyRange(oldBegin, oldValueStart, 0, 0, idx);
      	}
      	
      	if (lesserSize == size) {
      	  // shrink: skipping the old record at idx
      	  if (idx < size) {
      	  	_copyRange(oldBegin, oldValueStart, idx, idx+1, lesserSize-idx, !memChanged);
      	  }
      	} else {
      	  // expand: leaving new (uninitialized/unchanged) record at idx
      	  if (idx < oldSize) {
      	  	_copyRange(oldBegin, oldValueStart, idx+1, idx, lesserSize-idx, !memChanged);
      	  }
      	}
      }
    }
    
    if (memChanged && oldBegin != NULL) {
      // destruct old Values
      for(IndexType i = 0; i < oldCapacity; i++)
        (oldValueStart + i)->~Value();
      // free old memory
      std::free(oldBegin);
    }
  }
  
  void build(TwoBitEncoding *twoBit, SequenceLengthType length, bool leastComplement = false)
  {
  	setExclusiveLock();
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
    if (leastComplement) {
      TEMP_KMER(least);
      for(SequenceLengthType i=0; i < numKmers ; i++)
        if (! kmers[i].buildLeastComplement(least) )
          kmers[i] = least;
    }
    unsetExclusiveLock();
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

protected:
  IndexType _find(const Kmer &target) const {
    for(IndexType i=0; i<_size; i++)
      if (target.compare(get(i)) == 0) {
        return i;
      }
    return MAX_INDEX;
  }
public:
  IndexType find(const Kmer &target) const {
    setSharedLock();
    IndexType idx = _find(target);
    unsetSharedLock();
    return idx;
  }
  
protected:
  IndexType _findSorted(const Kmer &target, bool &targetIsFound) const {
  	// binary search
  	IndexType min = 0;
  	IndexType max = size();
 
  	if (max == 0)
    {
       targetIsFound = false;
       return 0;
    }
    max--; // never let mid == size()
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
  
public:
  IndexType findSorted(const Kmer &target, bool &targetIsFound) const {
 	setSharedLock();
 	IndexType idx = _findSorted(target, targetIsFound);
 	unsetSharedLock();
 	return idx;
  }

protected:
  void _insertAt(IndexType idx, const Kmer &target) {
  	if (idx > size())
  	  throw std::invalid_argument("attempt to access index greater than size in KmerArray insertAt");
  	resize(size() + 1, idx);
  	get(idx) = target;
  }
  void _insertAt(IndexType idx, const Kmer &target, const Value &value) {
  	_insertAt(idx, target);
  	valueAt(idx) = value;
  }

public:  
  void insertAt(IndexType idx, const Kmer &target) {
  	setExclusiveLock();
  	_insertAt(idx, target);
  	unsetExclusiveLock();
  }
  void insertAt(IndexType idx, const Kmer &target, const Value &value) {
  	setExclusiveLock();
  	_insertAt(idx,target,value);
  	unsetExclusiveLock();
  }
  
  IndexType append(const Kmer &target) {
  	setExclusiveLock();
   	IndexType idx = size();
   	_insertAt(idx, target);
   	unsetExclusiveLock();
   	return idx;
  }
  IndexType append(const Kmer &target, const Value &value) {
  	setExclusiveLock();
   	IndexType idx = size();
   	_insertAt(idx, target, value);
   	unsetExclusiveLock();
   	return idx;
  }

protected:
  IndexType _insertSorted(const Kmer &target) {
  	bool isFound;
  	IndexType idx = findSorted(target, isFound);
  	if (!isFound)
  	  _insertAt(idx, target);
  	return idx;
  }
  IndexType _insertSorted(const Kmer &target, const Value &value) {
  	IndexType idx = _insertSorted(target);
  	valueAt(idx) = value;
  	return idx;
  }
  
public:
  IndexType insertSorted(const Kmer &target) {
  	setExclusiveLock();
  	IndexType idx = _insertSorted(target);
  	unsetExclusiveLock();
  	return idx;
  }
  IndexType insertSorted(const Kmer &target, const Value &value) {
  	setExclusiveLock();
  	IndexType idx = _insertSorted(target,value);
  	unsetExclusiveLock();
  	return idx;
  }
  
  void remove(const Kmer &target) {
  	setExclusiveLock();
  	bool isFound;
  	IndexType idx = find(target, isFound);
  	if (isFound)
  	  remove(idx);
  	unsetExclusiveLock();
  }
  void remove(IndexType idx) {
  	setExclusiveLock();
    resize(size()-1,idx);
    unsetExclusiveLock();
  }
  
  void swap(IndexType idx1, IndexType idx2) {
  	if (idx1 == idx2) 
  	  return;
  	if (idx1 >= size() || idx2 >= size())
  	  throw std::invalid_argument("attempt to access index greater than size in KmerArray swap()");
  	  
  	setExclusiveLock();
  	get(idx1).swap(get(idx2));
  	if (sizeof(ValueType) > 0) {
  	  ValueType tmp = valueAt(idx1);
  	  valueAt(idx1) = valueAt(idx2);
  	  valueAt(idx2) = tmp; 
  	}
  	unsetExclusiveLock();
  }

  std::string toString() {
  	setSharedLock();
  	std::stringstream ss;
  	ss <<  "{";
  	IndexType idx=0;
  	for(idx=0; idx<size() && idx < 30; idx++) {
  		ss << get(idx).toFasta() << ":" << valueAt(idx) << ", ";
  	} 
  	if (idx < size())
  	    ss << " ... " << size() - idx << " more ";
  	ss << "}";
  	unsetSharedLock();
  	return ss.str();
  }

public:
  class Iterator : public std::iterator<std::forward_iterator_tag, KmerArray>
  {
    private:
      KmerArray *_tgt;
      IndexType _idx;
      ElementType thisElement;
    public:
      Iterator(KmerArray *target, IndexType idx = 0): _tgt(target), _idx(idx) {
      	setElement();
      }
      Iterator(const Iterator &copy) { *this = copy; }
      ~Iterator() {}
      Iterator& operator=(const Iterator& other) {
        _tgt = other._tgt;
        _idx = other._idx;
        thisElement = other.thisElement;
        return *this;
      }
      bool operator==(const Iterator& other) const { return _idx == other._idx && _tgt == other._tgt; }
      bool operator!=(const Iterator& other) const { return _idx != other._idx || _tgt != other._tgt; }
      Iterator& operator++() { if ( !isEnd() ) { ++_idx; setElement(); } return *this; }
      Iterator operator++(int unused) { Iterator tmp(*this) ; ++(*this); return tmp; }
      ElementType &operator*() { return thisElement; }
      Kmer &key() { return thisElement.key(); }
      Value &value() { return thisElement.value(); }
      ElementType *operator->() { return &thisElement; }
      bool isEnd() const { return _idx >= _tgt->size(); }
      void setElement() { if ( !isEnd() ) thisElement = _tgt->getElement(_idx); }
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
   typedef Kmer::NumberType NumberType;
   typedef KmerArray<Value> BucketType;
   typedef typename BucketType::Iterator BucketTypeIterator;
   typedef typename BucketType::ElementType ElementType;
   
   typedef std::vector< BucketType > BucketsVector;
   typedef typename BucketsVector::iterator BucketsVectorIterator;

private:
   BucketsVector _buckets;
   NumberType BUCKET_MASK;
   
public:
   KmerMap(IndexType bucketCount = 1024) {
   	
   	 // ensure buckets are a precicise power of two
   	 // with at least bucketCount buckets
   	 NumberType powerOf2 = bucketCount;
   	 if (powerOf2 == 0) {
   	    powerOf2 = 1;
     } else if ( (powerOf2 & (powerOf2 -1)) == 0 ) {
   	 	// argument is a power of 2
   	 } else {
   	 	powerOf2--;
   	 	for (unsigned int i = 1; i < sizeof(NumberType)*8 ; i<<=1)
   	 	    powerOf2 |= powerOf2 >> i;
   	 	powerOf2++;
   	 }  	 
   	 
   	 BUCKET_MASK = powerOf2 - 1;
     _buckets.resize(powerOf2);
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
   
   void setReadOnlyOptimization() {
   	for(unsigned int i = 0; i<_buckets.size(); i++)
   	  _buckets[i].setReadOnlyOptimization();
   }
   void unsetReadOnlyOptimization() {
   	for(unsigned int i = 0; i<_buckets.size(); i++)
   	  _buckets[i].unsetReadOnlyOptimization();
   }
   
   inline unsigned short getLocalThreadId(NumberType hash, unsigned short numThreads) const {
   	 // use the bottom bits of hash which are (used to sort by bucket)
   	 // partition by numThreads blocks
   	 return (hash & BUCKET_MASK) / (_buckets.size() / numThreads + 1);
   }
   inline unsigned short getLocalThreadId(const KeyType &key, unsigned short numThreads) const {
     return getLocalThreadId(key.hash(), numThreads);
   }
   inline unsigned short getDistributedThreadId(NumberType hash, unsigned short numThreads) const {
   	 // use top bits of hash (unused to sort by bucket)
   	 // partition roundrobin by numThreads
   	 return (hash >> 32 & BUCKET_MASK) % numThreads;
   }
   inline unsigned short getDistributedThreadId(const KeyType &key, unsigned short numThreads) const {
   	 return getDistributedThreadId(key.hash(), numThreads);
   }
   
   inline NumberType getBucketIdx(NumberType hash) const {
   	 return hash & BUCKET_MASK;
   }
   inline NumberType getBucketIdx(const KeyType &key) const {
   	 return getBucketIdx(key.hash());
   }
   
   inline BucketType &getBucket(NumberType hash) {
   	return _buckets[getBucketIdx(hash)];
   }
   inline const BucketType &getBucket(NumberType hash) const {
   	return _buckets[getBucketIdx(hash)];
   }
   
   inline BucketType &getBucket(const KeyType &key)  {
     return getBucket(getBucketIdx(key));
   }
   inline const BucketType &getBucket(const KeyType &key) const {
     return getBucket(getBucketIdx(key));
   }

   ElementType insert(const KeyType &key, const ValueType &value, BucketType &bucket) {
   	  bucket.setExclusiveLock();
   	  IndexType idx = bucket.insertSorted(key,value);
   	  ElementType element = bucket.getElement(idx);
   	  bucket.unsetExclusiveLock();
   	  return element;
   }
   ElementType insert(const KeyType &key, const ValueType &value) {
     return insert(key,value, getBucket(key));
   }
   

   bool remove(const KeyType &key, BucketType &bucket) {
   	  bool isFound;
   	  bucket.setExclusiveLock();
   	  IndexType idx = bucket.findSorted(key, isFound);
   	  if (isFound && idx != BucketType::MAX_INDEX)
   	    bucket.remove(idx);
   	  bucket.unsetExclusiveLock();
   	  return isFound;
   }
   bool remove(const KeyType &key) {
   	  return remove(key, getBucket(key));
   }
   
   bool _exists(const KeyType &key, IndexType &existingIdx, const BucketType &bucket) const {
   	 bool isFound;
   	 existingIdx = bucket.findSorted(key, isFound);
   	 return isFound;
   }
   bool _exists(const KeyType &key, IndexType &existingIdx) const {
   	 return _exists(key,existingIdx, getBucket(key));
   }
   
   bool exists(const KeyType &key, const BucketType &bucket) const {
   	 IndexType dummy;
   	 return _exists(key, dummy, bucket);
   }
   bool exists(const KeyType &key) const {
   	 return exists(key, getBucket(key));
   }
   
   ElementType getElementIfExists(const KeyType &key, BucketType &bucket) {
   	 IndexType existingIdx;
   	 ElementType element;
     bucket.setSharedLock();
   	 if (_exists(key, existingIdx, bucket))
   	 	element = bucket.getElement(existingIdx);
   	 bucket.unsetSharedLock();
   	 return element;
   }
   ElementType getElementIfExists(const KeyType &key) {
   	 return getElementIfExists(key, getBucket(key));
   }

   ElementType getElement(const KeyType &key, BucketType &bucket) {
     IndexType existingIdx;
   	 ElementType element;
   	 bucket.setExclusiveLock();
   	 if (_exists(key, existingIdx, bucket))
   	   element = bucket.getElement(existingIdx);
     else 
       element = insert(key, Value(), bucket);
     bucket.unsetExclusiveLock();
     return element;
   }
   ElementType getElement(const KeyType &key) {
   	 return getElement(key, getBucket(key));
   }
   
   // not thread safe!
   ValueType &operator[](const KeyType &key) {
     BucketType &bucket = getBucket(key);
     IndexType existingIdx;
   	 if (_exists(key, existingIdx, bucket))
   	   return bucket.valueAt(existingIdx);
     else 
       return insert(key, Value(), bucket).value();       
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
// Revision 1.64  2009-12-24 00:39:22  cfurman
// getAverageWeight() added
//
// Revision 1.63  2009-12-22 18:31:14  regan
// moved Open MP includes to a central config file
//
// Revision 1.62  2009-12-21 06:34:26  regan
// used openmp and clever partitioning to speed up building spectrum
//
// Revision 1.61  2009-12-18 19:05:09  regan
// added compile-time thread-safety option to kmer classes
//
// Revision 1.60  2009-12-14 05:31:35  regan
// optimized array resizing to malloc at logarithmic stepping
// fixed a bug in KmerArray<>::findSorted
//
// Revision 1.59  2009-11-29 19:04:45  regan
// optimized a bit, fixed a few bugs
//
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
