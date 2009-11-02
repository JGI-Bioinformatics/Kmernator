// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.1 2009-11-02 18:47:34 regan Exp $
//

#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <iostream>
#include <fstream>

#include "Kmer.h"
#include "Sequence.h"
#include "ReadSet.h"

template<typename Raw, typename Store>
class BucketedData
{
public:
  typedef Store StoreType;
  typedef Raw   RawType;
  
private:
  RawType _minValue, _maxValue;
  StoreType steps;
  
public:
  // TODO inclusive permutations: (), [), (], []
  //      presently it truncates fraction, so [)
  // TODO implement log scale
  BucketedData(RawType minValue, RawType maxValue): 
    _minValue(minValue), 
    _maxValue(maxValue)
  {
  	// assumes unsigned Store...
  	steps = ((StoreType) -1);
  }
  ~BucketedData() {}
  
  StoreType getStore(RawType value) {
  	if (value < _minValue || value > _maxValue)
  	  throw std::invalid_argument("getStore() out of range");
    return ((RawType) steps) * (value - _minValue) / (_maxValue - _minValue);
  }
  
  RawType getValue(StoreType store) {
  	return (_maxValue - _minValue) * ((RawType) store) / ((RawType) steps);
  }
};

typedef unsigned char   OneByte;
typedef unsigned short  TwoByte;
typedef unsigned int    FourByte;
typedef unsigned long   EightByte;


typedef BucketedData< double, OneByte   > DoubleToOneByte;
typedef BucketedData< double, TwoByte   > DoubleToTwoByte;
typedef BucketedData< double, FourByte  > DoubleToFourByte;
typedef BucketedData< double, EightByte > DoubleToEightByte;


unsigned long estimateWeakKmerBucketSize( ReadSet &store, unsigned long targetKmersPerBucket ) {
	unsigned long baseCount = store.getBaseCount();
	unsigned long avgSequenceLength = baseCount / store.getSize();
	unsigned long kmersPerRead = (avgSequenceLength - KmerSizer::getSequenceLength() + 1);
	if (kmersPerRead > avgSequenceLength + 1)
	  throw std::invalid_argument("Attempting to use a kmer greater than the avgerage sequence length");
	unsigned long rawKmers = kmersPerRead * store.getSize();
	unsigned long estimatedUniqueKmers =  rawKmers * ( std::max(pow(1.01, KmerSizer::getSequenceLength()),2.0) - 1.0 );
	unsigned long targetBuckets = estimatedUniqueKmers / targetKmersPerBucket;
	unsigned long maxBuckets = 128*1024*1024;
	unsigned long minBuckets = 128;
	return std::max( std::min(targetBuckets,maxBuckets), minBuckets );
}

class SolidTrackingData 
{
public:
  static const double minimumWeight = 0.1;
  
  unsigned short count;
  unsigned short directionBias;
  float weightedCount;
  SolidTrackingData(): count(0), directionBias(0), weightedCount(0.0) {}
  ~SolidTrackingData() {}
  
  bool track(double weight, bool forward) {
  	if (weight < minimumWeight)
  	  return false;
    if (count < -1) {
      count++;
      if (forward) 
        directionBias++;
      weightedCount += weight;
      return true;
    } else
      return false;
  }
};
std::ostream &operator<<(std::ostream &stream, SolidTrackingData ob)
{
  stream << ob.count << ':' << ((double)ob.directionBias / (double)ob.count) << ':' << ((double)ob.count / ob.weightedCount);
  return stream;
}


class SolidKmerTag : public KmerValue<SolidTrackingData> {};

typedef KmerArray<double> KmerWeights;

KmerWeights buildWeightedKmers(Read &read) {
  KmerWeights kmers(read.getTwoBitSequence(), read.getLength());
  std::string quals = read.getQuals();
  BaseLocationVectorType markups = read.getMarkups();
  for(int i=0; i< kmers.size(); i++) {
  	double weight = 1.0;
  	for(int j=0; j < KmerSizer::getSequenceLength(); j++) 
  	  weight *= Read::qualityToProbability[ quals[i+j] ];
  	for(int j=0; j < markups.size(); j++)
  	  if (markups[j].second >= i && markups[j].second < i+KmerSizer::getSequenceLength())
  	    weight = 0.0;
  	kmers.valueAt(i) = weight;
  }
  return kmers;
}


typedef KmerMap<SolidKmerTag>   KmerSolidMap;
typedef KmerMap<WeakKmerTag>    KmerWeakMap;
typedef KmerMap<unsigned short> KmerCountMap;

#endif

//
// $Log: Utils.h,v $
// Revision 1.1  2009-11-02 18:47:34  regan
// added some code without a permanent home
//
//
