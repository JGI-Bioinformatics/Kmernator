// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.3 2009-11-03 17:15:40 regan Exp $
//

#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <iostream>
#include <cstdlib>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
using namespace boost::accumulators;


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

class KmerSpectrum
{
public:
  KmerSolidMap weak, solid;
  KmerSpectrum(unsigned long buckets): weak(buckets), solid(buckets/64) {}
  ~KmerSpectrum() {}
};

void printStats(unsigned long pos, KmerSpectrum &stats) {
	KmerSolidMap::Iterator it = stats.solid.begin();
    std::cerr << pos << " reads, " << stats.solid.size() << " / " << stats.weak.size() << " kmers so far ";
    for(int i=0; i < 5 && it != stats.solid.end(); i++,it++)
        std::cerr << it.bucket().toString() << "; ";
     std::cerr << " minDepth: " << (unsigned long)SolidTrackingData::minimumDepth <<  std::endl;         
}

void buildKmerSpectrum( ReadSet &store, KmerSpectrum &spectrum ) 
{
	KmerSolidMap &weak  = spectrum.weak;
	KmerSolidMap &solid = spectrum.solid;
	weak.clear();
	solid.clear();
	   
    unsigned long numBuckets = estimateWeakKmerBucketSize( store, 64 );

    std::cerr << "targetting " << numBuckets << std::endl;
 
	KmerArray<> tmp;
	tmp.resize(2);
	KmerPtr least;
	
	for (int i=0 ; i < store.getSize(); i++)
    {
       KmerWeights kmers = buildWeightedKmers(store.getRead(i));
  
       for (int j=0; j < kmers.size(); j++)
       {
       	  
       	  TwoBitSequence::reverseComplement((TwoBitEncoding*)kmers[j].get(), (TwoBitEncoding*)tmp[0].get(), KmerSizer::getSequenceLength());
       	  
       	  bool keepDirection = kmers[j] < tmp[0];
       	  
       	  if (keepDirection)
       	     least = kmers[j].get();
       	  else
       	     least = tmp[0].get();
       	  
       	  if ( solid.exists( *least ) ) {
       	  	// track solid stats
       	  	solid[ *least ].value.track( kmers.valueAt(j), keepDirection );
       	  	
       	  } else {
       	  	// track weak stats
       	  	SolidTrackingData &data = weak[ *least ].value;
       	  	
       	  	data.track( kmers.valueAt(j), keepDirection );
       	  	if ( data.count > SolidTrackingData::minimumDepth ) {
          	  // track stats and pop out of weak hash
          	  solid[ *least ].value = data;          	
          	  weak.remove(*least);
       	    }
       	  }
       }
       if (i % 1000000 == 0) {
       	 printStats(i, spectrum);  
       }
    }
    printStats(store.getSize(), spectrum);    

};

#endif

//
// $Log: Utils.h,v $
// Revision 1.3  2009-11-03 17:15:40  regan
// minor refactor
//
// Revision 1.2  2009-11-02 21:19:25  regan
// fixed types and boundary tests
//
// Revision 1.1  2009-11-02 18:47:34  regan
// added some code without a permanent home
//
//
