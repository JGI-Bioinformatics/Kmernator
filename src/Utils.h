// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.4 2009-11-04 18:26:17 regan Exp $
//

#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/error_of.hpp>
#include <boost/accumulators/statistics/error_of_mean.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_density.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_p_square_quantile.hpp>


using namespace boost;
using namespace boost::accumulators;


#include "Kmer.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "MemoryUtils.h"

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
  KmerWeakMap  weak;
  KmerSolidMap solid;
  KmerSpectrum(unsigned long buckets): weak(buckets), solid(buckets/64) {}
  ~KmerSpectrum() {}
  
  void printHistograms(bool solidOnly = false) {

  	double maxCount = log10(TrackingData::maxCount);
  	double maxWeightedCount = log10(TrackingData::maxWeightedCount);
  	unsigned long numDataPoints = std::max(10ul, ((solidOnly ? 0 : weak.size()) + solid.size())/1000);
  
    unsigned long bins = 30;
    unsigned long cacheSize = numDataPoints;//std::max(10ul, numDataPoints / 100000);
    
    typedef 
  	accumulator_set<double, 
  	                stats< tag::variance, 
  	                       tag::density,
  	                       stats<tag::median(with_density)>
  	                       >
  	                > StdAccumulatorType;
  	StdAccumulatorType countAcc( tag::density::num_bins = bins, tag::density::cache_size = cacheSize         );
  	StdAccumulatorType weightedCountAcc( tag::density::num_bins = bins, tag::density::cache_size = cacheSize  );
  	StdAccumulatorType directionAcc( tag::density::num_bins = bins, tag::density::cache_size = cacheSize  );
  	                
  	typedef
   	accumulator_set<double, 
  	                stats< tag::weighted_variance,
  	                       tag::weighted_density 
  	                       >,
  	                double
  	                > WeightedAccumulatorType;
  	WeightedAccumulatorType countWAcc( tag::weighted_density::num_bins = bins, tag::weighted_density::cache_size = cacheSize         );
   	WeightedAccumulatorType weightedCountWAcc( tag::weighted_density::num_bins = bins, tag::weighted_density::cache_size = cacheSize  );
  	WeightedAccumulatorType directionWAcc( tag::weighted_density::num_bins = bins, tag::weighted_density::cache_size = cacheSize  );
   	  	
   	if (!solidOnly) {
  	 for(KmerWeakMap::Iterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
  	  TrackingData &data = it->value().value;
  	  countAcc( std::max(-1.0, (double)(data.getCount()) ) );
  	  weightedCountAcc( std::max(-1.0, (double)(data.getWeightedCount())) );
  	  directionAcc( data.getNormalizedDirectionBias() );

  	  countWAcc(  std::max(-1.0, (double)(data.getCount()) ),               weight = data.getCount() );
  	  weightedCountWAcc( std::max(-1.0, (double)(data.getWeightedCount())), weight = data.getCount() );
  	  directionWAcc( data.getNormalizedDirectionBias(),                     weight = data.getCount() );
  	 }
  	}
  	for(KmerSolidMap::Iterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++) {
  	  TrackingData &data = it->value().value;
  	  //std::cerr << data.getCount() << std::endl;
  	  countAcc( std::max(-1.0, (double)(data.getCount()) ) );
  	  weightedCountAcc( std::max(-1.0, (double)(data.getWeightedCount())) );
  	  directionAcc( data.getNormalizedDirectionBias() );

  	  countWAcc(  std::max(-1.0, (double)(data.getCount()) ),               weight = data.getCount() );
  	  weightedCountWAcc( std::max(-1.0, (double)(data.getWeightedCount())), weight = data.getCount() );
  	  directionWAcc( data.getNormalizedDirectionBias(),                     weight = data.getCount() );
  	}
  	
    typedef iterator_range< std::vector< std::pair<double, double> >::iterator > HistogramType;
  	HistogramType countHist          = density(countAcc);
  	HistogramType weightedCountHist  = density(weightedCountAcc);
  	HistogramType directionHist      = density(directionAcc);
   	
   	HistogramType countHistW         = weighted_density(countWAcc);
  	HistogramType weightedCountHistW = weighted_density(weightedCountWAcc);
  	HistogramType directionHistW     = weighted_density(directionWAcc);
   	
  	std::cerr << "Counts, Weights and Directions" << std::endl;
  	
  	std::cerr << "Counts:\t";
  	std::cerr << std::fixed << std::setprecision(2) << mean(countAcc) << "/" << median(countAcc) << " +-";
  	std::cerr << std::fixed << std::setprecision(2) << sqrt(variance(countAcc)) << "\t";
  	std::cerr << std::fixed << std::setprecision(2) << weighted_mean(countWAcc) << " +-";
  	std::cerr << std::fixed << std::setprecision(2) << sqrt(weighted_variance(countWAcc)) << "\t";
  
  	std::cerr << "WeightedCounts:\t";
  	std::cerr << std::fixed << std::setprecision(2) << mean(weightedCountAcc) << "/" << median(weightedCountAcc) << " +-";
  	std::cerr << std::fixed << std::setprecision(2) << sqrt(variance(weightedCountAcc)) << "\t";
  	std::cerr << std::fixed << std::setprecision(2) << weighted_mean(weightedCountWAcc) << " +-";
  	std::cerr << std::fixed << std::setprecision(2) << sqrt(weighted_variance(weightedCountWAcc)) << "\t";
  	
  	std::cerr << std::endl;
  	
  	for (int i=0; i < countHist.size(); i++) {
  		std::cerr << std::fixed << std::setprecision(2) << countHist[i].first << "\t";
  		std::cerr << std::fixed << std::setprecision(2) << countHist[i].second*100.0<< "\t";
  		//assert(countHist[i].first == countHistW[i].first);
  		std::cerr << std::fixed << std::setprecision(2) << countHistW[i].second *100.0<< "\t"; 
  		
  		std::cerr << std::fixed << std::setprecision(2) << weightedCountHist[i].first  << "\t";
  		std::cerr << std::fixed << std::setprecision(2) << weightedCountHist[i].second*100.0 << "\t";
  		//assert(weightedCountHist[i].first == weightedCountHistW[i].first);
  		std::cerr << std::fixed << std::setprecision(2) << weightedCountHistW[i].second*100.0 << "\t";
  		
  		std::cerr << std::fixed << std::setprecision(2) << directionHist[i].first  << "\t";
  		std::cerr << std::fixed << std::setprecision(2) << directionHist[i].second *100.0<< "\t"; 
  		//assert(directionHist[i].first == directionHistW[i].first);
   		std::cerr << std::fixed << std::setprecision(2) << directionHistW[i].second *100.0<< "\t"; 
   		
  		std::cerr << std::endl;
  	}
  }
  
  unsigned long promote(double probabilityQuantile) {
  	unsigned long promoted = 0;
  	double maxCount = log10(TrackingData::maxCount);
  	double maxWeightedCount = log10(TrackingData::maxWeightedCount);
  	unsigned long numDataPoints = std::max(10ul, (weak.size() + solid.size())/1000);
  
    unsigned long bins = 30;
    unsigned long cacheSize = numDataPoints;
    
  	typedef accumulator_set<
  	                unsigned long, 
  	                stats< tag::weighted_p_square_quantile         
  	                       >,
  	                double
  	                > AccumulatorType;
  	int count = 11;
  	double probabilities[] = { probabilityQuantile, 0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05 };
  	std::vector< AccumulatorType > accumulators;
  	for(int i=0; i < count; i++)
  	   accumulators.push_back( AccumulatorType( quantile_probability = probabilities[i] ) );
  	for(KmerWeakMap::Iterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++)
  	   for(int i=0; i< count; i++) {
  	   	  double count = it->value().value.getCount();
  	      accumulators[i]( count, weight = count );
  	   }
  	for(KmerSolidMap::Iterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
  	   for(int i=0; i< count; i++) {
  	   	  double count = it->value().value.getCount();
  	      accumulators[i]( count, weight = count );
  	   }
  	
  	std::cerr << "Quantiles: ";
    for(int i=0; i< count; i++) {
       std::cerr << std::fixed << std::setprecision(4) << probabilities[i] << ": ";
       std::cerr << std::fixed << std::setprecision(4) << weighted_p_square_quantile(accumulators[i]) << "\t";
    }
    std::cerr << std::endl;
  	
  	double minimumCount = weighted_p_square_quantile(accumulators[0]);
    for(KmerWeakMap::Iterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++)
      if (it->value().value.getCount() > minimumCount) {
      	solid[ it->key() ].value = it->value().value;
      	promoted++;
      }
  	for(KmerSolidMap::Iterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
  	  weak.remove( it->key() );

  	return promoted;
  }
  
  std::string contrastSpectrums(KmerSpectrum &reference, double minSolidDepth = -1.0) {
  	// ignores reference.weak and everything in reference.solid passed criteria
  	 
  	std::stringstream ss;
  	unsigned long thisSolidSize(0), thisWeakSize(0), referenceSolidSize(0), 
  	              matchingSolid(0), missingSolid(0), incorrectWeak(0), extraSolid(0), correctWeak(0);
  	
  	thisWeakSize = weak.size();
  	// compare this->solid to reference.solid 
  	for(KmerSolidMap::Iterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++) {
  	  if (minSolidDepth >= 0 && it->value().value.getCount() < minSolidDepth) {
  	  	thisWeakSize++;
  	    if ( reference.solid.exists( it->key() ) ) {
  	      // exists in reference but in this is weak: incorrectWeak
  	      incorrectWeak++;
  	      ss << "IW " << it->key().toFasta() << " " << it->value().value << std::endl;
  	    } else {
  	      // Weak and not in reference: correctWeak
  	      // count later
  	    }
  	  } else {
  	    thisSolidSize++;
  	    if ( reference.solid.exists( it->key() ) ) {
  	      // correctly matches solid in reference: matchingSolid
  	  	  matchingSolid++;
  	    } else {
  	      // exists in this but not reference: extraSolid
  	      extraSolid++;
  	      ss << "ES  " << it->key().toFasta() << " " << it->value().value << std::endl;
  	    }
  	  }
  	}
  	
  	// compare reference.solid to this->solid and this->weak
  	for(KmerSolidMap::Iterator it(reference.solid.begin()), itEnd(reference.solid.end()); it != itEnd; it++) {
  	  referenceSolidSize++;
  	  
  	  if ( solid.exists( it->key() ) ) {
  	  	// already accounted  	  	
  	  } else if ( weak.exists ( it->key() )) {
      	// exists in reference but in this is weak: incorrectWeak
      	incorrectWeak++;
      	ss << "IW " << it->key().toFasta() << " " << it->value().value << std::endl;
  	  } else {
  	  	// exists in reference but not in either solid nor weak: missingSolid
  	    ss << "MS " << it->key().toFasta() << " N/A " << std::endl;
  	    missingSolid++;
  	  }
  	}
  	correctWeak = thisWeakSize - incorrectWeak;

  	ss << "SolidSizes: " << thisSolidSize << " / " << referenceSolidSize << std::endl;
  	ss << "MatchingSolid: " << matchingSolid << " MissingSolid: " << missingSolid << " IncorrectWeak: " << incorrectWeak << " ExtraSolid: " << extraSolid << std::endl;
  	ss << "CorrectWeak: " << correctWeak << " / " << thisWeakSize << std::endl;
   	return ss.str();
  }
};

void printStats(unsigned long pos, KmerSpectrum &stats, bool solidOnly = false) {
	stats.printHistograms(solidOnly);
	KmerSolidMap::Iterator it = stats.solid.begin();
    std::cerr << pos << " reads, " << stats.solid.size() << " / " << stats.weak.size() << " kmers so far ";
    for(int i=0; i < 5 && it != stats.solid.end(); i++,it++)
        std::cerr << it.bucket().toString() << "; ";
    std::cerr << " minDepth: " << (unsigned long)TrackingData::minimumDepth <<  std::endl;
    std::cerr << MemoryUtils::getMemoryUsage() << std::endl;     
}

void buildKmerSpectrum( ReadSet &store, KmerSpectrum &spectrum ) 
{
	KmerWeakMap  &weak  = spectrum.weak;
	KmerSolidMap &solid = spectrum.solid;
	weak.clear();
	solid.clear();

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
       	  	TrackingData &data = weak[ *least ].value;
       	  	
       	  	data.track( kmers.valueAt(j), keepDirection );

       	  }
       }
       if (i % 1000000 == 0) {
       	 printStats(i, spectrum);  
       }
    }
    printStats(store.getSize(), spectrum);    
    unsigned long promoted = spectrum.promote(0.05);
    std::cerr << "Promoted " << promoted << " kmers" << std::endl;
    spectrum.printHistograms(true);

};

#endif

//
// $Log: Utils.h,v $
// Revision 1.4  2009-11-04 18:26:17  regan
// refactored
// added statistics calculations and histograms
//
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
