// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.8 2009-11-10 07:05:37 regan Exp $
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

#include "TwoBitSequence.h"
#include "Kmer.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "MemoryUtils.h"
#include "Options.h"

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
	if (baseCount == 0 || store.getSize() == 0)
	  return 128;
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



KmerWeights buildWeightedKmers(Read &read) {
  KmerWeights kmers(read.getTwoBitSequence(), read.getLength());
  std::string quals = read.getQuals();
  SequenceLengthType markupIdx = 0;
  
  BaseLocationVectorType markups = read.getMarkups();
  double weight = 0.0;
  double change = 0.0;
  
  for(SequenceLengthType i=0; i< kmers.size(); i++) {
  	if (i%100 == 0 || weight == 0.0) {
  	  weight = 1.0;
  	  for(SequenceLengthType j=0; j < KmerSizer::getSequenceLength(); j++) 
  	    weight *= Read::qualityToProbability[ quals[i+j] ];
  	} else {
  		change = Read::qualityToProbability[ quals[i+KmerSizer::getSequenceLength()-1] ] / 
  		         Read::qualityToProbability[ quals[i-1] ];
  		weight *= change;
  	}
  	while( markupIdx < markups.size() && markups[markupIdx].second < i)
  	  markupIdx++;
  	if (markupIdx < markups.size() && markups[markupIdx].second < i+KmerSizer::getSequenceLength()) {
  	    weight = 0.0;
//  	    std::cerr << "markupAt " << markups[markupIdx].first << " " << markups[markupIdx].second << "\t";
    }
  	kmers.valueAt(i) = weight;
//  	std::cerr << i << ":" << std::fixed << std::setprecision(3) << quals[i+KmerSizer::getSequenceLength()-1] << "-" << change << "/" << weight << "\t";
  	
  }
//  std::cerr << std::endl;
  return kmers;
}


class KmerSpectrum
{
public:
  KmerWeakMap  weak;
  KmerSolidMap solid;
  KmerSpectrum(unsigned long buckets): weak(buckets), solid(buckets/4) 
  {
  	// set the minimum weight that will be used to track kmers
  	// based on the given options
  	TrackingData::minimumWeight = Options::getMinKmerQuality();
  }
  ~KmerSpectrum() {}
  KmerSpectrum(const KmerSpectrum &copy) {
  	*this = copy;
  }
  KmerSpectrum &operator=(const KmerSpectrum &other) {
  	this->weak = other.weak;
  	this->solid = other.solid;
  }
  
  void printHistograms(bool solidOnly = false) {

  	double maxCount = (TrackingData::maxCount);
  	double maxWeightedCount = (TrackingData::maxWeightedCount);
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
    
  	typedef accumulator_set<
  	                unsigned long, 
  	                stats< tag::weighted_p_square_quantile         
  	                       >,
  	                double
  	                > AccumulatorType;
  	int count = 12;
  	double probabilities[] = { probabilityQuantile, 1.0-probabilityQuantile, 0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05 };
  	std::vector< AccumulatorType > countAcc, weightedCountAcc, directionAcc;
  	for(int i=0; i < count; i++) {
  	   countAcc.push_back( AccumulatorType( quantile_probability = probabilities[i] ) );
  	   weightedCountAcc.push_back( AccumulatorType( quantile_probability = probabilities[i] ) );
  	   directionAcc.push_back( AccumulatorType( quantile_probability = probabilities[i] ) );
  	}
  	for(KmerWeakMap::Iterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++)
  	   for(int i=0; i< count; i++) {
  	   	  TrackingData &data = it->value().value;
  	   	  double count = data.getCount();
  	      countAcc[i]( count, weight = count );
  	      weightedCountAcc[i]( data.getWeightedCount(), weight = count );
  	      directionAcc[i]( data.getNormalizedDirectionBias(), weight = count);
  	   }
  	for(KmerSolidMap::Iterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
  	   for(int i=0; i< count; i++) {
  	   	  TrackingData &data = it->value().value;
  	   	  double count = data.getCount();
  	      countAcc[i]( count, weight = count );
  	      weightedCountAcc[i]( data.getWeightedCount(), weight = count );
  	      directionAcc[i]( data.getNormalizedDirectionBias(), weight = count);
  	   }
  	
  	std::cerr << "Quantiles: ";
    for(int i=0; i< count; i++) {
       std::cerr << std::fixed << std::setprecision(4) << probabilities[i] << ": ";
       std::cerr << std::fixed << std::setprecision(4) << weighted_p_square_quantile(countAcc[i]) << "\t";
       std::cerr << std::fixed << std::setprecision(4) << weighted_p_square_quantile(weightedCountAcc[i]) << "\t";
       std::cerr << std::fixed << std::setprecision(4) << weighted_p_square_quantile(directionAcc[i]) << "\t";
       std::cerr << std::fixed << std::setprecision(4) << weighted_p_square_quantile(directionAcc[i]) << "\t";
       std::cerr << std::endl;
    }
  	double minimumCount = weighted_p_square_quantile(countAcc[0]);
  	double minimumWeightedCount = weighted_p_square_quantile(weightedCountAcc[0]);
  	double minimumDirBias = weighted_p_square_quantile(directionAcc[0]);
  	double maximumDirBias = weighted_p_square_quantile(directionAcc[1]);
  	
    for(KmerWeakMap::Iterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
      TrackingData &data = it->value().value;
      double dirBias = data.getNormalizedDirectionBias();
      if (data.getCount() > minimumCount && data.getWeightedCount() > minimumWeightedCount
          //&& dirBias < maximumDirBias && dirBias > minimumDirBias
          //&& dirBias < 0.8 && dirBias > 0.2
          ) {
        std::pair<double,double> scores = getPermutedScores( it->key(), Options::getFirstOrderWeight(), Options::getSecondOrderWeight() );
        double permuteScore = (scores.first+scores.second);
        if ( (double)data.getCount() / permuteScore > 1.0 ) {
      	  solid[ it->key() ].value = it->value().value;
      	  promoted++;
        }
      }
    }
  	for(KmerSolidMap::Iterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
  	  weak.remove( it->key() );

  	return promoted;
  }
  
  
  std::pair<double,double> getPermutedScores( KmerPtr kmer, double permutationWeight = 0.1, double secondWeight = 0.0 ) {
  	std::pair<double,double> score(0.0,0.0);
  	//std::cerr << "Permuting " << kmer->toFasta() << std::endl;
  	bool isSolid; double base = 0.0;
  	if ( solid.exists( *kmer ) ) {
  	  base = solid[ *kmer ].value.getCount();
  	  isSolid = true;
  	} else if (weak.exists( *kmer ) ) {
  	  base = weak[ *kmer ].value.getCount();
  	  isSolid = false;
  	}
  	KmerWeights permutations = KmerWeights::permuteBases(kmer, true);
  	for(int i=0; i<permutations.size(); i++) {
      //std::cerr << "Looking at " << permutations[i].toFasta() << std::endl;		
  	  if( solid.exists( permutations[i] ) ) {
  	  	score.first += solid[ permutations[i] ].value.getCount() * permutationWeight;
  	  } else if ( weak.exists( permutations[i] ) ) {
  	  	score.second += weak[ permutations[i] ].value.getCount() * permutationWeight;
  	  }
  	  if (secondWeight > 0.0) {
  	    std::pair<double,double> sScores = getPermutedScores( permutations[i], secondWeight );
  	    score.first += sScores.first;
  	    score.second += sScores.second;
  	    if (isSolid)
  	      score.first -= secondWeight * base;
  	    else 
  	      score.second -= secondWeight * base;
  	  }
  	}
  	return score;
  }
  
  int countGC(std::string fasta) {
     int count=0;
     for(int i=0; i<fasta.length(); i++)
       if (fasta[i] == 'G' || fasta[i] == 'C')
         count++;
     return count;
  }
  std::string pretty( KmerPtr kmer, std::string note ) {
  	std::stringstream ss;
  	
  	TwoBitEncoding _rev[ KmerSizer::getByteSize() ];
  	KmerPtr rev(&_rev);
  	kmer->buildReverseComplement( *rev );
  	
  	ss << kmer->toFasta() << " " << rev->toFasta() << " " << countGC(kmer->toFasta()) << "\t" << note << std::endl;
  	return ss.str();
  }
  std::string pretty( KmerPtr kmer, TrackingData &data ) {
  	std::stringstream ss;
  	ss << data << "\t";

  	std::pair<double,double> scores = getPermutedScores( kmer, 0.1, 0.01 );
  	ss << std::fixed << std::setprecision(2) << scores.first  << "\t";
  	ss << std::fixed << std::setprecision(2) << scores.second << "\t";
  	ss << std::fixed << std::setprecision(2) << (double)data.getCount() / (double)(scores.first+scores.second) << "\t";
  	
  	return pretty( kmer, ss.str() );
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
  	      ss << "IW " << pretty( it->key(), it->value().value);
  	    } else {
  	      // Weak and not in reference: correctWeak
  	      // count later
  	      if (Options::getVerbosity() > 1)
  	        ss << "CW " << pretty( it->key(), it->value().value);
  	    }
  	  } else {
  	    thisSolidSize++;
  	    if ( reference.solid.exists( it->key() ) ) {
  	      // correctly matches solid in reference: matchingSolid
  	  	  matchingSolid++;
  	  	  if (Options::getVerbosity() > 0)
  	  	    ss << "CS " << pretty( it->key(), it->value().value);
  	    } else {
  	      // exists in this but not reference: extraSolid
  	      extraSolid++;
  	      ss << "ES " << pretty( it->key(), it->value().value);
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
      	ss << "IW " << pretty( it->key(), weak[ it->key() ].value );
  	  } else {
  	  	// exists in reference but not in either solid nor weak: missingSolid
  	    ss << "MS " << pretty( it->key(), "N/A");
  	    missingSolid++;
  	  }
  	}
  	if (Options::getVerbosity() > 2) {
  	  for(KmerWeakMap::Iterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
  	    if (! reference.solid.exists( it->key() ))
  	      ss << "CW " << pretty( it->key(), it->value().value);
  	  }
  	}
  	correctWeak = thisWeakSize - incorrectWeak;

    std::stringstream header;
  	header << "SolidSizes: " << thisSolidSize << " / " << referenceSolidSize << std::endl;
  	header << "MatchingSolid: " << matchingSolid << " MissingSolid: " << missingSolid << " IncorrectWeak: " << incorrectWeak << " ExtraSolid: " << extraSolid << std::endl;
  	header << "CorrectWeak: " << correctWeak << " / " << thisWeakSize << std::endl;
   	return header.str() + ss.str();
  }
  
  void append( KmerWeights &kmers, bool isSolid = false ) {
      
     TwoBitEncoding _tmp[KmerSizer::getByteSize()];
     KmerPtr least(&_tmp);
     
     
     for (int j=0; j < kmers.size(); j++)
     {
     	bool keepDirection = kmers[j].buildLeastComplement(*least);
       	
       	if ( solid.exists( *least ) || isSolid ) {
       	  	// track solid stats
       	  	solid[ *least ].value.track( kmers.valueAt(j), keepDirection );
       	} else {
       	  	// track weak stats
       	  	TrackingData &data = weak[ *least ].value;
       	  	
       	  	data.track( kmers.valueAt(j), keepDirection );
       	}
     }  	
  }
  
};

void printStats(unsigned long pos, KmerSpectrum &stats, bool solidOnly = false) {
	stats.printHistograms(solidOnly);
	KmerSolidMap::Iterator it = stats.solid.begin();
    std::cerr << pos << " reads, " << stats.solid.size() << " / " << stats.weak.size() << " kmers so far ";
    if(it != stats.solid.end())
        std::cerr << it.bucket().toString() << "; ";
    std::cerr << " minDepth: " << (unsigned long)TrackingData::minimumDepth <<  std::endl;
    std::cerr << MemoryUtils::getMemoryUsage() << std::endl;     
}

void buildKmerSpectrum( ReadSet &store, KmerSpectrum &spectrum, bool isSolid = false ) 
{
	KmerWeakMap  &weak  = spectrum.weak;
	KmerSolidMap &solid = spectrum.solid;
	weak.reset();
	solid.reset();
	
	for (int i=0 ; i < store.getSize(); i++)
    {
       KmerWeights kmers = buildWeightedKmers(store.getRead(i));
  
       spectrum.append(kmers, isSolid);

       if (i % 1000000 == 0) {
       	 printStats(i, spectrum);  
       }
    }
    printStats(store.getSize(), spectrum);    
    
}

void experimentOnSpectrum( KmerSpectrum &spectrum ) {
    KmerSpectrum everything(spectrum);
    
    unsigned long promoted = spectrum.promote(0.05);
    std::cerr << "Promoted " << promoted << " kmers" << std::endl;
    spectrum.printHistograms(true);
    
    for (double prob = 0.10 ; prob > 0.0 ; prob -= 0.005) {
    	KmerSpectrum copy(everything);
    	std::cerr << "Testing " << prob << std::endl;
    	promoted = copy.promote(prob);
    	std::cerr << "Promoted " << promoted << " kmers" << std::endl;
    	copy.printHistograms(true);
        std::cerr << copy.contrastSpectrums(spectrum) << std::endl;;
    }
    
    for(KmerWeakMap::Iterator it(everything.weak.begin()), itEnd(everything.weak.end()); it != itEnd; it++) 
    	everything.solid[ it->key() ].value = it->value().value;
    everything.weak.clear();
    
    for (double minDepth = 1; minDepth < 10 ; minDepth++) {
    	std::cerr << "Comparing minDepth " << minDepth << std::endl;
    	std::cerr << everything.contrastSpectrums(spectrum, minDepth) << std::endl;
    }
};

#endif

//
// $Log: Utils.h,v $
// Revision 1.8  2009-11-10 07:05:37  regan
// changes for debugging
//
// Revision 1.7  2009-11-09 19:37:17  regan
// enhanced some debugging / analysis output
//
// Revision 1.6  2009-11-06 16:59:11  regan
// added base substitution/permutations table and build function
//
// Revision 1.5  2009-11-06 04:10:21  regan
// refactor of cmd line option handling
// added methods to evaluate spectrums
//
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
