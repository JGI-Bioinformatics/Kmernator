// $Header: /repository/PI_annex/robsandbox/KoMer/src/KmerSpectrum.h,v 1.7 2009-12-21 06:34:26 regan Exp $

#ifndef _KMER_SPECTRUM_H
#define _KMER_SPECTRUM_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

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

template< typename S, typename W>
class KmerSpectrum
{
public:

  typedef std::pair<double,double> DoublePairType;
  typedef DoublePairType           SolidWeakWeightType;
  typedef DoublePairType           MeanType;
  
  typedef std::vector< DoublePairType > DoublePairVectorType;
  typedef DoublePairVectorType          MeanVectorType;

  typedef S                              SolidDataType;
  typedef W                              WeakDataType;
  typedef TrackingDataSingleton          SingletonDataType;
  
  typedef typename WeakDataType::ReadPositionWeightVector           WeakReadPositionWeightVector;
  typedef typename WeakDataType::ReadPositionWeightVector::iterator WeakReadPositionWeightVectorIterator;

  typedef KmerMap<SolidDataType>      SolidMapType;
  typedef KmerMap<WeakDataType>       WeakMapType;
  typedef KmerMap<SingletonDataType>  SingletonMapType;
  
  typedef typename SolidMapType::Iterator        SolidIterator;
  typedef typename SolidMapType::ElementType     SolidElementType;
  typedef typename WeakMapType::Iterator         WeakIterator;
  typedef typename WeakMapType::ElementType      WeakElementType;
  typedef typename SingletonMapType::Iterator    SingletonIterator;
  typedef typename SingletonMapType::ElementType SingletonElementType;
  
  typedef accumulator_set<double, 
  	                stats< tag::variance(lazy), 
  	                       tag::mean
  	                       >
  	                > StdAccumulatorType;
  
public:
  SolidMapType                   solid;
  WeakMapType                    weak;
  SingletonMapType               singleton;
  bool                           hasSolids;

public:
  KmerSpectrum(unsigned long buckets): solid(buckets/8), weak(buckets/4), singleton(buckets), hasSolids(false)
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
  	this->singleton = other.singleton;
  	this->hasSolids  = other.hasSolids;
  }
  
  
  bool hasSolid( const Kmer &kmer ) const {
  	return hasSolids && solid.exists( kmer );
  }
  SolidElementType getIfExistsSolid( const Kmer &kmer ) {
  	if (hasSolids) 
  	  return solid.getElementIfExists( kmer );
  	else
  	  return SolidElementType();
  }
  SolidElementType getSolid( const Kmer &kmer ) {
  	hasSolids = true;
  	return solid.getElement(kmer);
  }
  
  bool hasWeak( const Kmer &kmer ) const {
  	return weak.exists( kmer );
  }
  WeakElementType getIfExistsWeak( const Kmer &kmer ) {
  	return weak.getElementIfExists( kmer );
  }
  WeakElementType getWeak( const Kmer &kmer ) {
  	return weak.getElement( kmer );
  }
  
  bool hasSingleton( const Kmer &kmer ) const {
  	return singleton.exists( kmer );
  }
  SingletonElementType getIfExistsSingleton( const Kmer &kmer )  {
  	return singleton.getElementIfExists( kmer );
  }
  SingletonElementType getSingleton( const Kmer &kmer ) {
  	return singleton.getElement( kmer );
  }

  class DataPointers
  {
  public:
    KmerSpectrum          *spectrum;
    SolidElementType      solid;
    WeakElementType       weak;
    SingletonElementType  singleton;
    DataPointers(KmerSpectrum &_spectrum) : spectrum(&_spectrum) {}
    DataPointers(KmerSpectrum &_spectrum, const Kmer &kmer): spectrum(&_spectrum) {
    	set(kmer);
    }
    void reset() {
    	solid = SolidElementType();
    	weak = WeakElementType();
    	singleton = SingletonElementType();
    }
    void set(const Kmer &kmer) {
    	reset();
    	if (spectrum->hasSolids)
    	  solid = spectrum->getIfExistsSolid(kmer);
    	  
    	if (!solid.isValid()) {
    	  weak = spectrum->getIfExistsWeak(kmer);
    	  if (!weak.isValid())
    	    singleton = spectrum->getIfExistsSingleton(kmer);
    	}
    }
  };
  
  void printHistograms(bool solidOnly = false) {

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
  	 for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
  	  WeakDataType &data = it->value();
  	  countAcc( std::max(-1.0, (double)(data.getCount()) ) );
  	  weightedCountAcc( std::max(-1.0, (double)(data.getWeightedCount())) );
  	  directionAcc( data.getNormalizedDirectionBias() );

  	  countWAcc(  std::max(-1.0, (double)(data.getCount()) ),               weight = data.getCount() );
  	  weightedCountWAcc( std::max(-1.0, (double)(data.getWeightedCount())), weight = data.getCount() );
  	  directionWAcc( data.getNormalizedDirectionBias(),                     weight = data.getCount() );
  	 }
  	}
  	for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++) {
  	  SolidDataType &data = it->value();
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
  	std::cerr << std::fixed << std::setprecision(2) << mean(countWAcc) << " +-";
  	std::cerr << std::fixed << std::setprecision(2) << sqrt(variance(countWAcc)) << "\t";
  
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
  
  unsigned long autoPromote(unsigned int minKmerCount = 3, double minKmerWeightedCount = 0.0,
                            double minWeakRatio = 0.50, double minSolidRatio = 0.05)
  {
  	unsigned long promoted = 0;
  	// build vector of WeakKmerMap
  	// screen out singletons 
  	// getCount() >= minKmerCount && getWeightedCount() >= minKmerWeightedCount

    typedef std::vector< WeakElementType > HeapType;
    HeapType weakHeap;
    
    WeakIterator mapIt  = weak.begin();
    WeakIterator mapEnd = weak.end();
    for(; mapIt != mapEnd; mapIt++ ) {
      WeakDataType &data = mapIt->value();
      if (data.getCount() >= minKmerCount && data.getWeightedCount() >= minKmerWeightedCount) {
        weakHeap.push_back( *mapIt );
      }
    }
    
    if (weakHeap.size() == 0) {
    	std::cerr << "There are no eligible kmers to promote" << std::endl;
    	return 0;
    }
    std::cerr << "Heap: " <<  weakHeap.size() << " " << weakHeap[0].value() << std::endl;
    //for(int i=0; i<weakHeap.size(); i++)
    //  std::cerr << i << ": " << weakHeap[i].value() << std::endl;
        
    //heapify kmer counts 
    std::make_heap( weakHeap.begin(), weakHeap.end() );
    //std::cerr << "Heap: " <<  weakHeap.size() << " " << weakHeap.front().value() << std::endl;
    //for(int i=0; i<weakHeap.size(); i++)
    //  std::cerr << i << ": " << weakHeap[i].value() << std::endl;
 
    // work in batches, stop when 0 new solids are added
    unsigned long count = 0;
    unsigned long promotedInBatch = 0;

    while ( weakHeap.begin() != weakHeap.end() ) {
    	WeakElementType &element = weakHeap.front();
    	unsigned long lastCount = element.value().getCount();
        if (shouldBeSolid( element, minWeakRatio, minSolidRatio )) {
        	SolidElementType solidElement = getSolid( element.key() );
        	solidElement.value() = element.value();
        	
        //	std::cerr << "Added " << prettyW(element.key(), element.value());
        	element.value().reset(); // to avoid double counting
        	promoted++;
        	promotedInBatch++;
        } else {
       // 	std::cerr << "Skipping " << prettyW(element.key(), element.value());
        }
        if (++count == 1000) {
          std::cerr << "Added " << promotedInBatch << ", Heap size: " <<  weakHeap.size() << " lastCount: " << lastCount << std::endl;
        	if (promotedInBatch == 0)
        	   break;
        	promotedInBatch = count = 0;
        	getErrorRates( solid );
        }
    	std::pop_heap(weakHeap.begin(), weakHeap.end());
    	weakHeap.pop_back();
    }

    for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
         weak.remove( it->key() );
    return promoted;
  }
  
  bool shouldBeSolid( WeakElementType &weakElement, double minWeakRatio, double minSolidRatio ) {
  	Kmer &kmer = weakElement.key();
  	WeakDataType &data = weakElement.value();
  	SolidWeakWeightType permutedScores = getPermutedScores(kmer, 1.0, 0.0);
  	if (   permutedScores.first  <= minSolidRatio  * data.getCount()
  	    && permutedScores.second <= minWeakRatio * data.getCount()
        )
  	    return true;
    else
        return false;  	
  }
  
    //repeat:
    // pop one
    // if solid (sum of all 1st order permutations < thresholdPermutationRatio * mycount)
    //    add to solid
    // else
    //    if rolling avg of tested kmer to added solid < thresholdSolidIncorporation 
    //         stop
    
    
    
    
    
    // find peak-getCount() (from maximum backwards) (above trivial at 1)
    // seeds Solid kmers
    
    // find all 1st order permutations of solids, knock out weaks
    // walk down sorted list of kmers
    // test 1st order permutation:
    //    if permute < threshold -> promote, knock out weaks
    //    else next
    
    // then check permutations of kmers with ratios > thresholdMax
    
    
    
  
  
  unsigned long promote(double probabilityQuantile) {
  	unsigned long promoted = 0;
    
  	typedef accumulator_set<
  	                unsigned long, 
  	                stats< tag::weighted_p_square_quantile         
  	                       >,
  	                double
  	                > AccumulatorType;
  	int count = 25;
  	double probabilities[] = { probabilityQuantile, 1.0-probabilityQuantile, 
  		                       0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,
  		                       0.06,  0.07,  0.08, 0.09,  0.10, 0.11,  0.12, 0.13,  0.14, 0.15,  0.16, 0.20 };
  	std::vector< AccumulatorType > countAcc, weightedCountAcc, directionAcc;
  	for(int i=0; i < count; i++) {
  	   countAcc.push_back( AccumulatorType( quantile_probability = probabilities[i] ) );
  	   weightedCountAcc.push_back( AccumulatorType( quantile_probability = probabilities[i] ) );
  	   directionAcc.push_back( AccumulatorType( quantile_probability = probabilities[i] ) );
  	}
  	for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++)
  	   for(int i=0; i< count; i++) {
  	   	  WeakDataType &data = it->value();
  	   	  double count = data.getCount();
  	      countAcc[i]( count, weight = count );
  	      weightedCountAcc[i]( data.getWeightedCount(), weight = count );
  	      directionAcc[i]( data.getNormalizedDirectionBias(), weight = count);
  	   }
  	for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
  	   for(int i=0; i< count; i++) {
  	   	  SolidDataType &data = it->value();
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
  	
    for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
      WeakDataType &data = it->value();
      double dirBias = data.getNormalizedDirectionBias();
      if (data.getCount() > minimumCount && data.getWeightedCount() > minimumWeightedCount
          //&& dirBias < maximumDirBias && dirBias > minimumDirBias
          //&& dirBias < 0.8 && dirBias > 0.2
          ) {
        SolidWeakWeightType scores = getPermutedScores( it->key(), Options::getFirstOrderWeight(), Options::getSecondOrderWeight() );
        double permuteScore = (scores.first+scores.second);
        if ( (double)data.getCount() / permuteScore > 1.0 ) {
      	  getSolid( it->key() ) = it->value();
      	  promoted++;
        }
      }
    }
  	for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
  	  weak.remove( it->key() );

  	return promoted;
  }
  
  MeanVectorType getErrorRates( SolidMapType &solidReference, bool useWeighted = false ) {
  	
  	std::vector< StdAccumulatorType > accumulators;
    for( SolidIterator it = solidReference.begin(); it != solidReference.end(); it++) {
      std::vector< double > errorRatios = getErrorRatios(it->key(), useWeighted);
      if (accumulators.size() < errorRatios.size())
        accumulators.resize( errorRatios.size() );
      for(unsigned int i = 0 ; i< errorRatios.size(); i++)
        accumulators[i](errorRatios[i]);
      //std::cerr << "Error Ratio for " << it->key().toFasta() << " "  << std::fixed << std::setprecision(6) << errorRatio << std::endl;
    }
    
    MeanVectorType rates;
    for(unsigned int i = 0 ; i< accumulators.size(); i++) {
    	rates.push_back( MeanType(mean(accumulators[i]), sqrt( variance(accumulators[i]) ) ) );
        std::cerr << "Error Rate for pos " << i << ", " << std::fixed << std::setprecision(6) << rates[i].first << " +/- " << std::fixed << std::setprecision(6) << rates[i].second << std::endl;
    }
    return rates;
  }
  
  std::vector< double > getErrorRatios(const Kmer &kmer, bool useWeighted = false) {
  	std::vector< double > errorRatios;
  	double baseValue;
  	DataPointers pointers( *this, kmer );

  	if ( pointers.solid.isValid() ) 
  	  baseValue = useWeighted ? pointers.solid.value().getWeightedCount() : pointers.solid.value().getCount();
  	else if ( pointers.weak.isValid() )
  	  baseValue = useWeighted ? pointers.weak.value().getWeightedCount() : pointers.weak.value().getCount();
    else if ( pointers.singleton.isValid() )
      baseValue = useWeighted ? pointers.singleton.value().getWeightedCount() : pointers.singleton.value().getCount();
    else
  	  return errorRatios;
  	  
  	Kmers permutations = Kmers::permuteBases(kmer,true);
  	std::vector<double> previouslyObservedSolids;
  	std::vector< std::vector<double> > weakCounts;
  	for(unsigned int i=0; i<permutations.size(); i++) {
  		pointers.set( permutations[i] );
  		if ( pointers.solid.isValid() ) {
  		  previouslyObservedSolids.push_back( useWeighted ? pointers.solid.value().getWeightedCount() : pointers.solid.value().getCount() );
  		} else if ( pointers.weak.isValid() ) {
  		  WeakDataType &data = pointers.weak.value();
  		  double count = useWeighted ? data.getWeightedCount() : data.getCount();
  		  WeakReadPositionWeightVector positions = data.getEachInstance();
  		  std::vector< double > positionSums;
  		  for(WeakReadPositionWeightVectorIterator it=positions.begin(); it!=positions.end(); it++) {
  		  	unsigned short pos = it->position;
  		  	if (positionSums.size() < pos + 1ul)
  		  	  positionSums.resize(pos + 1ul);
  		  	positionSums[pos] += useWeighted ? it->weight : 1.0;
  		  }
  		  for(unsigned int i=0; i< positionSums.size(); i++) {
  		  	if(positionSums[i] > 0) {
  		  		if (weakCounts.size() < i+1ul)
  		  		  weakCounts.resize(i+1ul);
  		  		weakCounts[i].push_back( count * positionSums[i] / positions.size() );
  		  	}
  		  }
  		} else if ( pointers.singleton.isValid() ) {
  		  TrackingDataSingleton &data = pointers.singleton.value();
  		  if (weakCounts.size() < data.getPosition()+1ul)
  		    weakCounts.resize(data.getPosition()+1);
  		  weakCounts[data.getPosition()].push_back( useWeighted ? pointers.singleton.value().getWeightedCount() : pointers.singleton.value().getCount() );
  		}
  	}
  	std::vector< double > median, nonZeroCount;
  	if (weakCounts.size() > 0) {
  	  median.resize( weakCounts.size() );
  	  nonZeroCount.resize( weakCounts.size() );
  	  errorRatios.resize( weakCounts.size() );
  	  for(unsigned int pos = 0 ; pos < weakCounts.size(); pos++) {
  	  	if (weakCounts[pos].size() > 0) {
  	  	  std::vector< double > &tmp = weakCounts[pos];
  		  sort(tmp.begin(), tmp.end());
  		  median[pos] = tmp[ tmp.size()/2 ];
  		  nonZeroCount[pos] = tmp.size();
  	  	  errorRatios[pos] = (median[pos] * nonZeroCount[pos]) / baseValue;
  	  	}
  	  }
  	  if (previouslyObservedSolids.size() == 0) 
  	    return errorRatios;  		
  		
      // some code to correct for previously Observed Solids
      for(unsigned int pos = 0 ; pos < weakCounts.size(); pos++) {
      	  double rawErrorRatio = 0.0;
  	      double sum = baseValue;
  	      for(std::vector<double>::iterator it = previouslyObservedSolids.begin(); it!=previouslyObservedSolids.end(); it++) {
  	         sum += *it;
  	         rawErrorRatio -= (median[pos] * nonZeroCount[pos]) / *it;
  	      }
  		
  	      rawErrorRatio += (median[pos] * nonZeroCount[pos]) / sum;
  		  
  	      errorRatios[pos] = rawErrorRatio;
      }
      return errorRatios;
  	} else {
  		return errorRatios;
  	}  	
  }
  
  SolidWeakWeightType getPermutedScores( const Kmer &kmer, double permutationWeight = 0.1, double secondWeight = 0.0 ) {
  	SolidWeakWeightType score(0.0,0.0);
  	if (permutationWeight == 0.0)
  	  return score;
  	DataPointers pointers(*this, kmer);

  	//std::cerr << "Permuting " << kmer->toFasta() << std::endl;
  	bool isSolid = false;
  	double base = 0.0;
  	if ( pointers.solid.isValid() ) {
  	  base = pointers.solid.value().getCount();
  	  isSolid = true;
  	} else if ( pointers.weak.isValid() ) {
  	  base = pointers.weak.value().getCount();
  	} else if ( pointers.singleton.isValid() ) {
  	  base = pointers.singleton.value().getCount();
  	}
  	KmerWeights permutations = KmerWeights::permuteBases(kmer, true);
  	for(unsigned int i=0; i<permutations.size(); i++) {
  	  pointers.set( permutations[i] );
      //std::cerr << "Looking at " << permutations[i].toFasta() << std::endl;		
  	  if( pointers.solid.isValid() ) {
  	  	score.first += pointers.solid.value().getCount() * permutationWeight;
  	  } else if ( pointers.weak.isValid() ) {
  	  	score.second += pointers.weak.value().getCount() * permutationWeight;
  	  } else if ( pointers.singleton.isValid() ) {
  	  	score.second += pointers.singleton.value().getCount() * permutationWeight;
  	  }
  	  if (secondWeight > 0.0) {
  	    SolidWeakWeightType sScores = getPermutedScores( permutations[i], secondWeight );
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
     unsigned int count=0;
     for(unsigned int i=0; i<fasta.length(); i++)
       if (fasta[i] == 'G' || fasta[i] == 'C')
         count++;
     return count;
  }
  std::string pretty( const Kmer &kmer, std::string note ) {
  	std::stringstream ss;
  	
  	TEMP_KMER(rev);
  	kmer.buildReverseComplement(  rev );
  	
  	ss << kmer.toFasta() << " " << rev.toFasta() << " " << countGC(kmer.toFasta()) << "\t" << note << std::endl;
  	return ss.str();
  }
  std::string prettyW( const Kmer &kmer, WeakDataType &data ) {
  	std::stringstream ss;
  	ss << data << "\t";

  	SolidWeakWeightType scores = getPermutedScores( kmer, 1.0, 0.0 );
  	ss << std::fixed << std::setprecision(2) << scores.first  << "\t";
  	ss << std::fixed << std::setprecision(2) << scores.second << "\t";
  	ss << std::fixed << std::setprecision(2) << (double)data.getCount() / (double)(scores.first+scores.second) << "\t";
  	
  	return pretty( kmer, ss.str() );
  }
  std::string prettyS( const Kmer &kmer, SolidDataType &data ) {
  	std::stringstream ss;
  	ss << data << "\t";

  	SolidWeakWeightType scores = getPermutedScores( kmer, 1.0, 0.0 );
  	ss << std::fixed << std::setprecision(2) << scores.first  << "\t";
  	ss << std::fixed << std::setprecision(2) << scores.second << "\t";
  	ss << std::fixed << std::setprecision(2) << (double)data.getCount() / (double)(scores.first+scores.second) << "\t";
  	
  	return pretty( kmer, ss.str() );
  }
  std::string contrastSpectrums(KmerSpectrum &reference, double minSolidDepth = -1.0) {
  	// ignores reference.weak and everything in reference.solid passed criteria
  	 
  	unsigned long thisSolidSize(0), thisWeakSize(0), referenceSolidSize(0), 
  	              trueSolid(0), missingSolid(0), falseWeak(0), falseSolid(0), trueWeak(0);
  	
  	thisWeakSize = weak.size();
  	// compare this->solid to reference.solid 
  	for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++) {
  	  if (minSolidDepth >= 0 && it->value().getCount() < minSolidDepth) {
  	  	thisWeakSize++;
  	    if ( reference.solid.exists( it->key() ) ) {
  	      // exists in reference but in this is weak: falseWeak
  	      falseWeak++;
  	      std::cerr << "FW " << prettyS( it->key(), it->value());
  	    } else {
  	      // Weak and not in reference: trueWeak
  	      // count later
  	      if (Options::getVerbosity() > 1)
  	        std::cerr << "TW " << prettyS( it->key(), it->value());
  	    }
  	  } else {
  	    thisSolidSize++;
  	    if ( reference.solid.exists( it->key() ) ) {
  	      // correctly matches solid in reference: trueSolid
  	  	  trueSolid++;
  	  	  if (Options::getVerbosity() > 0)
  	  	    std::cerr << "TS " << prettyS( it->key(), it->value());
  	    } else {
  	      // exists in this but not reference: falseSolid
  	      falseSolid++;
  	      std::cerr << "FS " << prettyS( it->key(), it->value());
  	    }
  	  }
  	}
  	
  	// compare reference.solid to this->solid and this->weak
  	for(SolidIterator it(reference.solid.begin()), itEnd(reference.solid.end()); it != itEnd; it++) {
  	  referenceSolidSize++;
  	  
  	  if ( solid.exists( it->key() ) ) {
  	  	// already accounted  	  	
  	  } else if ( weak.exists ( it->key() )) {
      	    // exists in reference but in this is weak: falseWeak
      	    falseWeak++;
      	    std::cerr << "FW " << prettyW( it->key(), weak[ it->key() ] );
  	  } else {
  	  	// exists in reference but not in either solid nor weak: missingSolid
  	    std::cerr << "MS " << pretty( it->key(), "N/A");
  	    missingSolid++;
  	  }
  	}
  	if (Options::getVerbosity() > 2) {
  	  for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
  	    if (! reference.solid.exists( it->key() ))
  	      std::cerr << "TW " << prettyW( it->key(), it->value());
  	  }
  	}
  	trueWeak = thisWeakSize - falseWeak;

    std::stringstream header;
  	header << "SolidSizes: " << thisSolidSize << " / " << referenceSolidSize << std::endl;
  	header << "trueSolid: " << trueSolid << " MissingSolid: " << missingSolid << " falseWeak: " << falseWeak << " falseSolid: " << falseSolid << std::endl;
  	header << "trueWeak: " << trueWeak << " / " << thisWeakSize << std::endl;
   	return header.str();
  }
  
  inline void append( KmerWeights &kmers, unsigned long readIdx, bool isSolid = false) {
  	return append(kmers, readIdx, isSolid, 0, kmers.size());
  }
  
  void append( KmerWeights &kmers, unsigned long readIdx, bool isSolid, unsigned long kmerIdx, unsigned long kmerLen ) {
      
     TEMP_KMER (least);
     
     DataPointers pointers(*this);
     
     //unsigned long j=0;
     //for (KmerWeights::Iterator it(kmers.begin()), itEnd(kmers.end()); it != itEnd ; it++,j++)
     for(unsigned long j=0; j < kmerLen; j++)
     {
     	bool keepDirection = kmers[kmerIdx+j].buildLeastComplement( least );
       	double weight = kmers.valueAt(kmerIdx+j);
       	
       	if ( isSolid ) {
       	  pointers.reset();
       	  SolidElementType elem = getSolid( least );
       	  elem.value().track( weight, keepDirection, readIdx, j );
       	} else if (! TrackingData::isDiscard( weight ) ) {
       	  pointers.set( least );
       	  
       	  if (pointers.solid.isValid()) {
       	  	// track solid stats
       	  	pointers.solid.value().track( weight, keepDirection, readIdx, j );
       	  } else if (pointers.weak.isValid()) {
       	  	// track weak stats       	  	
       	  	pointers.weak.value().track( weight, keepDirection, readIdx, j );
       	  } else {
       	  	  if ( pointers.singleton.isValid() ) {
       	  	  	
       	  	  	
       	  	  	// promote singleton to weak & track
       	  	  	SingletonDataType singleData = pointers.singleton.value();
       	  	  	pointers.singleton.value().reset();
                pointers.reset();
                //singleton.remove( least );
			    
                WeakElementType weakElem = getWeak( least );
       	  	  	weakElem.value() = singleData;
       	  	    weakElem.value().track( weight, keepDirection, readIdx, j);
			    
			    
       	  	  } else {
       	  	  	
       	  	  	// record this singleton
       	  	  	SingletonElementType elem = getSingleton(least);
       	  	  	elem.value().track( weight, keepDirection, readIdx, j );
			    
       	  	  }
       	  }
       	}
     }
  }
  

  void printStats(unsigned long pos, bool solidOnly = false) {
	//stats.printHistograms(solidOnly);
	SolidIterator it = solid.begin();
    std::cerr << pos << " reads, " << solid.size() << " solid / " << weak.size() << " weak / " << TrackingData::discarded << " discarded / " << TrackingData::singletonCount << " : " << singleton.size() << " singleton - kmers so far ";
    if(it != solid.end())
        std::cerr << it.bucket().toString() << "; ";
    std::cerr << " minDepth: " << (unsigned long)TrackingData::minimumDepth <<  std::endl;
    std::cerr << MemoryUtils::getMemoryUsage() << std::endl;     
  }

  inline unsigned int getSMPThread( Kmer &kmer, unsigned int numThreads ) {
  	// return the smallest KmerMap in the spectrum's getLocalThreadId()
  	TEMP_KMER(tmp)
  	kmer.buildLeastComplement(tmp);
  	return solid.getLocalThreadId(tmp, numThreads);
  }
  inline unsigned int getDMPThread( Kmer &kmer, unsigned int numThreads ) {
  	TEMP_KMER(tmp)
  	kmer.buildLeastComplement(tmp);
  	return solid.getDistributedThreadId(tmp, numThreads);
  }
  
  void buildKmerSpectrum( ReadSet &store, bool isSolid = false ) 
  {
	
	weak.reset();
	solid.reset();
	singleton.reset();
	
    #ifdef _USE_OPENMP
	
	  // allocate a square matrix of buffers: KmerWeights[writingThread][readingThread]
	  KmerWeights kmerBuffers[omp_get_num_procs()][omp_get_num_procs()];
	  KmerWeights kmers;
	  typedef std::pair< long, long > ReadPosType;
	  std::vector< ReadPosType > startReadIdx[ omp_get_num_procs() ][ omp_get_num_procs() ];

      int numProcs = omp_get_num_procs();

	  long batchIdx = 0;
	  long batch = 1000000;
	  while (batchIdx < (long) store.getSize()) 
	  {
	  	printStats(batchIdx, isSolid); 
	    #pragma omp parallel
	    {	    	
	    	for(int i = 0; i < numProcs; i++) {
	    	   kmerBuffers [ omp_get_thread_num() ][ i ].reset(false);
	    	   startReadIdx[ omp_get_thread_num() ][ i ].resize(0);
	    	   startReadIdx[ omp_get_thread_num() ][ i ].reserve(batch);
	    	}	    
	    }
	    //std::cerr << "Reallocated buffers" << std::endl;
	    
	    #pragma omp parallel for private(kmers)
	    for (long i=0 ; i < batch; i++)
	    {
	       long readIdx = batchIdx + i;
	       if (readIdx >= store.getSize() )
	         continue;
	  	   kmers = buildWeightedKmers(store.getRead( readIdx ));
	  	   for (long j = 0; j < omp_get_num_threads(); j++)
	  	     startReadIdx[ omp_get_thread_num() ][ j ].push_back( ReadPosType(readIdx, kmerBuffers[ omp_get_thread_num() ][j].size()) );
	  	   for (IndexType j = 0 ; j < kmers.size(); j++) {
	  	 	 kmerBuffers[ omp_get_thread_num() ][ getSMPThread( kmers[j], omp_get_num_threads() ) ].append(kmers[j], kmers.valueAt(j));
	  	   }
	    }
	    //std::cerr << "Loaded buffers" << std::endl;
	    //for(long i=0; i<numProcs; i++)
	    //  for (long j=0; j<numProcs; j++)
	    //    std::cerr << i << "," << j << "\t" << startReadIdx[i][j].size() << "\t" << kmerBuffers[i][j].size() << std::endl;
	        
	    #pragma omp parallel
	    {
	    	//#pragma omp critical
	    	//{ std::cerr << "running with " << omp_get_num_threads() << " this is " << omp_get_thread_num() << std::endl; }
	    	for(int threads = 0 ; threads < omp_get_num_threads(); threads++) 
	    	{
	    		for(unsigned long idx = 0; idx < startReadIdx[ threads ][ omp_get_thread_num() ].size(); idx++)
	    	    {
	    	    	ReadPosType readPos = startReadIdx[ threads ][ omp_get_thread_num() ][idx];
	    	    	unsigned long startIdx = readPos.second;
	    	    	unsigned long len = 0;
	    	    	if ( idx < startReadIdx[ threads ][ omp_get_thread_num() ].size() -1) {
	    	    		len = startReadIdx[ threads ][ omp_get_thread_num() ][idx+1].second - startIdx;
	    	    	} else {
	    	    		len = kmerBuffers[threads][ omp_get_thread_num() ].size() - startIdx;
	    	    	}
	    	    	if (len > 0)
	    	    	  append(kmerBuffers[threads][ omp_get_thread_num() ], readPos.first, isSolid, startIdx, len);
	    	    }
	    	}
	    }
	    //std::cerr << "appended buffers" << std::endl;
	    batchIdx += batch; 
	  }

    #else // not using OpenMP
	
	    for (long i=0 ; i < (long) store.getSize(); i++)
        {
         KmerWeights kmers = buildWeightedKmers(store.getRead(i));
       
         append(kmers, i, isSolid);

        if (i % 1000000 == 0) {
       	  printStats(i, isSolid);  
        }
      }
      
    #endif
	
    printStats(store.getSize(), isSolid);    
    if (!isSolid)
      printHistograms();
  }

  static void experimentOnSpectrum( KmerSpectrum &spectrum ) {
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
    
    for(WeakIterator it(everything.weak.begin()), itEnd(everything.weak.end()); it != itEnd; it++) {
    	everything.getSolid( it->key() ) = it->value();
    }
    everything.weak.clear();
    
    for (double minDepth = 1; minDepth < 10 ; minDepth++) {
    	std::cerr << "Comparing minDepth " << minDepth << std::endl;
    	std::cerr << everything.contrastSpectrums(spectrum, minDepth) << std::endl;
    }
  }

};

#endif

