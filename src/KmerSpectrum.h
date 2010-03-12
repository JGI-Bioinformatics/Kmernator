// $Header: /repository/PI_annex/robsandbox/KoMer/src/KmerSpectrum.h,v 1.26 2010-03-12 19:08:42 regan Exp $

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

using namespace boost::accumulators;

#include "config.h"
#include "TwoBitSequence.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerReadUtils.h"
#include "Options.h"

template<typename So, typename We, typename Si = TrackingDataSingleton>
class KmerSpectrum {
public:

	typedef std::pair<double, double> DoublePairType;
	typedef DoublePairType SolidWeakWeightType;
	typedef DoublePairType MeanType;

	typedef std::vector<DoublePairType> DoublePairVectorType;
	typedef DoublePairVectorType MeanVectorType;

	typedef So SolidDataType;
	typedef We WeakDataType;
	typedef Si SingletonDataType;

    typedef	typename WeakDataType::ReadPositionWeightVector WeakReadPositionWeightVector;
	typedef typename WeakDataType::ReadPositionWeightVector::iterator WeakReadPositionWeightVectorIterator;

	typedef KmerMap<SolidDataType> SolidMapType;
	typedef KmerMap<WeakDataType> WeakMapType;
	typedef KmerMap<SingletonDataType> SingletonMapType;

	typedef typename SolidMapType::Iterator SolidIterator;
	typedef typename SolidMapType::ElementType SolidElementType;
	typedef typename WeakMapType::Iterator WeakIterator;
	typedef typename WeakMapType::ElementType WeakElementType;
	typedef typename SingletonMapType::Iterator SingletonIterator;
	typedef typename SingletonMapType::ElementType SingletonElementType;

	typedef accumulator_set<double,
	stats< tag::variance(lazy),
	tag::mean
	>
	> StdAccumulatorType;

public:
	SolidMapType solid;
	WeakMapType weak;
	SingletonMapType singleton;
	bool hasSolids;
	unsigned long purgedSingletons;

public:
	KmerSpectrum(unsigned long buckets): solid(buckets/64), weak(buckets/8), singleton(buckets), hasSolids(false), purgedSingletons(0)
	{
		// set the minimum weight that will be used to track kmers
		// based on the given options
		TrackingData::minimumWeight = Options::getMinKmerQuality();
		TrackingData::minimumDepth = Options::getMinDepth();
		// apply the minimum quality automatically
		Read::setMinQualityScore( Options::getMinQuality() );
	}
	~KmerSpectrum() {}
	KmerSpectrum(const KmerSpectrum &copy) {
		*this = copy;
	}
	KmerSpectrum &operator=(const KmerSpectrum &other) {
		this->weak = other.weak;
		this->solid = other.solid;
		this->singleton = other.singleton;
		this->hasSolids = other.hasSolids;
		this->purgedSingletons = other.purgedSingletons;
	}

	static unsigned long estimateWeakKmerBucketSize( ReadSet &store, unsigned long targetKmersPerBucket = 64) {
		unsigned long baseCount = store.getBaseCount();
		if (baseCount == 0 || store.getSize() == 0)
		return 128;
		unsigned long avgSequenceLength = baseCount / store.getSize();
		unsigned long kmersPerRead = (avgSequenceLength - KmerSizer::getSequenceLength() + 1);

		unsigned long rawKmers = kmersPerRead * store.getSize();
		// assume 1% error rate
		unsigned long estimatedUniqueKmers = rawKmers * (unsigned long) ( std::max(pow(1.01, KmerSizer::getSequenceLength()),2.0) - 1.0 );
		unsigned long targetBuckets = estimatedUniqueKmers / targetKmersPerBucket;
		unsigned long maxBuckets = 128*1024*1024;
		unsigned long minBuckets = 128;
		return std::max( std::min(targetBuckets,maxBuckets), minBuckets );
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
	SingletonElementType getIfExistsSingleton( const Kmer &kmer ) {
		return singleton.getElementIfExists( kmer );
	}
	SingletonElementType getSingleton( const Kmer &kmer ) {
		return singleton.getElement( kmer );
	}

	class DataPointers
	{
	public:
		KmerSpectrum *spectrum;
		SolidElementType solidElem;
		WeakElementType weakElem;
		SingletonElementType singletonElem;
		DataPointers(KmerSpectrum &_spectrum) : spectrum(&_spectrum) {}
		DataPointers(KmerSpectrum &_spectrum, const Kmer &kmer): spectrum(&_spectrum) {
			set(kmer);
		}
		void reset() {
			solidElem = SolidElementType();
			weakElem = WeakElementType();
			singletonElem = SingletonElementType();
		}
		void set(const Kmer &kmer) {
			reset();
			if (spectrum->hasSolids)
			solidElem = spectrum->getIfExistsSolid(kmer);

			if (!solidElem.isValid()) {
				weakElem = spectrum->getIfExistsWeak(kmer);
				if (!weakElem.isValid())
				singletonElem = spectrum->getIfExistsSingleton(kmer);
			}
		}
	};

	class Histogram {
	public:

		class HistogramElement {
		public:
			unsigned long visits;
			unsigned long visitedCount;
			double visitedWeight;
			HistogramElement() : visits(0), visitedCount(0), visitedWeight(0.0) {}
		};

		typedef std::vector< HistogramElement > BucketsType;

	private:
		double logBase;
		double logFactor;
		unsigned int zoomMax;
		unsigned int zoomLogSkip;
		unsigned int maxIdx;
		unsigned int lastBucket;
		BucketsType buckets;
		unsigned long count;
		double totalCount, totalWeightedCount;

		unsigned int getIdx(unsigned int count) {
			return count <= zoomMax ? count : (unsigned int) ( log((double)count)/logFactor - zoomLogSkip + zoomMax );
		}
		unsigned int getBucketValue(unsigned int idx) {
			return idx <= zoomMax ? idx : (unsigned int) ( pow(logBase, (double) (idx + zoomLogSkip - zoomMax)) );
		}
	public:
		Histogram(unsigned int _zoomMax, double _logBase = 2.0) :
		logBase(_logBase), logFactor(log(logBase)),
		zoomMax(_zoomMax), zoomLogSkip(0), maxIdx(0), lastBucket(0), buckets(),
		count(0), totalCount(0.0), totalWeightedCount(0.0) {
			zoomLogSkip = (unsigned int) ( log((double) zoomMax+1.0)/logFactor - 1.0 );
			int maxLog2 = 16;
			maxIdx = (1<<maxLog2) + 1 + zoomMax;
			buckets.resize(maxIdx+1);
		}
		inline unsigned long getCount() {return count;}
		inline double getTotalCount() {return totalCount;}
		inline double getTotalWeight() {return totalWeightedCount;}
		inline BucketsType getBuckets() {return buckets;}

		void addRecord(unsigned long count, double weight) {
			if (count == 0)
				return;
			unsigned int idx = getIdx(count);
			HistogramElement &elem = buckets[ idx ];
			elem.visits++;
			elem.visitedCount += count;
			elem.visitedWeight += weight;
			if (count < weight)
				throw std::invalid_argument("count < weight!");
		}
		void resetTotals() {
			count = 0;
			lastBucket = 0;
			totalCount = totalWeightedCount = 0.0;
		}
		void finish() {
			resetTotals();
			for(unsigned int i = 0; i<buckets.size(); i++) {
				HistogramElement &elem = buckets[ i ];
				if (elem.visits > 0) {
					count += elem.visits;
					totalCount += elem.visitedCount;
					totalWeightedCount += elem.visitedWeight;
					lastBucket = i;
				}
			}
		}
		std::string toString() {
			std::stringstream ss;

			finish();
			ss << std::fixed << std::setprecision(3);
			ss << "Counts, Weights and Directions" << std::endl;
			ss << "Counts:\t" << count
			<< "\t" << totalCount
			<< "\t" << (totalCount/count)
			<< "\t" << std::endl;

			//std::cerr << mean(countAcc) << "/" << median(countAcc) << " +-";
			//std::cerr << sqrt(variance(countAcc)) << "\t";
			//std::cerr << mean(countWAcc) << " +-";
			//std::cerr << sqrt(variance(countWAcc)) << "\t";

			ss << "Weights:\t" << totalWeightedCount
			<< "\t" << (totalWeightedCount/count)
			<< "\t" << (totalWeightedCount/totalCount)
			<< std::endl;
			//std::cerr << mean(weightedCountAcc) << "/" << median(weightedCountAcc) << " +-";
			//std::cerr << sqrt(variance(weightedCountAcc)) << "\t";
			//std::cerr << weighted_mean(weightedCountWAcc) << " +-";
			//std::cerr << sqrt(weighted_variance(weightedCountWAcc)) << "\t";

			ss << std::endl;

			ss << "Bucket\tUnique\t%Unique\tCount\t%Count\tWeight\tQualProb\t%Weight" << std::endl;
			for (unsigned int i=1; i < lastBucket+1; i++) {
				ss << getBucketValue(i) << "\t";
				ss << buckets[i].visits << "\t";
				ss << 100.0 * buckets[i].visits / count << "\t";
				ss << buckets[i].visitedCount << "\t";
				ss << 100.0 * buckets[i].visitedCount / totalCount << "\t\t";
				ss << buckets[i].visitedWeight << "\t";
				ss << buckets[i].visitedWeight / buckets[i].visitedCount << "\t";
				ss << 100.0 * buckets[i].visitedWeight / totalWeightedCount << "\t";
				ss << std::endl;
			}
			return ss.str();
		}
	};

	void printHistograms(bool solidOnly = false) {

		Histogram histogram(63);

		if (!solidOnly) {
			setWeakHistogram(histogram);
			setSingletonHistogram(histogram);
		}
		setSolidHistogram(histogram);

		std::cerr << histogram.toString();
	}

	void setWeakHistogram(Histogram &histogram) {
		for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
			WeakDataType &data = it->value();
			histogram.addRecord( data.getCount(), data.getWeightedCount() );
		}
	}
	void setSolidHistogram(Histogram &histogram) {
		for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++) {
			SolidDataType &data = it->value();
			histogram.addRecord( data.getCount(), data.getWeightedCount() );
		}
	}
	void setSingletonHistogram(Histogram &histogram) {
		for(SingletonIterator it(singleton.begin()), itEnd(singleton.end()); it != itEnd; it++) {
			SingletonDataType &data = it->value();
			histogram.addRecord( data.getCount(), data.getWeightedCount() );
		}
		for(unsigned long i; i < purgedSingletons ; i++) {
			histogram.addRecord( 1, 1.0);
		}
	}

	void printHistogramsOld(bool solidOnly = false) {

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
		StdAccumulatorType countAcc( tag::density::num_bins = bins, tag::density::cache_size = cacheSize );
		StdAccumulatorType weightedCountAcc( tag::density::num_bins = bins, tag::density::cache_size = cacheSize );
		StdAccumulatorType directionAcc( tag::density::num_bins = bins, tag::density::cache_size = cacheSize );

		typedef
		accumulator_set<double,
		stats< tag::weighted_variance,
		tag::weighted_density
		>,
		double
		> WeightedAccumulatorType;
		WeightedAccumulatorType countWAcc( tag::weighted_density::num_bins = bins, tag::weighted_density::cache_size = cacheSize );
		WeightedAccumulatorType weightedCountWAcc( tag::weighted_density::num_bins = bins, tag::weighted_density::cache_size = cacheSize );
		WeightedAccumulatorType directionWAcc( tag::weighted_density::num_bins = bins, tag::weighted_density::cache_size = cacheSize );

		if (!solidOnly) {
			for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
				WeakDataType &data = it->value();
				countAcc( std::max(-1.0, (double)(data.getCount()) ) );
				weightedCountAcc( std::max(-1.0, (double)(data.getWeightedCount())) );
				directionAcc( data.getNormalizedDirectionBias() );

				countWAcc( std::max(-1.0, (double)(data.getCount()) ), weight = data.getCount() );
				weightedCountWAcc( std::max(-1.0, (double)(data.getWeightedCount())), weight = data.getCount() );
				directionWAcc( data.getNormalizedDirectionBias(), weight = data.getCount() );
			}
		}
		for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++) {
			SolidDataType &data = it->value();
			//std::cerr << data.getCount() << std::endl;
			countAcc( std::max(-1.0, (double)(data.getCount()) ) );
			weightedCountAcc( std::max(-1.0, (double)(data.getWeightedCount())) );
			directionAcc( data.getNormalizedDirectionBias() );

			countWAcc( std::max(-1.0, (double)(data.getCount()) ), weight = data.getCount() );
			weightedCountWAcc( std::max(-1.0, (double)(data.getWeightedCount())), weight = data.getCount() );
			directionWAcc( data.getNormalizedDirectionBias(), weight = data.getCount() );
		}

		typedef boost::iterator_range< std::vector< std::pair<double, double> >::iterator > HistogramType;
		HistogramType countHist = density(countAcc);
		HistogramType weightedCountHist = density(weightedCountAcc);
		HistogramType directionHist = density(directionAcc);

		HistogramType countHistW = weighted_density(countWAcc);
		HistogramType weightedCountHistW = weighted_density(weightedCountWAcc);
		HistogramType directionHistW = weighted_density(directionWAcc);

		std::cerr << std::fixed << std::setprecision(3);

		std::cerr << "Counts, Weights and Directions" << std::endl;

		std::cerr << "Counts:\t";
		std::cerr << mean(countAcc) << "/" << median(countAcc) << " +-";
		std::cerr << sqrt(variance(countAcc)) << "\t";
		std::cerr << mean(countWAcc) << " +-";
		std::cerr << sqrt(variance(countWAcc)) << "\t";

		std::cerr << "Weights:\t";
		std::cerr << mean(weightedCountAcc) << "/" << median(weightedCountAcc) << " +-";
		std::cerr << sqrt(variance(weightedCountAcc)) << "\t";
		std::cerr << weighted_mean(weightedCountWAcc) << " +-";
		std::cerr << sqrt(weighted_variance(weightedCountWAcc)) << "\t";

		std::cerr << std::endl;

		for (int i=0; i < countHist.size(); i++) {
			std::cerr << countHist[i].first << "\t";
			std::cerr << countHist[i].second*100.0<< "\t";
			//assert(countHist[i].first == countHistW[i].first);
			std::cerr << countHistW[i].second *100.0<< "\t|\t";

			std::cerr << weightedCountHist[i].first << "\t";
			std::cerr << weightedCountHist[i].second*100.0 << "\t";
			//assert(weightedCountHist[i].first == weightedCountHistW[i].first);
			std::cerr << weightedCountHistW[i].second*100.0 << "\t";

			std::cerr << directionHist[i].first << "\t";
			std::cerr << directionHist[i].second *100.0<< "\t";
			//assert(directionHist[i].first == directionHistW[i].first);
			std::cerr << directionHistW[i].second *100.0<< "\t";

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

		WeakIterator mapIt = weak.begin();
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
		std::cerr << "Heap: " << weakHeap.size() << " " << weakHeap[0].value() << std::endl;
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
				std::cerr << "Added " << promotedInBatch << ", Heap size: " << weakHeap.size() << " lastCount: " << lastCount << std::endl;
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
		if ( permutedScores.first <= minSolidRatio * data.getCount()
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
		double probabilities[] = {probabilityQuantile, 1.0-probabilityQuantile,
			0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05,
			0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.20};
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
			for(unsigned int i = 0; i< errorRatios.size(); i++)
			accumulators[i](errorRatios[i]);
			//std::cerr << "Error Ratio for " << it->key().toFasta() << " "  << std::fixed << std::setprecision(6) << errorRatio << std::endl;
		}

		MeanVectorType rates;
		for(unsigned int i = 0; i< accumulators.size(); i++) {
			rates.push_back( MeanType(mean(accumulators[i]), sqrt( variance(accumulators[i]) ) ) );
			std::cerr << "Error Rate for pos " << i << ", " << std::fixed << std::setprecision(6) << rates[i].first << " +/- " << std::fixed << std::setprecision(6) << rates[i].second << std::endl;
		}
		return rates;
	}

	std::vector< double > getErrorRatios(const Kmer &kmer, bool useWeighted = false) {
		std::vector< double > errorRatios;
		double baseValue;
		DataPointers pointers( *this, kmer );

		if ( pointers.solidElem.isValid() )
		baseValue = useWeighted ? pointers.solidElem.value().getWeightedCount() : pointers.solidElem.value().getCount();
		else if ( pointers.weakElem.isValid() )
		baseValue = useWeighted ? pointers.weakElem.value().getWeightedCount() : pointers.weakElem.value().getCount();
		else if ( pointers.singletonElem.isValid() )
		baseValue = useWeighted ? pointers.singletonElem.value().getWeightedCount() : pointers.singletonElem.value().getCount();
		else
		return errorRatios;

		Kmers permutations = Kmers::permuteBases(kmer,true);
		std::vector<double> previouslyObservedSolids;
		std::vector< std::vector<double> > weakCounts;
		for(unsigned int i=0; i<permutations.size(); i++) {
			pointers.set( permutations[i] );
			if ( pointers.solidElem.isValid() ) {
				previouslyObservedSolids.push_back( useWeighted ? pointers.solidElem.value().getWeightedCount() : pointers.solidElem.value().getCount() );
			} else if ( pointers.weakElem.isValid() ) {
				WeakDataType &data = pointers.weakElem.value();
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
			//} else if ( pointers.singletonElem.isValid() ) {
			//	TrackingDataSingleton &data = pointers.singletonElem.value();
			//	//if (weakCounts.size() < data.getPosition()+1ul)
			//	//weakCounts.resize(data.getPosition()+1);
			//	//weakCounts[data.getPosition()].push_back( useWeighted ? pointers.singletonElem.value().getWeightedCount() : pointers.singletonElem.value().getCount() );
			}
		}
		std::vector< double > median, nonZeroCount;
		if (weakCounts.size() > 0) {
			median.resize( weakCounts.size() );
			nonZeroCount.resize( weakCounts.size() );
			errorRatios.resize( weakCounts.size() );
			for(unsigned int pos = 0; pos < weakCounts.size(); pos++) {
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
			for(unsigned int pos = 0; pos < weakCounts.size(); pos++) {
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
		if ( pointers.solidElem.isValid() ) {
			base = pointers.solidElem.value().getCount();
			isSolid = true;
		} else if ( pointers.weakElem.isValid() ) {
			base = pointers.weakElem.value().getCount();
		} else if ( pointers.singletonElem.isValid() ) {
			base = pointers.singletonElem.value().getCount();
		}
		KmerWeights permutations = KmerWeights::permuteBases(kmer, true);
		for(unsigned int i=0; i<permutations.size(); i++) {
			pointers.set( permutations[i] );
			//std::cerr << "Looking at " << permutations[i].toFasta() << std::endl;
			if( pointers.solidElem.isValid() ) {
				score.first += pointers.solidElem.value().getCount() * permutationWeight;
			} else if ( pointers.weakElem.isValid() ) {
				score.second += pointers.weakElem.value().getCount() * permutationWeight;
			} else if ( pointers.singletonElem.isValid() ) {
				score.second += pointers.singletonElem.value().getCount() * permutationWeight;
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
		kmer.buildReverseComplement( rev );

		ss << kmer.toFasta() << " " << rev.toFasta() << " " << countGC(kmer.toFasta()) << "\t" << note << std::endl;
		return ss.str();
	}
	std::string prettyW( const Kmer &kmer, WeakDataType &data ) {
		std::stringstream ss;
		ss << data << "\t";

		SolidWeakWeightType scores = getPermutedScores( kmer, 1.0, 0.0 );
		ss << std::fixed << std::setprecision(2) << scores.first << "\t";
		ss << std::fixed << std::setprecision(2) << scores.second << "\t";
		ss << std::fixed << std::setprecision(2) << (double)data.getCount() / (double)(scores.first+scores.second) << "\t";

		return pretty( kmer, ss.str() );
	}
	std::string prettyS( const Kmer &kmer, SolidDataType &data ) {
		std::stringstream ss;
		ss << data << "\t";

		SolidWeakWeightType scores = getPermutedScores( kmer, 1.0, 0.0 );
		ss << std::fixed << std::setprecision(2) << scores.first << "\t";
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

				if (pointers.solidElem.isValid()) {
					// track solid stats
					pointers.solidElem.value().track( weight, keepDirection, readIdx, j );
				} else if (pointers.weakElem.isValid()) {
					// track weak stats
					pointers.weakElem.value().track( weight, keepDirection, readIdx, j );
				} else {
					if ( pointers.singletonElem.isValid() ) {

						// promote singleton to weak & track
						SingletonDataType singleData = pointers.singletonElem.value();
						pointers.reset();

						WeakElementType weakElem = getWeak( least );
						weakElem.value() = singleData;
						weakElem.value().track( weight, keepDirection, readIdx, j);

						singleton.remove( least );

					} else {
						// record this new singleton
						SingletonElementType elem = getSingleton(least);
						elem.value().track( weight, keepDirection, readIdx, j );
					}
				}
			} else if (Options::getDebug()){
				std::cerr << "discarded kmer " << readIdx << "@" << j << " " << weight << " " << least.toFasta() << std::endl;
			}
		}
	}

	void printStats(unsigned long pos, bool solidOnly = false, bool fullStats = false) {
		//stats.printHistograms(solidOnly);
		std::cerr << pos << " reads" << "\t";
		if (fullStats) {
			std::cerr << ", " << solid.size() << " solid / " << weak.size() << " weak / "
			<< TrackingData::discarded << " discarded / "
			<< singleton.size() + purgedSingletons << " singleton - kmers so far " << std::endl;
		}
		std::cerr << MemoryUtils::getMemoryUsage() << std::endl;
	}

	void printBucketStats(unsigned long pos, bool solidOnly = false, bool fullStats = false) {
		//stats.printHistograms(solidOnly);
		std::cerr << pos << ": ";

		unsigned long ss = solid.size();
		unsigned long buckets = solid.getNumBuckets();
		static unsigned long prev_ss = 0;
		std::cerr << ss << " added : " << ss - prev_ss << ",  average:"
		<< ss/buckets << " Max:" << solid.maxBucket() << " " << MemoryUtils::getMemoryUsage() << std::endl;
		prev_ss = ss;

	}

	inline unsigned int getSMPThread( Kmer &kmer, unsigned int numThreads ) {
		// return the smallest KmerMap in the spectrum's getLocalThreadId()
		TEMP_KMER(tmp);
		kmer.buildLeastComplement(tmp);
		return solid.getLocalThreadId(tmp, numThreads);
	}
	inline unsigned int getDMPThread( Kmer &kmer, unsigned int numThreads ) {
		TEMP_KMER(tmp);
		kmer.buildLeastComplement(tmp);
		return solid.getDistributedThreadId(tmp, numThreads);
	}

	unsigned long resetSingletons() {
		unsigned long singletonCount = singleton.size();
		purgedSingletons += singletonCount;
		std::cerr << "Purging Singletons: " << singletonCount << " (" << purgedSingletons << " total so far)" << std::endl;
		singleton.reset();
		return singletonCount;
	}
	void buildKmerSpectrum( ReadSet &store, bool isSolid = false )
	{

		weak.reset();
		solid.reset();
		singleton.reset();

		long purgeEvery = Options::getPeriodicSingletonPurge();
		long purgedSingletons = 0;
        long purgeCount = 0;
		long batchIdx = 0;
		long batch = 1000000;

#ifdef _USE_OPENMP

		int numThreads = omp_get_max_threads();

		// allocate a square matrix of buffers: KmerWeights[writingThread][readingThread]
		KmerWeights kmerBuffers[ numThreads ][ numThreads ];
		KmerWeights kmers;
		typedef std::pair< long, long > ReadPosType;
		std::vector< ReadPosType > startReadIdx[ numThreads ][ numThreads ];

		if (store.getSize() < batch)
		batch = store.getSize() / 10 + 1;

#pragma omp parallel num_threads(numThreads)
		{
			if (numThreads != omp_get_num_threads())
			throw "OMP thread count mis-match";

#pragma omp single
			{
				std::cerr << "Executing parallel buildKmerSpectrum with " << numThreads << std::endl;
			}
		}

		while (batchIdx < (long) store.getSize())
		{
			printStats(batchIdx, isSolid);
			for(int tmpThread1=0; tmpThread1 < numThreads; tmpThread1++)
			{
				for(int tmpThread2 = 0; tmpThread2 < numThreads; tmpThread2++) {
					kmerBuffers [ tmpThread1 ][ tmpThread2 ].reset(false);
					startReadIdx[ tmpThread1 ][ tmpThread2 ].resize(0);
					startReadIdx[ tmpThread1 ][ tmpThread2 ].reserve(batch);
				}
			}
			//std::cerr << "Reallocated buffers" << std::endl;

#pragma omp parallel for private(kmers) schedule(dynamic) num_threads(numThreads)
			for (long i=0; i < batch; i++)
			{
				long readIdx = batchIdx + i;
				if (readIdx >= store.getSize() )
				continue;
				kmers = KmerReadUtils::buildWeightedKmers(store.getRead( readIdx ));
				for (long j = 0; j < numThreads; j++)
				startReadIdx[ omp_get_thread_num() ][ j ].push_back( ReadPosType(readIdx, kmerBuffers[ omp_get_thread_num() ][j].size()) );
				for (IndexType j = 0; j < kmers.size(); j++) {
					kmerBuffers[ omp_get_thread_num() ][ getSMPThread( kmers[j], numThreads ) ].append(kmers[j], kmers.valueAt(j));
				}
			}
			//std::cerr << "Loaded buffers" << std::endl;
			//for(long i=0; i<numThreads; i++) {
			//   std::cerr << i ;
			//   for (long j=0; j<numThreads; j++) {
			//     std::cerr << "\t" << kmerBuffers[i][j].size();
			//   }
			//   std::cerr << std::endl;
			// }

#pragma omp parallel num_threads(numThreads)
			{
				//#pragma omp critical
				//{ std::cerr << "running with " << omp_get_num_threads() << " this is " << omp_get_thread_num() << std::endl; }
				if (numThreads != omp_get_num_threads())
				throw "OMP thread count mis-match";

				for(int threads = 0; threads < numThreads; threads++)
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
				//std::cerr << omp_get_thread_num() << ", ";
			}
			//std::cerr << " appended buffers" << std::endl;
			batchIdx += batch;
			if (purgeEvery > 0 && batchIdx >= (purgeCount+1)*purgeEvery) {
				purgedSingletons += resetSingletons();
				purgeCount++;
			}
		}

#else // not using OpenMP
		for (long i=0; i < (long) store.getSize(); i++)
		{
			KmerWeights kmers = KmerReadUtils::buildWeightedKmers(store.getRead(i));

			append(kmers, i, isSolid);

			if (i % batch == 0) {
				printStats(i, isSolid);
				batchIdx += batch;
				if (purgeEvery > 0 && batchIdx >= (purgeCount+1)*purgeEvery) {
					purgedSingletons += resetSingletons();
					purgeCount++;
				}
			}
		}

#endif

		printStats(store.getSize(), isSolid, true);
		if (!isSolid) {
			printHistograms();
		}
	}

	static void experimentOnSpectrum( KmerSpectrum &spectrum ) {
		KmerSpectrum everything(spectrum);

		unsigned long promoted = spectrum.promote(0.05);
		std::cerr << "Promoted " << promoted << " kmers" << std::endl;
		spectrum.printHistograms(true);

		for (double prob = 0.10; prob > 0.0; prob -= 0.005) {
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

		for (double minDepth = 1; minDepth < 10; minDepth++) {
			std::cerr << "Comparing minDepth " << minDepth << std::endl;
			std::cerr << everything.contrastSpectrums(spectrum, minDepth) << std::endl;
		}
	}

	void analyseSingletons()
	{
		int found = 0;
		for(WeakIterator it( weak.begin()), itEnd( weak.end()); it != itEnd; it++) {
			Kmers permutations = Kmers::permuteBases(it->key(), true);
			for(int i=0; i<permutations.size(); i++) {
				if (singleton.exists( permutations[i] )) {

					it->value().track(1.0, true);

					found++;
				}
			}
		}
		std::cerr << "Singleton permutations found in weak : " << found << std::endl;
	}

	void findSecondOrderPermutations(const Kmer &kmer)
	{
		Kmers permutations = Kmers::permuteBases(kmer, true);
		for(int i=0; i<permutations.size(); i++) {

			Kmers permutations2 = Kmers::permuteBases(permutations[i],true);
			for (int j = i; j<permutations2.size(); j++) {
				if (permutations2[j] != kmer)
				// if  (SolidDataType *s = solid.getIfExists(permutations2[j])  )
				{
					std::cerr << "   " << permutations2[j].toFasta() /*<< " : " << s->getCount() */<< std::endl;
				}
			}
		}
		std::cerr << "\n";
	}

	bool firstOrderMatch(const WeakElementType &element)
	{
		const int multiple = 5;
		int count = element.value().getCount();

		count *= multiple;
		const Kmer &kmer = element.key();
		Kmers permutations = Kmers::permuteBases(kmer, true);
		for(int i=0; i<permutations.size(); i++) {
			/*
			 SolidDataType *s = solid.getIfExists(permutations[i]);
			 if (s &&  s->getCount() > count)
			 return false;
			 */
			if (weak.exists(permutations[i])) {
				WeakElementType w = weak.getElement(permutations[i]);
				if (w.value().getCount() > count)
				return true;
			}
		}
		return false;
	}

	bool noSecondOrderMatch(const WeakElementType &element)
	{

		int count = element.value().getCount();

		count *= 10;
		const Kmer &kmer = element.key();
		Kmers permutations = Kmers::permuteBases(kmer, true);
		for(int i=0; i<permutations.size(); i++) {
			SolidDataType *s = solid.getIfExists(permutations[i]);
			Kmers permutations2 = Kmers::permuteBases(permutations[i],true);
			for (int j = i; j<permutations2.size(); j++) {
				if (permutations2[j] != kmer) {
					WeakDataType *w = weak.getIfExists(permutations2[j]);
					if (w && w->getCount() > count)
					return false;
				}
			}
		}
		return true;
	}
	/*
	 bool noSecondOrderMatch(const WeakElementType &element)
	 {

	 Kmers unique((32*3)*(32*3));

	 const Kmer &kmer = element.key();
	 Kmers permutations = Kmers::permuteBases(kmer, true);
	 for(int i=0; i<permutations.size(); i++) {
	 SolidDataType *s = solid.getIfExists(permutations[i]);
	 Kmers permutations2 = Kmers::permuteBases(permutations[i],true);
	 for (int j = i ; j<permutations2.size(); j++) {
	 unique.insertSorted(permutations2[j]);
	 }
	 }
	 int count = element.value().getCount();

	 for (int i =0; i < unique.size(); i++) {
	 if (unique[i] != kmer) {
	 SolidDataType *s = solid.getIfExists(unique[i]);
	 if (s &&  s->getCount() > count)
	 return false;
	 }
	 }
	 return true;
	 }*/

	void calcWeightStats(double &mean, double &deviation)
	{

		StdAccumulatorType weightAcc;

		for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
			WeakDataType &data = it->value();

			weightAcc(data.getAverageWeight() );

		}

		mean = boost::accumulators::mean(weightAcc);
		deviation = sqrt( boost::accumulators::variance(weightAcc));

	}

	unsigned long filter(KmerSpectrum &reference)
	{

		typedef std::vector< WeakElementType > HeapType;
		HeapType weakHeap;

		unsigned long promoted = 0;

		const unsigned long minKmerCount = 3;
		unsigned long count = 0;
		unsigned long falseSolids =0;

		unsigned long secondOrderTrue = 0;
		unsigned long secondOrderFalse = 0;
		unsigned long pq_count = 0;
		unsigned long pq_correct = 0;
		unsigned long pq_wrong = 0;

		unsigned long gq_count = 0;
		unsigned long gq_correct = 0;
		unsigned long gq_wrong = 0;

		unsigned long fo_count = 0;
		unsigned long fo_correct = 0;
		unsigned long fo_wrong = 0;

		double mean,deviation;
		calcWeightStats(mean,deviation);

		double poorQualityCutoff = 0.779;// mean - 3 * deviation;
		double goodQualityCutoff = mean + (1.00 - mean) * 0.9;

		std::cerr << mean << " +-" << deviation << " => cutoffs: " << poorQualityCutoff << " " << goodQualityCutoff << "\n";

		WeakIterator mapIt = weak.begin();
		WeakIterator mapEnd = weak.end();
		for(; mapIt != mapEnd; mapIt++ )
		{

			WeakElementType &element = *mapIt;
			unsigned long lastCount = element.value().getCount();

			if (lastCount < minKmerCount)
			continue;

			//         bool noBias = true;
			//         if (lastCount > 20)
			//         {
			//            double dirBias  = element.value().getNormalizedDirectionBias();
			//             noBias =  (dirBias > 0.02 && dirBias < 0.98);
			// //             if (!noBias) {
			// //
			// //
			// //                weakHeap.push_back( element);
			// //             }
			//         }


			count++;
			bool inRef = reference.solid.exists( element.key());

			bool poorQuality = element.value().getAverageWeight() < poorQualityCutoff;
			bool goodQuality = element.value().getAverageWeight() > goodQualityCutoff;

			/*
			 if (poorQuality)
			 {
			 pq_count++;
			 if (!inRef)
			 pq_correct++;
			 else
			 {
			 pq_wrong++;
			 //std::cerr << "FP "   << prettyS( element.key(), element.value());
			 }
			 }else*/
			if (firstOrderMatch(element))
			{
				fo_count++;
				if (!inRef)
				fo_correct++;
				else
				{
					fo_wrong++;
					//std::cerr << "FP "   << prettyS( element.key(), element.value());
				}

			}
			/*
			 if (goodQuality)
			 {
			 gq_count++;
			 if (inRef)
			 gq_correct++;
			 else
			 {
			 gq_wrong++;
			 //std::cerr << "FP "   << prettyS( element.key(), element.value());
			 }
			 }*/

			//         if ( //   noBias &&
			//         //      enoughWeight &&
			//           //    noFirstOrderMatch(element)&&
			//         //      noSecondOrderMatch(element)&&
			//            1)
			//         {
			//           SolidElementType solidElement = getSolid( element.key() );
			//             solidElement.value() = element.value();
			//            if (!inRef) {
			//                  // std::cerr << "FS " << prettyS( element.key(), element.value());
			//                  falseSolids++;
			//            }
			//            promoted++;
			//         }
			//


		}

		std::cerr << "Poor Quality: " << pq_count << " (" << pq_count*100.0/count << "%) Correct:  " << pq_correct << " Wrong: " << pq_wrong << std::endl;
		//    std::cerr << "Good Quality: " << gq_count <<  " ("  << gq_count*100.0/count <<  "%) Correct:  " << gq_correct << " Wrong: " << gq_wrong   << std::endl;
		std::cerr << "1st Order  : " << fo_count << " (" << fo_count*100.0/count << "%) Correct:  " << fo_correct << " Wrong: " << fo_wrong << std::endl;

		//     std::make_heap( weakHeap.begin(), weakHeap.end() );
		//
		//     std::cerr << " Heap size: " <<  weakHeap.size() << std::endl;
		//
		//     int i=0;
		//     while ( weakHeap.begin() != weakHeap.end() ) {
		//         WeakElementType &element = weakHeap.front();
		//         bool inRef = reference.solid.exists( element.key());
		//         const char * st =  (inRef) ? "*" : " ";
		//          std::cerr << "DB " << st  << prettyS( element.key(), element.value());
		//
		//         std::cout << "> Kmer" << ++i << std::endl;
		//         std::cout << element.key().toFasta() << std::endl;
		//
		//         std::pop_heap(weakHeap.begin(), weakHeap.end());
		//         weakHeap.pop_back();
		//     }

		for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
		weak.remove( it->key() );

		return promoted;
	}
};

#endif

