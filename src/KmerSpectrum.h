//
// Kmernator/src/KmerSpectrum.h
//
// Author: Rob Egan, Craig Furman
//
// Copyright 2010 The Regents of the University of California.
// All rights reserved.
//
// The United States Government has rights in this work pursuant
// to contracts DE-AC03-76SF00098, W-7405-ENG-36 and/or
// W-7405-ENG-48 between the United States Department of Energy
// and the University of California.
//
// Redistribution and use in source and binary forms are permitted
// provided that: (1) source distributions retain this entire
// copyright notice and comment, and (2) distributions including
// binaries display the following acknowledgement:  "This product
// includes software developed by the University of California,
// JGI-PSF and its contributors" in the documentation or other
// materials provided with the distribution and in all advertising
// materials mentioning features or use of this software.  Neither the
// name of the University nor the names of its contributors may be
// used to endorse or promote products derived from this software
// without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE.
//

#ifndef _KMER_SPECTRUM_H
#define _KMER_SPECTRUM_H

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <sys/mman.h>

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
#include "KmerTrackingData.h"
#include "Options.h"
#include "Log.h"

template<typename So, typename We, typename Si = TrackingDataSingleton>
class KmerSpectrum {
public:

	typedef Kmernator::KmerIndexType IndexType;

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

	typedef TrackingData::WeightType WeightType;
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	typedef TrackingData::PositionType PositionType;

	typedef typename SolidMapType::Iterator SolidIterator;
	typedef typename SolidMapType::ElementType SolidElementType;
	typedef typename SolidMapType::BucketType SolieBucketType;
	typedef typename WeakMapType::Iterator WeakIterator;
	typedef typename WeakMapType::ElementType WeakElementType;
	typedef typename WeakMapType::BucketType WeakBucketType;
	typedef typename SingletonMapType::Iterator SingletonIterator;
	typedef typename SingletonMapType::ElementType SingletonElementType;
	typedef typename SingletonMapType::BucketType SingletonBucketType;

	typedef std::vector< KmerSpectrum > Vector;

	typedef Kmer::NumberType NumberType;

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
	bool hasSingletons;
	unsigned long purgedSingletons;

public:
	KmerSpectrum(unsigned long buckets = 0, bool separateSingletons = true):
		solid(buckets/64), weak(separateSingletons ? buckets/8 : buckets), singleton(separateSingletons ? buckets : 1),
		hasSolids(false), hasSingletons(separateSingletons), purgedSingletons(0)
	{
		// set the minimum weight that will be used to track kmers
		// based on the given options
		TrackingData::minimumWeight = Options::getMinKmerQuality();
		TrackingData::minimumDepth = Options::getMinDepth();
		// apply the minimum quality automatically
		Read::setMinQualityScore( Options::getMinQuality(), Read::FASTQ_START_CHAR );
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
		this->hasSingletons = other.hasSingletons;
		this->purgedSingletons = other.purgedSingletons;
		return *this;
	}

	void setSolidOnly() {
		weak.clear(true);
		singleton.clear(true);
		hasSolids = true;
	}
	void reset(bool releaseMemory = true) {
		solid.clear(releaseMemory);
		weak.clear(releaseMemory);
		singleton.clear(releaseMemory);
		hasSolids = false;
		purgedSingletons = 0;
	}

	static unsigned long estimateWeakKmerBucketSize( ReadSet &store, unsigned long targetKmersPerBucket = 128) {
		unsigned long baseCount = store.getBaseCount();
		if (baseCount == 0 || store.getSize() == 0)
		return 128;
		unsigned long avgSequenceLength = baseCount / store.getSize();
		unsigned long kmersPerRead = (avgSequenceLength - KmerSizer::getSequenceLength() + 1);
		if (KmerSizer::getSequenceLength() > avgSequenceLength)
			kmersPerRead = 1;

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
		if (hasSingletons)
		    return singleton.exists( kmer );
		else
			return false;
	}
	SingletonElementType getIfExistsSingleton( const Kmer &kmer ) {
		if (hasSingletons)
		    return singleton.getElementIfExists( kmer );
		else
			return SingletonElementType();
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
				if (spectrum->hasSingletons && !weakElem.isValid())
				singletonElem = spectrum->getIfExistsSingleton(kmer);
			}
		}
		double getCount(bool useWeights = false) const {
			if (useWeights) {
				if (spectrum->hasSolids && solidElem.isValid())
					return solidElem.value().getWeightedCount();
				else if (weakElem.isValid())
					return weakElem.value().getWeightedCount();
				else if (singletonElem.isValid())
					return singletonElem.value().getWeightedCount();
				else
					return 0.0;
			} else {
				if (spectrum->hasSolids && solidElem.isValid())
					return solidElem.value().getCount();
				else if (weakElem.isValid())
					return weakElem.value().getCount();
				else if (singletonElem.isValid())
					return singletonElem.value().getCount();
				else
					return 0.0;
			}
		}
	};

	// returns the count for a single kmer
	double getCount(const Kmer &kmer, bool useWeights = false) {
		DataPointers pointer(*this, kmer);
		return getCount(pointer, useWeights);
	}
	double getCount(DataPointers &pointer, bool useWeights = false) {
		return pointer.getCount(useWeights);
	}

	// sets the counts for each kmer
	void getCounts(KmerWeights &weights, bool useWeights = false) {
		DataPointers pointers(*this);
		for(SequenceLengthType idx = 0; idx < weights.size(); idx++) {
			pointers.set(weights[idx]);
			weights.valueAt(idx) = pointers.getCount(useWeights);
		}
	}

	void consolidate(const Kmer &myKmer, bool useWeights = false) {
		// find all matches to this spectrum
		KmerWeights weights = KmerWeights::permuteBases(myKmer, useWeights);
		getCounts(weights, useWeights);
		double myCount = getCount(myKmer, useWeights);
		consolidate(myKmer, myCount, weights);
	}

	// test and optionally shuffle lesser kmer counts from 1st degree permutations into myKmer
	void consolidate(const Kmer &myKmer, double myCount, const KmerWeights &weights, bool useWeights = false) {
		SequenceLengthType size = weights.size();
		for(SequenceLengthType idx = 0; idx < size; idx++) {
			double testCount = weights.valueAt(idx);
			if (testCount > 0.0 && myCount >= testCount) {
				#pragma omp critical (KS_consolidate)
				{
				  // refresh myCount within critical block
				  myCount = getCount(myKmer, useWeights);
				  if (myCount >= testCount && ( myCount > testCount || myKmer > weights[idx]) ) {
					migrateKmerData(weights[idx], myKmer);
				  }
				}
			}
		}
	}

	// adds/moves tracking from one kmer to another
	void migrateKmerData(const Kmer &srcKmer, const Kmer &dstKmer) {
		bool trans = false;
		{
			LOG_DEBUG(6, "Merging kmers " << srcKmer.toFasta() << " to " << dstKmer.toFasta());

			if (hasSolids) {
			    SolidElementType src = getIfExistsSolid(srcKmer);
			    if (src.isValid()) {
			    	SolidElementType dst = getSolid(dstKmer);
			    	dst.value().add(src.value());
			    	trans = true;
			    	src.value().reset();
			    }
		    }
		    if (!trans) {
		    	WeakElementType src = getIfExistsWeak(srcKmer);
		    	if (src.isValid()) {
		    		WeakElementType dst = getWeak(dstKmer);
		    		dst.value().add(src.value());
		    		trans = true;
		    		src.value().reset();
		    	}
		    }
		    if (!trans && hasSingletons) {
		    	SingletonElementType src = getIfExistsSingleton(srcKmer);
		    	if (src.isValid()) {
		    		WeakElementType dst = getIfExistsWeak(dstKmer);
		    		if (dst.isValid()) {
		    			// dst Weak exists, so add this src singlton
		    			dst.value().add(src.value());
		    		} else {
		    			// dst Weak does not exist
		    			SingletonElementType dst2 = getSingleton(dstKmer);
		    			if (dst2.isValid()) {
		    				// dst Singleton exists, promote both to new Weak entry
		    				dst = getWeak(dstKmer);
		    				dst.value() = dst2.value();
		    				dst.value().add(src.value());
		    				dst2.value().reset();
		    				singleton.remove(dstKmer);
		    			} else {
		    				// migrate singleton to dst
		    				dst2.value() = src.value();
		    			}
		    		}
		    		trans = true;
		    		src.value().reset();
		    	}
		    }
		}
		assert(trans);
	}


	class Histogram {
	public:

		class HistogramElement {
		public:
			unsigned long visits;
			unsigned long visitedCount;
			double visitedWeight;
			HistogramElement() : visits(0), visitedCount(0), visitedWeight(0.0) {}
			HistogramElement &operator+(const HistogramElement rh) {
				visits += rh.visits;
				visitedCount += rh.visitedCount;
				visitedWeight += rh.visitedWeight;
				return *this;
			}
		};

		typedef std::vector< HistogramElement > BucketsType;

	protected:
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
		Histogram operator+(const Histogram rh) {
			for(int i = 0 ; i < buckets.size(); i++)
				buckets[i] += rh.buckets[i];
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
		void reset() {
			int s = buckets.size();
			resetTotals();
			buckets.clear();
			buckets.resize(s);
		}
		void finish() {
			resetTotals();
			for(size_t i = 0; i<buckets.size(); i++) {
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

			ss << "Weights:\t" << totalWeightedCount
			<< "\t" << (totalWeightedCount/count)
			<< "\t" << (totalWeightedCount/totalCount)
			<< std::endl;

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
		void set(KmerSpectrum &ks, bool solidOnly = false, int rank = 0, int size = 1) {
			reset();
			if ( ! solidOnly ) {
				for(WeakIterator it(ks.weak.begin(rank, size)), itEnd(ks.weak.end()); it != itEnd; it++) {
					WeakDataType &data = it->value();
					addRecord( data.getCount(), data.getWeightedCount() );
				}
				for(SingletonIterator it(ks.singleton.begin(rank, size)), itEnd(ks.singleton.end()); it != itEnd; it++) {
					SingletonDataType &data = it->value();
					addRecord( data.getCount(), data.getWeightedCount() );
				}
				for(unsigned long i = 0; i < ks.purgedSingletons ; i++) {
					addRecord( 1, 1.0);
				}
			}
			for(SolidIterator it(ks.solid.begin(rank, size)), itEnd(ks.solid.end()); it != itEnd; it++) {
				SolidDataType &data = it->value();
				addRecord( data.getCount(), data.getWeightedCount() );
			}

		}
	};

	void printHistograms(bool printSolidOnly = false) {
		if (Log::isVerbose(1))
			printHistograms(Log::Verbose(), printSolidOnly);
	}
	void printHistograms(std::ostream &os, bool printSolidOnly = false) {
		os << getHistogram(printSolidOnly);
	}
	virtual std::string getHistogram(bool solidOnly = false) {
		Histogram histogram(63);

		histogram.set(*this, solidOnly);
		return histogram.toString();
	}

	class GCCoverageHeatMap {
	private:
		unsigned long maxCover;
		std::vector< double > zeroes;
		std::vector< std::vector< double > > gcCoverCounts;
	public:
		GCCoverageHeatMap() : maxCover(2), zeroes(KmerSizer::getSequenceLength()+1, 0), gcCoverCounts(maxCover, zeroes) {}

		void addRecord(const Kmer& kmer, unsigned long count, double weight) {
			if (count >= maxCover) {
				maxCover = count + 1;
				gcCoverCounts.resize(maxCover, zeroes);

			}
			gcCoverCounts[count][kmer.getGC()] += weight;
		}

		string toString() const {
			stringstream ss;

			for(SequenceLengthType gc = 0 ; gc <= KmerSizer::getSequenceLength(); gc++) {
				ss << "depth\t" << 100.0 * (double) gc / (double) KmerSizer::getSequenceLength();
			}
			ss << endl;
			for(unsigned long cover = 0 ; cover < maxCover; cover++) {
				ss << cover;
				for(SequenceLengthType gc = 0 ; gc < KmerSizer::getSequenceLength() + 1; gc++) {
					ss << "\t";
					double counts = gcCoverCounts[cover][gc];
					if (counts != 0.0)
						ss << counts;
				}
				ss << endl;
			}
			return ss.str();
		}
	};

	void printGC(bool printSolidOnly = false) {
		if (Log::isVerbose(1))
			printGC(Log::Verbose(), printSolidOnly);
	}
	void printGC(std::ostream &os, bool printSolidOnly = false) {
		GCCoverageHeatMap hm;

		if (!printSolidOnly) {
			setWeakGCHM(hm);
			setSingletonGCHM(hm);
		}
		setSolidGCHM(hm);

		os << hm.toString();
	}

	void setWeakGCHM(GCCoverageHeatMap &hm) {
		for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
			WeakDataType &data = it->value();
			hm.addRecord( it->key(), data.getCount(), data.getWeightedCount() );
		}
	}
	void setSolidGCHM(GCCoverageHeatMap &hm) {
		for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++) {
			SolidDataType &data = it->value();
			hm.addRecord( it->key(), data.getCount(), data.getWeightedCount() );
		}
	}
	void setSingletonGCHM(GCCoverageHeatMap &hm) {
		for(SingletonIterator it(singleton.begin()), itEnd(singleton.end()); it != itEnd; it++) {
			SingletonDataType &data = it->value();
			hm.addRecord( it->key(), data.getCount(), data.getWeightedCount() );
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
			LOG_DEBUG(2, "There are no eligible kmers to promote");
			return 0;
		}
		LOG_DEBUG(4, "Heap: " << weakHeap.size() << " " << weakHeap[0].value());

		//heapify kmer counts
		std::make_heap( weakHeap.begin(), weakHeap.end() );

		// work in batches, stop when 0 new solids are added
		unsigned long count = 0;
		unsigned long promotedInBatch = 0;

		while ( weakHeap.begin() != weakHeap.end() ) {
			WeakElementType &element = weakHeap.front();
			unsigned long lastCount = element.value().getCount();
			if (shouldBeSolid( element, minWeakRatio, minSolidRatio )) {
				SolidElementType solidElement = getSolid( element.key() );
				solidElement.value() = element.value();

				element.value().reset(); // to avoid double counting
				promoted++;
				promotedInBatch++;
			}
			if (++count == 1000) {
				LOG_DEBUG(3, "Added " << promotedInBatch << ", Heap size: " << weakHeap.size() << " lastCount: " << lastCount);
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


		if (Log::isDebug(1)) {
		  ostream &debug = Log::Debug() << "Quantiles: ";

		  for(int i=0; i< count; i++) {
			  debug << std::fixed << std::setprecision(4) << probabilities[i] << ": ";
			  debug << std::fixed << std::setprecision(4) << weighted_p_square_quantile(countAcc[i]) << "\t";
			  debug << std::fixed << std::setprecision(4) << weighted_p_square_quantile(weightedCountAcc[i]) << "\t";
			  debug << std::fixed << std::setprecision(4) << weighted_p_square_quantile(directionAcc[i]) << "\t";
			  debug << std::fixed << std::setprecision(4) << weighted_p_square_quantile(directionAcc[i]) << "\t";
			  debug << std::endl;
		  }
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
				SolidWeakWeightType scores = getPermutedScores( it->key(), 1.0, 0.01 );
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
			for(size_t i = 0; i< errorRatios.size(); i++)
			accumulators[i](errorRatios[i]);
		}

		MeanVectorType rates;
		for(size_t i = 0; i< accumulators.size(); i++) {
			rates.push_back( MeanType(mean(accumulators[i]), sqrt( variance(accumulators[i]) ) ) );
			LOG_DEBUG(2, "Error Rate for pos " << i << ", " << std::fixed << std::setprecision(6) << rates[i].first << " +/- " << std::fixed << std::setprecision(6) << rates[i].second );
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
		for(Kmer::IndexType i=0; i<permutations.size(); i++) {
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
				for(size_t i=0; i< positionSums.size(); i++) {
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
			for(size_t pos = 0; pos < weakCounts.size(); pos++) {
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
			for(size_t pos = 0; pos < weakCounts.size(); pos++) {
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
		for(KmerWeights::IndexType i=0; i<permutations.size(); i++) {
			pointers.set( permutations[i] );

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

	SequenceLengthType countGC(std::string fasta) {
		SequenceLengthType count=0;
		for(size_t i=0; i<fasta.length(); i++)
		if (fasta[i] == 'G' || fasta[i] == 'C') {
			count++;
		}
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
	std::string contrastSpectrums(ostream &os, KmerSpectrum &reference, double minSolidDepth = -1.0, int verbosity = 0) {
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
					os << "FW " << prettyS( it->key(), it->value());
				} else {
					// Weak and not in reference: trueWeak
					// count later
					if (verbosity > 1)
						os << "TW " << prettyS( it->key(), it->value());
				}
			} else {
				thisSolidSize++;
				if ( reference.solid.exists( it->key() ) ) {
					// correctly matches solid in reference: trueSolid
					trueSolid++;
					if (verbosity > 0)
						os << "TS " << prettyS( it->key(), it->value());
				} else {
					// exists in this but not reference: falseSolid
					falseSolid++;
					os << "FS " << prettyS( it->key(), it->value());
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
				os << "FW " << prettyW( it->key(), weak[ it->key() ] );
			} else {
				// exists in reference but not in either solid nor weak: missingSolid
				os << "MS " << pretty( it->key(), "N/A");
				missingSolid++;
			}
		}
		if (Options::getVerbosity() > 2) {
			for(WeakIterator it(weak.begin()), itEnd(weak.end()); it != itEnd; it++) {
				if (! reference.solid.exists( it->key() ))
				os << "TW " << prettyW( it->key(), it->value());
			}
		}
		trueWeak = thisWeakSize - falseWeak;

		std::stringstream header;
		header << "SolidSizes: " << thisSolidSize << " / " << referenceSolidSize << std::endl;
		header << "trueSolid: " << trueSolid << " MissingSolid: " << missingSolid << " falseWeak: " << falseWeak << " falseSolid: " << falseSolid << std::endl;
		header << "trueWeak: " << trueWeak << " / " << thisWeakSize << std::endl;
		return header.str();
	}

	inline void append( KmerWeights &kmers, unsigned long readIdx, bool isSolid = false, NumberType partIdx = 0, NumberType partBitMask = 0) {
		return append(kmers, readIdx, isSolid, 0, kmers.size(), partIdx, partBitMask);
	}

	inline void append(DataPointers &pointers, Kmer &least, WeightType &weight, ReadSetSizeType &readIdx, PositionType &readPos, bool isSolid = false) {
		bool keepDirection = true;
		if (weight < 0.0) {
			keepDirection = false;
			weight = 0.0-weight;
		}
		if ( TrackingData::isDiscard( weight ) )  {
			LOG_DEBUG(4, "discarded kmer " << readIdx << "@" << readPos << " " << weight << " " << least.toFasta());
			return;
		}

		if ( isSolid ) {
			pointers.reset();
			SolidElementType elem = getSolid( least );
			elem.value().track( weight, keepDirection, readIdx, readPos );
		} else {
			pointers.set( least );

			if (pointers.solidElem.isValid()) {
				// track solid stats
				pointers.solidElem.value().track( weight, keepDirection, readIdx, readPos );
			} else if (pointers.weakElem.isValid()) {
				// track weak stats
				pointers.weakElem.value().track( weight, keepDirection, readIdx, readPos );
			} else {
				if ( pointers.singletonElem.isValid() ) {

					// promote singleton to weak & track
					SingletonDataType singleData = pointers.singletonElem.value();
					pointers.reset();

					WeakElementType weakElem = getWeak( least );
					weakElem.value() = singleData;
					weakElem.value().track( weight, keepDirection, readIdx, readPos);

					singleton.remove( least );

				} else {
					if (hasSingletons) {
					    // record this new singleton
					    SingletonElementType elem = getSingleton(least);
					    elem.value().track( weight, keepDirection, readIdx, readPos);
					} else {
						// record in the weak spectrum
						WeakElementType weakElem = getWeak( least );
						weakElem.value().track( weight, keepDirection, readIdx, readPos);
					}
				}
			}
		}

	}
	void append( KmerWeights &kmers, ReadSetSizeType readIdx, bool isSolid, ReadSetSizeType kmerIdx, PositionType kmerLen, int partIdx, NumberType partBitMask ) {

		DataPointers pointers(*this);

		for(PositionType readPos=0; readPos < kmerLen; readPos++)
		{
			Kmer &least = kmers[kmerIdx+readPos];
			if (partBitMask > 1 && getDMPThread(least, partBitMask, true) != partIdx)
				continue;
			WeightType weight = kmers.valueAt(kmerIdx+readPos);
			append(pointers, least, weight, readIdx, readPos, isSolid);
		}
	}

	ostream &printStats(ostream &os, unsigned long pos, bool printSolidOnly = false, bool fullStats = false) {
		//stats.printHistograms(printSolidOnly);
		os << pos << " reads" << "\t";
		if (fullStats) {
			os << ", " << solid.size() << " solid / " << weak.size() << " weak / "
			<< TrackingData::discarded << " discarded / "
			<< singleton.size() + purgedSingletons << " singleton (" << purgedSingletons << " purged) - kmers so far " << std::endl;
		}
		os << MemoryUtils::getMemoryUsage() << std::endl;
		return os;
	}

	ostream &printBucketStats(ostream &os, unsigned long pos, bool printSolidOnly = false, bool fullStats = false) {
		os << pos << ": ";

		unsigned long ss = solid.size();
		unsigned long buckets = solid.getNumBuckets();
		static unsigned long prev_ss = 0;
		os << ss << " added : " << ss - prev_ss << ",  average:"
		<< ss/buckets << " Max:" << solid.maxBucket() << " " << MemoryUtils::getMemoryUsage() << std::endl;
		prev_ss = ss;
		return os;
	}

	inline bool getSMPThread( Kmer &kmer, int &smpThreadId, int numSMPThreads, int dmpThreadId, int threadBitMask, bool isLeastComplement = false) {
		// return the smallest KmerMap in the spectrum's getLocalThreadId()
		if (isLeastComplement) {
			return solid.getLocalThreadId(kmer, smpThreadId, numSMPThreads, dmpThreadId, threadBitMask);
		} else {
		    TEMP_KMER(tmp);
		    kmer.buildLeastComplement(tmp);
		    return solid.getLocalThreadId(tmp, smpThreadId, numSMPThreads, dmpThreadId, threadBitMask);
		}
	}

	inline int getSMPThread( Kmer &kmer, int numThreads, bool isLeastComplement = false) {
		// return the smallest KmerMap in the spectrum's getLocalThreadId()
		if (isLeastComplement) {
			return solid.getLocalThreadId(kmer, numThreads);
		} else {
		    TEMP_KMER(tmp);
		    kmer.buildLeastComplement(tmp);
		    return solid.getLocalThreadId(tmp, numThreads);
		}
	}
	inline int getDMPThread( Kmer &kmer, NumberType threadBitMask, bool isLeastComplement = false ) {
		if (isLeastComplement) {
			return solid.getDistributedThreadId(kmer, threadBitMask);
		} else {
		    TEMP_KMER(tmp);
		    kmer.buildLeastComplement(tmp);
		    return solid.getDistributedThreadId(tmp, threadBitMask);
		}
	}

	unsigned long resetSingletons() {
		if (Options::getMinDepth() <= 1) {
			return 0;
		}
		unsigned long singletonCount = singleton.size();
		purgedSingletons += singletonCount;
		LOG_DEBUG(2, "Purging Singletons: " << singletonCount << " (" << purgedSingletons << " total so far)");
		singleton.reset();
		return singletonCount;
	}

	void purgeMinDepth(long minimumCount, bool purgeSolidsToo = false) {
		if (hasSolids && purgeSolidsToo)
			solid.purgeMinCount(minimumCount);
		weak.purgeMinCount(minimumCount);
		if (hasSingletons && minimumCount > 1)
			singleton.clear(false);
	}

	// important! returned memory maps must remain in scope!
	Kmernator::MmapFileVector buildKmerSpectrumInParts(ReadSet &store, NumberType numParts) {
		bool isSolid = false; // not supported for references...
		if (numParts == 0) {
			buildKmerSpectrum(store, isSolid);

			// purge
			purgeMinDepth(Options::getMinDepth());
			return Kmernator::MmapFileVector();
		}

		assert((numParts & (numParts-1)) == 0); // numParts must be a power of 2
		Kmernator::MmapFileVector mmaps(numParts*2);

		// build each part of the spectrum
		for (NumberType partIdx = 0; partIdx < numParts; partIdx++) {
			LOG_VERBOSE(2, "Building part of spectrum: " << (partIdx+1) << " of " << numParts << std::endl << MemoryUtils::getMemoryUsage());

			// build
			buildKmerSpectrum(store, isSolid, partIdx, numParts);

			// purge
			purgeMinDepth(Options::getMinDepth());

			// store
			mmaps[partIdx] = weak.store("WeakKmer");
			if (Options::getMinDepth() <= 1)
				mmaps[partIdx + numParts] = singleton.store("SingletonKmer");
			else if (Options::getVerbosity())
				LOG_VERBOSE(1, "Not storing singletons which would have been this size: " << singleton.getSizeToStore() );

			// optimize memory preallocations
			weak.rotateDMPBuffers(numParts);
			singleton.rotateDMPBuffers(numParts);
		}


		// first free up memory
		if (Options::getMinDepth() <= 1) {
			LOG_DEBUG(2, "Clearing memory from singletons" << std::endl << MemoryUtils::getMemoryUsage());
		    singleton.clear();
		}

		LOG_VERBOSE(2, "Merging partial spectrums" );
		LOG_DEBUG(2, MemoryUtils::getMemoryUsage() );

		const WeakMapType constWeak(weak.getNumBuckets());
		const SingletonMapType constSingleton(singleton.getNumBuckets());

		// restore and merge
		for (NumberType partIdx = 0; partIdx < numParts; partIdx++) {
			const WeakMapType tmpWMap = WeakMapType::restore( mmaps[partIdx].data() ) ;
			constWeak.merge( tmpWMap );
			if (mmaps[partIdx+numParts].is_open()) {
				const SingletonMapType tmpSMap = SingletonMapType::restore( mmaps[partIdx+numParts].data() );
				constSingleton.merge( tmpSMap );
			}
		}
		weak.swap( const_cast<WeakMapType&>(constWeak) );
		if (Options::getMinDepth() <= 1)
			singleton.swap( const_cast<SingletonMapType&>( constSingleton) );

		LOG_VERBOSE(2, "Finished merging partial spectrums" << std::endl << MemoryUtils::getMemoryUsage());

		if (Log::isVerbose(1)) {
			printStats(Log::Verbose("Final Stats"), store.getSize(), isSolid, true);
			if (!isSolid) {
				printHistograms(Log::Verbose("Final Histogram"));
			}
		}

		for(Kmernator::MmapFileVector::iterator it = mmaps.begin(); it != mmaps.end(); it++) {
			if (it->is_open())
				madvise(const_cast<char*>(it->data()), it->size(), MADV_RANDOM);
		}
		return mmaps;
	}

	void _evaluateBatch(bool isSolid, long batchIdx, long purgeEvery, long purgeCount) {
		if (Log::isVerbose(2)) {
			printStats(Log::Verbose("Batch Stats"), batchIdx, isSolid);
		}
		if (purgeEvery > 0 && batchIdx >= (purgeCount+1)*purgeEvery) {
			purgedSingletons += resetSingletons();
			purgeCount++;
		}

	}
	void _buildKmerSpectrumSerial(ReadSet &store, bool isSolid, NumberType partIdx, NumberType partBitMask, long batch, long purgeEvery, long &purgeCount) {
		for (long batchIdx=0; batchIdx < (long) store.getSize(); batchIdx++)
		{
			const Read &read = store.getRead( batchIdx );
			if ( !read.isDiscarded() ) {
				KmerWeights kmers = KmerReadUtils::buildWeightedKmers(read, true, true);

				append(kmers, batchIdx, isSolid, partIdx, partBitMask);
			}

			if (batchIdx % batch == 0) {
				_evaluateBatch(isSolid, batchIdx, purgeEvery, purgeCount);
			}
		}

	}
	void _buildKmerSpectrumParallel(ReadSet &store, bool isSolid, NumberType partIdx, NumberType partBitMask, long batch, long purgeEvery, long &purgeCount) {

		long numThreads = omp_get_max_threads();
		long batchIdx = 0;

		// allocate a square matrix of buffers: KmerWeights[writingThread][readingThread]
		KmerWeights kmerBuffers[ numThreads ][ numThreads ];
		typedef std::pair< ReadSet::ReadSetSizeType, Sequence::SequenceLengthType > ReadPosType;
		std::vector< ReadPosType > startReadIdx[ numThreads ][ numThreads ];

		long size = store.getSize();
		if (size < batch) {
		    batch = store.getSize() / 10 + 1;
		}

		long reservation;
		if (KmerSizer::getSequenceLength() > store.getMaxSequenceLength()) {
			LOG_WARN(1, "KmerSize " << KmerSizer::getSequenceLength() << " is larger than your data set: " << store.getMaxSequenceLength());
			reservation = batch * 2;
		} else {
			reservation = batch * (store.getMaxSequenceLength() - KmerSizer::getSequenceLength() + 1);
		}
		if (numThreads > 1)
			reservation /= numThreads;

		#pragma omp parallel num_threads(numThreads)
		{
			if (numThreads != omp_get_num_threads())
				throw "OMP thread count mis-match";
		}
		LOG_DEBUG(1, "Executing parallel buildKmerSpectrum with " << numThreads);

		while (batchIdx < (long) store.getSize())
		{
			for(int tmpThread1=0; tmpThread1 < numThreads; tmpThread1++)
			{
				for(int tmpThread2 = 0; tmpThread2 < numThreads; tmpThread2++) {
					kmerBuffers [ tmpThread1 ][ tmpThread2 ].reset(false);
					kmerBuffers [ tmpThread1 ][ tmpThread2 ].reserve(reservation);
					startReadIdx[ tmpThread1 ][ tmpThread2 ].resize(0);
					startReadIdx[ tmpThread1 ][ tmpThread2 ].reserve(reservation);
				}
			}


            #pragma omp parallel for schedule(dynamic) num_threads(numThreads)
			for (long i=0; i < batch; i++)
			{
				long readIdx = batchIdx + i;

				if (readIdx >= size ) {
					continue;
				}
				const Read &read = store.getRead( readIdx );
				if (read.isDiscarded())
					continue;

				KmerWeights kmers = KmerReadUtils::buildWeightedKmers(read, true, true);
				for (long j = 0; j < numThreads; j++)
				startReadIdx[ omp_get_thread_num() ][ j ].push_back( ReadPosType(readIdx, kmerBuffers[ omp_get_thread_num() ][j].size()) );
				for (IndexType j = 0; j < kmers.size(); j++) {
					int smpThreadId;
					if ( getSMPThread( kmers[j], smpThreadId, numThreads, partIdx, partBitMask, true) )
					    kmerBuffers[ omp_get_thread_num() ][ smpThreadId ].append(kmers[j], kmers.valueAt(j));
				}
			}

			#pragma omp parallel num_threads(numThreads)
			{
				if (numThreads != omp_get_num_threads()) {
					throw "OMP thread count mis-match";
				}

				for(int threads = 0; threads < numThreads; threads++)
				{
					int threadId = omp_get_thread_num();
					for(unsigned long idx = 0; idx < startReadIdx[ threads ][ threadId ].size(); idx++)
					{
						ReadPosType readPos = startReadIdx[ threads ][ threadId ][idx];
						unsigned long startIdx = readPos.second;
						unsigned long len = 0;
						if ( idx < startReadIdx[ threads ][ threadId ].size() -1) {
							len = startReadIdx[ threads ][ threadId ][idx+1].second - startIdx;
						} else {
							len = kmerBuffers[threads][ threadId ].size() - startIdx;
						}
						if (len > 0) {
							// do not include partIdx or partBitmMask , as getSMPThread() already accounted for partitioning
						    append(kmerBuffers[threads][ threadId ], readPos.first, isSolid, startIdx, len, 0, 0);
						}
					}
				}

			}

			batchIdx += batch;
			_evaluateBatch(isSolid, batchIdx, purgeEvery, purgeCount);
		}

	}

	virtual void buildKmerSpectrum( ReadSet &store ) {
		return this->buildKmerSpectrum(store, false);
	}
	virtual void buildKmerSpectrum( ReadSet &store, bool isSolid ) {
		return this->buildKmerSpectrum(store, isSolid, 0, 1);
	}
	void buildKmerSpectrum( ReadSet &store, bool isSolid, NumberType partIdx, NumberType numParts)
	{
		assert( numParts != 0 && (numParts & (numParts-1)) == 0); // numParts must be a power of 2
		assert(partIdx < numParts);
		NumberType partBitMask = numParts-1;

		solid.reset(false);
		if (isSolid) {
			// no need to keep this memory around then...
			weak.clear();
			singleton.clear();
		} else {
		    weak.reset(false);
		    singleton.reset(false);
		}

		long purgeEvery = Options::getPeriodicSingletonPurge();
		long purgeCount = 0;
		long batch = Options::getBatchSize();

		if (omp_get_max_threads() > 1) {
#ifdef _USE_OPENMP
			_buildKmerSpectrumParallel(store, isSolid, partIdx, partBitMask, batch, purgeEvery, purgeCount);
#endif
		} else {
			_buildKmerSpectrumSerial  (store, isSolid, partIdx, partBitMask, batch, purgeEvery, purgeCount);
		}

		if (Log::isVerbose(2)) {
		    printStats(Log::Verbose("Stats"), store.getSize(), isSolid, true);
		}
		if (Log::isVerbose(1)) {
		    if (!isSolid) {
			    printHistograms(Log::Verbose("Histogram"));
		    }
		}
	}

	bool _setPurgeVariant(DataPointers &pointers, const Kmer &kmer, double threshold, double &v) {
		bool returnValue = false;
		pointers.set( kmer );
		WeakElementType &w = pointers.weakElem;
		if (w.isValid()) {
			v = w.value().getWeightedCount();
			if (v > 0.0 && v < threshold) {
				w.value().reset();
				LOG_DEBUG(4, "Purged Variant " << kmer.toFasta() << " " << v);
				returnValue = true;
			} else if (v < threshold) {
				// already been visited
				v = -1.0;
			}
		} else {
			v = 0.0;
		}
		return returnValue;
	}

	// recursively purge kmers within edit distance and below threshold
	virtual long _purgeVariants(DataPointers &pointers, const Kmer &kmer, WeakBucketType &variants, double threshold, short editDistance) {
		long purgedKmers = 0;

		if (editDistance == 0)
			return purgedKmers;

		WeakBucketType::permuteBases(kmer, variants, editDistance, true);

		for(SequenceLengthType i = 0 ; i < variants.size(); i++) {
			Kmer &varKmer = variants[i];
			double dummy;
			if (this->_setPurgeVariant(pointers, varKmer, threshold, dummy))
				purgedKmers++;
		}
		return purgedKmers;
	}
	virtual void _preVariants(double variantSigmas, double minDepth) {}
	virtual long _postVariants(long purgedKmers) {
		LOG_VERBOSE(1, "Removed " << purgedKmers << " kmer-variants");
		return purgedKmers;
	}
	virtual void _variantThreadSync(long processed, long remaining, double maxDepth) {
		LOG_DEBUG(2, "Batch processed " << processed << " remaining: " << remaining << " threshold: " << maxDepth);
	}
	virtual long _variantBatchSync(long remaining, long purgedKmers, double maxDepth, double threshold) {
		LOG_DEBUG(2, "Purged " << purgedKmers << " variants below: " << maxDepth << " / " <<  threshold << " with " << remaining << " remaining");
		return remaining;
	}

	double getVariantThreshold(double count, double variantSigmas) {
		return count - sqrt(count)*variantSigmas;
	}

	double getNewMaxDepth(double maxDepth, double variantSigmas) {
		// quadratic equation solving maxDepth = newMaxDepth - variantSigmas * sqrt(newMaxDepth)
		// 0 == newMaxDepth*newMaxDepth + (variantSigmas*variantSigmas - 2*maxDepth) + maxDepth*maxDepth
		// a=1 b=(-variantsigmas**2 - 2*maxDepth) c=maxDepth**2
		double a=1.0;
		double b= 0.0 - variantSigmas*variantSigmas - 2.0 * maxDepth;
		double c=maxDepth*maxDepth;
		maxDepth = (-b + sqrt(b*b - 4*a*c)) / (2.0*a);
		return maxDepth;
	}
	long purgeVariants(double variantSigmas = Options::getVariantSigmas(), short editDistance = 2) {
		if (variantSigmas < 0.0)
			return 0;
		long purgedKmers = 0;
		LOG_VERBOSE(1, "Purging kmer variants within " << editDistance << " edit distance which are >= " << variantSigmas << " sigmas less abundant than a more abundant version");

		if (hasSolids) {
			throw; // TODO unsupported
		}
		double minDepth = Options::getMinDepth();
		int numThreads = omp_get_max_threads();
		this->_preVariants(variantSigmas, minDepth);
		LOG_DEBUG(1, "Purging with " << numThreads << " threads");

		std::vector<DataPointers> pointers(numThreads, DataPointers(*this));
		std::vector<WeakBucketType> variants(numThreads);

		// sweep through bottom up.  Start at first possible variant to compare at
		double maxDepth = getNewMaxDepth(minDepth, variantSigmas);
		long remaining = 1;
		long processed = 0;
		while (remaining > 0) {
			double lastDepth = maxDepth;
			// increment maxDepth by a step
			maxDepth = getNewMaxDepth(maxDepth, variantSigmas);

			remaining = 0;
			processed = 0;
			LOG_DEBUG_OPTIONAL(2, true, "Starting threaded kmer scan");
			#pragma omp parallel num_threads(numThreads) reduction(+: purgedKmers) reduction(+: processed) reduction(+: remaining)
			{
				int threadId = omp_get_thread_num();

				for(WeakIterator it = weak.beginThreaded(); it != weak.endThreaded(); it++) {
					double count = it->value().getWeightedCount();
					if (count <= lastDepth)
						continue;
					if (count > maxDepth) {
						remaining++;
						continue;
					}
					double threshold = getVariantThreshold(count, variantSigmas);
					if (threshold > minDepth) {
						LOG_DEBUG(3, "Purging Variants of " << it->key().toFasta() << " below " << threshold << " (" << count << ")");
						purgedKmers += this->_purgeVariants(pointers[threadId], it->key(), variants[threadId], threshold, editDistance);
					}
					if (++processed % 10000 == 0)
						LOG_DEBUG(2, "progress processed " << processed);
				}
				this->_variantThreadSync(processed, remaining, maxDepth);
			}
			remaining = this->_variantBatchSync(remaining, purgedKmers, maxDepth, getVariantThreshold(maxDepth, variantSigmas));
		}

		LOG_DEBUG(3, "Finished processing variants: " << purgedKmers << " waiting for _postVariants");
		purgedKmers = this->_postVariants(purgedKmers);

		LOG_DEBUG(1, "Purging to min depth: " << Options::getMinDepth());
		this->purgeMinDepth(Options::getMinDepth());

		return purgedKmers;
	}

	static void experimentOnSpectrum( ostream &os, KmerSpectrum &spectrum ) {
		KmerSpectrum everything(spectrum);

		unsigned long promoted = spectrum.promote(0.05);
		os << "Promoted " << promoted << " kmers" << std::endl;
		spectrum.printHistograms(true);

		for (double prob = 0.10; prob > 0.0; prob -= 0.005) {
			KmerSpectrum copy(everything);
			os << "Testing " << prob << std::endl;
			promoted = copy.promote(prob);
			os << "Promoted " << promoted << " kmers" << std::endl;
			copy.printHistograms(true);
			copy.contrastSpectrums(os, spectrum);
		}

		for(WeakIterator it(everything.weak.begin()), itEnd(everything.weak.end()); it != itEnd; it++) {
			everything.getSolid( it->key() ) = it->value();
		}
		everything.weak.clear();

		for (double minDepth = 1; minDepth < 10; minDepth++) {
			os << "Comparing minDepth " << minDepth << std::endl;
			everything.contrastSpectrums(os, spectrum, minDepth);
		}
	}

	void analyseSingletons(ostream &os)
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
		os << "Singleton permutations found in weak : " << found << std::endl;
	}

	void findSecondOrderPermutations(ostream &os, const Kmer &kmer)
	{
		Kmers permutations = Kmers::permuteBases(kmer, true);
		for(int i=0; i<permutations.size(); i++) {

			Kmers permutations2 = Kmers::permuteBases(permutations[i],true);
			for (int j = i; j<permutations2.size(); j++) {
				if (permutations2[j] != kmer)
				// if  (SolidDataType *s = solid.getIfExists(permutations2[j])  )
				{
					os << "   " << permutations2[j].toFasta() /*<< " : " << s->getCount() */<< std::endl;
				}
			}
		}
		os << "\n";
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

	unsigned long filter(ostream &os, KmerSpectrum &reference)
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

		os << mean << " +-" << deviation << " => cutoffs: " << poorQualityCutoff << " " << goodQualityCutoff << "\n";

		WeakIterator mapIt = weak.begin();
		WeakIterator mapEnd = weak.end();
		for(; mapIt != mapEnd; mapIt++ )
		{

			WeakElementType &element = *mapIt;
			unsigned long lastCount = element.value().getCount();

			if (lastCount < minKmerCount)
			continue;

			count++;
			bool inRef = reference.solid.exists( element.key());

			bool poorQuality = element.value().getAverageWeight() < poorQualityCutoff;
			bool goodQuality = element.value().getAverageWeight() > goodQualityCutoff;

			if (firstOrderMatch(element))
			{
				fo_count++;
				if (!inRef)
				fo_correct++;
				else
				{
					fo_wrong++;
				}

			}
        }
        
		os << "Poor Quality: " << pq_count << " (" << pq_count*100.0/count << "%) Correct:  " << pq_correct << " Wrong: " << pq_wrong << std::endl;
		os << "1st Order  : " << fo_count << " (" << fo_count*100.0/count << "%) Correct:  " << fo_correct << " Wrong: " << fo_wrong << std::endl;

		for(SolidIterator it(solid.begin()), itEnd(solid.end()); it != itEnd; it++)
		weak.remove( it->key() );

		return promoted;
	}

	void static mergeVector(Vector &vec, int minimumCount = 2) {
		// solid
		if (vec.size() <= 1)
			return;
		if (vec[0].hasSolids) {
		  for (size_t i = 1; i < vec.size(); i++) {
			vec[0].solid.mergeAdd(vec[i].solid);
		  }
		}
		// weak
		for (size_t i = 1; i < vec.size(); i++) {
			vec[0].weak.mergeAdd(vec[i].weak);
		}
		if (vec[0].hasSingletons) {
		  // singleton
		  for (size_t i = 1; i < vec.size(); i++) {
			vec[0].singleton.mergePromote(vec[i].singleton, vec[0].weak);
		  }
		  if (minimumCount > 1)
			  vec[0].singleton.clear();
		}
		if (minimumCount > 1)
		    vec[0].purgeMinDepth(minimumCount);
	}
};

#endif


// $Log: KmerSpectrum.h,v $
// Revision 1.39  2010-08-18 17:50:39  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.38.4.1  2010-07-20 20:02:56  regan
// autodetect fastq quality range
//
// Revision 1.38  2010-06-22 23:06:31  regan
// merged changes in CorruptionBugfix-20100622 branch
//
// Revision 1.37.4.1  2010-06-22 23:00:03  regan
// named all critical sections
//
// Revision 1.37  2010-05-24 21:48:46  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.36.2.4  2010-05-20 18:26:43  regan
// attempt to fix a race condition when consolidating/merging edit-distance spectrums
//
// Revision 1.36.2.3  2010-05-19 23:41:06  regan
// re-added memory optimizations
//
// Revision 1.36.2.2  2010-05-19 23:38:18  regan
// bug and performance fixes
//
// Revision 1.36.2.1  2010-05-19 21:53:20  regan
// bugfixes
//
// Revision 1.36  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.35.2.8  2010-05-18 16:43:47  regan
// added GC heatmap output .. still refining
//
// Revision 1.35.2.7  2010-05-17 17:52:43  regan
// working through some optimizations
//
// Revision 1.35.2.6  2010-05-13 22:37:10  regan
// minor performance opt to skip discarded reads when building spectrum
//
// Revision 1.35.2.5  2010-05-13 20:29:30  regan
// new methods to support changes to CompareSpectrum
//
// Revision 1.35.2.4  2010-05-12 22:44:00  regan
// bugfix
//
// Revision 1.35.2.3  2010-05-12 19:52:23  regan
// refactored options
// removed obsolete parameters
// added new ones
//
// Revision 1.35.2.2  2010-05-10 17:58:08  regan
// fixing types
//
// Revision 1.35.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.35  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.34  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
//
