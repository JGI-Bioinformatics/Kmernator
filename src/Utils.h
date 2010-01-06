// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.22 2010-01-06 15:20:24 regan Exp $
//

#ifndef _UTILS_H
#define _UTILS_H

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
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/unordered_map.hpp>


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
	// assume 1% error rate
	unsigned long estimatedUniqueKmers = rawKmers * ( std::max(pow(1.01, KmerSizer::getSequenceLength()),2.0) - 1.0 );
	unsigned long targetBuckets = estimatedUniqueKmers / targetKmersPerBucket;
	unsigned long maxBuckets = 128*1024*1024;
	unsigned long minBuckets = 128;
	return std::max( std::min(targetBuckets,maxBuckets), minBuckets );
}



KmerWeights buildWeightedKmers(Read &read, bool leastComplement = false) {
  KmerWeights kmers(read.getTwoBitSequence(), read.getLength(), leastComplement);
  std::string quals = read.getQuals();
  SequenceLengthType markupIdx = 0;
  
  BaseLocationVectorType markups = read.getMarkups();
  double weight = 0.0;
  double change = 0.0;
  
  for(SequenceLengthType i=0; i< kmers.size(); i++) {
  	if (i%100 == 0 || weight == 0.0) {
  	  weight = 1.0;
  	  for(SequenceLengthType j=0; j < KmerSizer::getSequenceLength(); j++) 
  	    weight *= Read::qualityToProbability[ (unsigned char) quals[i+j] ];
  	} else {
  		change = Read::qualityToProbability[ (unsigned char) quals[i+KmerSizer::getSequenceLength()-1] ] / 
  		         Read::qualityToProbability[ (unsigned char) quals[i-1] ];
  		weight *= change;
  	}
  	while( markupIdx < markups.size() && markups[markupIdx].second < i)
  	  markupIdx++;
  	if (markupIdx < markups.size() && markups[markupIdx].second < i+KmerSizer::getSequenceLength()) {
  	    weight = 0.0;
//  	    std::cerr << "markupAt " << markups[markupIdx].first << " " << markups[markupIdx].second << "\t";
    }
    if (kmers[i].isTrivial())
        weight = 0.0;
  	kmers.valueAt(i) = weight;
//  	std::cerr << i << ":" << std::fixed << std::setprecision(3) << quals[i+KmerSizer::getSequenceLength()-1] << "-" << change << "/" << weight << "\t";
  	
  }
//  std::cerr << std::endl;
  return kmers;
}

template< typename M >
class ReadSelector
{
public:
  struct ReadTrimType { 
  	ReadSet::ReadSetSizeType readIdx;
    Sequence::SequenceLengthType trimOffset; 
    std::string label;
         };
  typedef M DataType;
  typedef typename DataType::ReadPositionWeightVector ReadPositionWeightVector;
  typedef KmerMap<DataType> KMType;
  typedef typename KMType::ConstIterator KMIterator;
  typedef KmerMap<unsigned short> KMCacheType;
  typedef std::vector< ReadTrimType > ReadTrimVector;

private:
  const ReadSet &_reads;
  const KMType &_map;
  ReadTrimVector _picks;

public:
  ReadSelector(const ReadSet &reads, const KMType &map): _reads(reads), _map(map), _picks() { }
  
  const ReadTrimVector &getPicks() const {
  	return _picks;
  }
  
  std::ostream &writePicks(std::ostream &os) const {
  	foreach( ReadTrimType readTrim, _picks)
  	  os << _reads.getRead( readTrim.readIdx ).toFastq( readTrim.trimOffset, readTrim.label );
  	return os;
  }
  
  void pickLeastCoveringSubset() {
  	KMCacheType visits = KMCacheType(_map.getNumBuckets());
  	
  	std::vector< SequenceLengthType > readHits;
  	readHits.resize( _reads.getSize() );
  	KMIterator it( _map.begin() ), end( _map.end() );
  	for(; it != end; it++) {
  		const DataType &data = it->value();
  		ReadPositionWeightVector rpos = data.getEachInstance();
  		for(unsigned int i=0; i< rpos.size(); i++) {
  			readHits[ rpos[i].readId ]++;
  		}
  	}
  	#ifdef _USE_OPENMP

    #else
	
	#endif  	
  }
  
  void pickAllCovering() {
  	
  }
  
  void pickCoverageNormalizedSubset() {
  	
  }
  
};

class FilterKnownOddities
{
private:
  ReadSet sequences;
  boost::unordered_map< unsigned long, std::string > descriptions;
  typedef boost::unordered_map< unsigned long, unsigned long > CountMap;
  CountMap counts;
  unsigned long mask;
  unsigned short length;
public:
  FilterKnownOddities(int _length = 16) : length(_length) {
  	std::string fasta = getFasta();
  	sequences.appendFastaFile( fasta );
  	mask = ((unsigned long) -1) & ~( 0x01<<(65-length/2)- 1 );
  	prepareMaps();
  }
  
  void prepareMaps() {
  	unsigned long oldKmerLength = KmerSizer::getSequenceLength();
  	KmerSizer::set(length);
  	
  	TEMP_KMER(reverse);
  	unsigned long key;
  	for(unsigned int i = 0 ; i < sequences.getSize() ; i++) {
  		// TODO fix bit-shift boundaries, make 4 masks or make peace with sloppy frameshifts
  		
  		Read &read = sequences.getRead(i);
  		KmerWeights kmers = buildWeightedKmers(read);
  		for(unsigned int j = 0; j < kmers.size(); j++) {
  			kmers[j].buildReverseComplement(reverse);  	        
  			
  	        key = *( (unsigned long *) kmers[j].get() ) & mask;
  	        descriptions[ key ] = std::string(read.getName());
  	        counts[ key ] = 0;  	        
  	        
  	        key = *( (unsigned long *) reverse.get() ) & mask;
  	        descriptions[ key ] = std::string(read.getName() + "-reverse");
  	        counts[ key ] = 0;
  		}
  	}
  	
  	KmerSizer::set(oldKmerLength);
  }
  
  unsigned long applyFilter(ReadSet &reads) {
  	unsigned long affectedCount = 0;
  	for(unsigned long i = 0; i < reads.getSize(); i++) {
  		Read &read = reads.getRead(i);
  		bool wasAffected = false;
  		
  		TwoBitEncoding *ptr = read.getTwoBitSequence();
  		CountMap::iterator it;
  		unsigned long lastByte = read.getTwoBitEncodingSequenceLength() - TwoBitSequence::fastaLengthToTwoBitLength(length);
  		for(unsigned long bytes = 0; bytes < lastByte ; bytes++) {
  			unsigned long key = *((unsigned long *) ptr) & mask;
  			it = counts.find( key );
  			if (it != counts.end()) {
  				counts[key]++;
  				read.zeroQuals( bytes * 4, length);
  				if (!wasAffected) {
  				  wasAffected = true;
  				  affectedCount++;
  				}
  			}
  		}
  		
  	}
  	return affectedCount;
  }
  
  
  const ReadSet &getSequences() { return sequences; }
  std::string getFasta() {
  	// TODO logic to check Options
  	std::stringstream ss;
  	ss << ">PE-Adaptor-P" << std::endl;
  	ss << "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG" << std::endl;
  	ss << ">PE-PCR-PRIMER1" << std::endl;
  	ss << "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" << std::endl;
  	ss << ">PE-PCR-PRIMER2" << std::endl;
  	ss << "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT" << std::endl;
  	ss << ">SE-Adaptor-P" << std::endl;
  	ss << "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG" << std::endl;
   	ss << ">SE-Primer2" << std::endl;
  	ss << "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT" << std::endl;
  	ss << ">Homopolymer-A" << std::endl;
  	ss << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;
  	ss << ">Homopolymer-C" << std::endl;
  	ss << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << std::endl;
  	ss << ">Homopolymer-G" << std::endl;
  	ss << "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" << std::endl;
  	ss << ">Homopolymer-T" << std::endl;
  	ss << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" << std::endl;
    return ss.str();
  }
  
};


#endif

//
// $Log: Utils.h,v $
// Revision 1.22  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.21  2010-01-05 06:43:47  regan
// bugfix
//
// Revision 1.20  2009-12-24 00:56:11  regan
// started class to output picked reads
//
// Revision 1.19  2009-12-21 06:34:26  regan
// used openmp and clever partitioning to speed up building spectrum
//
// Revision 1.18  2009-11-28 01:00:07  regan
// fixed bugs and warnings
//
// Revision 1.17  2009-11-26 09:03:29  regan
// refactored and stuff
//
// Revision 1.16  2009-11-24 13:35:29  cfurman
// removed KmerPtr class.
//
// Revision 1.15  2009-11-22 08:16:41  regan
// some fixes some bugs... optimized vs debug vs deb4/5 give different results
//
// Revision 1.14  2009-11-21 18:46:53  regan
// added bugs
//
// Revision 1.13  2009-11-21 15:58:29  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.12  2009-11-12 17:01:51  regan
// checkpoint
//
// Revision 1.11  2009-11-12 01:29:22  cfurman
// Solid picking logic bug fixed, params tweaked.
//
// Revision 1.10  2009-11-11 17:23:23  regan
// fixed bugs in heap generation
// solid picking logic needs work
//
// Revision 1.9  2009-11-11 07:57:23  regan
// built framework for autoPromote (not working) - make_heap is broken
//
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
