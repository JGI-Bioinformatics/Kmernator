// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.27 2010-01-13 07:20:08 regan Exp $
//

#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "config.h"
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


class KmerReadUtils 
{
public:
  static KmerWeights buildWeightedKmers(Read &read, bool leastComplement = false) {
	  KmerWeights kmers(read.getTwoBitSequence(), read.getLength(), leastComplement);
	  std::string quals = read.getQuals();
	  SequenceLengthType markupIdx = 0;
	  
	  BaseLocationVectorType markups = read.getMarkups();
	  double weight = 0.0;
	  double change = 0.0;
	  bool isRef = false;
	  if (quals.length()>0 && quals[0] == Read::REF_QUAL)
	    isRef = true;
	  
	  for(SequenceLengthType i=0; i< kmers.size(); i++) {
	  	if (isRef) {
	  	  weight = 1.0;
	  	} else if (i%100 == 0 || weight == 0.0) {
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
	  	kmers.valueAt(i) = weight;
	//  	std::cerr << i << ":" << std::fixed << std::setprecision(3) << quals[i+KmerSizer::getSequenceLength()-1] << "-" << change << "/" << weight << "\t";
	  	
	  }
	//  std::cerr << std::endl;
	  return kmers;
  }
};

template< typename M >
class ReadSelector
{
public:
  typedef ReadSet::ReadSetSizeType ReadSetSizeType;
  typedef Sequence::SequenceLengthType SequenceLengthType;
  class ReadTrimType {
    SequenceLengthType trimOffset;
    float score;
    std::string label;
    ReadTrimType() : trimOffset(0), score(0.0), label() {}
  };
  typedef M DataType;
  typedef typename DataType::ReadPositionWeightVector ReadPositionWeightVector;
  typedef KmerMap<DataType> KMType;
  typedef typename KMType::ConstIterator KMIterator;
  typedef typename KMType::ElementType ElementType;
  typedef KmerMap<unsigned short> KMCacheType;
  typedef std::vector< ReadTrimType > ReadTrimVector;
  typedef std::vector< ReadSetSizeType > PicksVector;

private:
  const ReadSet &_reads;
  const KMType &_map;
  ReadTrimVector _trims;
  PicksVector _picks;

public:
  ReadSelector(const ReadSet &reads, const KMType &map): _reads(reads), _map(map), _trims(), _picks() 
  { 
  	scoreAndTrimReads();
  }

  void pickAllPassingReads(float minimumScore = 0.0, SequenceLengthType minimumLength = KmerSizer::getSequenceLength()) {
  	if (_reads.getPairSize() > 0)
  	  // pick pairs of reads
  	  for(long i = 0; i < (long) _reads.getPairSize(); i++) {
  		ReadSet::Pair &pair = _reads.getPair(i);
  		if ( _reads.isValidRead(pair.read1) )
  		  if ( _trims[pair.read1].score < minimumScore || _trims[pair.read1].trimOffset < minimumLength-1)
  		     continue;
        if ( _reads.isValidRead(pair.read2) )
  		  if ( _trims[pair.read2].score < minimumScore || _trims[pair.read2].trimOffset < minimumLength-1)
  		     continue;  
  		_picks.push_back(i);		
  	  }
  	else
  	  // pick individual reads
  	  for(long i = 0; i < (long) _reads.getSize(); i++) {
  	  	if (_trims[i].score < minimumScore || _trims[i].trimOffset < minimumLength-1)
  		  continue;
  		_picks.push_back(i);
  	  }
  }
    
  void scoreAndTrimReads() {
  	_trims.resize(_reads.getSize());
  	#ifdef _USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
	#endif
  	for(long i = 0; i < (long) _reads.getSize(); i++) {
  		Read &read = _reads.getRead(i);
  		ReadTrimType &trim = _trims[i];
  		
  		KmerArray<char> kmers(read.getTwoBitSequence(), read.getLength(), true);
  		for(unsigned long j = 0; j < kmers.size(); j++) {
  			ElementType elem = _map.getElementIfExists(kmers[j]);
  			if (elem.isValid()) {
  				trim.trimOffset++;
  				trim.score += elem.getWeightedCount();
  			} else
  			    break;
  		}
  		if (trim.trimOffset > 0) {
  		  trim.trimOffset += KmerSizer::getSequenceLength() - 1;
  		  // TODO label the read
  		}
  	}
  }
  
  void pickLeastCoveringSubset() {
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

  const ReadTrimVector &getPairedPicks() const {
  	return _picks;
  }
  
  std::ostream &writePicks(std::ostream &os) const {
  	foreach( ReadSetSizeType pairIdx, _picks) {
  	  ReadSet::Pair &pair = _reads.getPair(pairIdx);
  	  if (_reads.isValidRead(pair.read1)) {
  	  	ReadTrimType &readTrim = _trims[pair.read1];
  	    os << _reads.getRead( pair.read1 ).toFastq( readTrim.trimOffset, readTrim.label );
  	  }
  	  if (_reads.isValidRead(pair.read2)) {
  	  	ReadTrimType &readTrim = _trims[pair.read2];
  	    os << _reads.getRead( pair.read2 ).toFastq( readTrim.trimOffset, readTrim.label );
  	  }  	  
  	}
  	return os;
  }
  
};



#endif

//
// $Log: Utils.h,v $
// Revision 1.27  2010-01-13 07:20:08  regan
// refactored filter
// checkpoint on read picker
//
// Revision 1.26  2010-01-13 00:25:43  regan
// use less memory for reference sequences and those without quality
//
// Revision 1.25  2010-01-11 19:14:10  regan
// minor performance enhancements
//
// Revision 1.24  2010-01-11 05:40:09  regan
// believe that FilterKnownOddities is working
//
// Revision 1.23  2010-01-08 06:25:18  regan
// refactored some code
// still working on FilterKnownOddities
//
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
