// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSelector.h,v 1.1 2010-01-13 23:48:51 regan Exp $
//

#ifndef _READ_SELECTOR_H
#define _READ_SELECTOR_H

#include <iostream>
#include <cstdlib>
#include <cstring>

#include <vector>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "Sequence.h"
#include "ReadSet.h"
#include "Kmer.h"

template< typename M >
class ReadSelector
{
public:
  typedef Sequence::SequenceLengthType SequenceLengthType;
  typedef ReadSet::ReadSetSizeType ReadSetSizeType;
  class ReadTrimType {
  public:
    SequenceLengthType trimOffset;
    float score;
    std::string label;
    bool isAvailable;
    ReadTrimType() : trimOffset(0), score(0.0), label(), isAvailable(true){}
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

  bool pickIfNew(ReadSetSizeType readIdx) {
  	if ( _reads.isValidRead(readIdx) && _trims[readIdx].isAvailable ) {
  	  _picks.push_back(readIdx);
  	  _trims[readIdx].isAvailable = false;
  	  // TODO set label
  	  return true;
  	} else
  	  return false;	  
  }
  
  bool isPassingRead(ReadSetSizeType readIdx) {
  	if (! _reads.isValidRead(readIdx) )
  	  return false;
  }
  bool isPassingRead(ReadSetSizeType readIdx, float minimumScore, SequenceLengthType minimumLength) {
  	if (! isPassingRead(readIdx) )
  	  return false;
  	ReadTrimType &trim = _trims[readIdx];
  	return trim.isAvailable && trim.score >= minimumScore && trim.trimOffset > minimumLength;
  }
  
  int pickAllPassingReads(float minimumScore = 0.0, SequenceLengthType minimumLength = KmerSizer::getSequenceLength()) {
  	int picked = 0;
  	for(long i = 0; i < (long) _reads.getSize(); i++) {
  	  	if (isPassingRead(i,minimumScore,minimumLength) && pickIfNew(i))
  		  picked++;
  	  }
  	return picked;
  }
  int pickAllPassingPairs(float minimumScore = 0.0, SequenceLengthType minimumLength = KmerSizer::getSequenceLength()) {
  	  int picked = 0;
  	  for(long i = 0; i < (long) _reads.getPairSize(); i++) {
  		const ReadSet::Pair &pair = _reads.getPair(i);
  		if ( isPassingRead(pair.read1,minimumScore,minimumLength) && isPassingRead(pair.read2, minimumScore, minimumLength) )
  		  if (pickIfNew(pair.read1))
  		    picked++;
  		  if (pickIfNew(pair.read2))
  		    picked++;		
  	  }
  	  return picked;  	  
  }
    
  void scoreAndTrimReads() {
  	_trims.resize(_reads.getSize());
  	#ifdef _USE_OPENMP
    #pragma omp parallel for schedule(dynamic)
	#endif
  	for(long i = 0; i < (long) _reads.getSize(); i++) {
  		const Read &read = _reads.getRead(i);
  		ReadTrimType &trim = _trims[i];
  		
  		KmerArray<char> kmers(read.getTwoBitSequence(), read.getLength(), true);
  		for(unsigned long j = 0; j < kmers.size(); j++) {
  			const ElementType elem = _map.getElementIfExists(kmers[j]);
  			if (elem.isValid()) {
  				trim.trimOffset++;
  				trim.score += elem.value().getWeightedCount();
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

  const PicksVector &getPickedReads() const {
  	return _picks;
  }
  
  std::ostream &writePicks(std::ostream &os) const {
  	foreach( ReadSetSizeType readIdx, _picks) {
  	  ReadTrimType &trim = _trims[readIdx];
  	  os << _reads.getRead( readIdx ).toFastq( trim.trimOffset, trim.label );  	  	  
  	}
  	return os;
  }
};

#endif

