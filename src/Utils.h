// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.25 2010-01-11 19:14:10 regan Exp $
//

#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <iostream>
#include <iomanip>
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
#include <boost/lexical_cast.hpp>


using namespace boost::accumulators;

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
};

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
public:
  typedef boost::unordered_map< Kmer::NumberType, std::string > DescriptionMap;
  typedef boost::unordered_map< Kmer::NumberType, unsigned long > CountMap;
  typedef std::vector< unsigned long > PatternVector;

private:
  ReadSet sequences;
  DescriptionMap descriptions;
  CountMap counts;
  Kmer::NumberType mask;
  unsigned short length;
  unsigned short maskBytes;
  
public:
  FilterKnownOddities(int _length = 24, int numErrors = 2) : length(_length) {
  	std::string fasta = getFasta();
  	sequences.appendFastaFile( fasta );
  	// mask is low bits
  	maskBytes = TwoBitSequence::fastaLengthToTwoBitLength(length);
  	mask = 1;
  	mask <<= (length*2);
  	mask--;
  	prepareMaps(numErrors);
  }
  
  void prepareMaps(int numErrors) {
  	unsigned long oldKmerLength = KmerSizer::getSequenceLength();
  	KmerSizer::set(length);
  	
  	TEMP_KMER(reverse);
  	for(unsigned int i = 0 ; i < sequences.getSize() ; i++) {
  		// TODO refine search around byte-boundary problem
  		
  		Read &read = sequences.getRead(i);
  		KmerWeights kmers = KmerReadUtils::buildWeightedKmers(read);
  		for(unsigned int j = 0; j < kmers.size(); j++) {
  			
  			_setMaps( kmers[j].toNumber() , read.getName() + "@" + boost::lexical_cast<std::string>(j));  			
   			kmers[j].buildReverseComplement(reverse);  	        
  		  	_setMaps( reverse.toNumber()  , read.getName() + "-reverse@" + boost::lexical_cast<std::string>(j));
   		}
  		
  	}
  	for (int error = 0 ; error < numErrors; error++) {
  		// build a vector of all patterns
  		
  		PatternVector patterns;
  		patterns.reserve(descriptions.size());
  		for(DescriptionMap::iterator it = descriptions.begin(); it != descriptions.end(); it++)
  		  patterns.push_back(it->first);
  		
  		for(PatternVector::iterator it = patterns.begin(); it != patterns.end(); it++) {
  		  KmerWeights kmers = KmerWeights::permuteBases( *( (Kmer *) &(*it) ) );
  		  for(unsigned int j = 0; j < kmers.size(); j++) {
  		  	kmers[j].buildReverseComplement(reverse);
  		  	
  		  	_setMaps( kmers[j].toNumber(), descriptions[*it] + "-error" + boost::lexical_cast<std::string>(error) + "/" + boost::lexical_cast<std::string>(j));
  			
  	        _setMaps( reverse.toNumber() , descriptions[*it] + "-reverror" + boost::lexical_cast<std::string>(error) + "/" + boost::lexical_cast<std::string>(j));
  	        
  		  }
  		}
  		
  		// permute bases
  		// store new patterns (if new)
  	}
  	Kmer::NumberType tmp[10];
  	for(int i=0;i<10;i++)
  	  tmp[i] = 0;
  	std::cerr << "Filter size " << descriptions.size() << " " << std::setbase(16) << mask << std::setbase(10) << std::endl;
  	//foreach( DescriptionMap::value_type _desc, descriptions) {
    //  tmp[0] = _desc.first;
  	//  std::cerr << "Filter "  << ((Kmer *)&tmp[0])->toFasta() << " " << _desc.first << " " << _desc.second << std::endl;
  	//}
  	
  	KmerSizer::set(oldKmerLength);
  }
  
  void _setMapsBitShifted(Kmer::NumberType chunk, std::string description) {
  	for(int shift=0; shift<4; shift++) {
      Kmer::NumberType key = (chunk<<shift) & mask;
      _setMaps(key, description);
    }  
  }
  void _setMaps( Kmer::NumberType key, std::string description) {
  	if (descriptions.find(key) == descriptions.end()) {
      	descriptions[key] = description;
      	counts[key] = 0;
    }
  }
  unsigned long applyFilter(ReadSet &reads) {
  	unsigned long affectedCount = 0;
  	#ifdef _USE_OPENMP
    #pragma omp parallel for schedule(dynamic) reduction(+:affectedCount)
	#endif
  	for(long i = 0; i < (long) reads.getSize(); i++) {
  		Read &read = reads.getRead(i);
  		bool wasAffected = false;
  		
  		//std::cerr << "Checking " << read.getName() << "\t" << read.getFasta() << std::endl;
  		
  		unsigned long lastByte = read.getTwoBitEncodingSequenceLength() - 1;
  		if (lastByte < maskBytes) 
  		  continue;
  		  
  		TwoBitEncoding *ptr = read.getTwoBitSequence() + lastByte;
  		CountMap::iterator it;

  		Kmer::NumberType chunk = 0;
  		for(unsigned long loop = lastByte + maskBytes; loop >= maskBytes ; loop--) {
  			// shift one byte
  			chunk <<= 8;
  			chunk |= *(ptr--);
  			unsigned long bytes = loop - maskBytes;
  			//{
  			//	unsigned long maskedChunk = chunk & mask;
  			//    std::cerr << read.getName() << " " << bytes << " " 
  			//        <<  TwoBitSequence::getFasta( (TwoBitEncoding *) &chunk, 32) << "\t" 
  			//        << TwoBitSequence::getFasta( (TwoBitEncoding *) &maskedChunk, 32)<< std::endl;
  			//}		
  			if (loop <= lastByte + 1) {  
  			  Kmer::NumberType key;
  			  for(int baseShift = 0 ; baseShift <4 ; baseShift++) {
                if (loop > lastByte && baseShift != 0)
                  break;
  			  	key = (chunk>> (baseShift*2)) & mask;
  
   			    //std::cerr << bytes << "\t" << std::setbase(16) << key << "\t" << chunk << std::setbase(10) << std::endl;
  
  			    it = counts.find( key );
  			    if (it != counts.end()) {
  				  counts[key]++;
  				  unsigned long offset = bytes * 4 + baseShift;
  				  read.zeroQuals( offset, length);
  				  if (!wasAffected) {
  				    wasAffected = true;
  			  	    affectedCount++;
  				    if (0)
  				      std::cerr << "FilterMatch to " << read.getName() 
  				        << " against " << descriptions[key] << " at " << offset << "\t" 
  				        << read.getFasta() << "\t" << std::setbase(16) << key << "\t" << chunk << std::setbase(10) << std::endl;
  				  }
  			    }
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
  	ss << ">PrimerDimer" << std::endl;
  	ss << "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG" << std::endl;
  	//ss << ">PE-Adaptor-P" << std::endl;
  	//ss << "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG" << std::endl;
  	//ss << ">PE-PCR-PRIMER1" << std::endl;
  	//ss << "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" << std::endl;
  	//ss << ">PE-PCR-PRIMER2" << std::endl;
  	//ss << "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT" << std::endl;
  	//ss << ">SE-Adaptor-P" << std::endl;
  	//ss << "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG" << std::endl;
   	//ss << ">SE-Primer2" << std::endl;
  	//ss << "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT" << std::endl;
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
