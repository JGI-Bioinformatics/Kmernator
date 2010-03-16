// $Header: /repository/PI_annex/robsandbox/KoMer/src/FilterKnownOddities.h,v 1.14 2010-03-16 18:57:55 regan Exp $

#ifndef _FILTER_H
#define _FILTER_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>

#include "config.h"
#include "Kmer.h"
#include "ReadSet.h"
#include "KmerReadUtils.h"
#include "KmerSpectrum.h"

class FilterKnownOddities {
public:
	typedef Kmer::NumberType NumberType;
	typedef boost::unordered_map<NumberType, std::string> DescriptionMap;
	typedef boost::array< DescriptionMap, 4 > DescriptionMapArray;
	typedef boost::unordered_map<NumberType, unsigned long> CountMap;
	typedef boost::array< CountMap, 4 > CountMapArray;
	typedef boost::array< NumberType, 4 > BitShiftNumberArray;
	typedef std::vector<unsigned long> PatternVector;


private:
	ReadSet sequences;
	BitShiftNumberArray masks;
	DescriptionMapArray descriptions;
	CountMapArray counts;
	unsigned short length;
	unsigned short twoBitLength;
	unsigned short maskBytes;
	std::string maskBases;
	std::string maskAs;

public:
	FilterKnownOddities(int _length = 26, int numErrors = 2) :
		length(_length) {
		if (length > 28) {
			throw std::invalid_argument("FilterKnownOddities must use 7 bytes or less");
		}
		twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(length);
		// T is 11, A is 00, so mask is all T's surrounded by A's
		maskBases = std::string("TTTTTTTTTTTTTTTTTTTTTTTTTTTT").substr(0,length);
		maskAs = std::string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		std::string fasta = getFasta();
		sequences.appendFastaFile(fasta);

		maskBytes = TwoBitSequence::fastaLengthToTwoBitLength(length);
		NumberType mask;
		TwoBitSequence::compressSequence(        maskBases + maskAs.substr(0,32-length), (TwoBitEncoding *) &mask);
		setBitShiftNumbers( masks, mask);

		if (Options::getDebug() > 3) {
			std::cerr
			<< TwoBitSequence::getFasta( (TwoBitEncoding *) &masks[0], 32) << "\t"
			<< TwoBitSequence::getFasta( (TwoBitEncoding *) &masks[1], 32) << "\t"
			<< TwoBitSequence::getFasta( (TwoBitEncoding *) &masks[2], 32) << "\t"
			<< TwoBitSequence::getFasta( (TwoBitEncoding *) &masks[3], 32) << "\t" << std::endl;
		}
		prepareMaps(numErrors);
	}

	void setBitShiftNumbers(BitShiftNumberArray &nums, std::string fasta, int len) {
		TwoBitSequence::compressSequence(        fasta.substr(0,len) + maskAs.substr(0,32-len), (TwoBitEncoding *) &nums[0]);
		TwoBitSequence::compressSequence("A"   + fasta.substr(0,len) + maskAs.substr(0,31-len), (TwoBitEncoding *) &nums[1]);
		TwoBitSequence::compressSequence("AA"  + fasta.substr(0,len) + maskAs.substr(0,30-len), (TwoBitEncoding *) &nums[2]);
		TwoBitSequence::compressSequence("AAA" + fasta.substr(0,len) + maskAs.substr(0,29-len), (TwoBitEncoding *) &nums[3]);
	}
	void setBitShiftNumbers(BitShiftNumberArray &nums, NumberType num) {
		NumberType right4 = num << 8;
		NumberType tmp;
		nums[0] = num;
		for(int i = 1; i < 4 ; i++) {
			TwoBitSequence::shiftLeft(&right4, &tmp, twoBitLength+1, i);
			nums[4-i] = tmp;
		}
	}

	void prepareMaps(int numErrors) {
		unsigned long oldKmerLength = KmerSizer::getSequenceLength();
		KmerSizer::set(length);

		TEMP_KMER(reverse);
		for (unsigned int i = 0; i < sequences.getSize(); i++) {

			const Read read = sequences.getRead(i);
			KmerWeights kmers = KmerReadUtils::buildWeightedKmers(read);
			for (unsigned int j = 0; j < kmers.size(); j++) {

				_setMaps(kmers[j].toNumber(), read.getName() + "@"
						+ boost::lexical_cast<std::string>(j));
				kmers[j].buildReverseComplement(reverse);
				_setMaps(reverse.toNumber(), read.getName() + "-reverse@"
						+ boost::lexical_cast<std::string>(j));
			}

		}
		for (int error = 0; error < numErrors; error++) {
			// build a vector of all patterns

			PatternVector patterns;
			patterns.reserve(descriptions[0].size());
			for (DescriptionMap::iterator it = descriptions[0].begin(); it
					!= descriptions[0].end(); it++)
				patterns.push_back(it->first);

			for (PatternVector::iterator it = patterns.begin(); it
					!= patterns.end(); it++) {
				KmerWeights kmers = KmerWeights::permuteBases(
						*((Kmer *) &(*it)));
				for (unsigned int j = 0; j < kmers.size(); j++) {
					kmers[j].buildReverseComplement(reverse);

					_setMaps(kmers[j].toNumber(), descriptions[0][*it] + "-error"
							+ boost::lexical_cast<std::string>(error) + "/"
							+ boost::lexical_cast<std::string>(j));

					_setMaps(reverse.toNumber(), descriptions[0][*it]
							+ "-reverror" + boost::lexical_cast<std::string>(
							error) + "/" + boost::lexical_cast<std::string>(j));

				}
			}

			// permute bases
			// store new patterns (if new)
		}
		Kmer::NumberType tmp[10];
		for (int i = 0; i < 10; i++)
			tmp[i] = 0;
		std::cerr << "Filter size " << descriptions[0].size()+descriptions[1].size()+descriptions[2].size()+descriptions[3].size() << " "
				<< std::setbase(16) << masks[0] << std::setbase(10) << std::endl;
		//foreach( DescriptionMap::value_type _desc, descriptions) {
		//  tmp[0] = _desc.first;
		//  std::cerr << "Filter "  << ((Kmer *)&tmp[0])->toFasta() << " " << _desc.first << " " << _desc.second << std::endl;
		//}

		KmerSizer::set(oldKmerLength);
	}

	// NOT USED!
	void _setMapsBitShifted(Kmer::NumberType chunk, std::string description) {
		for (int shift = 0; shift < 4; shift++) {
			Kmer::NumberType key = (chunk << shift) & masks[0];
			_setMaps(key, description);
		}
	}
	void _setMaps(Kmer::NumberType key, std::string description) {
		BitShiftNumberArray numbers;
		setBitShiftNumbers(numbers, key);
		for(int i = 0; i < 4 ; i++) {
		  if (descriptions[i].find(numbers[i]) == descriptions[i].end()) {
			descriptions[i][numbers[i]] = description;
			counts[i][numbers[i]] = 0;
		  }
		}
	}

	unsigned long applyFilter(ReadSet &reads) {
		unsigned long affectedCount = 0;
#ifdef _USE_OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:affectedCount)
#endif
		for (long readIdx = 0; readIdx < (long) reads.getSize(); readIdx++) {
			Read &read = reads.getRead(readIdx);
			bool wasAffected = false;

			//std::cerr << "Checking " << read.getName() << "\t" << read.getFasta() << std::endl;

			unsigned long lastByte = read.getTwoBitEncodingSequenceLength() - 1;
			if (lastByte < maskBytes)
				continue;


			CountMap::iterator it;

			Kmer::NumberType chunk = 0;
			for (unsigned long loop = lastByte + maskBytes; loop >= maskBytes; loop--) {
				// shift one byte
				chunk <<= 8;
				TwoBitEncoding *ptr = read.getTwoBitSequence() + loop - maskBytes;
				TwoBitEncoding newByte = *ptr;
				chunk |= newByte;
				unsigned long bytes = loop - maskBytes;
				if (Options::getDebug() >= 3){
					unsigned long maskedChunk = chunk & masks[0];
			        std::cerr << read.getName() << " " << bytes << " "
				        << TwoBitSequence::getFasta( (TwoBitEncoding *) &chunk, 32) << "\t"
				        << TwoBitSequence::getFasta( (TwoBitEncoding *) &maskedChunk, 32) << "\t"
				        << TwoBitSequence::getFasta( (TwoBitEncoding *) &newByte, 4) << "\t"
				        << (size_t) &(*ptr) << "\t" << (size_t) &( *(read.getTwoBitSequence()) ) << "\t"
				        << std::endl;
				}
				if (loop <= lastByte + 1) {
					Kmer::NumberType key;
					for (int baseShift = 0; baseShift < 4; baseShift++) {
						if (loop > lastByte && baseShift != 0)
							break;
						key = chunk & masks[baseShift];

						if (Options::getDebug() >= 3){
					      std::cerr << bytes << "\t" << baseShift << "\t" << TwoBitSequence::getFasta( (TwoBitEncoding *) &key, 32)
					      << "\t" << TwoBitSequence::getFasta( (TwoBitEncoding *) &chunk, 32)
					      << "\t" << std::setbase(16) << key << "\t" << chunk << std::setbase(10) << std::endl;
						}

						it = counts[baseShift].find(key);
						if (it != counts[baseShift].end()) {
							counts[baseShift][key]++;
							unsigned long offset = bytes * 4 + baseShift;
							read.markupBases(offset, length, 'X');
							if (Options::getDebug()>=2) {
									std::cerr << "FilterMatch to "
									        << readIdx << " "
											<< read.getName() << " against "
											<< TwoBitSequence::getFasta( (TwoBitEncoding *) &key, 32) << ":"
											<< offset << " (" << bytes*4 << "+" << baseShift << ")\t"
											<< descriptions[baseShift][key] << " at "
											<< read.getFasta() << "\t"
											<< std::setbase(16) << key << "\t"
											<< chunk << std::setbase(10)
											<< std::endl;
							}
							if (!wasAffected) {
						        wasAffected = true;
								affectedCount++;
						    }
						}
					}
				}

			}
		}
		return affectedCount;
	}

	const ReadSet &getSequences() {
		return sequences;
	}
	std::string getFasta() {
		// TODO logic to check Options
		std::stringstream ss;
		ss << ">PrimerDimer" << std::endl;
		ss
				<< "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"
				<< std::endl;
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

	typedef KmerSpectrum<TrackingDataWithAllReads, TrackingDataWithAllReads, TrackingDataWithAllReads> KS;
	typedef TrackingData::ReadPositionWeightVector RPW;
	typedef RPW::iterator RPWIterator;
	typedef ReadSet::Pair Pair;
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;

	static unsigned long filterIdenticalFragmentPairs(ReadSet &reads, unsigned char sequenceLength = 16, unsigned int cutoffThreshold = 2) {
	  unsigned long affectedCount = 0;

	  // select the number of bytes from each pair to scan
	  unsigned char bytes = sequenceLength / 4;
	  if (bytes == 0) {
			bytes = 1;
	  }
	  SequenceLengthType oldKmerSize = KmerSizer::getSequenceLength();
	  KmerSizer::set(bytes * 4 * 2);
	  {

		// TODO parallelize - build one KS per thread, then merge, skipping singletons
        // no need to include quality scores

		// build the kmer spectrum with the concatenated prefixes
		// from each read in the pair
		int numThreads = 1;
#ifdef _USE_OPENMP
		numThreads = omp_get_max_threads();
#endif
		KS::Vector ksv(numThreads);
		KmerWeights::Vector tmpKmerv(numThreads);
		for(int i = 0; i < numThreads; i++) {
		  ksv[i] = KS(reads.getPairSize() / 64 / numThreads, false);
		  tmpKmerv[i].resize(1);
		  tmpKmerv[i].valueAt(0) = 1.0;
		}
		ReadSetSizeType pairSize =  reads.getPairSize();
#pragma omp parallel for
		for(long pairIdx = 0; pairIdx < pairSize; pairIdx++) {
			Pair &pair = reads.getPair(pairIdx);
			int threadNum = 0;
#ifdef _USE_OPENMP
			threadNum = omp_get_thread_num();
#endif
			if (reads.isValidRead(pair.read1) && reads.isValidRead(pair.read2)) {
			  const Read read1 = reads.getRead(pair.read1);
			  const Read read2 = reads.getRead(pair.read2);
			  if ((read1.getMarkupBasesCount() == 0 || read1.getMarkups()[0].second > sequenceLength)
			     && (read2.getMarkupBasesCount() == 0 || read2.getMarkups()[0].second > sequenceLength)) {

				  // create read1 + the reverse complement of read2
				  // when it is represented as a kmer, the leastcomplement will be stored
				  // and properly account for identical fragment pairs

				  memcpy(tmpKmerv[threadNum][0].getTwoBitSequence()        , read1.getTwoBitSequence(), bytes);
			      TwoBitSequence::reverseComplement( read2.getTwoBitSequence(), tmpKmerv[threadNum][0].getTwoBitSequence() + bytes, bytes*4);

				  // store the pairIdx (not readIdx)
			      ksv[threadNum].append(tmpKmerv[threadNum], pairIdx);

			  }
			}
		}
		if (Options::getDebug() > 3) {
		  for (int i = 0; i < numThreads; i++) {
			  std::cerr << "spectrum " << i << std::endl;
			  ksv[i].printHistograms();
		  }
		} else {
			std::cerr << "merging duplicate fragment spectrums" << std::endl;
		}

		KS::mergeVector(ksv, 2);
		KS &ks = ksv[0];

		// analyze the spectrum
		ks.printHistograms();
		// TODO parallelize (partition by bucket range)

		ReadSet newReads;
		for(KS::WeakIterator it = ks.weak.begin(); it != ks.weak.end(); it++) {
		    if (it->value().getCount() >= cutoffThreshold) {
		    	RPW rpw = it->value().getEachInstance();

		    	ReadSet tmpReadSet1;
		    	ReadSet tmpReadSet2;

		    	// std::stringstream out1, out2, out1a, out2a;

		    	for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {

		    		// iterator readId is actually the pairIdx built above
		    		ReadSetSizeType pairIdx = rpwit->readId;

		    		Pair &pair = reads.getPair(pairIdx);
		    		const Read &read1 = reads.getRead(pair.read1);
		    		const Read &read2 = reads.getRead(pair.read2);

		    		tmpReadSet1.append( read1 );
		    		tmpReadSet2.append( read2 );
		    		//out1 << read1.getFasta() << std::endl;
		    		//out1a << read1.getQuals() << std::endl;
		    		//out2 << read2.getFasta() << std::endl;
		    		//out2a << read2.getQuals() << std::endl;
		    	}

		    	Read consensus1 = tmpReadSet1.getConsensusRead();
		    	Read consensus2 = tmpReadSet2.getConsensusRead();
		    	newReads.append(consensus1);
		    	newReads.append(consensus2);


		    	//	std::cerr << consensus1.getName() << std::endl
		    	//	<< out1.str() << consensus1.getFasta() << std::endl << out1a.str()
		    	//	<< consensus1.getQuals() << std::endl;
		    	//	std::cerr << consensus2.getName() << std::endl
		    	//	<< out2.str() << consensus2.getFasta() << std::endl << out2a.str()
		    	//	<< consensus2.getQuals() << std::endl;


		    	//ReadSetSizeType readIdx = tmpReadSet.getCentroidRead();
		    	//ReadSetSizeType count = 0;

		    	for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {
		    	    //if (count++ == readIdx)
		    	    //	continue;
		    	    ReadSetSizeType pairIdx = rpwit->readId;

		    	    Pair &pair = reads.getPair(pairIdx);
		    	    Read &read1 = reads.getRead(pair.read1);
		    	    Read &read2 = reads.getRead(pair.read2);
		    		read1.markupBases(0, -1, 'X');
		    		read2.markupBases(0, -1, 'X');
		    		affectedCount += 2;
		    	}

		    }
		}
		std::cerr << "Built " << newReads.getSize() << " new consensus reads" << std::endl;
		newReads.identifyPairs();
		reads.append(newReads);
	  }
      KmerSizer::set(oldKmerSize);
	  return affectedCount;
	}

};

#endif
