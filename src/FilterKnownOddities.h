// $Header: /repository/PI_annex/robsandbox/KoMer/src/FilterKnownOddities.h,v 1.3 2010-02-26 13:01:17 regan Exp $

#ifndef _FILTER_H
#define _FILTER_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>

#include "config.h"
#include "Kmer.h"
#include "ReadSet.h"
#include "KmerReadUtils.h"

class FilterKnownOddities {
public:
	typedef boost::unordered_map<Kmer::NumberType, std::string> DescriptionMap;
	typedef boost::unordered_map<Kmer::NumberType, unsigned long> CountMap;
	typedef std::vector<unsigned long> PatternVector;

private:
	ReadSet sequences;
	DescriptionMap descriptions;
	CountMap counts;
	Kmer::NumberType mask;
	unsigned short length;
	unsigned short maskBytes;

public:
	FilterKnownOddities(int _length = 24, int numErrors = 2) :
		length(_length) {
		std::string fasta = getFasta();
		sequences.appendFastaFile(fasta);
		// mask is low bits
		maskBytes = TwoBitSequence::fastaLengthToTwoBitLength(length);
		mask = 1;
		mask <<= (length * 2);
		mask--;
		prepareMaps(numErrors);
	}

	void prepareMaps(int numErrors) {
		unsigned long oldKmerLength = KmerSizer::getSequenceLength();
		KmerSizer::set(length);

		TEMP_KMER(reverse);
		for (unsigned int i = 0; i < sequences.getSize(); i++) {
			// TODO refine search around byte-boundary problem

			Read &read = sequences.getRead(i);
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
			patterns.reserve(descriptions.size());
			for (DescriptionMap::iterator it = descriptions.begin(); it
					!= descriptions.end(); it++)
				patterns.push_back(it->first);

			for (PatternVector::iterator it = patterns.begin(); it
					!= patterns.end(); it++) {
				KmerWeights kmers = KmerWeights::permuteBases(
						*((Kmer *) &(*it)));
				for (unsigned int j = 0; j < kmers.size(); j++) {
					kmers[j].buildReverseComplement(reverse);

					_setMaps(kmers[j].toNumber(), descriptions[*it] + "-error"
							+ boost::lexical_cast<std::string>(error) + "/"
							+ boost::lexical_cast<std::string>(j));

					_setMaps(reverse.toNumber(), descriptions[*it]
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
		std::cerr << "Filter size " << descriptions.size() << " "
				<< std::setbase(16) << mask << std::setbase(10) << std::endl;
		//foreach( DescriptionMap::value_type _desc, descriptions) {
		//  tmp[0] = _desc.first;
		//  std::cerr << "Filter "  << ((Kmer *)&tmp[0])->toFasta() << " " << _desc.first << " " << _desc.second << std::endl;
		//}

		KmerSizer::set(oldKmerLength);
	}

	void _setMapsBitShifted(Kmer::NumberType chunk, std::string description) {
		for (int shift = 0; shift < 4; shift++) {
			Kmer::NumberType key = (chunk << shift) & mask;
			_setMaps(key, description);
		}
	}
	void _setMaps(Kmer::NumberType key, std::string description) {
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
		for (long i = 0; i < (long) reads.getSize(); i++) {
			Read &read = reads.getRead(i);
			bool wasAffected = false;

			//std::cerr << "Checking " << read.getName() << "\t" << read.getFasta() << std::endl;

			unsigned long lastByte = read.getTwoBitEncodingSequenceLength() - 1;
			if (lastByte < maskBytes)
				continue;

			TwoBitEncoding *ptr = read.getTwoBitSequence() + lastByte;
			CountMap::iterator it;

			Kmer::NumberType chunk = 0;
			for (unsigned long loop = lastByte + maskBytes; loop >= maskBytes; loop--) {
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
					for (int baseShift = 0; baseShift < 4; baseShift++) {
						if (loop > lastByte && baseShift != 0)
							break;
						key = (chunk >> (baseShift * 2)) & mask;

						//std::cerr << bytes << "\t" << std::setbase(16) << key << "\t" << chunk << std::setbase(10) << std::endl;

						it = counts.find(key);
						if (it != counts.end()) {
							counts[key]++;
							unsigned long offset = bytes * 4 + baseShift;
							read.zeroQuals(offset, length);
							if (!wasAffected) {
								wasAffected = true;
								affectedCount++;
								if (0)
									std::cerr << "FilterMatch to "
											<< read.getName() << " against "
											<< descriptions[key] << " at "
											<< offset << "\t"
											<< read.getFasta() << "\t"
											<< std::setbase(16) << key << "\t"
											<< chunk << std::setbase(10)
											<< std::endl;
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

};

#endif
