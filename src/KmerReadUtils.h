//
// Kmernator/src/KmerReadUtils.h
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

#ifndef _KMER_READ_UTILS_H
#define _KMER_READ_UTILS_H

#include "config.h"
#include "Sequence.h"
#include "KmerTrackingData.h"
#include "Kmer.h"

class KmerReadUtils {
public:
	static KmerWeightedExtensions buildWeightedKmers(const Read &read, bool leastComplement = false, bool leastComplementForNegativeWeight = false) {
		SequenceLengthType readLength = read.getLength();
		bool needMalloc = readLength > MAX_STACK_SIZE;
		bool _bools[ needMalloc ? 0 : readLength ];
		bool *bools = _bools;
		if (needMalloc) {
			bools = new bool[readLength];
		}
		std::string fasta = read.getFasta();
		int kmerLen = KmerSizer::getSequenceLength();

		KmerWeightedExtensions kmers(read.getTwoBitSequence(), readLength, leastComplement, bools);
		std::string quals = read.getQuals();
		size_t markupIdx = 0;

		BaseLocationVectorType markups = read.getMarkups();
		double weight = 0.0;
		double change = 0.0;
		bool isRef = false;
		if (quals.length() > 0 && quals[0] == Read::REF_QUAL)
			isRef = true;

		SequenceLengthType size = (SequenceLengthType) kmers.size();
		assert(size < readLength);
		Extension left = Extension('X', ExtensionTracking::getMinQuality()), right;
		for (SequenceLengthType i = 0; i < size; i++) {
			if (isRef) {
				weight = 1.0;
			} else if (i % 1024 == 0 || weight == 0.0) {
				weight = 1.0;
				for (SequenceLengthType j = 0; j
						< KmerSizer::getSequenceLength(); j++)
					weight
							*= Read::qualityToProbability[(unsigned char) quals[i
									+ j]];
			} else {
				change = Read::qualityToProbability[(unsigned char) quals[i
						+ KmerSizer::getSequenceLength() - 1]]
						/ Read::qualityToProbability[(unsigned char) quals[i
								- 1]];
				weight *= change;
			}
			while (markupIdx < markups.size() && markups[markupIdx].second < i)
				markupIdx++;
			if (markupIdx < markups.size() && markups[markupIdx].second < i
					+ KmerSizer::getSequenceLength()) {
				weight = 0.0;
			}

			// set the weight
			kmers.valueAt(i).setWeight( leastComplementForNegativeWeight && bools[i] ? weight : (0.0-weight) );

			SequenceLengthType rightBase = i + kmerLen;
			if (rightBase < fasta.length())
				right = Extension(fasta[rightBase], quals[rightBase] - Read::FASTQ_START_CHAR);
			else
				right = Extension('X', ExtensionTracking::getMinQuality());

			// set the extensions
			if (bools[i])
				kmers.valueAt(i).setExtensions(left, right);
			else // kmer is reverse complement, so reverse extensions
				kmers.valueAt(i).setExtensions(right.getReverseComplement(), left.getReverseComplement());

			left = Extension(fasta[i], quals[i] - Read::FASTQ_START_CHAR);
		}
		if (Log::isDebug(5)) {
			ostream &debug = Log::Debug() << "KmerWeights: idx valueAt toFasta" << std::endl;;
		    for(Kmer::IndexType i = 0 ; i < kmers.size(); i++) {
			  debug << i << " " << kmers.valueAt(i).getWeight() << " " << kmers[i].toFasta() << std::endl;
		    }
		}
		if (needMalloc) {
			delete [] bools;
		}
		return kmers;
	}
};

#endif


// $Log: KmerReadUtils.h,v $
// Revision 1.10  2010-08-18 17:50:40  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.9.10.1  2010-07-19 18:43:41  regan
// fixed stack vs malloc problem
//
// Revision 1.9  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.8.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.8  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
//
