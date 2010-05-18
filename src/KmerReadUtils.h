// $Header: /repository/PI_annex/robsandbox/KoMer/src/KmerReadUtils.h,v 1.9 2010-05-18 20:50:24 regan Exp $
//

#ifndef _KMER_READ_UTILS_H
#define _KMER_READ_UTILS_H

#include "config.h"

class KmerReadUtils {
public:
	static KmerWeights buildWeightedKmers(const Read &read, bool leastComplement = false, bool leastComplementForNegativeWeight = false) {
		SequenceLengthType readLength = read.getLength();
		bool bools[readLength];
		KmerWeights kmers(read.getTwoBitSequence(), readLength, leastComplement, bools);
		std::string quals = read.getQuals();
		size_t markupIdx = 0;

		BaseLocationVectorType markups = read.getMarkups();
		double weight = 0.0;
		double change = 0.0;
		bool isRef = false;
		if (quals.length() > 0 && quals[0] == Read::REF_QUAL)
			isRef = true;

		SequenceLengthType size = (SequenceLengthType) kmers.size();
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
				//  	    std::cerr << "markupAt " << markups[markupIdx].first << " " << markups[markupIdx].second << "\t";
			}
			kmers.valueAt(i) = leastComplementForNegativeWeight && bools[i] ? weight : (0.0-weight);
			//  	std::cerr << i << ":" << std::fixed << std::setprecision(3) << quals[i+KmerSizer::getSequenceLength()-1] << "-" << change << "/" << weight << "\t";

		}
		//  std::cerr << std::endl;
		if (Options::getDebug() > 3) {
		  for(Kmer::IndexType i = 0 ; i < kmers.size(); i++) {
			std::cerr << i << " " << kmers.valueAt(i) << " " << kmers[i].toFasta() << std::endl;
		  }
		}
		return kmers;
	}
};

#endif


// $Log: KmerReadUtils.h,v $
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
