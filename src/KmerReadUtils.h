// $Header: /repository/PI_annex/robsandbox/KoMer/src/KmerReadUtils.h,v 1.5 2010-03-14 16:55:13 regan Exp $
//

#ifndef _KMER_READ_UTILS_H
#define _KMER_READ_UTILS_H

class KmerReadUtils {
public:
	static KmerWeights buildWeightedKmers(const Read &read, bool leastComplement =
			false) {
		KmerWeights kmers(read.getTwoBitSequence(), read.getLength(),
				leastComplement);
		std::string quals = read.getQuals();
		SequenceLengthType markupIdx = 0;

		BaseLocationVectorType markups = read.getMarkups();
		double weight = 0.0;
		double change = 0.0;
		bool isRef = false;
		if (quals.length() > 0 && quals[0] == Read::REF_QUAL)
			isRef = true;

		for (SequenceLengthType i = 0; i < kmers.size(); i++) {
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
			kmers.valueAt(i) = weight;
			//  	std::cerr << i << ":" << std::fixed << std::setprecision(3) << quals[i+KmerSizer::getSequenceLength()-1] << "-" << change << "/" << weight << "\t";

		}
		//  std::cerr << std::endl;
		if (Options::getDebug() > 3) {
		  for(unsigned int i = 0 ; i < kmers.size(); i++) {
			std::cerr << i << " " << kmers.valueAt(i) << " " << kmers[i].toFasta() << std::endl;
		  }
		}
		return kmers;
	}
};

#endif

