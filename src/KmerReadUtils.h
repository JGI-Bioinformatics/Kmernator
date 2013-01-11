//
// Kmernator/src/KmerReadUtils.h
//
// Author: Rob Egan
//
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/

#ifndef _KMER_READ_UTILS_H
#define _KMER_READ_UTILS_H

#include "config.h"
#include "Sequence.h"
#include "KmerTrackingData.h"
#include "Kmer.h"

class KmerReadUtils {
private:
	KmerWeightedExtensions kmers;
public:
	KmerReadUtils() {
		LOG_DEBUG_OPTIONAL(2, true, "KmerReadUtils()" << &kmers);
	}
	~KmerReadUtils() {}
	KmerWeightedExtensions &buildWeightedKmers(const Read &read, bool leastComplement = false, bool leastComplementForNegativeWeight = false) {
		if (read.isDiscarded()) {
			kmers.resize(0);
			return kmers;
		}

		SequenceLengthType readLength = read.getLength();
		STACK_ALLOC(bool, bools, readLength);
		std::string fasta = read.getFastaNoMarkup();
		int kmerLen = KmerSizer::getSequenceLength();

		kmers.build(read.getTwoBitSequence(), readLength, leastComplement, bools);
		std::string quals = read.getQuals();
		size_t markupIdx = 0;

		BaseLocationVectorType markups = read.getMarkups();
		double weight = 0.0;
		double change = 0.0;
		bool isRef = false;
		if (quals.length() > 0 && quals[0] == Read::REF_QUAL)
			isRef = true;

		SequenceLengthType size = (SequenceLengthType) kmers.size();
		assert(size == 0 || size < readLength);
		Extension left = Extension('X', ExtensionTracking::getMinQuality()), right;
		for (SequenceLengthType i = 0; i < size; i++) {
			if (isRef) {
				weight = 1.0;
			} else if (i % 1024 == 0 || weight == 0.0) {
				weight = 1.0;
				for (SequenceLengthType j = 0; j
				< KmerSizer::getSequenceLength(); j++)
					weight
					*= Read::qualityToProbability[(unsigned char) quals[i + j]];
			} else {
				change = Read::qualityToProbability[(unsigned char) quals[i + KmerSizer::getSequenceLength() - 1]] / Read::qualityToProbability[(unsigned char) quals[i - 1]];
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

		STACK_DEALLOC(bools);

		return kmers;
	}
};

#endif

