//
// Kmernator/src/TwoBitSequence.h
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

#ifndef _TWO_BIT_SEQUENCE_H
#define _TWO_BIT_SEQUENCE_H

#include <cstring>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "config.h"

namespace TwoBitSequenceBase {
class _TwoBitEncodingPtr;
};

namespace boost {
void intrusive_ptr_add_ref(TwoBitSequenceBase::_TwoBitEncodingPtr* r);
void intrusive_ptr_release(TwoBitSequenceBase::_TwoBitEncodingPtr* r);
};

#include <boost/intrusive_ptr.hpp>

namespace TwoBitSequenceBase {

typedef unsigned char TwoBitEncoding;

typedef Kmernator::SequenceLengthType SequenceLengthType;
typedef Kmernator::SequenceLengthType2 SequenceLengthType2;
typedef Kmernator::SequenceLengthType1 SequenceLengthType1;

typedef std::pair<char, SequenceLengthType> BaseLocationType;
typedef std::vector<BaseLocationType> BaseLocationVectorType;
// alternate forms of base location type to save space
typedef std::pair<char, SequenceLengthType2> BaseLocationType2;
typedef std::pair<char, SequenceLengthType1> BaseLocationType1;
typedef std::pair<SequenceLengthType1, SequenceLengthType1> MarkupElementSizeType;

class _TwoBitEncodingPtr {
public:
	typedef TwoBitEncoding T;
	typedef Kmernator::UI8 C;

public:
	static const T MAX_COUNTER = MAX_UI8;

	_TwoBitEncodingPtr(); // never explicitly construct
	~_TwoBitEncodingPtr(); // never explicitly destruct

	static _TwoBitEncodingPtr *allocate(unsigned long size);
	static void release(_TwoBitEncodingPtr *ptr);

	inline C *getCounter() const { return ((C*) this) - 1; }
	inline C &count() const { return *getCounter(); }

	void increment() const;

	bool decrement() const;

	inline T *get() const {
		return (T*) (( (C*)this ) + 1);
	}
	inline T &operator*() const {
		return *get();
	}
	inline T *operator->() const {
		return get();
	}
	inline T *operator&() const {
		return get();
	}
};

typedef boost::intrusive_ptr< _TwoBitEncodingPtr > TwoBitEncodingPtr;

};


using namespace TwoBitSequenceBase;

class TwoBitSequence // yee hah!
{
private:
	TwoBitSequence();
	static TwoBitSequence singleton;
	static void initReverseComplementTable();
	static void initPermutationsTable();
	static void initGCTable();
	static void initShiftLeftMatrix();
	static void initUncompressSequenceLookupTable();

public:
	static TwoBitEncoding bitMasks[8];
	static TwoBitEncoding reverseComplementTable[256];
	static TwoBitEncoding permutations[256 * 12];
	static SequenceLengthType  gcCount[256];
	static TwoBitEncoding shiftLeftMatrix[3][65536];
	static char uncompressSequenceLookupTable[256][4];

public:
	static char uncompressBase(unsigned int v);
	static unsigned char compressBase(char base);

	// NULL for out is okay to just get markups
	static BaseLocationVectorType compressSequence(const char *bases,TwoBitEncoding *out);
	static BaseLocationVectorType compressSequence(const std::string &bases, TwoBitEncoding *out) {
		return compressSequence(bases.c_str(), out);
	}
	static void uncompressSequence(const TwoBitEncoding *in , int num_bases, char *bases);
	static void uncompressSequence(const TwoBitEncoding *in , int num_bases, std::string &bases);

	static void applyMarkup(std::string &bases, const BaseLocationVectorType &markupBases);
	static void applyMarkup(std::string &bases, SequenceLengthType markupBasesSize, const BaseLocationType *markupBases);

	// returns 0 if no markup otherwise the base position (1+)
	static SequenceLengthType firstMarkup(const BaseLocationVectorType &markups);
	// returns 0 if no markup otherwise the base position (1+)
	static SequenceLengthType firstMarkupN(const BaseLocationVectorType &markups);
	// returns 0 if no markup otherwise the base position (1+)
	static SequenceLengthType firstMarkupX(const BaseLocationVectorType &markups);
	// returns 0 if no markup otherwise the base position (1+)
	static SequenceLengthType firstMarkupNorX(const BaseLocationVectorType &markups);

	static std::string getFasta(const TwoBitEncoding *in, SequenceLengthType offset, SequenceLengthType length);

	static void reverseComplement(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType length);

	static std::string getReverseComplementFasta(const TwoBitEncoding *in, SequenceLengthType length);

	static void shiftLeft(const void *in, void *out, SequenceLengthType twoBitLength, unsigned char shiftAmountInBases, bool hasExtraByte = false);

	static void extendBase(std::string _fasta, char base, void *twoBitOut, bool toRight);

	static MarkupElementSizeType getMarkupElementSize(const BaseLocationVectorType &markups);
	static MarkupElementSizeType getMarkupElementSize(const BaseLocationVectorType &markups, long &totalMarkupSize);

	static void permuteBase(const TwoBitEncoding *in, TwoBitEncoding *out1, TwoBitEncoding *out2, TwoBitEncoding *out3,
			SequenceLengthType sequenceLength, SequenceLengthType permuteBaseIdx);

	static SequenceLengthType getGC(const TwoBitEncoding *in, SequenceLengthType length);

	static inline SequenceLengthType fastaLengthToTwoBitLength(SequenceLengthType fastaLength)
	{
		return (fastaLength + 3)/4;
	}
};

#endif
