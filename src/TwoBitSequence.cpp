//
// Kmernator/src/TwoBitSequence.cpp
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

#include "TwoBitSequence.h"

// _TwoBitEncodingPtr static methods
TwoBitSequenceBase::_TwoBitEncodingPtr *TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(unsigned long size) {
	C *counterPtr = (C*) malloc( size+sizeof(C) );
	if (counterPtr == NULL)
		throw std::bad_alloc();
	// construct / initialize count
	*counterPtr = 0;

	// return address one C into the allocation
	_TwoBitEncodingPtr *ptr = ((_TwoBitEncodingPtr*) (counterPtr+1));

	return ptr;
}

void TwoBitSequenceBase::_TwoBitEncodingPtr::release(TwoBitSequenceBase::_TwoBitEncodingPtr *ptr) {
	C *counterPtr = (C*) ptr;
	// free the original allocation, one C before this reference
	free( counterPtr-1 );
}

// _TwoBitEncodingPtr instance methods
void TwoBitSequenceBase::_TwoBitEncodingPtr::increment() const {
	C &counter = count();
	if (counter == MAX_COUNTER)
		return;

#pragma omp atomic
	++counter;

	if (counter == 0)
		counter = MAX_COUNTER;
}

bool TwoBitSequenceBase::_TwoBitEncodingPtr::decrement() const {
	C &counter = count();
	if (counter == MAX_COUNTER)
		return false;

#pragma omp atomic
	--counter;

	return counter == 0;
}

// intrusive pointer methods for _TwoBitEncodingPtr / TwoBitEncodingPtr
void boost::intrusive_ptr_add_ref(TwoBitSequenceBase::_TwoBitEncodingPtr* r)
{
	r->increment();
}

void boost::intrusive_ptr_release(TwoBitSequenceBase::_TwoBitEncodingPtr* r)
{
	if (r->decrement())
		TwoBitSequenceBase::_TwoBitEncodingPtr::release(r);
}

char TwoBitSequence::uncompressBase(unsigned int v) {
	static const char bases[] = {'A', 'C', 'G', 'T'};
	assert(v <= 3);
	return bases[v];
}

unsigned char TwoBitSequence::compressBase(char base) {
	switch (base) {
	case 'A':
	case 'a':
		return 0;
	case 'C':
	case 'c':
		return 1;
	case 'G':
	case 'g':
		return 2;
	case 'T':
	case 't':
		return 3;
	case '\0':
		return 254;
	default:
		return 255;
	}
}

// TwoBitSequence static methods

// initialize the singleton
TwoBitSequence TwoBitSequence::singleton = TwoBitSequence();
TwoBitSequence::TwoBitSequence() {
	TwoBitSequence::initReverseComplementTable();
	TwoBitSequence::initPermutationsTable();
	TwoBitSequence::initGCTable();
	TwoBitSequence::initShiftLeftMatrix();
	TwoBitSequence::initUncompressSequenceLookupTable();
}

/* initialize reverse complement table
 00000000 -> 11111111
 10000000 -> 11111110
 01000000 -> 11111101
 11000000 -> 11111100
 ...
 11111110 -> 10000000
 11111111 -> 00000000
 */
TwoBitEncoding TwoBitSequence::reverseComplementTable[256];
void TwoBitSequence::initReverseComplementTable() {
	for (int c = 0; c < 256; c++) {
		TwoBitEncoding i = c;
		TwoBitEncoding complement = ~i;
		reverseComplementTable[i] = ((complement << 6) & (0x03 << 6))
						| ((complement << 2) & (0x03 << 4)) | ((complement >> 2)
								& (0x03 << 2)) | ((complement >> 6) & (0x03 << 0));
	}
}

/*
 initialize permutations table - [256][12]
 AAAA -> CAAA GAAA TAAA  ACAA AGAA ATAA  AACA AAGA AATA  AAAC AAAG AAAT
 AAAC -> CAAC GAAC TAAC  ACAC ...
 ...
 TTTT -> ATTT CTTT GTTT  TATT ...
 */
TwoBitEncoding TwoBitSequence::permutations[256* 12 ];
void TwoBitSequence::initPermutationsTable() {
	for (int i = 0; i < 256; i++) {
		TwoBitEncoding tb = i;
		int j = 0;
		for (int baseIdx = 0; baseIdx < 4; baseIdx++) {
			for (int k = 0; k < 4; k++) {
				TwoBitEncoding newBase = k << (6 - baseIdx * 2);
				TwoBitEncoding baseMask = 0x03 << (6 - baseIdx * 2);
				TwoBitEncoding test = (tb & (~baseMask)) | newBase;
				if (test != i)
					permutations[i * 12 + j++] = test;
			}
		}
	}
}

/*
 * initialize GC count table
 */
SequenceLengthType TwoBitSequence::gcCount[256];
// C is 01, G is 10 (1,2 in decimal)
void TwoBitSequence::initGCTable() {
	for (int i = 0; i < 256; i++) {
		TwoBitEncoding tb = i;
		TwoBitEncoding test;
		test = tb & 0x03;
		gcCount[i] = (test == 0x01 || test == 0x02 ? 1 : 0);
		test = tb & 0x0c;
		gcCount[i] += (test == 0x04 || test == 0x08 ? 1 : 0);
		test = tb & 0x30;
		gcCount[i] += (test == 0x10 || test == 0x20 ? 1 : 0);
		test = tb & 0xC0;
		gcCount[i] += (test == 0x40 || test == 0x80 ? 1 : 0);
	}
}

SequenceLengthType TwoBitSequence::getGC(const TwoBitEncoding *in, SequenceLengthType length) {
	SequenceLengthType count = 0;
	SequenceLengthType twoBitLen = fastaLengthToTwoBitLength(length);
	SequenceLengthType len = twoBitLen - 1;
	for(SequenceLengthType i = 0 ; i < len; i++) {
		count += gcCount[ *(in+i) ];
	}
	SequenceLengthType remainder = length & 0x03;
	TwoBitEncoding mask = 0xff;
	switch(remainder) {
	case 1 : mask = 0xc0; break;
	case 2 : mask = 0xf0; break;
	case 3 : mask = 0xfc; break;
	}
	count += gcCount [ *(in+len) & mask ];
	return count;
}

// NULL for out is okay to just get markups
BaseLocationVectorType TwoBitSequence::compressSequence(const char *bases,
		TwoBitEncoding *out) {
	BaseLocationVectorType otherBases;
	SequenceLengthType offset = 0;
	while (bases[offset] != '\0') {
		TwoBitEncoding c = 0;
		for (int i = 6; i >= 0 && *bases; i -= 2) {
			TwoBitEncoding cbase = compressBase(bases[offset]);
			if (cbase == 254)
				break;

			if (cbase == 255) {
				char base = bases[offset];
				// translate . to N
				if (base == '.')
					base = 'N';
				otherBases.push_back(BaseLocationType(base, offset));
				cbase = 0;
			}
			offset++;
			c |= cbase << i;
		}
		if (out != NULL)
			*out++ = c;
	}

	return otherBases;
}

char TwoBitSequence::uncompressSequenceLookupTable[256][4];
void TwoBitSequence::initUncompressSequenceLookupTable() {
	char btable[4] = { 'A', 'C', 'G', 'T' };
	for(int twoBitEnc = 0; twoBitEnc < 256; twoBitEnc++) {
		char *bases = uncompressSequenceLookupTable[twoBitEnc];
		int num_bases = 4;
		while (num_bases) {
			for (int i = 6; i >= 0 && num_bases; i -= 2) {
				char base = btable[(twoBitEnc >> i) & 3];
				*bases++ = base;
				num_bases--;
			}
		}
	}
}
void TwoBitSequence::uncompressSequence(const TwoBitEncoding *in, int num_bases, char *bases) {
	assert(num_bases >= 0);
	for(int i = 0; i < num_bases; i+=4) {
		char *seq = uncompressSequenceLookupTable[*in++];
		memcpy(bases+i, seq, 4);
	}
	*(bases+num_bases) = '\0';	
	if (false) {
	static char btable[4] = { 'A', 'C', 'G', 'T' };
	while (num_bases) {
		for (int i = 6; i >= 0 && num_bases; i -= 2) {
			char base = btable[(*in >> i) & 3];
			*bases++ = base;
			num_bases--;
		}
		in++;
	}
	*bases = '\0';
	}
}

void TwoBitSequence::uncompressSequence(const TwoBitEncoding *in , int num_bases, std::string &bases) {
	STACK_ALLOC(char, tmp, num_bases + 4);

	uncompressSequence(in, num_bases, tmp);
	bases = std::string(tmp,num_bases);

	STACK_DEALLOC(tmp);
}


void TwoBitSequence::applyMarkup(std::string &bases,
		const BaseLocationVectorType &markupBases) {
	SequenceLengthType len = (SequenceLengthType) bases.length();
	for (BaseLocationVectorType::const_iterator ptr = markupBases.begin();
			ptr != markupBases.end();
			ptr++) {
		if (ptr->second >= len)
			break;
		bases[ptr->second] = ptr->first;
	}
}

void TwoBitSequence::applyMarkup(std::string &bases,
		SequenceLengthType markupBasesSize, const BaseLocationType *markupBases) {
	if (markupBasesSize > 0) {
		SequenceLengthType len = (SequenceLengthType) bases.length();
		for (const BaseLocationType *ptr = markupBases;
				ptr < markupBases + markupBasesSize && ptr->second < len;
				ptr++) {
			if (ptr->second >= len)
				break;
			bases[ptr->second] = ptr->first;
		}
	}
}

SequenceLengthType TwoBitSequence::firstMarkup(const BaseLocationVectorType &markups) {
	if (markups.empty()) {
		return 0;
	} else {
		return markups[0].second + 1;
	}
}
SequenceLengthType TwoBitSequence::firstMarkupN(const BaseLocationVectorType &markups) {
	if (markups.empty()) {
		return 0;
	} else {
		for(BaseLocationVectorType::const_iterator it = markups.begin(); it != markups.end(); it++)
			if (it->first == 'N')
				return it->second + 1;
		return 0;
	}
}
SequenceLengthType TwoBitSequence::firstMarkupX(const BaseLocationVectorType &markups) {
	if (markups.empty()) {
		return 0;
	} else {
		for(BaseLocationVectorType::const_iterator it = markups.begin(); it != markups.end(); it++)
			if (it->first == 'X')
				return it->second + 1;
		return 0;
	}
}
SequenceLengthType TwoBitSequence::firstMarkupNorX(const BaseLocationVectorType &markups) {
	if (markups.empty()) {
		return 0;
	} else {
		for(BaseLocationVectorType::const_iterator it = markups.begin(); it != markups.end(); it++)
			if (it->first == 'N' || it->first == 'X')
				return it->second + 1;
		return 0;
	}
}

std::string TwoBitSequence::getFasta(const TwoBitEncoding *in,
		SequenceLengthType offset, SequenceLengthType length) {

	STACK_ALLOC(char, buffer, offset+length+4);

	uncompressSequence(in, offset+length, buffer);

	std::string str(buffer, offset+length);

	STACK_DEALLOC(buffer);

	return str.substr(offset, length);
}

void TwoBitSequence::reverseComplement(const TwoBitEncoding *in,
		TwoBitEncoding *out, SequenceLengthType length) {
	SequenceLengthType twoBitLength = fastaLengthToTwoBitLength(length);

	TwoBitEncoding *tmpOut = out + twoBitLength;

	unsigned long bitShift = length & 0x03;

	while (tmpOut != out)
		*(--tmpOut) = reverseComplementTable[*(in++)];

	if (bitShift > 0) {
		shiftLeft(out, out, twoBitLength, 4 - bitShift);
	}
}

std::string TwoBitSequence::getReverseComplementFasta(const TwoBitEncoding *in, SequenceLengthType length) {
	SequenceLengthType twoBitLength = fastaLengthToTwoBitLength(length);
	TwoBitEncoding rev[twoBitLength];
	reverseComplement(in, rev, length);
	return getFasta(rev, 0, length);
}

void TwoBitSequence::shiftLeft(const void *twoBitIn, void *twoBitOut,
		SequenceLengthType twoBitLength, unsigned char shiftAmountInBases,
		bool hasExtraByte) {
	assert(shiftAmountInBases <= 3);
	assert(twoBitLength > 0);

	TwoBitEncoding *in = (TwoBitEncoding*) twoBitIn;
	TwoBitEncoding *out = (TwoBitEncoding*) twoBitOut;

	if (shiftAmountInBases == 0 && (in != out)) {
		memcpy(out, in, twoBitLength);
		return;
	}

	in += twoBitLength;
	out += twoBitLength;

	unsigned short buffer;
	if (hasExtraByte) {
		buffer = *((unsigned short *)--in);
	} else {
		buffer = *(--in);
	}

	if (false) {
	const int shift = (8-shiftAmountInBases*2);
	for (SequenceLengthType i = 0; i < twoBitLength; i++) {
		TwoBitEncoding byte = ((buffer >> 8) | (buffer << 8)) >> shift;
		if (i < twoBitLength -1)
			buffer = *((unsigned short *)--in);
		*--out = byte;
	}
	} else {
		TwoBitEncoding *shiftLookup = shiftLeftMatrix[shiftAmountInBases-1];
		bool cont = true;
		while (cont) {
			TwoBitEncoding byte = *(shiftLookup + buffer);

			if (in != twoBitIn)
				buffer = *((unsigned short *)--in);
			else
				cont = false;

			*--out = byte;
		}
	}

}

TwoBitEncoding TwoBitSequence::shiftLeftMatrix[3][65536];
void TwoBitSequence::initShiftLeftMatrix() {
	for (int shiftAmountInBases = 1; shiftAmountInBases <= 3; shiftAmountInBases++) {
		const int shift = (8-shiftAmountInBases*2);
		for (int i = 0; i < 65536 ; i++) {
			unsigned short buffer = i;
			TwoBitEncoding byte = ((buffer >> 8) | (buffer << 8)) >> shift;
			shiftLeftMatrix[shiftAmountInBases-1][buffer] = byte;
		}
	}
}

void TwoBitSequence::extendBase(std::string _fasta, char base, void *twoBitOut, bool toRight) {
	char bases[_fasta.length()+2];
	unsigned int idx = 0;
	if ( toRight )
		bases[idx++] = base;
	const char *fasta = _fasta.c_str();
	for (unsigned int i = 0 ; i < _fasta.length(); i++)
		bases[idx++] = fasta[i];
	if (toRight)
		bases[idx++] = base;
	bases[idx++] = '\0';
	assert(idx == _fasta.length() + 2);
	compressSequence(bases,(TwoBitEncoding*) twoBitOut);
}


TwoBitSequenceBase::MarkupElementSizeType TwoBitSequence::getMarkupElementSize(const BaseLocationVectorType &markups) {
	long totalMarkupSize;
	return getMarkupElementSize(markups, totalMarkupSize);
}
TwoBitSequenceBase::MarkupElementSizeType TwoBitSequence::getMarkupElementSize(const BaseLocationVectorType &markups, long &totalMarkupSize) {
	MarkupElementSizeType markupElementSize(0,0);
	SequenceLengthType size = markups.size();
	if (size > 0) {
		if (markups[0].first == 'X' && markups[0].second == 0) {
			// if the first base is masked, the entire read is discarded... no markups will be stored
			// and size should be zero
		} else if (size < 255 && markups[size-1].second < 255) {
			markupElementSize.first = sizeof(SequenceLengthType1);
			markupElementSize.second = sizeof(BaseLocationType1);
		} else if (size < 65535 && markups[size-1].second < 65535) {
			markupElementSize.first = sizeof(SequenceLengthType2);
			markupElementSize.second = sizeof(BaseLocationType2);
		} else {
			markupElementSize.first = sizeof(SequenceLengthType);
			markupElementSize.second = sizeof(BaseLocationType);
		}
		totalMarkupSize = markupElementSize.first + markupElementSize.second * size;
	} else {
		totalMarkupSize = 0;
	}
	return markupElementSize;
}

void TwoBitSequence::permuteBase(const TwoBitEncoding *in, TwoBitEncoding *out1, TwoBitEncoding *out2, TwoBitEncoding *out3,
		SequenceLengthType sequenceLength, SequenceLengthType permuteBaseIdx) {

	assert(permuteBaseIdx < sequenceLength);

	SequenceLengthType bytes = fastaLengthToTwoBitLength(sequenceLength);
	SequenceLengthType permuteByteIdx = fastaLengthToTwoBitLength(permuteBaseIdx+1) - 1;

	// start off identical
	memcpy(out1, in, bytes);
	memcpy(out2, in, bytes);
	memcpy(out3, in, bytes);

	SequenceLengthType j = (permuteBaseIdx & 0x03) * 3;

	unsigned int arrayIdx = *(in + permuteByteIdx) * 12;

	TwoBitEncoding *ptr;

	ptr = out1 + permuteByteIdx;
	*ptr = permutations[ arrayIdx + j ];

	ptr = out2 + permuteByteIdx;
	*ptr = permutations[ arrayIdx + j + 1 ];

	ptr = out3 + permuteByteIdx;
	*ptr = permutations[ arrayIdx + j + 2 ];

}

