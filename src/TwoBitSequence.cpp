//
// Kmernator/src/TwoBitSequence.cpp
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

#include "TwoBitSequence.h"

// _TwoBitEncodingPtr static methods
TwoBitSequenceBase::_TwoBitEncodingPtr *TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(unsigned long size) {
	C *counterPtr = (C*) malloc( size+sizeof(C) );
	if (counterPtr == NULL)
		throw;
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
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
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

void TwoBitSequence::uncompressSequence(const TwoBitEncoding *in,
		int num_bases, char *bases) {
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

void TwoBitSequence::uncompressSequence(const TwoBitEncoding *in , int num_bases, std::string &bases) {
        if (num_bases < 1024) {
            char tmp[num_bases];
            uncompressSequence(in, num_bases, tmp);
            bases = std::string(tmp,num_bases);
        } else {
            char *tmp = (char*) malloc(num_bases);
            if (tmp == NULL) throw std::bad_alloc();
            uncompressSequence(in, num_bases, tmp);
            bases = std::string(tmp,num_bases);
            free(tmp);
       }
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

	bool needMalloc = (offset+length) < MAX_STACK_SIZE - 1;

	char _buffer[needMalloc ? 0 : offset+length+1];
	char *buffer = _buffer;
	if (needMalloc) {
		buffer = new char[offset+length+1];
	}
	uncompressSequence(in, offset+length, buffer);

	std::string str(buffer, offset+length);
	if (needMalloc)
		delete [] buffer;

	return str.substr(offset, length);
}

void TwoBitSequence::reverseComplement(const TwoBitEncoding *in,
		TwoBitEncoding *out, SequenceLengthType length) {
	SequenceLengthType twoBitLength = fastaLengthToTwoBitLength(length);

	out += twoBitLength;
	unsigned long bitShift = length & 0x03;

	for (SequenceLengthType i = 0; i < twoBitLength; i++)
		*(--out) = reverseComplementTable[*in++];

	if (bitShift > 0) {
		shiftLeft(out, out, twoBitLength, 4 - bitShift);
	}
}

void TwoBitSequence::shiftLeft(const void *twoBitIn, void *twoBitOut,
		SequenceLengthType twoBitLength, unsigned char shiftAmountInBases,
		bool hasExtraByte) {
	TwoBitEncoding *in = (TwoBitEncoding*) twoBitIn;
	TwoBitEncoding *out = (TwoBitEncoding*) twoBitOut;

	if (shiftAmountInBases == 0 && (in != out)) {
		memcpy(out, in, twoBitLength);
		return;
	}

	if (shiftAmountInBases > 3) {
		throw ;
	}

	in+= twoBitLength;
	out += twoBitLength;

	unsigned short buffer;
	if (hasExtraByte) {
		buffer = *((unsigned short *)--in);
	} else {
		buffer = *(--in);
	}

	const int shift = (8-shiftAmountInBases*2);
	for (SequenceLengthType i = 0; i < twoBitLength; i++) {
		TwoBitEncoding byte = ((buffer >> 8) | (buffer << 8)) >> shift;
		if (i < twoBitLength -1)
		buffer = *((unsigned short *)--in);
		*--out = byte;
	}

}

void TwoBitSequence::extendBase(std::string _fasta, char base, void *twoBitOut, bool toRight) {
	char bases[_fasta.length()+2];
	int idx = 0;
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


//
// $Log: TwoBitSequence.cpp,v $
// Revision 1.26  2010-06-22 23:06:31  regan
// merged changes in CorruptionBugfix-20100622 branch
//
// Revision 1.25.6.1  2010-06-22 22:59:12  regan
// modified markups to actually save memory when discarded
//
// Revision 1.25  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.24.2.3  2010-05-18 16:43:31  regan
// added count gc methods and lookup tables
//
// Revision 1.24.2.2  2010-05-10 19:41:33  regan
// minor refactor moved code into cpp
//
// Revision 1.24.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.24  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
// Revision 1.23.2.1  2010-05-04 19:49:51  regan
// minor rework on include headers
//
// Revision 1.23  2010-05-01 21:57:53  regan
// merged head with serial threaded build partitioning
//
// Revision 1.21.2.4  2010-04-29 20:33:23  regan
// minor optimization
//
// Revision 1.21.2.3  2010-04-28 23:45:29  regan
// made '.' recognized as 'N'
//
// Revision 1.21.2.2  2010-04-28 22:59:19  regan
// made '.' recognized as 'N'
//
// Revision 1.21.2.1  2010-04-23 17:46:20  regan
// merged bugfixes from head
//
// Revision 1.22  2010-04-22 23:41:32  regan
// fixed a few bugs
//
// Revision 1.21  2010-04-21 00:33:20  regan
// merged with branch to detect duplicated fragment pairs with edit distance
//
// Revision 1.20.2.1  2010-04-19 18:20:51  regan
// refactored base permutation
//
// Revision 1.20  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.19.2.3.2.2  2010-04-16 17:43:19  regan
// changed allocation / deallocation rules to hide counter
//
// Revision 1.19.2.3.2.1  2010-04-16 05:30:00  regan
// checkpoint.. broke it
//
// Revision 1.19.2.3  2010-04-14 03:51:19  regan
// checkpoint. compiles but segfaults
//
// Revision 1.19.2.2  2010-04-04 16:22:31  regan
// bugfix
//
// Revision 1.19.2.1  2010-04-04 15:29:28  regan
// migrated markup types and methods
//
// Revision 1.19  2010-03-02 14:27:00  regan
// reformatted
//
// Revision 1.18  2010-02-26 13:01:16  regan
// reformatted
//
// Revision 1.17  2010-01-13 23:34:59  regan
// made const class modifications
//
// Revision 1.16  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.15  2009-11-06 16:59:11  regan
// added base substitution/permutations table and build function
//
// Revision 1.14  2009-11-02 18:24:29  regan
// *** empty log message ***
//
// Revision 1.13  2009-10-31 00:16:35  regan
// minor changes and optimizations
//
// Revision 1.12  2009-10-27 22:13:41  cfurman
// removed bit shift table
//
// Revision 1.11  2009-10-26 17:38:43  regan
// moved KmerSizer to Kmer.h
//
// Revision 1.10  2009-10-23 23:22:41  regan
// checkpoint
//
// Revision 1.9  2009-10-23 01:24:53  cfurman
// ReadSet test created
//
// Revision 1.8  2009-10-23 00:13:54  cfurman
// reverse complement now works
//
// Revision 1.7  2009-10-22 21:46:49  regan
// fixed ushort to ulong conversion problems
//
// Revision 1.6  2009-10-22 20:49:15  cfurman
// tests added
//
// Revision 1.5  2009-10-22 07:04:06  regan
// added a few unit tests
// minor refactor
//
// Revision 1.4  2009-10-22 01:39:43  cfurman
// bug fix in kmer.h
//
// Revision 1.3  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
// Revision 1.2  2009-10-21 18:44:20  regan
// checkpoint
//
// Revision 1.1  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
//
