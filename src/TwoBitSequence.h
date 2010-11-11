//
// Kmernator/src/TwoBitSequence.h
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

#ifndef _TWO_BIT_SEQUENCE_H
#define _TWO_BIT_SEQUENCE_H

#include <cstring>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include <boost/intrusive_ptr.hpp>

#include "config.h"

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

}

namespace boost {
void intrusive_ptr_add_ref(TwoBitSequenceBase::_TwoBitEncodingPtr* r);
void intrusive_ptr_release(TwoBitSequenceBase::_TwoBitEncodingPtr* r);
}

using namespace TwoBitSequenceBase;

class TwoBitSequence // yee hah!
{
private:
	TwoBitSequence();
	static TwoBitSequence singleton;
	static void initReverseComplementTable();
	static void initPermutationsTable();
	static void initGCTable();

public:
	static TwoBitEncoding bitMasks[8];
	static TwoBitEncoding reverseComplementTable[256];
	static TwoBitEncoding permutations[256 * 12];
	static SequenceLengthType  gcCount[256];

public:
	// NULL for out is okay to just get markups
	static BaseLocationVectorType compressSequence(const char *bases,TwoBitEncoding *out);
	static BaseLocationVectorType compressSequence(const std::string &bases, TwoBitEncoding *out) {
		return compressSequence(bases.c_str(), out);
	}
    static void uncompressSequence(const TwoBitEncoding *in , int num_bases, char *bases);
    static void uncompressSequence(const TwoBitEncoding *in , int num_bases, std::string &bases);

    static void applyMarkup(std::string &bases, const BaseLocationVectorType &markupBases);
    static void applyMarkup(std::string &bases, SequenceLengthType markupBasesSize, const BaseLocationType *markupBases);

    static SequenceLengthType firstMarkup(const BaseLocationVectorType &markups);
    static SequenceLengthType firstMarkupN(const BaseLocationVectorType &markups);
    static SequenceLengthType firstMarkupX(const BaseLocationVectorType &markups);
    static SequenceLengthType firstMarkupNorX(const BaseLocationVectorType &markups);

    static std::string getFasta(const TwoBitEncoding *in, SequenceLengthType offset, SequenceLengthType length);
    
    static void reverseComplement(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType length);

    static void shiftLeft(const void *in, void *out, SequenceLengthType twoBitLength, unsigned char shiftAmountInBases, bool hasExtraByte = false);

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

//
// $Log: TwoBitSequence.h,v $
// Revision 1.26  2010-05-24 21:48:46  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.25.2.1  2010-05-19 21:36:54  regan
// refactored duplicate fragment filter code
// added duplicate fragment on single ended reads
//
// Revision 1.25  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.24.2.4  2010-05-18 16:43:31  regan
// added count gc methods and lookup tables
//
// Revision 1.24.2.3  2010-05-13 00:19:49  regan
// bugfix in types that affected some reference markups
//
// Revision 1.24.2.2  2010-05-10 19:41:33  regan
// minor refactor moved code into cpp
//
// Revision 1.24.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.24  2010-05-06 21:46:53  regan
// merged changes from PerformanceTuning-20100501
//
// Revision 1.23.2.1  2010-05-04 19:49:51  regan
// minor rework on include headers
//
// Revision 1.23  2010-05-01 21:57:53  regan
// merged head with serial threaded build partitioning
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
// Revision 1.19.2.6.2.3  2010-04-16 17:49:42  regan
// fixed comments
//
// Revision 1.19.2.6.2.2  2010-04-16 17:43:19  regan
// changed allocation / deallocation rules to hide counter
//
// Revision 1.19.2.6.2.1  2010-04-16 05:30:00  regan
// checkpoint.. broke it
//
// Revision 1.19.2.6  2010-04-14 22:36:06  regan
// round of bugfixes
//
// Revision 1.19.2.5  2010-04-14 03:51:19  regan
// checkpoint. compiles but segfaults
//
// Revision 1.19.2.4  2010-04-04 16:22:31  regan
// bugfix
//
// Revision 1.19.2.3  2010-04-04 15:58:02  regan
// fixed assertion code to obey debug rules
//
// Revision 1.19.2.2  2010-04-04 15:29:28  regan
// migrated markup types and methods
//
// Revision 1.19.2.1  2010-04-04 12:00:13  regan
// added special smaller markup types
//
// Revision 1.19  2010-03-02 15:03:18  regan
// added a method
//
// Revision 1.18  2010-03-02 13:51:38  regan
// reformatted
//
// Revision 1.17  2010-02-26 13:01:16  regan
// reformatted
//
// Revision 1.16  2010-01-13 23:34:59  regan
// made const class modifications
//
// Revision 1.15  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.14  2009-11-06 16:59:11  regan
// added base substitution/permutations table and build function
//
// Revision 1.13  2009-11-02 18:24:29  regan
// *** empty log message ***
//
// Revision 1.12  2009-10-27 22:13:41  cfurman
// removed bit shift table
//
// Revision 1.11  2009-10-26 23:03:05  regan
// checkpoint
//
// Revision 1.10  2009-10-26 17:38:42  regan
// moved KmerSizer to Kmer.h
//
// Revision 1.9  2009-10-23 23:22:41  regan
// checkpoint
//
// Revision 1.8  2009-10-23 01:24:53  cfurman
// ReadSet test created
//
// Revision 1.7  2009-10-23 00:13:54  cfurman
// reverse complement now works
//
// Revision 1.6  2009-10-22 20:49:16  cfurman
// tests added
//
// Revision 1.5  2009-10-22 07:04:06  regan
// added a few unit tests
// minor refactor
//
// Revision 1.4  2009-10-22 01:39:44  cfurman
// bug fix in kmer.h
//
// Revision 1.3  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
// Revision 1.2  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
//
