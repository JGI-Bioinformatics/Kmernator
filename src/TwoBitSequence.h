// $Header: /repository/PI_annex/robsandbox/KoMer/src/TwoBitSequence.h,v 1.21 2010-04-21 00:33:20 regan Exp $
//

#ifndef _TWO_BIT_SEQUENCE_H
#define _TWO_BIT_SEQUENCE_H

#include <string>
#include <vector>
#include <cstdlib>

#include <boost/intrusive_ptr.hpp>

#include "config.h"

namespace TwoBitSequenceBase {
typedef unsigned int SequenceLengthType;
typedef unsigned short SequenceLengthType2;
typedef unsigned char SequenceLengthType1;

typedef unsigned char TwoBitEncoding;

typedef std::pair<char, SequenceLengthType> BaseLocationType;
typedef std::vector<BaseLocationType> BaseLocationVectorType;
// alternate forms of base location type to save space
typedef std::pair<char, SequenceLengthType2> BaseLocationType2;
typedef std::pair<char, SequenceLengthType1> BaseLocationType1;
typedef std::pair<SequenceLengthType1, SequenceLengthType1> MarkupElementSizeType;

class _TwoBitEncodingPtr {
public:
	typedef TwoBitEncoding T;
	typedef unsigned char C;

public:
    static const T MAX_COUNTER = -1;

    _TwoBitEncodingPtr(); // never explicitly construct
    ~_TwoBitEncodingPtr(); // never explicitly destruct

    static _TwoBitEncodingPtr *allocate(unsigned long size) {
    	C *counterPtr = (C*) malloc( size+sizeof(C) );
    	if (counterPtr == NULL)
    		throw;
    	// construct / initialize count
        *counterPtr = 0;

        // return address one C into the allocation
    	_TwoBitEncodingPtr *ptr = ((_TwoBitEncodingPtr*) (counterPtr+1));

    	return ptr;
    }
    static void release(_TwoBitEncodingPtr *ptr) {
    	C *counterPtr = (C*) ptr;
    	// free the original allocation, one C before this reference
    	free( counterPtr-1 );
    }

    C *getCounter() const { return ((C*) this) - 1; }
    C &count() const { return *getCounter(); }

    void increment() const {
    	C &counter = count();
    	if (counter == MAX_COUNTER)
    		return;
#pragma omp atomic
        ++counter;
    	if (counter == 0)
    		counter = MAX_COUNTER;
    }

    bool decrement() const {
    	C &counter = count();
    	if (counter == MAX_COUNTER)
    		return false;
#pragma omp atomic
    	--counter;
        return counter == 0;
    }

    T *get() const {
      return (T*) (( (C*)this ) + 1);
    }
    T &operator*() const {
      return *get();
    }
    T *operator->() const {
      return get();
    }
    T *operator&() const {
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

public:
	static TwoBitEncoding bitMasks[8];
	static TwoBitEncoding reverseComplementTable[256];
	static TwoBitEncoding permutations[256* 12 ];

public:
	// NULL for out is okay to just get markups
	static BaseLocationVectorType compressSequence(const char *bases,TwoBitEncoding *out);
	static BaseLocationVectorType compressSequence(const std::string &bases, TwoBitEncoding *out) {
		return compressSequence(bases.c_str(), out);
	}
    static void uncompressSequence(const TwoBitEncoding *in , int num_bases, char *bases);
    static void uncompressSequence(const TwoBitEncoding *in , int num_bases, std::string &bases) {
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
    static void applyMarkup(std::string &bases, const BaseLocationVectorType &markupBases);
    static void applyMarkup(std::string &bases, SequenceLengthType markupBasesSize, const BaseLocationType *markupBases);

    static std::string getFasta(const TwoBitEncoding *in, SequenceLengthType length);
    static void reverseComplement(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType length);
    static void shiftLeft(const void *in, void *out, SequenceLengthType twoBitLength, unsigned char shiftAmountInBases, bool hasExtraByte = false);

	static MarkupElementSizeType getMarkupElementSize(const BaseLocationVectorType &markups);
	static MarkupElementSizeType getMarkupElementSize(const BaseLocationVectorType &markups, long &totalMarkupSize);

	static void permuteBase(const TwoBitEncoding *in, TwoBitEncoding *out1, TwoBitEncoding *out2, TwoBitEncoding *out3,
			SequenceLengthType sequenceLength, SequenceLengthType permuteBaseIdx);

    static SequenceLengthType fastaLengthToTwoBitLength(SequenceLengthType fastaLength)
    {
       return (fastaLength + 3)/4;
    }
};

#endif

//
// $Log: TwoBitSequence.h,v $
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
