// $Header: /repository/PI_annex/robsandbox/KoMer/src/TwoBitSequence.h,v 1.8 2009-10-23 01:24:53 cfurman Exp $
//

#ifndef _TWO_BIT_SEQUENCE_H
#define _TWO_BIT_SEQUENCE_H

#include <string>
#include <vector>

typedef unsigned int SequenceLengthType;
typedef unsigned char TwoBitEncoding;

typedef std::pair<char,SequenceLengthType>  BaseLocationType;
typedef std::vector<BaseLocationType> BaseLocationVectorType;

class TwoBitSequence // yee hah!
{
private:
   TwoBitSequence();
   static TwoBitSequence singleton;
   
   static void initBitMasks();
   static void initReverseComplementTable();
   static void initBitShiftTable();

public:
   static TwoBitEncoding bitMasks[8];
   static TwoBitEncoding reverseComplementTable[256];
   static TwoBitEncoding bitShiftTable[256*256*3];
   
public:
   static BaseLocationVectorType compressSequence(const char *bases,  TwoBitEncoding *out);
   static void uncompressSequence(const TwoBitEncoding *in , int num_bases, char *bases);
   static void applyMarkup(char *bases, BaseLocationVectorType markupBases);
   static void applyMarkup(std::string &bases, SequenceLengthType markupBasesSize, BaseLocationType *markupBases);
   
   static std::string getFasta(const TwoBitEncoding *in, SequenceLengthType length);
   static void reverseComplement(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType length);
   static void shiftLeft(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType twoBitLength, unsigned char shiftAmountInBases);


    static SequenceLengthType fastaLengthToTwoBitLength(SequenceLengthType fastaLength)
    {
       return (fastaLength + 3)/4;
    }
};

class KmerSizer
{
private:
    static KmerSizer singleton;

private:
   
   KmerSizer(SequenceLengthType sequenceLength, unsigned long extraBytes) {
     set(sequenceLength,extraBytes);
  }


  SequenceLengthType _sequenceLength;
  unsigned long _extraBytes;
  
  SequenceLengthType _twoBitLength;
  unsigned long _totalSize;
public:

  static void set(SequenceLengthType sequenceLength, unsigned long extraBytes=0)
  {
    singleton._sequenceLength = sequenceLength;
    singleton._extraBytes = extraBytes;
    singleton._twoBitLength =  TwoBitSequence::fastaLengthToTwoBitLength(singleton._sequenceLength);
    singleton._totalSize = singleton._twoBitLength + extraBytes;
  }

  static SequenceLengthType getSequenceLength()  {
    return singleton._sequenceLength;
  }
  static unsigned long getExtraBytes()  {
    return singleton._extraBytes;
  }
  static SequenceLengthType getTwoBitLength()   { 
    return singleton._twoBitLength;
  }
  static unsigned long getTotalSize()  {
    return singleton._totalSize;
  }
//   int compare(const void *l, const void *r) const {
//     return memcmp(l, r, _twoBitLength);
//   }
//   void *nextKmer(const void *kmer) const {
//     return (TwoBitEncoding*)kmer + _totalSize;
//   }
//   void *kmerAt(const void *kmer, unsigned long idx) const {
//     return (TwoBitEncoding*)kmer + _totalSize * idx;
//  }
};
#endif

//
// $Log: TwoBitSequence.h,v $
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
