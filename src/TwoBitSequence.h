// $Header: /repository/PI_annex/robsandbox/KoMer/src/TwoBitSequence.h,v 1.4 2009-10-22 01:39:44 cfurman Exp $
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
   static BaseLocationVectorType compressSequence(const char *bases,  unsigned char *out);
   static void uncompressSequence(const unsigned char *in , int num_bases, char *bases);

   static std::string getFasta(const TwoBitEncoding *NCBI2NA, SequenceLengthType length);
   static void reverseComplement(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType length);
   static void shiftLeft(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType length, unsigned char bases)
   { throw; }
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
    singleton._twoBitLength =  (singleton._sequenceLength+3)/4;
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
