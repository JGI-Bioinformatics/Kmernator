// $Header: /repository/PI_annex/robsandbox/KoMer/src/TwoBitSequence.h,v 1.2 2009-10-21 06:51:34 regan Exp $
//

#ifndef _TWO_BIT_SEQUENCE_H
#define _TWO_BIT_SEQUENCE_H

#include <string>
#include <vector>

typedef unsigned int SequenceLengthType;
typedef unsigned char NCBI2NA_Type;

typedef std::pair<char,SequenceLengthType>  BaseLocationType;
typedef std::vector<BaseLocationType> BaseLocationVectorType;

class TwoBitSequence // yee hah!
{
private:
   TwoBitSequence();
   static TwoBitSequence singleton;
   
   static NCBI2NA_Type bitMasks[8];
   static NCBI2NA_Type reverseComplementTable[256];
   static NCBI2NA_Type bitShiftTable[256*256*3];
   
   static void initBitMasks();
   static void initReverseComplementTable();
   static void initBitShiftTable();
public:
   static BaseLocationVectorType compressSequence(const char *bases,  unsigned char *out);
   static void uncompressSequence(const unsigned char *in , int num_bases, char *bases);

   static std::string getFasta(const NCBI2NA_Type *NCBI2NA, SequenceLengthType length);
   static void reverseComplement(const NCBI2NA_Type *in, NCBI2NA_Type *out, SequenceLengthType length);
   static void shiftLeft(const NCBI2NA_Type *in, NCBI2NA_Type *out, SequenceLengthType length, unsigned char bases)
   { throw; }
};


#endif

//
// $Log: TwoBitSequence.h,v $
// Revision 1.2  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
//
