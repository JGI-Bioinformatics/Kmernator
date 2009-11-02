// $Header: /repository/PI_annex/robsandbox/KoMer/src/TwoBitSequence.h,v 1.13 2009-11-02 18:24:29 regan Exp $
//

#ifndef _TWO_BIT_SEQUENCE_H
#define _TWO_BIT_SEQUENCE_H

#include <string>
#include <vector>

typedef unsigned int SequenceLengthType;
typedef unsigned char TwoBitEncoding;
typedef TwoBitEncoding *TwoBitEncodingPtr;


typedef std::pair<char,SequenceLengthType>  BaseLocationType;
typedef std::vector<BaseLocationType> BaseLocationVectorType;

class TwoBitSequence // yee hah!
{
private:
   TwoBitSequence();
   static TwoBitSequence singleton;  
   static void initReverseComplementTable();
 

public:
   static TwoBitEncoding bitMasks[8];
   static TwoBitEncoding reverseComplementTable[256];
 
   
public:
   static BaseLocationVectorType compressSequence(const char *bases,  TwoBitEncoding *out);
   static void uncompressSequence(const TwoBitEncoding *in , int num_bases, char *bases);
   static void applyMarkup(char *bases, BaseLocationVectorType markupBases);
   static void applyMarkup(std::string &bases, SequenceLengthType markupBasesSize, BaseLocationType *markupBases);
   
   static std::string getFasta(const TwoBitEncoding *in, SequenceLengthType length);
   static void reverseComplement(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType length);
   static void shiftLeft(const void *in, void *out, SequenceLengthType twoBitLength, unsigned char shiftAmountInBases, bool hasExtraByte = false);


    static SequenceLengthType fastaLengthToTwoBitLength(SequenceLengthType fastaLength)
    {
       return (fastaLength + 3)/4;
    }
};

#endif

//
// $Log: TwoBitSequence.h,v $
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
