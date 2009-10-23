// $Header: /repository/PI_annex/robsandbox/KoMer/src/TwoBitSequence.cpp,v 1.9 2009-10-23 01:24:53 cfurman Exp $
//

#include <cstring>
#include "TwoBitSequence.h"

static unsigned char compressBase(char base)
{
  switch (base) {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    case '\0': return 254;
    default  :  return 255;
  }
}

// initialize the singleton
TwoBitSequence TwoBitSequence::singleton = TwoBitSequence(); 
TwoBitSequence::TwoBitSequence()
{
	TwoBitSequence::initBitMasks();
	TwoBitSequence::initReverseComplementTable();
	TwoBitSequence::initBitShiftTable();
}

/* initialize bitMasks[]
  0: 00000001
  1: 00000010
  2: 00000100
  3: 00001000
  4: 00010000
  5: 00100000
  6: 01000000
  7: 10000000
  */
TwoBitEncoding TwoBitSequence::bitMasks[8];
void TwoBitSequence::initBitMasks()
{
  bitMasks[0] = 0x01;
  for(int i=1; i<8; i++)
    bitMasks[i] = bitMasks[i-1]<<1;
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
void TwoBitSequence::initReverseComplementTable() 
{
  for(TwoBitEncoding i=0; i<255; i++) {
  	TwoBitEncoding complement = ~i;
    reverseComplementTable[i] = ((complement << 6) & (0x03 << 6)) |
                                ((complement << 2) & (0x03 << 4)) |
                                ((complement >> 2) & (0x03 << 2)) |
                                ((complement >> 6) & (0x03 << 0)) ;
  }
}


/* initialize bit shift table
   TwoBitSequence::bitShiftTable
   
   every two-byte possibility (65536) x 3 entries (excluding trivial non-bit-shift)
   */
TwoBitEncoding TwoBitSequence::bitShiftTable[256*256*3];
void TwoBitSequence::initBitShiftTable()   
{
  unsigned short i=0;
  do {
  	const TwoBitEncoding *ptr = (TwoBitEncoding *) &i;
 
//     TwoBitSequence::bitShiftTable[(unsigned long)i*3ul+0ul] = (i >> 6)&0xff;
//     TwoBitSequence::bitShiftTable[(unsigned long)i*3ul+1ul] = (i >> 4)&0xff;
//     TwoBitSequence::bitShiftTable[(unsigned long)i*3ul+2ul] = (i >> 2)&0xff;

  	TwoBitSequence::bitShiftTable[(unsigned long)i*3ul+0ul] = ((*ptr)<<2 & 0xfc) | ((*(ptr+1))>>6 & 0x03);
  	TwoBitSequence::bitShiftTable[(unsigned long)i*3ul+1ul] = ((*ptr)<<4 & 0xf0) | ((*(ptr+1))>>4 & 0x0f);
  	TwoBitSequence::bitShiftTable[(unsigned long)i*3ul+2ul] = ((*ptr)<<6 & 0xc0) | ((*(ptr+1))>>2 & 0x3f); 
 
  } while(++i != 0);
}

BaseLocationVectorType TwoBitSequence::compressSequence(const char *bases,  TwoBitEncoding *out)
{
  BaseLocationVectorType otherBases;
  SequenceLengthType offset = 0;
  while (bases[offset] != '\0') {
    TwoBitEncoding c = 0;
    for (int i = 6; i >= 0  && *bases; i-= 2) {
      TwoBitEncoding cbase = compressBase(bases[offset]);
      if (cbase == 254)
      	break;
      
      if (cbase == 255)
      {
         otherBases.push_back(BaseLocationType(bases[offset],offset));
         cbase = 0;
      }
      offset++;
      c |= cbase <<  i;
    }
    *out++ = c;
  }

  return otherBases;
}

void TwoBitSequence::uncompressSequence(const TwoBitEncoding *in , int num_bases, char *bases)
{
  static char btable[4] = { 'A','C','G','T'};
  while(num_bases) {
    for (int i = 6; i >= 0  && num_bases; i-= 2) {
      char base = btable[(*in >> i) & 3];
      *bases++ = base;
      num_bases--;
    }
    in++;
  }
  *bases = '\0';
}

void TwoBitSequence::applyMarkup(char *bases, BaseLocationVectorType markupBases) {
  for(BaseLocationVectorType::iterator ptr = markupBases.begin(); ptr != markupBases.end() ; ptr++)
    bases[ptr->second] = ptr->first;
}
void TwoBitSequence::applyMarkup(std::string &bases, SequenceLengthType markupBasesSize, BaseLocationType *markupBases)
{
  if (markupBasesSize > 0) {
    for (BaseLocationType *ptr = markupBases; ptr < markupBases+markupBasesSize; ptr++)
        bases[ptr->second] = ptr->first;
  }
}

std::string TwoBitSequence::getFasta(const TwoBitEncoding *in, SequenceLengthType length)
{
  char buffer[length+1];
  uncompressSequence(in,length, buffer);

  return std::string(buffer);
}

void TwoBitSequence::reverseComplement(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType length)
{
   
 
  SequenceLengthType twoBitLength = fastaLengthToTwoBitLength(length);
   
  out+=twoBitLength;
  unsigned long bitShift = length %4;
 
  for(SequenceLengthType i = 0; i<twoBitLength; i++)
    *(--out) = reverseComplementTable[*in++];

  if (bitShift > 0)
  {
     shiftLeft(out,out,twoBitLength,4-bitShift);
  }
  
}

#define BIT_SHIFT(basesIn, targetBases, baseShift)\
    TwoBitSequence::compressSequence(basesIn, in);\
    out[0] = TwoBitSequence::bitShiftTable[((unsigned long)*((unsigned short*)in))*3ul+baseShift-1];\
    TwoBitSequence::uncompressSequence(out,4,fasta);\
    BOOST_CHECK_EQUAL(targetBases,fasta);


void TwoBitSequence::shiftLeft(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType twoBitLength, unsigned char shiftAmountInBases)
{
    if (shiftAmountInBases == 0 && (in != out)) {
       memcpy(out,in,twoBitLength);
       return;
    }
    
    if (shiftAmountInBases > 3)
      throw;
      
    in+= twoBitLength;
    out += twoBitLength;
    
    unsigned short buffer = *(--in);
     
    for (SequenceLengthType i = 0; i < twoBitLength; i++) {
       TwoBitEncoding byte = bitShiftTable[(unsigned long)buffer*3ul+shiftAmountInBases-1];
       if (i < twoBitLength -1)
          buffer = *((unsigned short *)--in);
       *--out = byte;
    }

}

KmerSizer KmerSizer::singleton = KmerSizer(21,0);

//
// $Log: TwoBitSequence.cpp,v $
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
