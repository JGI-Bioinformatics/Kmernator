#ifndef _TWO_BIT_SEQUENCE_H
#define _TWO_BIT_SEQUENCE_H

#include <string>

typedef unsigned int SequenceLengthType;
typedef unsigned char NCBI2NA_Type;

typedef std::pair<char,SequenceLengthType>  BaseLocationType;
typedef std::vector<BaseLocationType> BaseLocationVectorType;

namespace TwoBitSequence // yee hah!
{
   static BaseLocationVectorType compressSequence(const char *bases,  unsigned char *out);
   static void uncompressSequence(const unsigned char *in , int num_bases, char *bases);

   static std::string getFasta(const NCBI2NA_Type *NCBI2NA, SequenceLengthType length);
   static void reverseComplement(const NCBI2NA_Type *in, const NCBI2NA_Type *out, SequenceLengthType length);
};

static unsigned char compressBase(char base)
{
  switch (base) {
    case 'A' : return 0;
    case 'C' : return 1;
    case 'G' : return 2;
    case 'T' : return 3;
    default :  return 255;
  }
}

BaseLocationVectorType TwoBitSequence::compressSequence(const char *bases,  unsigned char *out)
{
  BaseLocationVectorType otherBases;
  SequenceLengthType offset = 0;
  while (bases[offset] != '\0') {
    unsigned char c = 0;
    for (int i = 6; i >= 0  && *bases; i-= 2) {
      unsigned char cbase = compressBase(bases[offset]);
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

void TwoBitSequence::uncompressSequence(const unsigned char *in , int num_bases, char *bases)
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


std::string TwoBitSequence::getFasta(const NCBI2NA_Type *NCBI2NA, SequenceLengthType length)
{
  char buffer[length+1];
  uncompressSequence(NCBI2NA,length, buffer);

  return std::string(buffer);
}

 void TwoBitSequence::reverseComplement(const NCBI2NA_Type *in, const NCBI2NA_Type *out, SequenceLengthType length)
{
  throw;
}



#endif
