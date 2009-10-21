// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.5 2009-10-21 06:51:34 regan Exp $
//

#ifndef _KMER_H
#define _KMER_H

#include "TwoBitSequence.h"


#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 32
#endif

class Kmer
{
  const static SequenceLengthType BYTES = MAX_KMER_SIZE/4;
  typedef unsigned char MerArrayType[BYTES];
  MerArrayType _mer;
public:

   Kmer(NCBI2NA_Type *twoBit, SequenceLengthType offset )
   {
     memcpy(_mer,twoBit+offset/4,sizeof (_mer));
     throw;
     //memm
   }

   int compare(const Kmer &other)
   {
     return memcmp(_mer,other._mer,sizeof(_mer));
   }

   bool operator ==(const Kmer &other)
   {
      return compare(other) == 0;
   }

   Kmer &operator <<(const Kmer &other)
   {
   	 throw;
   }
   
   Kmer reverseComplement()
   {
     throw;
   }

  std::vector<Kmer> generateKmers(NCBI2NA_Type *twoBit, SequenceLengthType length)
  {
    std::vector<Kmer> mers;
    mers.reserve(length - BYTES + 1);
   
    memcpy(&mers[0], twoBit, BYTES);

    for(int i = 1 ; i <  mers.size(); i++)
    {
      //mers[i] <<= mers[i-1] ;
      //mers[i] |= 
    }
    return mers;
  }

};
#endif

//
// $Log: Kmer.h,v $
// Revision 1.5  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
//
