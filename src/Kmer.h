#ifndef _KMER_H
#define _KMER_H

#include "TwoBitSequence.h"


#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 32
#endif

class Kmer
{
  typedef unsigned char MerArrayType[MAX_KMER_SIZE/4];
  MerArrayType _mer;
public:

   Kmer(NCBI2NA_Type *twoBit, SequenceLengthType offset )
   {
     memcpy(_mer,twoBit+offset/4,sizeof (_mer));
     
     memm
   }

   int compare(const Kmer &other)
   {
     return memcmp(_mer,other._mer,sizeof(_mer));
   };

   bool operator ==(const Kmer &other)
   {
      return compare(other) == 0;
   }

   Kmer reverseComplement()
   {

   }

//    static MerArrayType reverseComplement(MerArrayType &mer)
//    {
//      MerArrayType rc = mer;
//      return rc;
//    }
};

std::vector<Kmer> generateKmers(NCBI2NA_Type *twoBit,SequenceLengthType length)
{
   std::vector<Kmer> mers;
   mers.reserve(length-MAX_KMER_SIZE/4 + 1)
   
   memcpy(&mers[0],twoBit, MAX_KMER_SIZE/4);

   for(int i = 1 ; i <  mers.size(); i++)
   {
     mers[i] <<= mers[i-1] ;
     mers[i] |= 
   }

}

#endif
