#include "Sequence.h"


template <unsigned char BYTES>
class _kmer  
{
  typedef unsigned char MerArrayType[BYTES];
  MerArrayType _mer;
public:

   _kmer(NCBI2NA_Type *seq)
   {
     memcpy(_mer,seq,BYTES);
   }

   int compare(const _kmer<BYTES> &other)
   {
     return memcmp(_mer,other._mer,BYTES);
   };
   
   bool operator ==(const _kmer<BYTES> &other)
   {
      return compare(other) == 0;
   }

//    static MerArrayType reverseComplement(MerArrayType &mer)
//    {
//      MerArrayType rc = mer;
//      return rc;
//    }
};


typedef _kmer<4> Kmer16;
typedef _kmer<5> Kmer20;
typedef _kmer<6> Kmer24;
typedef _kmer<7> Kmer28;