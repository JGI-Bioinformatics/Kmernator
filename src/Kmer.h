// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.2 2009-10-20 17:25:50 regan Exp $
//

#ifndef _KMER_H
#define _KMER_H

template <unsigned char BYTES>
class _kmer
{
  typedef unsigned char MerArrayType[BYTES];
  MerArrayType _mer;
public:

   int compare(const _kmer<BYTES> &other)
   {
     return memcmp(_mer,other._mer,BYTES);
   };

   bool operator ==(const _kmer<BYTES> &other)
   {
      return compare(other) == 0;
   }

   static MerArrayType reverseComplement(MerArrayType &mer)
   {
     throw;
   }
};





typedef _kmer<4> Kmer16;

#endif

//
// $Log: Kmer.h,v $
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//
