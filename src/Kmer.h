



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
     M
     return 
   }
};





typedef _kmer<4> Kmer16;