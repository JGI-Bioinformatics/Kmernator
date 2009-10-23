#include "Kmer.h"


KmerPtr Kmer::operator&()
{
   return KmerPtr(_data());
}
