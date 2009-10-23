#include "Kmer.h"


KmerPtr Kmer::operator&()
{
   return KmerPtr(_data());
}

void KmerArray::build(TwoBitEncoding *twoBit, SequenceLengthType length)
{
  SequenceLengthType numKmers = length - Kmer::getLength() + 1;
  if (_size != numKmers)
    throw;
    
  KmerArray &kmers = *this;  
  for(SequenceLengthType i=0; i < numKmers ; i+=4) {
    TwoBitEncoding *ref = twoBit+i/4;
    for (int bitShift=0; bitShift < 4 && i+bitShift < numKmers; bitShift++) 
      TwoBitSequence::shiftLeft(ref, &kmers[i+bitShift], Kmer::getTwoBitLength(), bitShift, bitShift != 0);
  }
}
