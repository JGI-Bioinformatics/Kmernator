// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.7 2009-10-21 18:58:44 regan Exp $
//

#ifndef _KMER_H
#define _KMER_H
#include <tr1/memory>

#include "TwoBitSequence.h"

class KmerContainer
{
private:
  SequenceLengthType _sequenceLength;
  unsigned long _extraBytes;
  
  SequenceLengthType _twoBitLength;
  unsigned long _totalSize;
public:
  KmerContainer(SequenceLengthType sequenceLength, unsigned long extraBytes) {
  	_sequenceLength = sequenceLength;
  	_extraBytes = extraBytes;
  	_twoBitLength =  (_sequenceLength+3)/4;
  	_totalSize = _twoBitLength + extraBytes;
  }
  SequenceLengthType getSequenceLength() const {
  	return _sequenceLength;
  }
  unsigned long getExtraBytes() const {
  	return _extraBytes;
  }
  SequenceLengthType getTwoBitLength() const { 
  	return _twoBitLength;
  }
  unsigned long getTotalSize() const {
  	return _totalSize;
  }
  int compare(const void *l, const void *r) const {
  	return memcmp(l, r, _twoBitLength);
  }
  void *nextKmer(const void *kmer) const {
  	return (NCBI2NA_Type*)kmer + _totalSize;
  }
  void *kmerAt(const void *kmer, unsigned long idx) const {
  	return (NCBI2NA_Type*)kmer + _totalSize * idx;
  }
}

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 1024
#endif

typedef std::tr1::shared_ptr<NCBI2NA_Type> KmerSharedPtr;

class Kmer
{

private:
   Kmer() { throw; } // never construct, just use as cast
   
   static KmerContainer container;
   NCBI2NA_Type _data[MAX_KMER_SIZE]; // need somedata to hold a pointer and a large amount to avoid memory warnings
   
public:

   static void setKmerContainer(const KmerContainer &kmerContainer) {
   	  if (kmerContainer.getTwoBitLength() > MAX_KMER_SIZE) {
   	  	throw;
   	  }
   	  container = KmerContainer(kmerContainer.getSequenceLength(), kmerContainer.getExtraBytes());
   }
   static unsigned long getByteSize() {
   	  return container.getTotalSize();
   }
   
   int compare(const Kmer &other) const
   {
     return memcmp(&_data, &other._data, container.getTwoBitLength());
   }
   bool operator ==(const Kmer &other)
   {
      return compare(other) == 0;
   }
   bool operator <(const Kmer &other)
   {
   	  return compare(other) < 0;
   }
   bool operator <=(const Kmer &other)
   {
   	  return compare(other) <= 0;
   }
   bool operator >(const Kmer &other)
   {
   	  return compare(other) > 0;
   }
   bool operator >=(const Kmer &other)
   {
   	  return compare(other) >= 0;
   }
   // override [], ptr++, ptr--, = 
   
   void reverseComplement(Kmer &output) const
   {
     TwoBitSequence::reverseComplement((NCBI2NA_Type*)&_data, (NCBI2NA_Type*)&output._data, container.getSequenceLength());
   }
   

   KmerSharedPtr generateKmers(NCBI2NA_Type *twoBit, SequenceLengthType length) const
   {
   	  SequenceLengthType numKmers = length - container.getSequenceLength() + 1;
      KmerSharedPtr ptr(calloc(getByteSize(), numKmers));
      
      if (ptr == NULL)
        throw;
      
      for(SequenceLengthType i=0; i < numKmers ; i+=4) {
      	// 0 base shift
      	memcpy(ptr+i, twoBit + i/4, container.getTwoBitLength());
      	
      	NCBI2NA_Type *ref = twoBit+i/4;
      	if (i+1 < numKmers) {
      	  unsigned short twoByte = *( (unsigned short*)ref );
      	
      	  for (int bitShift=1; bitShift < 4; bitShift++) {
            if (i+bitShift < numKmers) {
      	      for(SequenceLengthType bytes=0; bytes < container.getTwoBitLength(), bytes++) {
      	      	NCBI2NA_Type *kmer = (NCBI2NA_Type*) ptr+i+bitShift;
      	      	NCBI2NA_Type *kmerByte = kmer+bytes;
      	        *(kmerByte) = TwoBitSequence::bitShiftTable[ twoByte + (bitShift-1)];
      	      }
      	    }
          }
      	}
      }
      return ptr;
  }

};
#endif

//
// $Log: Kmer.h,v $
// Revision 1.7  2009-10-21 18:58:44  regan
// checkpoint
//
// Revision 1.6  2009-10-21 18:44:20  regan
// checkpoint
//
// Revision 1.5  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
//
