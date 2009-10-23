// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.h,v 1.13 2009-10-23 20:32:50 cfurman Exp $
//

#ifndef _KMER_H
#define _KMER_H
#include <tr1/memory>
#include <cstring>
#include <cstdlib>
#include "TwoBitSequence.h"



#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 1024
#endif

typedef std::tr1::shared_ptr<TwoBitEncoding> KmerSharedPtr;

class KmerPtr;
class KmerArray;

class Kmer
{

private:
   Kmer(); // never construct, just use as cast

#ifdef STRICT_MEM_CHECK
   TwoBitEncoding _someData[MAX_KMER_SIZE]; // need somedata to hold a pointer and a large amount to avoid memory warnings
  const void *_data() const { return _someData;}
   void *_data()  { return _someData;}
#else
   // No data for you!!!
   const void *_data() const { return this;}
   void *_data()  { return this;}
#endif

public:
   
   int compare(const Kmer &other) const
   {
     return memcmp(_data(), other._data(), KmerSizer::getTwoBitLength());
   }

   Kmer &operator=(const Kmer &other)
   {
      memcpy(_data(), other._data(), KmerSizer::getTwoBitLength());
      return *this;
   }
//    KmerPtr operator&() const;
   KmerPtr operator&() ;    

   bool operator ==(const Kmer &other) const
   {
      return compare(other) == 0;
   }
   bool operator !=(const Kmer &other) const
   {
   	  return compare(other) != 0;
   }
   bool operator <(const Kmer &other) const
   {
   	  return compare(other) < 0;
   }
   bool operator <=(const Kmer &other) const
   {
   	  return compare(other) <= 0;
   }
   bool operator >(const Kmer &other) const
   {
   	  return compare(other) > 0;
   }
   bool operator >=(const Kmer &other) const
   {
   	  return compare(other) >= 0;
   }

   void swap(Kmer &other)
   {
      TwoBitEncoding buffer[getByteSize()];
      Kmer &temp = (Kmer &)buffer;
      other = *this;
      *this = temp;
   }

   TwoBitEncoding *getTwoBitSequence() const
   {
     return (TwoBitEncoding *)_data();
   }

   void reverseComplement(Kmer &output) const
   {
     TwoBitSequence::reverseComplement((TwoBitEncoding*)_data(), (TwoBitEncoding*)output._data(), KmerSizer::getSequenceLength());
   }

   static SequenceLengthType getLength() {
     return KmerSizer::getSequenceLength();
   }

   static SequenceLengthType getTwoBitLength() {
     return KmerSizer::getTwoBitLength();
   }
   static unsigned long getByteSize() {
      return KmerSizer::getTotalSize();
   }

   std::string toFasta() const
   {
      return TwoBitSequence::getFasta(getTwoBitSequence(), getLength());
   }
   
};

typedef TwoBitEncoding *TwoBitEncodingPtr;
typedef Kmer *RawKmerPtr;
typedef void *VoidPtr;

class KmerPtr 
{
private:
   Kmer *_me;
    //Kmer  *operator->() const { return _me;  }
public:
   KmerPtr():
   _me(NULL)
   {  }
   
   KmerPtr(Kmer &in):
   _me((&in)._me)
   { }

   KmerPtr( void *in):
   _me((Kmer *)in)
   { }

   void *get() const  { return _me; };
   Kmer & operator*() const  { return *_me; }
  
      KmerPtr operator->() const { return KmerPtr(_me);  }

   KmerPtr &operator=(void *right)          { _me = (Kmer *)right; return *this; }
   KmerPtr &operator=(const KmerPtr &right) { _me = right._me ;    return *this; }
   
 //  bool operator==(const KmerPtr &right) const { return _me == right._me; }
 //  bool operator!=(const KmerPtr &right) const { return _me != right._me; }

   KmerPtr  operator+ (unsigned long right) const { return KmerPtr((Kmer *)((char *)_me + right * Kmer::getByteSize())); }
   KmerPtr  operator- (unsigned long right) const { return *this + (-right); }
   
   KmerPtr &operator+=(unsigned long right)       { *this = *this + right; return *this; }
   KmerPtr &operator-=(unsigned long right)       { *this = *this - right; return *this; }

   KmerPtr &operator++()           { return *this += 1;}
   KmerPtr operator++(int unused)  { KmerPtr saved = *this; ++(*this); return saved; }

   KmerPtr &operator--()           { return *this -= 1;}
   KmerPtr operator--(int unused)  { KmerPtr saved = *this; --(*this); return saved; }

   Kmer &operator[](unsigned long index) const { return *(*this + index); }

   // cast operator
     //  operator VoidPtr() { return (VoidPtr)_me ; }
   //  operator RawKmerPtr() { return _me; }
   // operator TwoBitEncodingPtr() { return (TwoBitEncodingPtr)_me; }
};

 

 
class KmerArray
{

private:
  KmerPtr _begin;
  unsigned long _size;

  KmerArray();
public:
  
  
  KmerArray(unsigned long size):
   _size(size)
  {
   if (_begin.get() == NULL)
       throw;
  }

  ~KmerArray()
  {
     reset();
  }

  KmerArray &operator=(const KmerArray &other)
  {
    reset();
    resize(other.size());
    memcpy(_begin.get(),other._begin.get(),_size*Kmer::getByteSize());
  }


  Kmer &operator[](unsigned long index) const
  {
    if (index >= _size)
       throw; 
    return _begin[index];
  }

  Kmer &get(unsigned long index) const
  {
    return (*this)[index];
  }

  unsigned int size() const { return _size; }

  void reset()
  {
    void *test = (void *)_begin.get();
    if (test != NULL)
     free(test);
    _begin =  NULL;
    _size = 0;
  }

  void resize(unsigned long size)
  {
    void *temp = _begin.get();
    _begin = KmerPtr(calloc(Kmer::getByteSize(),size));
    if (temp != NULL && _size > 0) {
       memcpy(_begin.get(),temp,_size*Kmer::getByteSize());
       free(temp);
    }
    _size = size;
  }
    

  KmerArray(TwoBitEncoding *twoBit, SequenceLengthType length):
  _size(0),
  _begin(NULL)
  {
    SequenceLengthType numKmers = length - Kmer::getLength() + 1;
    resize(numKmers);
    build(twoBit,length);
  }


   void build(TwoBitEncoding *twoBit, SequenceLengthType length)
   {
      SequenceLengthType numKmers = length - Kmer::getLength() + 1;
      if (_size != numKmers)
        throw;
        
      KmerArray &kmers = *this;  
      for(SequenceLengthType i=0; i < numKmers ; i+=4) {
        // 0 base shift
   //     memcpy(&kmers[i], twoBit + i/4, Kmer::getTwoBitLength());
        
        TwoBitEncoding *ref = twoBit+i/4;
        if (i+1 < numKmers) {
          unsigned short twoByte = *( (unsigned short*)ref );
        
          for (int bitShift=1; bitShift < 4; bitShift++) {
            if (i+bitShift < numKmers) {
              for(SequenceLengthType bytes=0; bytes < Kmer::getTwoBitLength(); bytes++) {
            //    TwoBitEncoding *kmer = (TwoBitEncoding*)((void *) &kmers[i + bitShift]);
            //    TwoBitEncoding *kmerByte = kmer+bytes;
            //    *(kmerByte) = TwoBitSequence::bitShiftTable[ twoByte + (bitShift-1)];
              }
            }
          }
        }
      }
    }
    
};




#endif


 

//
// $Log: Kmer.h,v $
// Revision 1.13  2009-10-23 20:32:50  cfurman
// more kmer changes
//
// Revision 1.12  2009-10-23 17:22:39  regan
// added more tests
//
// Revision 1.11  2009-10-23 07:06:59  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
// Revision 1.10  2009-10-23 01:24:53  cfurman
// ReadSet test created
//
// Revision 1.9  2009-10-22 01:39:43  cfurman
// bug fix in kmer.h
//
// Revision 1.8  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
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
