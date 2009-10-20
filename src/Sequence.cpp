// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.cpp,v 1.3 2009-10-20 20:56:27 cfurman Exp $
//

#include <cstring>
#include <iostream>
#include <stdexcept>

#include "Sequence.h"
#include <cstdlib>

using namespace std;


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

static BaseLocationVectorType compressSequence(const char *bases,  unsigned char *out)
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

static void uncompressSequence(const unsigned char *in , int num_bases, char *bases)
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

  return string(buffer);
}

 void TwoBitSequence::reverseComplement(const NCBI2NA_Type *in, const NCBI2NA_Type *out, SequenceLengthType length)
{
  throw;
}



/*----------------------------- SEQUENCE -------------------------------------------*/
static const std::tr1::shared_ptr<unsigned char> nullSequence(new unsigned char [0]);




 Sequence::Sequence()
 {
   reset();
 };

Sequence::Sequence(std::string name):
 _length(0)
{
  setSequence(name);
}

Sequence::~Sequence()
{
  reset();
}

void Sequence::setSequence(std::string fasta)
{
   setSequence(fasta,0);
}

void Sequence::setSequence(std::string fasta, unsigned int extraBytes)
{
   reset();
   _length = fasta.length();

   unsigned char buffer[get2NASequenceLength()];
   BaseLocationVectorType markupBases = compressSequence(fasta.c_str(),buffer);
   SequenceLengthType markupBasesSize = markupBases.size();
   try {
     unsigned int size =  get2NASequenceLength() + sizeof(SequenceLengthType)+
                          markupBasesSize*sizeof(BaseLocationVectorType) + extraBytes;
     _data = std::tr1::shared_ptr<unsigned char>(new unsigned char[size]);
   } catch (...) {
     throw new std::runtime_error("Cannot allocate memory in Sequence::setSequence()");
   }

   memcpy(get2NASequence(), buffer, get2NASequenceLength());

   memcpy(_getMarkupBaseCount(),&markupBasesSize,sizeof(markupBasesSize));
   BaseLocationType *ptr = _getMarkupBases();
   for (SequenceLengthType i = 0 ; i < markupBasesSize;i++) {
     memcpy(ptr++,&markupBases[i],sizeof(BaseLocationType));
   }

}


void Sequence::reset()
{
  _length = 0;
  _data = nullSequence;
}

string Sequence::getFasta()
{

  string fasta = TwoBitSequence::getFasta(get2NASequence() ,getLength());

  SequenceLengthType markupBasesSize  = *_getMarkupBaseCount();
  if (markupBasesSize > 0) {
    BaseLocationType *markupBases = _getMarkupBases();
    for (BaseLocationType *ptr = markupBases; ptr < markupBases+markupBasesSize; ptr++)
        fasta[ptr->second] = ptr->first;
  }

  return fasta;
}
NCBI2NA_Type *Sequence::get2NASequence()
{
   return _data.get();

}
SequenceLengthType *Sequence::_getMarkupBaseCount()
{
  return  (SequenceLengthType *)(get2NASequence() + get2NASequenceLength());
}

BaseLocationType *Sequence::_getMarkupBases()
{
   return (BaseLocationType *)(_getMarkupBaseCount()+1);
}

SequenceLengthType Sequence::getLength()
{
  return _length;
}
SequenceLengthType Sequence::get2NASequenceLength()
{
  return (getLength() + 3)/4;
}




/*------------------------------------ READ ----------------------------------------*/



Read::Read(std::string name, std::string fasta, std::string qualBytes)
{
  setRead(name,fasta,qualBytes);
}



char * Read::_getQual()
{
  return (char *)(_getMarkupBases() + *_getMarkupBaseCount());
}

char * Read::_getName()
{
  return (char *)(_getQual() + _length);
}



void Read::setRead(std::string name, std::string fasta, std::string qualBytes )
{
   if ( fasta.length() != qualBytes.length())
      throw new std::invalid_argument("fasta length != qual length for name = " + name);

   Sequence::setSequence(fasta,qualBytes.length() + (name.length() + 1) );

   memcpy(_getQual(), qualBytes.c_str(), _length);
   strcpy(_getName(), name.c_str());
}


string Read::getQuals()
{
  return  string(_getQual(), _length);
}

string Read::getName()
{
  return string(_getName());
}


string Read::toFastq()
{
  return  string ('@' + getName() + "\n" + getFasta() + "\n+\n" + getQuals() + "\n" )  ;
}

//
// $Log: Sequence.cpp,v $
// Revision 1.3  2009-10-20 20:56:27  cfurman
// Got it to compile!
//
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//


