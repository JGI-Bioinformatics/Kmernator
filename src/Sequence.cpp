// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.cpp,v 1.15 2009-11-28 01:00:07 regan Exp $
//

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include "Sequence.h"
#include <cstdlib>
#include <cmath>

using namespace std;





/*----------------------------- SEQUENCE -------------------------------------------*/
static const std::tr1::shared_ptr<TwoBitEncoding> nullSequence(new TwoBitEncoding [0]);




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

   unsigned long buffSize = getTwoBitEncodingSequenceLength();
   TwoBitEncoding *buffer = (TwoBitEncoding*) malloc(buffSize);
   if ( buffer == NULL )
     throw std::runtime_error("Could not allocate buffer memory in Sequence::setSequence");
     
   const char *f_str = fasta.c_str();
   BaseLocationVectorType markupBases = TwoBitSequence::compressSequence(f_str,buffer);
   SequenceLengthType markupBasesSize = markupBases.size();
   try {
     unsigned int size =  getTwoBitEncodingSequenceLength() + sizeof(SequenceLengthType)+
                          markupBasesSize*sizeof(BaseLocationVectorType) + extraBytes;
     _data = std::tr1::shared_ptr<unsigned char>(new unsigned char[size]);
   } catch (...) {
     throw new std::runtime_error("Cannot allocate memory in Sequence::setSequence()");
   }

   memcpy(getTwoBitSequence(), buffer, getTwoBitEncodingSequenceLength());
   free(buffer);
   
   memcpy(_getMarkupBasesCount(),&markupBasesSize,sizeof(markupBasesSize));
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
  if(_data == nullSequence)
    return string("");
  string fasta = TwoBitSequence::getFasta(getTwoBitSequence() ,getLength());
  TwoBitSequence::applyMarkup(fasta, *_getMarkupBasesCount(), _getMarkupBases());

  return fasta;
}
TwoBitEncoding *Sequence::getTwoBitSequence()
{
   return _data.get();

}
SequenceLengthType *Sequence::_getMarkupBasesCount()
{
  return  (SequenceLengthType *)(getTwoBitSequence() + getTwoBitEncodingSequenceLength());
}

BaseLocationType *Sequence::_getMarkupBases()
{
   return (BaseLocationType *)(_getMarkupBasesCount()+1);
}

SequenceLengthType Sequence::getLength()
{
  return _length;
}

SequenceLengthType Sequence::getTwoBitEncodingSequenceLength()
{
  return TwoBitSequence::fastaLengthToTwoBitLength(getLength());;
}

BaseLocationVectorType Sequence::getMarkups()
{
  SequenceLengthType size =  *_getMarkupBasesCount();
  BaseLocationVectorType markups( size );
  if (size > 0) {
  	BaseLocationType *ptr = _getMarkupBases();
  	for(unsigned int i=0; i<size; i++)
  	  markups.push_back( *(ptr++) );
  }
  return markups;
}

/*------------------------------------ READ ----------------------------------------*/

double Read::qualityToProbability[256];
int Read::initializeQualityToProbability() 
{
  for (int i=0; i<256; i++) {
  	qualityToProbability[i] = 0;
  }
  int start = 64;
  for (int i = start ; i < 164 ; i++ )
    qualityToProbability[i] = 1.0 - pow(10.0,( ( start - i ) / 10.0 ));

  qualityToProbability[255] = 1.0; // for reads with no quality data
  return 1;
}
int Read::qualityToProbabilityInitialized = Read::initializeQualityToProbability();


Read::Read(std::string name, std::string fasta, std::string qualBytes)
{
  setRead(name,fasta,qualBytes);
}



char * Read::_getQual()
{
  return (char *)(_getMarkupBases() + *_getMarkupBasesCount());
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


string Read::getName()
{
  if(_data == nullSequence)
    return string("");
  else
    return string(_getName());    
}

string Read::getQuals()
{
  return  string(_getQual(), _length);
}


string Read::toFastq()
{
  return  string ('@' + getName() + "\n" + getFasta() + "\n+\n" + getQuals() + "\n" )  ;
}
string Read::getFormattedQuals()
{
  string quals = getQuals();
  stringstream ss;
  for(unsigned int i =0; i < quals.length(); i++)
  {
     ss << (int)quals[i]-64 << ' ';  
  }
  return ss.str();
}
//
// $Log: Sequence.cpp,v $
// Revision 1.15  2009-11-28 01:00:07  regan
// fixed bugs and warnings
//
// Revision 1.14  2009-11-21 15:58:29  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.13  2009-11-10 07:04:59  regan
// bugfix in quality lookup table -- use float!
//
// Revision 1.12  2009-11-07 00:28:41  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
// Revision 1.11  2009-11-06 04:07:13  regan
// bugfix when stack size is limited
//
// Revision 1.10  2009-11-04 19:32:03  cfurman
// now reads in fasta (with optional qual) files
//
// Revision 1.9  2009-11-02 18:27:00  regan
// added getMarkups()
// added quality to probability lookup table
//
// Revision 1.8  2009-10-30 00:51:40  regan
// bug fix and working on executable
//
// Revision 1.7  2009-10-22 20:49:15  cfurman
// tests added
//
// Revision 1.6  2009-10-22 07:04:06  regan
// added a few unit tests
// minor refactor
//
// Revision 1.5  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
// Revision 1.4  2009-10-21 00:00:58  cfurman
// working on kmers....
//
// Revision 1.3  2009-10-20 20:56:27  cfurman
// Got it to compile!
//
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//


