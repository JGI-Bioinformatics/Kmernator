// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.h,v 1.4 2009-10-21 06:51:34 regan Exp $
//
#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <string>
#include <vector>

#include <tr1/memory>

#include "TwoBitSequence.h"

class Sequence
{

protected:
  SequenceLengthType _length;

/*
   _data contains a composite of 1 fields:
      +0                 : the sequence as NCBI 2NA (2 bits per base ACGT)
      += (length +3)/4   : non-ACGT bases: count followed by array of markups
*/
  std::tr1::shared_ptr<unsigned char> _data;

  SequenceLengthType *_getMarkupBaseCount();
  BaseLocationType   *_getMarkupBases();

  void reset();

  void setSequence( std::string fasta,unsigned int extraBytes);
public:

  Sequence();
  Sequence(std::string fasta);

  ~Sequence();

  void  setSequence(std::string fasta);

  SequenceLengthType getLength();
  std::string getFasta();

  SequenceLengthType get2NASequenceLength();
  NCBI2NA_Type *get2NASequence();

};


class Read : public Sequence
{
private:

/*
   _data is inherited and now contains a composite of 4 fields:
      +0                    : the sequence as NCBI 2NA (2 bits per base ACGT)
      += (length +3)/4      :  non-ACGT bases: count followed by array of markups
      += getMarkupLength()  :qualities as 1 byte per base, 0 = N 1..255 Phred Quality Score.
      += length             : null terminated name.
 */


  char * _getQual();
  char * _getName();

public:

  Read():Sequence(){ };
  Read(std::string name, std::string fasta, std::string qualBytes);

  void  setRead(std::string name, std::string fasta, std::string qualBytes);

  std::string getName();
  std::string getQuals();

  std::string toFastq();
};

#endif

//
// $Log: Sequence.h,v $
// Revision 1.4  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
// Revision 1.3  2009-10-21 00:00:58  cfurman
// working on kmers....
//
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//
