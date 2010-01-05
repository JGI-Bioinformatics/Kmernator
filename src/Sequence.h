// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.h,v 1.10 2010-01-05 06:44:39 regan Exp $
//
#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <string>
#include <vector>

#include <tr1/memory>

#include "TwoBitSequence.h"

class Sequence
{
public:
  typedef TwoBitSequenceBase::SequenceLengthType SequenceLengthType;
  const static SequenceLengthType MAX_SEQUENCE_LENGTH = (SequenceLengthType) -1;
  
protected:
  SequenceLengthType _length;

/*
   _data contains a composite of 1 fields:
      +0                 : the sequence as NCBI 2NA (2 bits per base ACGT)
      += (length +3)/4   : non-ACGT bases: count followed by array of markups
*/
  std::tr1::shared_ptr<unsigned char> _data;

  SequenceLengthType *_getMarkupBasesCount();
  BaseLocationType   *_getMarkupBases();

  void reset();

  void setSequence( std::string fasta,unsigned int extraBytes);
public:

  Sequence();
  Sequence(std::string fasta);

  ~Sequence();

  void  setSequence(std::string fasta);

  SequenceLengthType getLength();
  std::string getFasta(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH);
  BaseLocationVectorType getMarkups();

  SequenceLengthType getTwoBitEncodingSequenceLength();
  TwoBitEncoding *getTwoBitSequence();

};


class Read : public Sequence
{
private:

	/*
   _data is inherited and now contains a composite of 4 fields:
      +0                    : the sequence as NCBI 2NA (2 bits per base ACGT)
      += (length +3)/4      :  non-ACGT bases: count followed by array of markups
      += getMarkupLength()  : qualities as 1 byte per base, 0 = N 1..255 Phred Quality Score.
      += length             : null terminated name.
 */


  char * _getQual();
  char * _getName();

  static int qualityToProbabilityInitialized;
  static int initializeQualityToProbability();
  
public:

  Read():Sequence(){ };
  Read(std::string name, std::string fasta, std::string qualBytes);

  void  setRead(std::string name, std::string fasta, std::string qualBytes);

  std::string getName();
  std::string getQuals(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH);

  std::string toFastq(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH, std::string label = "");
  std::string getFormattedQuals(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH);

  
  static double qualityToProbability[256];
  
};

#endif

//
// $Log: Sequence.h,v $
// Revision 1.10  2010-01-05 06:44:39  regan
// fixed warnings
//
// Revision 1.9  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.8  2009-11-07 00:28:41  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
// Revision 1.7  2009-11-02 18:27:00  regan
// added getMarkups()
// added quality to probability lookup table
//
// Revision 1.6  2009-10-22 07:04:06  regan
// added a few unit tests
// minor refactor
//
// Revision 1.5  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
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
