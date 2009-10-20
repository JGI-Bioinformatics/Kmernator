#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <string>
#include <vector>

#include <tr1/memory>

#define MAX_SEQUENCE_LENGTH 1024*1024

typedef unsigned int SequenceLengthType;
typedef unsigned char NCBI2NA_Type;
class TwoBitSequence // yee hah!
{
private:
  TwoBitSequence();
  
public:
 
   static std::string getFasta(const NCBI2NA_Type *NCBI2NA, SequenceLengthType length);
   static void reverseComplement(const NCBI2NA_Type *in, const NCBI2NA_Type *out, SequenceLengthType length);
};



typedef std::pair<char,SequenceLengthType>  BaseLocationType;
typedef std::vector<BaseLocationType> BaseLocationVectorType;

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
  
  void  setSequence( std::string fasta);
  
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
