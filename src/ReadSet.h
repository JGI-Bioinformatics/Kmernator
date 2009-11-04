// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.h,v 1.7 2009-11-04 19:32:03 cfurman Exp $
//

#ifndef _READ_SET_H
#define _READ_SET_H
#include <string>
#include "Sequence.h"


typedef unsigned int ReadSetSizeType;

class ReadSet {
    std::vector<Read> _reads;
    unsigned long _baseCount;
    
private:
    void push_back(Read &read);
    
public:
    ReadSet();
    ~ReadSet();

    void appendFastq(std::string fastqFilePath);
    void appendFasta(std::string fastaFilePath, std::string qualFilePath = "");

    ReadSetSizeType getSize();
    unsigned long getBaseCount() {return _baseCount;}

    Read &getRead(ReadSetSizeType index);

};


class ReadIndexScore
{
public:
   ReadSetSizeType readIndex;
   float score; 
};


class KmerReadSetStats
{
public:
   float  kmerScore;
   std::vector<ReadIndexScore> linkedReads;
};


#endif

//
// $Log: ReadSet.h,v $
// Revision 1.7  2009-11-04 19:32:03  cfurman
// now reads in fasta (with optional qual) files
//
// Revision 1.6  2009-10-31 00:16:35  regan
// minor changes and optimizations
//
// Revision 1.5  2009-10-26 23:02:49  regan
// checkpoint
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
