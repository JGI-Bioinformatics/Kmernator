// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.h,v 1.4 2009-10-21 06:51:34 regan Exp $
//

#ifndef _READ_SET_H
#define _READ_SET_H
#include <string>
#include <tr1/unordered_map>
#include "Kmer.h"
#include "Sequence.h"


typedef unsigned int ReadSetSizeType;

class ReadSet {
    class body;
    body &_my;
public:
    ReadSet();
    ~ReadSet();

    void appendFastq(std::string fastqFilePath);

    ReadSetSizeType getSize();

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

typedef std::tr1::unordered_map<Kmer, KmerReadSetStats> KmerReadSetStatsMap;

#endif

//
// $Log: ReadSet.h,v $
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
