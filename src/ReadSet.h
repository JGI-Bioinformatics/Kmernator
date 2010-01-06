// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.h,v 1.12 2010-01-06 15:20:24 regan Exp $
//

#ifndef _READ_SET_H
#define _READ_SET_H
#include <string>

#include "config.h"
#include "Options.h"
#include "Sequence.h"


class ReadFileReader;

class ReadSet {
public:
    typedef unsigned int ReadSetSizeType;
    
private:
    std::vector<Read> _reads;
    unsigned long _baseCount;
    
private:
    void addRead(Read &read);
    
public:
    ReadSet();
    ~ReadSet();

    void appendAnyFile(std::string filePath, std::string filePath2 = "");
    void appendAllFiles(Options::FileListType &files);
    void appendFastaFile(std::string &is);
    
    void append(ReadSet &reads);

    ReadSetSizeType getSize();
    unsigned long getBaseCount() {return _baseCount;}

    Read &getRead(ReadSetSizeType index);

protected:
    void appendFasta(std::string fastaFilePath, std::string qualFilePath = "");
    void appendFasta(ReadFileReader &reader);
    
    void appendFastq(std::string fastqFilePath);
    void appendFastqBlockedOMP(std::string fastaFilePath, std::string qualFilePath = "");
    void appendFastqBatchedOMP(std::string fastaFilePath, std::string qualFilePath = "");

};


class ReadIndexScore
{
public:
   ReadSet::ReadSetSizeType readIndex;
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
// Revision 1.12  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.11  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.10  2009-12-23 07:16:52  regan
// fixed reading of fasta files
// parallelized reading of multiple files
//
// Revision 1.9  2009-12-22 18:31:41  regan
// parallelized reading fastq if openmp is enabled
//
// Revision 1.8  2009-11-07 00:28:41  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
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
