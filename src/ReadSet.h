// $Header: /repository/PI_annex/robsandbox/KoMer/src/ReadSet.h,v 1.15 2010-01-14 00:50:07 regan Exp $
//

#ifndef _READ_SET_H
#define _READ_SET_H
#include <string>
#include <boost/unordered_map.hpp>

#include "config.h"
#include "Options.h"
#include "Sequence.h"


class ReadFileReader;

class ReadSet {
public:
    typedef unsigned int ReadSetSizeType;
    static const ReadSetSizeType MAX_READ_IDX = (unsigned int) -1;
    class Pair {
    public:
    	ReadSetSizeType read1;
    	ReadSetSizeType read2;
    	
    	Pair() : read1(MAX_READ_IDX), read2(MAX_READ_IDX) {}
    	Pair(ReadSetSizeType _read1) : read1(_read1), read2(MAX_READ_IDX) {}
    	Pair(ReadSetSizeType _read1, ReadSetSizeType _read2) : read1(_read1), read2(_read2) {}
    	Pair(const Pair &copy) : read1(copy.read1), read2(copy.read2) {} 
    	bool operator==(const Pair &other) const { return (read1 == other.read1) && (read2 == other.read2); }
    	bool operator<(const Pair &other) const { return lesser() < other.lesser(); }
    	inline ReadSetSizeType lesser() const { return (read1 < read2 ? read1 : read2); }
    };
    
    typedef std::vector< Pair > PairedIndexType;
    
private:
    std::vector<Read> _reads;
    unsigned long _baseCount;
    PairedIndexType _pairs;
    
private:
    void addRead(Read &read);
    
public:
    ReadSet():_baseCount(0) {}
    ~ReadSet() {}

    void appendAnyFile(std::string filePath, std::string filePath2 = "");
    void appendAllFiles(Options::FileListType &files);
    void appendFastaFile(std::string &is);
    
    void append(ReadSet &reads);

    inline ReadSetSizeType getSize() const { return _reads.size(); }
    unsigned long getBaseCount() const {return _baseCount;}

    inline ReadSetSizeType getPairSize() const { return _pairs.size(); }
    
    inline bool isValidRead(ReadSetSizeType index) const { return index < getSize(); }
    
    inline Read &getRead(ReadSetSizeType index) { return _reads[index]; }
    inline const Read &getRead(ReadSetSizeType index) const { return _reads[index]; }
    

    // by default no pairs are identified
    ReadSetSizeType identifyPairs();
    bool hasPairs() { return getPairSize() != 0 && getPairSize() < getSize(); }
    
    // may return either as MAX_READ_IDX
    inline Pair &getPair(ReadSetSizeType pairIndex) { return _pairs[pairIndex]; }
    inline const Pair &getPair(ReadSetSizeType pairIndex) const { return _pairs[pairIndex]; }

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
// Revision 1.15  2010-01-14 00:50:07  regan
// fixes
//
// Revision 1.14  2010-01-13 23:47:44  regan
// made const class modifications
// fixed identify pairs
//
// Revision 1.13  2010-01-13 00:26:49  regan
// fixed some parallelism
// started pair indentification
//
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
