// $Header: /repository/PI_annex/robsandbox/KoMer/test/ReadSetTest.cpp,v 1.7 2009-12-23 07:16:50 regan Exp $
//
 

#include "ReadSet.h"
#define BOOST_TEST_MODULE ReadSetTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;


void testZeroReads()
{
	Read s;
	
	BOOST_CHECK_EQUAL( 0ul, s.getLength() );
	BOOST_CHECK_EQUAL( "", s.getQuals() );
	BOOST_CHECK_EQUAL( "", s.getFasta() );
	BOOST_CHECK_EQUAL( "", s.getName()  );
}


static string getFileContents(string filename)
{
    ifstream ifs(filename.c_str());
    BOOST_REQUIRE( ! ifs.fail() );
    string contents;
 
    getline(ifs,contents,'\0');
    return contents;
}

void testFastQFile(string filename)
{
    string fileContents = getFileContents(filename);
    
    ReadSet store;
    store.appendAnyFile(filename);

    string fastq;
    for (unsigned int i=0 ; i < store.getSize(); i++)
    {
       fastq += store.getRead(i).toFastq();
    }
 
    BOOST_CHECK_EQUAL(fileContents, fastq);
} 
  
void testFastaWithQualFile(string f,string q)
{
    string fFile = getFileContents(f);
    string qFile = getFileContents(q);
         
    ReadSet store;
    store.appendAnyFile(f,q);

    string fasta,qual;
    for (unsigned int i=0 ; i < store.getSize(); i++)
    {
       Read &s = store.getRead(i);

       string nameLine('>' + s.getName()  + "\n");
       fasta += nameLine;
       fasta += s.getFasta() + "\n";

       qual += nameLine;
       qual += s.getFormattedQuals() + "\n";

    }
 
    BOOST_CHECK_EQUAL(fFile, fasta);
    BOOST_CHECK_EQUAL(qFile,qual);
}


BOOST_AUTO_TEST_CASE( ReadSetTest )
{
  testZeroReads();
  testFastQFile("10.fastq");
  testFastaWithQualFile("10.fasta","10.qual");
}

//
// $Log: ReadSetTest.cpp,v $
// Revision 1.7  2009-12-23 07:16:50  regan
// fixed reading of fasta files
// parallelized reading of multiple files
//
// Revision 1.6  2009-11-28 01:00:10  regan
// fixed bugs and warnings
//
// Revision 1.5  2009-11-07 00:28:38  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
// Revision 1.4  2009-10-30 00:51:37  regan
// bug fix and working on executable
//
// Revision 1.3  2009-10-23 20:32:52  cfurman
// more kmer changes
//
// Revision 1.2  2009-10-23 07:06:57  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
// Revision 1.1  2009-10-23 01:24:56  cfurman
// ReadSet test created
//
 
