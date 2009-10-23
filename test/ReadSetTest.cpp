// $Header: /repository/PI_annex/robsandbox/KoMer/test/ReadSetTest.cpp,v 1.2 2009-10-23 07:06:57 regan Exp $
//
 

#include "ReadSet.h"
#define BOOST_TEST_MODULE ReadSetTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;



void testReadWriteFile(string filename)
{
    ifstream ifs(filename.c_str());
    BOOST_REQUIRE( ! ifs.fail() );
    string fileContents;
    char buffer[1024*1024];
    while (!ifs.eof()) {
    	ifs.getline(buffer, sizeof(buffer));
    	fileContents += string(buffer) + '\n';
    }
    
    string fastq;
    ReadSet store;
    store.appendFastq(filename.c_str());
    
    for (int i=0 ; i < store.getSize(); i++)
    {
       fastq += store.getRead(i).toFastq();
    }
    fastq += '\n'; // not sure why this works...
    
    BOOST_CHECK_EQUAL(fileContents, fastq);
} 

BOOST_AUTO_TEST_CASE( ReadSetTest )
{

  testReadWriteFile(string("10.fastq"));
  
}

//
// $Log: ReadSetTest.cpp,v $
// Revision 1.2  2009-10-23 07:06:57  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
// Revision 1.1  2009-10-23 01:24:56  cfurman
// ReadSet test created
//
 
