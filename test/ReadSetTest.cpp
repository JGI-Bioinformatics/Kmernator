// $Header: /repository/PI_annex/robsandbox/KoMer/test/ReadSetTest.cpp,v 1.1 2009-10-23 01:24:56 cfurman Exp $
//
 

#include "ReadSet.h"
#define BOOST_TEST_MODULE ReadSetTest
#include <boost/test/unit_test.hpp>

using namespace std;




BOOST_AUTO_TEST_CASE( ReadSetTest )
{

    string fastq;
    ReadSet store;
    store.appendFastq("10.fastq");

 
    for (int i=0 ; i < store.getSize(); i++)
    {
       fastq += store.getRead(i).toFastq();
    }



   BOOST_CHECK_EQUAL("hello", "there");
} 

//
// $Log: ReadSetTest.cpp,v $
// Revision 1.1  2009-10-23 01:24:56  cfurman
// ReadSet test created
//
 
