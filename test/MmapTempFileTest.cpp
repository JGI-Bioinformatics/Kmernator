// $Header: /repository/PI_annex/robsandbox/KoMer/test/MmapTempFileTest.cpp,v 1.2 2010-05-18 20:50:21 regan Exp $
//

#include "MmapTempFile.h"
#define BOOST_TEST_MODULE MmapTempfileTest
#include <boost/test/unit_test.hpp>


typedef MmapTempFile::MmapFile MmapFile;

#define INIT_AND_TEST(ptr,j) \
	for (long i = 0; i < 100; i++) \
		*(ptr+i) = i+j; \
	for (long i = 0; i < 100; i++) \
		BOOST_CHECK_EQUAL(*(ptr+i) , i+j);

void testMmapTempFile()
{
	MmapFile f = MmapTempFile::buildNewMmap(100 * sizeof(long));
	long *ptr = (long*) f.data();

	INIT_AND_TEST(ptr,0);
	INIT_AND_TEST(ptr,1);
	INIT_AND_TEST(ptr,15);

}

void testMmapAllocator()
{
	typedef MmapTempFile::MmapAllocator< sizeof(long) > Allocator;
	long *ptr = (long*) Allocator::malloc( 100 * sizeof(long) );
	long *ptr2 = (long*) Allocator::malloc( 100 * sizeof(long) );
	long *ptr3 = (long*) Allocator::malloc( 100 * sizeof(long) );

	BOOST_CHECK( ptr < ptr2 );
	BOOST_CHECK( ptr < ptr3 );
	BOOST_CHECK( ptr2 < ptr3 );

	INIT_AND_TEST(ptr,0);
	INIT_AND_TEST(ptr2,2);
	INIT_AND_TEST(ptr3,3);

}

//TODO
void testBoostPoolWithMmapAllocator()
{

}

BOOST_AUTO_TEST_CASE( MmapTempfileTest )
{
	testMmapTempFile();
	testMmapAllocator();
	testBoostPoolWithMmapAllocator();

}

//
// $Log: MmapTempFileTest.cpp,v $
// Revision 1.2  2010-05-18 20:50:21  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.1.2.1  2010-05-12 18:25:27  regan
// added tests
//
