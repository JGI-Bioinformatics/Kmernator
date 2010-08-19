//
// Kmernator/test/MmapTempFileTest.cpp
//
// Author: Rob Egan, Craig Furman
//
// Copyright 2010 The Regents of the University of California.
// All rights reserved.
//
// The United States Government has rights in this work pursuant
// to contracts DE-AC03-76SF00098, W-7405-ENG-36 and/or
// W-7405-ENG-48 between the United States Department of Energy
// and the University of California.
//
// Redistribution and use in source and binary forms are permitted
// provided that: (1) source distributions retain this entire
// copyright notice and comment, and (2) distributions including
// binaries display the following acknowledgement:  "This product
// includes software developed by the University of California,
// JGI-PSF and its contributors" in the documentation or other
// materials provided with the distribution and in all advertising
// materials mentioning features or use of this software.  Neither the
// name of the University nor the names of its contributors may be
// used to endorse or promote products derived from this software
// without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE.
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
