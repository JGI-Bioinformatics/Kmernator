//
// Kmernator/test/MmapTempFileTest.cpp
//
// Author: Rob Egan
//
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/

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
