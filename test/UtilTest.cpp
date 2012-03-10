//
// Kmernator/test/UtilTest.cpp
//
// Author: Rob Egan
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

#include "Utils.h"
#define BOOST_TEST_MODULE UtilTest
#include <boost/test/unit_test.hpp>

typedef Statistics::MeanStdCount MSC;

void testBimodalPartition()
{
	MSC f,s;
	double a[10] = { 30, 34, 31, 34, 33, 29, 4, 3, 2, 4 };
	double *p;

	// Test 0, 1 and 2 sized data sets
	p = Statistics::findBimodalPartition(0, f, s, a, a);
	BOOST_CHECK_EQUAL( p, a);
	p = Statistics::findBimodalPartition(0, f, s, a, a+1);
	BOOST_CHECK_EQUAL( p, a+1);
	p = Statistics::findBimodalPartition(0, f, s, a, a+2);
	BOOST_CHECK_EQUAL( p, a+2);
	p = Statistics::findBimodalPartition(0, f, s, a+4, a+6);
	BOOST_CHECK_EQUAL( p, a+6);
	p = Statistics::findBimodalPartition(0, f, s, a+5, a+6);
	BOOST_CHECK_EQUAL( p, a+6);
	p = Statistics::findBimodalPartition(0, f, s, a+5, a+7);
	BOOST_CHECK_EQUAL( p, a+7);

	// No bimodal partition in first 6 values with 2 sigma
	for(int i = 0; i < 6 ; i++) {
		p = Statistics::findBimodalPartition(2, f, s, a, a+i);
		if (p != a+i)
			std::cout << i << ":" << f.mean << " " << f.stdDev << " " << f.count << " : " << s.mean << " " << s.stdDev << " " << s.count << std::endl;
		BOOST_CHECK_EQUAL( p, a+i);
	}
	// No bimodal partition in last 4 values with 2 sigma
	for(int i = 6; i <= 10 ; i++) {
		p = Statistics::findBimodalPartition(2, f, s, a+6, a+i);
		if (p != a+i)
			std::cout << i << ":" << f.mean << " " << f.stdDev << " " << f.count << " : " << s.mean << " " << s.stdDev << " " << s.count << std::endl;
		BOOST_CHECK_EQUAL( p, a+i);
	}

	// detect bimodal partition with 2 sigma but not 50 sigma
	for(int i = 0; i < 6 ; i++) {
		for(int j = 7; j <= 10; j++) {
			if (j-i >= 3) {
				p = Statistics::findBimodalPartition(2, f, s, a+i, a+j);
				if (p != a+6)
					std::cout << i << ":" << j << " " << f.mean << " " << f.stdDev << " " << f.count << " : " << s.mean << " " << s.stdDev << " " << s.count << std::endl;
				BOOST_CHECK_EQUAL( p, a+6);
				p = Statistics::findBimodalPartition(50, f, s, a+i, a+j);
				if (p != a+j)
					std::cout << i << ":" << j << " " << f.mean << " " << f.stdDev << " " << f.count << " : " << s.mean << " " << s.stdDev << " " << s.count << std::endl;
				BOOST_CHECK_EQUAL( p, a+j);
			}
		}
	}

	p = Statistics::findBimodalPartition(0, f, s, a, a+10);
	BOOST_CHECK_EQUAL( p, a+6);
	BOOST_CHECK_EQUAL( f.count, 6);
	BOOST_CHECK_EQUAL( f.mean, (a[0]+a[1]+a[2]+a[3]+a[4]+a[5]) / 6.0);
	BOOST_CHECK_EQUAL( s.count, 4);
	BOOST_CHECK_EQUAL( s.mean, (a[6]+a[7]+a[8]+a[9]) / 4.0);

	p = Statistics::findBimodalPartition(1, f, s, a, a+10);
	BOOST_CHECK_EQUAL( p, a+6);
	BOOST_CHECK_EQUAL( f.count, 6);
	BOOST_CHECK_EQUAL( f.mean, (a[0]+a[1]+a[2]+a[3]+a[4]+a[5]) / 6.0);
	BOOST_CHECK_EQUAL( s.count, 4);
	BOOST_CHECK_EQUAL( s.mean, (a[6]+a[7]+a[8]+a[9]) / 4.0);

	p = Statistics::findBimodalPartition(2, f, s, a, a+10);
	BOOST_CHECK_EQUAL( p, a+6);
	BOOST_CHECK_EQUAL( f.count, 6);
	BOOST_CHECK_EQUAL( f.mean, (a[0]+a[1]+a[2]+a[3]+a[4]+a[5]) / 6.0);
	BOOST_CHECK_EQUAL( s.count, 4);
	BOOST_CHECK_EQUAL( s.mean, (a[6]+a[7]+a[8]+a[9]) / 4.0);

	p = Statistics::findBimodalPartition(3, f, s, a, a+10);
	BOOST_CHECK_EQUAL( p, a+6);
	BOOST_CHECK_EQUAL( f.count, 6);
	BOOST_CHECK_EQUAL( f.mean, (a[0]+a[1]+a[2]+a[3]+a[4]+a[5]) / 6.0);
	BOOST_CHECK_EQUAL( s.count, 4);
	BOOST_CHECK_EQUAL( s.mean, (a[6]+a[7]+a[8]+a[9]) / 4.0);


}

void testIpipestream()
{

	for(int i = 1 ; i < 4000 ; i++) {
		IPipestream p("head -" + boost::lexical_cast<std::string>(i) + " 1000.fastq");
		int c = 0;
		std::string line;
		while (p.good() && !p.eof()) {
			getline(p, line);
			if (line.length() > 0)
				c++;
		}
		p.close();
		BOOST_CHECK_EQUAL(c, i);
	}
}


BOOST_AUTO_TEST_CASE( MmapTempfileTest )
{
	testBimodalPartition();
	testIpipestream();
}
