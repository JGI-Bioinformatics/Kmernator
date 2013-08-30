//
// Kmernator/test/UtilTest.cpp
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

typedef RankVector<int> RV;
void testRankVector(int size = 100) {
	std::vector<int> v,u;
	v.reserve(size);
	u.reserve(size);

	// put 1..50 with duplicate entries
	for(int i = 0; i < size; i++) {
		v.push_back(i/2);
		u.push_back((size-i-1)/2);
	}

	RV rv(v), ru(u);
	BOOST_CHECK_EQUAL(rv.size(), v.size());
	BOOST_CHECK_EQUAL(ru.size(), u.size());
	for(int i = 0; i < size; i++) {
		// odd sized batches have last rank untied
		if ((size & 0x01) == 0x01) {
			if (i == 0) {
				BOOST_CHECK_EQUAL(rv[i], 1.5);
				BOOST_CHECK_EQUAL(ru[i], size);
			} else if (i == size - 1) {
				// last rank does not tie rv
				BOOST_CHECK_EQUAL(rv[i], size);
				BOOST_CHECK_EQUAL(ru[i], 1.5);
			} else {
				// all other ranks tie...
				BOOST_CHECK_EQUAL(rv[i], (2*(i/2)) + 1.5);
				BOOST_CHECK_EQUAL(ru[i], (2*((size-i-1)/2)) + 1.5);
			}
		} else {
			// all ranks tie
			BOOST_CHECK_EQUAL(rv[i], (2*(i/2)) + 1.5);
			BOOST_CHECK_EQUAL(ru[i], (2*((size-i-1)/2)) + 1.5);
		}
	}
	BOOST_CHECK_EQUAL( rv.getSpearmanDistance(rv), 0.0 );
	BOOST_CHECK_EQUAL( ru.getSpearmanDistance(ru), 0.0 );
	BOOST_CHECK( rv.getSpearmanDistance(ru) >= 0.99 );

	BOOST_CHECK(abs(rv.getSpearmanDistance(ru) - ru.getSpearmanDistance(rv)) < 0.01);
	BOOST_CHECK(rv.getSpearmanDistance(ru) >= 0.0);
	BOOST_CHECK(ru.getSpearmanDistance(rv) >= 0.0);
}


BOOST_AUTO_TEST_CASE( UtilTest )
{
	for(int i = 20; i < 200; i += 17)
		testRankVector(i);
	testBimodalPartition();
	testIpipestream();
}
