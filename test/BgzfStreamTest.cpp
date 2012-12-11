//
// Kmernator/test/BgzfStreamTest.cpp
//
// Author: Rob Egan
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

#include "BgzfStream.h"
#define BOOST_TEST_MODULE BgzfStreamTest
#include <boost/test/unit_test.hpp>
#include <string>
#include <ios>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


void testGzip(std::string test)
{
	std::ostringstream comp_oss;
	boost::iostreams::gzip_compressor comp;
	boost::iostreams::gzip_decompressor decomp, decomp2;

	boost::iostreams::filtering_ostream foss;
	foss.push(comp);
	foss.push(comp_oss);
	foss << test;
	foss.reset();

	//std::cerr << "compressed: " << comp_oss.tellp() << " vs " << test.length() << std::endl;
	std::string s = comp_oss.str();


	//std::cerr << s << std::endl;

	//std::cerr << "compressed string " << s.length() << std::endl;

	std::istringstream decomp_iss(s);
	boost::iostreams::filtering_istream fiss_decomp;
	fiss_decomp.push(decomp);
	fiss_decomp.push(decomp_iss);

	std::ostringstream ss;
	boost::iostreams::copy(fiss_decomp, ss);
	//std::cerr << "decompressed: " << ss.tellp() << " vs " << s.length() << std::endl;
	std::string s_decomp = ss.str();

	BOOST_CHECK_EQUAL( strcmp(test.c_str(), s_decomp.c_str()),0 );

	// pack it on again, making two concatenated blocks...
	std::string s2 = s + s;
	std::string test2 = test + test;

	std::istringstream decomp_iss2(s2);
	boost::iostreams::filtering_istream fiss_decomp2;
	fiss_decomp2.push(decomp2);
	fiss_decomp2.push(decomp_iss2);

	std::ostringstream ss2;
	boost::iostreams::copy(fiss_decomp2, ss2);
	//std::cerr << "decompressed: " << ss.tellp() << " vs " << s.length() << std::endl;
	std::string s2_decomp = ss2.str();

	BOOST_CHECK_EQUAL( strcmp(test2.c_str(), s2_decomp.c_str()),0 );

}

void testBGZFIO(std::string test, bool withEOF) {
	std::ostringstream oss, oss2;
	{
		bgzf_ostream bgzf_o(oss, withEOF);
		bgzf_o << test;
	}
	std::string s = oss.str();
	{
		std::istringstream iss(s);
		bgzf_istream bgzf_i(iss);
		boost::iostreams::copy(bgzf_i, oss2);
	}
	std::string s2 = oss2.str();
	BOOST_CHECK_EQUAL(test.length(), s2.length());
	BOOST_CHECK_EQUAL(strcmp(test.c_str(), s2.c_str()), 0);
}
void testBGZF(std::string test, bool withEOF)
{
	std::ostringstream bgzf_oss, decomp_oss;
	bgzf_compressor bgzfout_comp(withEOF);
	//bgzfout_comp.setAddEOFBlock(withEOF);
	boost::iostreams::filtering_ostream foss;
	foss.push(bgzfout_comp);
	foss.push(bgzf_oss);
	foss << test;
	foss.reset();

	//std::cerr << bgzf_oss.tellp() << std::endl;
	std::string s = bgzf_oss.str();
	//std::cerr << s << std::endl;
	if (false) {
		std::ofstream x("/tmp/test.bgzf");
		x << s;
		x.close();
	}

	// decompress the BGZF stream
	bgzf_decompressor decomp;
	std::istringstream decomp_iss(s);
	boost::iostreams::filtering_istream fiss_decomp;
	fiss_decomp.push(decomp);
	fiss_decomp.push(decomp_iss);


	std::ostringstream ss;
	boost::iostreams::copy(fiss_decomp, ss);
	//std::cerr << "decompressed: " << ss.tellp() << " vs " << s.length() << std::endl;
	std::string s2 = ss.str();

	//std::cerr << s2 << std::endl;
	BOOST_CHECK_EQUAL( test.length(), s2.length() );
	BOOST_CHECK_EQUAL( strcmp(test.c_str(), s2.c_str()),0 );

	/*
	// BGZF is fully gzip compliant
	boost::iostreams::gzip_decompressor decomp;
	std::istringstream decomp_iss(s);
	boost::iostreams::filtering_istream fiss_decomp;
	fiss_decomp.push(decomp);
	fiss_decomp.push(decomp_iss);


	std::ostringstream ss;
	boost::iostreams::copy(fiss_decomp, ss);
	std::cerr << "decompressed: " << ss.tellp() << " vs " << s.length() << std::endl;
	std::string s2 = ss.str();

	std::cerr << s2 << std::endl;
	BOOST_CHECK_EQUAL( strcmp(test.c_str(), s2.c_str()),0 );

	 */

}


std::string generateSomething(int length) {
	std::ostringstream oss;
	for(int i = 0; i < length; i++)
		oss << (char) ((i%64) + 33);
	return oss.str();
}

BOOST_AUTO_TEST_CASE( BgzfStreamTest )
{
	testGzip("Testing");
	testBGZF("Testing", true);
	testBGZF("Testing", false);
	testBGZFIO("Testing", true);
	testBGZFIO("Testing", false);
	int j = 1023;
	while (j < 256*1024) {
		std::string str = generateSomething(j);
		testGzip(str);
		testBGZF(str, true);
		testBGZF(str, false);
		testBGZFIO(str, true);
		testBGZFIO(str, false);
		j = j + 8192;
	}

}
