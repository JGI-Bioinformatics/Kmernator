//
// Kmernator/test/BgzfStreamTest.cpp
//
// Author: Rob Egan
//
// Copyright 2012 The Regents of the University of California.
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

void testBGZF(std::string test, bool withEOF)
{
	std::ostringstream bgzf_oss, decomp_oss;
	bgzf_compressor bgzfout_comp;
	bgzfout_comp.setAddEOFBlock(withEOF);
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
	int j = 1023;
	while (j < 256*1024) {
		std::string str = generateSomething(j);
		testGzip(str);
		testBGZF(str, true);
		testBGZF(str, false);
		j = j + 8192;
	}

}
