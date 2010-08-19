//
// Kmernator/test/ReadSetTest.cpp
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

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "ReadFileReader.h"
#include "Utils.h"
#define BOOST_TEST_MODULE ReadSetTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

void testZeroReads() {
	Read s;

	BOOST_CHECK_EQUAL(0ul, s.getLength());
	BOOST_CHECK_EQUAL("", s.getQuals());
	BOOST_CHECK_EQUAL("", s.getFasta());
	BOOST_CHECK_EQUAL("", s.getName());
}

static string getFileContents(string filename) {
	ifstream ifs(filename.c_str());
	BOOST_REQUIRE(!ifs.fail());
	string contents;

	getline(ifs, contents, '\0');
	return contents;
}

void testFastQFile(string filename) {
	string fileContents = getFileContents(filename);

	ReadSet store;
	store.appendAnyFile(filename);
	store.identifyPairs();

	BOOST_CHECK_EQUAL(store.getSize() / 2, store.getPairSize());

	string fastq;
	for (unsigned int i = 0; i < store.getSize(); i++) {
		fastq += store.getRead(i).toFastq();
	}

	BOOST_CHECK_EQUAL(fileContents, fastq);

	for (unsigned int i = 0; i < store.getPairSize(); i++) {
		ReadSet::Pair &pair = store.getPair(i);
		BOOST_CHECK(store.isValidRead(pair.read1) == true);
		BOOST_CHECK(store.isValidRead(pair.read2) == true);
	}
}

void testFastaWithQualFile(string f, string q) {
	string fFile = getFileContents(f);
	string qFile = getFileContents(q);

	ReadSet store;
	store.appendAnyFile(f, q);

	store.identifyPairs();

	BOOST_CHECK_EQUAL(store.getSize() / 2, store.getPairSize());

	string fasta, qual;
	for (unsigned int i = 0; i < store.getSize(); i++) {
		const Read s = store.getRead(i);

		string nameLine('>' + s.getName() + "\n");
		fasta += nameLine;
		fasta += s.getFasta() + "\n";

		qual += nameLine;
		qual += s.getFormattedQuals() + "\n";

	}

	BOOST_CHECK_EQUAL(fFile, fasta);
	BOOST_CHECK_EQUAL(qFile, qual);

	for (unsigned int i = 0; i < store.getPairSize(); i++) {
		ReadSet::Pair &pair = store.getPair(i);
		BOOST_CHECK(store.isValidRead(pair.read1) == true);
		BOOST_CHECK(store.isValidRead(pair.read2) == true);
	}

}

#define TEST_TRIM_NAME_PARSER(_target,_test) \
	target = _target; test = _test; \
	SequenceRecordParser::trimName(test); \
	BOOST_CHECK_EQUAL(target,test);

void testParser()
{
	std::string target, test;

	TEST_TRIM_NAME_PARSER("asdf1234",    ">asdf1234");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  ">asdf1234/1");
	TEST_TRIM_NAME_PARSER("asdf1234/A",  ">asdf1234/A");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  ">asdf1234/1 ");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  ">asdf1234/1 blah=blah2");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  ">asdf1234/1 blah=blah2\n");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  ">asdf1234/1\tblah=blah2\n");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  ">asdf1234/1  blah=blah2\n");

	TEST_TRIM_NAME_PARSER("asdf1234",    "@asdf1234");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  "@asdf1234/1");
	TEST_TRIM_NAME_PARSER("asdf1234/A",  "@asdf1234/A");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  "@asdf1234/1 ");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  "@asdf1234/1 blah=blah2");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  "@asdf1234/1 blah=blah2\n");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  "@asdf1234/1\tblah=blah2\n");
	TEST_TRIM_NAME_PARSER("asdf1234/1",  "@asdf1234/1  blah=blah2\n");

}

void testConsensus(string filename)
{
	ReadSet store;
	store.appendAnyFile(filename);
	store.identifyPairs();
	BOOST_CHECK_EQUAL(store.getSize() / 2, store.getPairSize());

	ReadSet reads1, reads2;
	for(unsigned int i = 0 ; i < store.getSize(); i+=2) {
		reads1.append(store.getRead(i).clone());
		reads2.append(store.getRead(i+1).clone());
	}
	BOOST_CHECK_EQUAL( reads1.getCentroidRead(), 0);
	BOOST_CHECK_EQUAL( reads2.getCentroidRead(), 0);

	Read consensus1 = reads1.getConsensusRead();
	Read consensus2 = reads2.getConsensusRead();

	BOOST_CHECK_EQUAL( consensus1.getFasta(), store.getRead(0).getFasta() );
	BOOST_CHECK_EQUAL( consensus2.getFasta(), store.getRead(1).getFasta() );

}

BOOST_AUTO_TEST_CASE( ReadSetTest )
{
	testParser();

	Sequence::clearCaches();
	testZeroReads();

	Sequence::clearCaches();
	testFastQFile("10.fastq");

	Sequence::clearCaches();
	testFastaWithQualFile("10.fasta","10.qual");

	Sequence::clearCaches();
	testConsensus("consensus1.fastq");
	Sequence::clearCaches();
	testConsensus("consensus2.fastq");
	Sequence::clearCaches();
	testConsensus("consensus3.fastq");
	Sequence::clearCaches();
	testConsensus("consensus2-diff.fastq");

	Sequence::clearCaches();

}

//
// $Log: ReadSetTest.cpp,v $
// Revision 1.15  2010-05-06 21:46:51  regan
// merged changes from PerformanceTuning-20100501
//
// Revision 1.14  2010-05-05 06:28:38  regan
// merged changes from FixPairOutput-20100504
//
// Revision 1.13.4.1  2010-05-05 05:57:56  regan
// fixed pairing
// fixed name to exclude labels and comments after whitespace
// applied some performance optimizations from other branch
// created FixPair application
//
// Revision 1.13.2.1  2010-05-03 21:34:37  regan
// added consensus unit tests
//
// Revision 1.13  2010-05-01 21:57:51  regan
// merged head with serial threaded build partitioning
//
// Revision 1.12.8.2  2010-04-30 23:53:17  regan
// attempt to fix a bug.  clearing Sequence caches when it makes sense
//
// Revision 1.12.8.1  2010-04-30 21:53:59  regan
// reuse memory efficiently for cache lookups
//
// Revision 1.12  2010-03-14 16:57:09  regan
// minor refactor
//
// Revision 1.11  2010-03-04 06:37:13  regan
// bugfix
//
// Revision 1.10  2010-02-26 13:01:21  regan
// reformatted
//
// Revision 1.9  2010-02-22 14:41:18  regan
// minor bug fixes
//
// Revision 1.8  2010-01-13 07:20:10  regan
// refactored filter
// checkpoint on read picker
//
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

