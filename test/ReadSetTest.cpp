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
#include "KmerReadUtils.h"

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

void testFastQFile(string filename, bool paired = true) {
	string fileContents = getFileContents(filename);

	ReadSet store;
	store.appendAnyFile(filename);
	store.identifyPairs();

	BOOST_CHECK_EQUAL(store.getSize() / (paired?2:1), store.getPairSize());

	string fastq;
	for (unsigned int i = 0; i < store.getSize(); i++) {
		fastq += store.getRead(i).toFastq();
	}

	if (GlobalOptions::isCommentStored())
		BOOST_CHECK_EQUAL(fileContents, fastq);

	for (unsigned int i = 0; paired && i < store.getPairSize(); i++) {
		ReadSet::Pair &pair = store.getPair(i);
		BOOST_CHECK(store.isValidRead(pair.read1) == true);
		BOOST_CHECK(store.isValidRead(pair.read2) == true);
	}
}

void testSplitFastQFile(string filename1, string filename2) {
	string fileContents1 = getFileContents(filename1);
	string fileContents2 = getFileContents(filename2);

	ReadSet store;
	store.appendAnyFile(filename1);
	store.appendAnyFile(filename2);
	store.identifyPairs();

	BOOST_CHECK(store.getSize() > 0);
	BOOST_CHECK_EQUAL(store.getSize() / 2, store.getPairSize());

	string fastq1;
	for (unsigned int i = 0; i < store.getSize() / 2; i++) {
		fastq1 += store.getRead(i).toFastq();
	}
	string fastq2;
	for (unsigned int i = store.getSize() / 2; i < store.getSize(); i++) {
		fastq2 += store.getRead(i).toFastq();
	}

	if (GlobalOptions::isCommentStored()) {
		BOOST_CHECK_EQUAL(fileContents1, fastq1);
		BOOST_CHECK_EQUAL(fileContents2, fastq2);
	}

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
		{ std::string target = _target, test = _test; \
		  std::string comment; \
		  SequenceRecordParser::trimName(test, comment); \
		  BOOST_CHECK_EQUAL(target,test); \
		}

void testParser()
{


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

	bool oldOpt = GlobalOptions::isCommentStored();
	GlobalOptions::isCommentStored() = false;

	TEST_TRIM_NAME_PARSER("asdf1234/1",  "@asdf1234 1:Y:0:A\n");
	TEST_TRIM_NAME_PARSER("asdf1234/2",  "@asdf1234 2:Y:0:A\n");

	GlobalOptions::isCommentStored() = true;

	TEST_TRIM_NAME_PARSER("asdf1234",  "@asdf1234 1:Y:0:A\n");
	TEST_TRIM_NAME_PARSER("asdf1234",  "@asdf1234 2:Y:0:A\n");

	GlobalOptions::isCommentStored() = oldOpt;
}

void testConsensus(string filename)
{
	ReadSet store;
	store.appendAnyFile(filename);
	store.identifyPairs();
	BOOST_CHECK_EQUAL(store.getSize() / 2, store.getPairSize());

	int minq = Options::getOptions().getMinQuality();
	Options::getOptions().getMinQuality() = 0;
	ReadSet reads1, reads2;
	for(unsigned int i = 0 ; i < store.getSize(); i+=2) {
		reads1.append(store.getRead(i).clone());
		reads2.append(store.getRead(i+1).clone());
	}
	BOOST_CHECK_EQUAL( reads1.getCentroidRead(), 0ul);
	BOOST_CHECK_EQUAL( reads2.getCentroidRead(), 0ul);

	Read consensus1 = reads1.getConsensusRead();
	Read consensus2 = reads2.getConsensusRead();

	BOOST_CHECK_EQUAL( consensus1.getFasta(), store.getRead(0).getFasta() );
	BOOST_CHECK_EQUAL( consensus2.getFasta(), store.getRead(1).getFasta() );
	Options::getOptions().getMinQuality() = minq;
}

void testStore(string filename) {
	ReadSet a, b, c;
	a.appendAnyFile(filename);
	long bufSize = a.getStoreSize();
	char *buf = new char[bufSize];
	a.store(buf);
	b.restore(buf);

	BOOST_CHECK_EQUAL( a.getSize(), b.getSize() );
	for(unsigned int i = 0 ; i < a.getSize(); i++) {
		const Read &reada = a.getRead(i);
		const Read &readb = a.getRead(i);

		BOOST_CHECK_EQUAL(reada.getName(), readb.getName());
		BOOST_CHECK_EQUAL(reada.getFastaNoMarkup(), readb.getFastaNoMarkup());
		BOOST_CHECK_EQUAL(reada.getQuals(), readb.getQuals());
	}

	delete [] buf;
}


void testKmerMap(SequenceLengthType size) {
	std::string
	A("ACGTCGTAACGTCGTA"),
	B("TACGACGTTACGACGT"),
	C("AAAACCCCGGGGTTTTTACGTCGTAGTACTACGAAAACCCCGGGGTTTTACGTCGTAGTACTACG");
	int oldSize = KmerSizer::getSequenceLength();
	KmerSizer::set(size);

	Read readA("A", A, "", ""), readAq("Aq", A, std::string("h",A.length()), "");
	Read readB("B", B, "", ""), readBq("Bq", B, std::string("h",B.length()), "");
	Read readC("C", C, "", ""), readCq("Cq", C, std::string("h",C.length()), "");

	KmerWeightedExtensions weights;
	weights = KmerReadUtils::buildWeightedKmers(readA);
	weights = KmerReadUtils::buildWeightedKmers(readAq);
	weights = KmerReadUtils::buildWeightedKmers(readB);
	weights = KmerReadUtils::buildWeightedKmers(readBq);
	weights = KmerReadUtils::buildWeightedKmers(readC);
	weights = KmerReadUtils::buildWeightedKmers(readCq);


	KmerSizer::set(oldSize);
}


BOOST_AUTO_TEST_CASE( ReadSetTest )
{
	testParser();

	Sequence::clearCaches();
	testZeroReads();



	{
	Sequence::clearCaches();
	bool oldOpt = GlobalOptions::isCommentStored();
	GeneralOptions::getOptions().getMmapInput() = false;
	testFastQFile("10.fastq");
	testFastQFile("10-cs18.fastq");
	testFastQFile("10-cs18.1.fastq", false);
	testFastQFile("10-cs18.2.fastq", false);
	testSplitFastQFile("10-cs18.1.fastq","10-cs18.2.fastq");
	testSplitFastQFile("10-cs18.2.fastq","10-cs18.1.fastq");
	GlobalOptions::isCommentStored() = !oldOpt;
	testFastQFile("10.fastq");
	testFastQFile("10-cs18.fastq");
	testFastQFile("10-cs18.1.fastq", false);
	testFastQFile("10-cs18.2.fastq", false);
	testSplitFastQFile("10-cs18.1.fastq","10-cs18.2.fastq");
	testSplitFastQFile("10-cs18.2.fastq","10-cs18.1.fastq");
	GlobalOptions::isCommentStored() = oldOpt;
	}

	{
	Sequence::clearCaches();
	GeneralOptions::getOptions().getMmapInput() = true;
	bool oldOpt = GlobalOptions::isCommentStored();
	testFastQFile("10.fastq");
	testFastQFile("10-cs18.fastq");
	testFastQFile("10-cs18.1.fastq", false);
	testFastQFile("10-cs18.2.fastq", false);
	testSplitFastQFile("10-cs18.1.fastq","10-cs18.2.fastq");
	testSplitFastQFile("10-cs18.2.fastq","10-cs18.1.fastq");
	GlobalOptions::isCommentStored() = !oldOpt;
	testFastQFile("10.fastq");
	testFastQFile("10-cs18.fastq");
	testFastQFile("10-cs18.1.fastq", false);
	testFastQFile("10-cs18.2.fastq", false);
	testSplitFastQFile("10-cs18.1.fastq","10-cs18.2.fastq");
	testSplitFastQFile("10-cs18.2.fastq","10-cs18.1.fastq");
	GlobalOptions::isCommentStored() = oldOpt;
	}

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
	testStore("10.fasta");
	Sequence::clearCaches();
	testStore("consensus1.fastq");
	Sequence::clearCaches();
	testStore("consensus2.fastq");
	Sequence::clearCaches();

}
