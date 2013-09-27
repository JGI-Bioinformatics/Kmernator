//
// Kmernator/test/ReadSetTest.cpp
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

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "ReadFileReader.h"
#include "Utils.h"
#include "KmerReadUtils.h"
#include "Options.h"

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
	KmerReadUtils kru;
	weights = kru.buildWeightedKmers(readA);
	weights = kru.buildWeightedKmers(readAq);
	weights = kru.buildWeightedKmers(readB);
	weights = kru.buildWeightedKmers(readBq);
	weights = kru.buildWeightedKmers(readC);
	weights = kru.buildWeightedKmers(readCq);


	KmerSizer::set(oldSize);
}


void _testKmerBuild(int kmerSize, string fasta) {
	int oldSize = KmerSizer::getSequenceLength();
	KmerSizer::set(kmerSize);
	Sequence seq(fasta);
	int len = fasta.length();
	BOOST_CHECK_EQUAL( len, (int) fasta.size() );
	BOOST_CHECK_EQUAL( len, (int) seq.getLength() );
	KmerWeightedExtensions kmers(seq.getTwoBitSequence(), len);
	int kmersLen = len - kmerSize + 1;
	BOOST_CHECK_EQUAL((int) kmers.size(), kmersLen);
	for(int i = 0; i < kmersLen ; i++)
		BOOST_CHECK_EQUAL( fasta.substr(i, kmerSize), kmers[i].toFasta());

	KmerSizer::set(oldSize);
}
void testKmerBuild() {
	string test = "ACGTCGTAGTATTACGTTTTCCCCAAAAGGGG";
	for (int p = 0; p < (int) test.length(); p++) {
		for (int l = 2; l < (int) test.length() - p - 1; l++) {
			for (int k = 1; k < l - 1; k++) {
				_testKmerBuild(k, test.substr(p,l));
			}
		}
	}
}



BOOST_AUTO_TEST_CASE( ReadSetTest )
{
	GeneralOptions::getOptions().getOutputFastqBaseQuality() = 33;
	GeneralOptions::getOptions().getFastqBaseQuality() = 64;

	testKmerBuild();

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
