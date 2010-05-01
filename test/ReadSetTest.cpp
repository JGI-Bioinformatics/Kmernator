// $Header: /repository/PI_annex/robsandbox/KoMer/test/ReadSetTest.cpp,v 1.13 2010-05-01 21:57:51 regan Exp $
//

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
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

BOOST_AUTO_TEST_CASE( ReadSetTest )
{
	Sequence::clearCaches();
	testZeroReads();

	Sequence::clearCaches();
	testFastQFile("10.fastq");

	Sequence::clearCaches();
	testFastaWithQualFile("10.fasta","10.qual");

	Sequence::clearCaches();
}

//
// $Log: ReadSetTest.cpp,v $
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

