//
// Kmernator/test/TwoBitSequenceTest.cpp
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

#include "TwoBitSequence.h"
#define BOOST_TEST_MODULE TwoBitSequenceTest
#include <boost/test/unit_test.hpp>

const char fastaA[] = "AAAA";
const char fastaC[] = "CCCC";
const char fastaG[] = "GGGG";
const char fastaT[] = "TTTT";

const char fasta1[] = "ACGTCGTAGTACTACG";
const char rev1[]   = "CGTAGTACTACGACGT";
const char fasta2[] = "ACGTCGTAGTACTACGA";
const char rev2[]   = "TCGTAGTACTACGACGT";
const char fasta3[] = "CACGTCGTAGTACTACGA";
const char rev3[]   = "TCGTAGTACTACGACGTG";
const char fasta4[] = "GCACGTCGTAGTACTACGA";
const char rev4[]   = "TCGTAGTACTACGACGTGC";
const char n1[] = "NCGTAGTACTACGACGTGCC";
const char n2[] = "NCGTANTACTACGACGTGCCA";
const char n3[] = "NCGTANTACTNCGACGTGCCAG";
const char n4[] = "NCGTANTACTNCGACNTGCCAGT";
const char n5[] = "NCGTANTACTNCGACNTGCCAGTN";

TwoBitEncoding in[1024], out[1024], out1[1024], out2[1024], out3[1024], test[1024];
char fasta[1024];
SequenceLengthType twoBitLength, sequenceLength;


//static void permuteBase(const TwoBitEncoding *in, TwoBitEncoding *out1, TwoBitEncoding *out2, TwoBitEncoding *out3, SequenceLengthType sequenceLength, SequenceLengthType permuteBaseIdx);

#define PERMUTE_BASE(inFasta, targetFasta1, targetFasta2, targetFasta3, baseIdx) \
		sequenceLength = std::strlen(inFasta);\
		twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(sequenceLength);\
		TwoBitSequence::compressSequence(inFasta, in);\
		TwoBitSequence::permuteBase(in, out1, out2, out3, sequenceLength, baseIdx);\
		TwoBitSequence::uncompressSequence(out1, sequenceLength, fasta);\
		BOOST_CHECK_EQUAL(targetFasta1, fasta);\
		TwoBitSequence::uncompressSequence(out2, sequenceLength, fasta);\
		BOOST_CHECK_EQUAL(targetFasta2, fasta);\
		TwoBitSequence::uncompressSequence(out3, sequenceLength, fasta);\
		BOOST_CHECK_EQUAL(targetFasta3, fasta);

void testPermuteBase() {
	PERMUTE_BASE("A", "C", "G", "T", 0);
	PERMUTE_BASE("AA", "CA", "GA", "TA", 0);
	PERMUTE_BASE("AC", "CC", "GC", "TC", 0);
	PERMUTE_BASE("AG", "CG", "GG", "TG", 0);
	PERMUTE_BASE("AT", "CT", "GT", "TT", 0);

	PERMUTE_BASE("AA", "AC", "AG", "AT", 1);
	PERMUTE_BASE("AC", "AA", "AG", "AT", 1);
	PERMUTE_BASE("AG", "AA", "AC", "AT", 1);
	PERMUTE_BASE("AT", "AA", "AC", "AG", 1);

	PERMUTE_BASE("ACGT", "CCGT", "GCGT", "TCGT", 0);
	PERMUTE_BASE("ACGT", "AAGT", "AGGT", "ATGT", 1);
	PERMUTE_BASE("ACGT", "ACAT", "ACCT", "ACTT", 2);
	PERMUTE_BASE("ACGT", "ACGA", "ACGC", "ACGG", 3);

	PERMUTE_BASE("ACGTT", "CCGTT", "GCGTT", "TCGTT", 0);
	PERMUTE_BASE("ACGTG", "AAGTG", "AGGTG", "ATGTG", 1);
	PERMUTE_BASE("ACGTC", "ACATC", "ACCTC", "ACTTC", 2);
	PERMUTE_BASE("ACGTA", "ACGAA", "ACGCA", "ACGGA", 3);
	PERMUTE_BASE("ACGTT", "ACGTA", "ACGTC", "ACGTG", 4);

}

//static void shiftLeft(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType twoBitLength, unsigned char shiftAmountInBases);

#define LEFT_SHIFT(inFasta,targetFasta,shiftAmount)\
		sequenceLength = std::strlen(inFasta);\
		twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(sequenceLength);\
		TwoBitSequence::compressSequence(inFasta, in);\
		TwoBitSequence::shiftLeft(in,out,twoBitLength,shiftAmount);\
		TwoBitSequence::uncompressSequence(out, sequenceLength, fasta);\
		BOOST_CHECK_EQUAL(targetFasta,fasta);

void testLeftShift() {
	LEFT_SHIFT("ACGT","ACGT",0);
	LEFT_SHIFT("ACGT", "CGTA",1);
	LEFT_SHIFT("ACGT", "GTAA",2);
	LEFT_SHIFT("ACGT", "TAAA",3);

	LEFT_SHIFT("AAAAAAAA", "AAAAAAAA",0);
	LEFT_SHIFT("AAAAAAAA", "AAAAAAAA",1);
	LEFT_SHIFT("AAAAAAAA", "AAAAAAAA",2);
	LEFT_SHIFT("AAAAAAAA", "AAAAAAAA",3);

	LEFT_SHIFT("CCCCCCCC", "CCCCCCCA",1);
	LEFT_SHIFT("CCCCCCCC", "CCCCCCAA",2);
	LEFT_SHIFT("CCCCCCCC", "CCCCCAAA",3);

	LEFT_SHIFT("CCCCCCC",  "CCCCCCA",1);
	LEFT_SHIFT("CCCCCCC",  "CCCCCAA",2);
	LEFT_SHIFT("CCCCCCC",  "CCCCAAA",3);

	LEFT_SHIFT("CCCCCC",  "CCCCCA",1);
	LEFT_SHIFT("CCCCCC",  "CCCCAA",2);
	LEFT_SHIFT("CCCCCC",  "CCCAAA",3);

	LEFT_SHIFT("CCCCC",  "CCCCA",1);
	LEFT_SHIFT("CCCCC",  "CCCAA",2);
	LEFT_SHIFT("CCCCC",  "CCAAA",3);

	LEFT_SHIFT("ACGTACGT","CGTACGTA",1);
	LEFT_SHIFT("ACGTACGT", "GTACGTAA",2);
	LEFT_SHIFT("ACGTACGT", "TACGTAAA",3);

	LEFT_SHIFT("TACGTACGT","ACGTACGTA",1);
	LEFT_SHIFT("TACGTACGT", "CGTACGTAA",2);
	LEFT_SHIFT("TACGTACGT", "GTACGTAAA",3);

	LEFT_SHIFT("GTACGTACGT","TACGTACGTA",1);
	LEFT_SHIFT("GTACGTACGT", "ACGTACGTAA",2);
	LEFT_SHIFT("GTACGTACGT", "CGTACGTAAA",3);

	LEFT_SHIFT("CGTACGTACGT","GTACGTACGTA",1);
	LEFT_SHIFT("CGTACGTACGT", "TACGTACGTAA",2);
	LEFT_SHIFT("CGTACGTACGT", "ACGTACGTAAA",3);

}

#define REV_COMP(fwd,rev) \
		TwoBitSequence::compressSequence(fwd, in);\
		TwoBitSequence::compressSequence(rev, test); \
		sequenceLength = std::strlen(fwd);\
		twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(sequenceLength);\
		TwoBitSequence::reverseComplement(in,out,sequenceLength);\
		TwoBitSequence::reverseComplement(out,out2,sequenceLength); \
		TwoBitSequence::reverseComplement(test,out3,sequenceLength); \
		TwoBitSequence::uncompressSequence(out, sequenceLength, fasta);\
		BOOST_CHECK_EQUAL(rev,fasta); \
		TwoBitSequence::uncompressSequence(out2, sequenceLength, fasta);\
		BOOST_CHECK_EQUAL(fwd, fasta); \
		TwoBitSequence::uncompressSequence(out3, sequenceLength, fasta);\
		BOOST_CHECK_EQUAL(fwd, fasta);

void testReverseComplement() {

	REV_COMP("A", "T");
	REV_COMP("C", "G");
	REV_COMP("AA", "TT");
	REV_COMP("CC", "GG");
	REV_COMP("AAA", "TTT");
	REV_COMP("CCC", "GGG");
	REV_COMP("AAAA", "TTTT");
	REV_COMP("CCCC", "GGGG");
	REV_COMP("AAAAA", "TTTTT");
	REV_COMP("CCCCC", "GGGGG");
	REV_COMP("AAAAAA", "TTTTTT");
	REV_COMP("CCCCCC", "GGGGGG");
	REV_COMP("AAAAAAA", "TTTTTTT");
	REV_COMP("CCCCCCC", "GGGGGGG");
	REV_COMP("AAAAAAAA", "TTTTTTTT");
	REV_COMP("CCCCCCCC", "GGGGGGGG");
	REV_COMP("AAAAAAAAA", "TTTTTTTTT");
	REV_COMP("CCCCCCCCC", "GGGGGGGGG");
	REV_COMP("AAAAAAAAAA", "TTTTTTTTTT");
	REV_COMP("CCCCCCCCCC", "GGGGGGGGGG");
	REV_COMP("AAAAAAAAAAA", "TTTTTTTTTTT");
	REV_COMP("CCCCCCCCCCC", "GGGGGGGGGGG");
	REV_COMP("AAAAAAAAAAAA", "TTTTTTTTTTTT");
	REV_COMP("CCCCCCCCCCCC", "GGGGGGGGGGGG");

	REV_COMP(fastaA,fastaT);
	REV_COMP(fastaC,fastaG);
	REV_COMP(fastaG,fastaC);
	REV_COMP(fastaT,fastaA);

	REV_COMP(fasta1,rev1);
	REV_COMP(fasta2,rev2);
	REV_COMP(fasta3,rev3);
	REV_COMP(fasta4,rev4);

	REV_COMP("ACGT", "ACGT");
	REV_COMP("AGCT", "AGCT");
	REV_COMP("GCTA", "TAGC");
	REV_COMP("TTAA", "TTAA");
	REV_COMP("ATAT", "ATAT");
	REV_COMP("TACC", "GGTA");
	REV_COMP("A", "T");
	REV_COMP("T", "A");
	REV_COMP("C", "G");
	REV_COMP("G", "C");
	REV_COMP("AA", "TT");
	REV_COMP("AC", "GT");
	REV_COMP("AG", "CT");
	REV_COMP("AT", "AT");
	REV_COMP("CA", "TG");
	REV_COMP("CC", "GG");
	REV_COMP("CG", "CG");
	REV_COMP("CT", "AG");
	REV_COMP("GA", "TC");
	REV_COMP("GC", "GC");
	REV_COMP("GG", "CC");
	REV_COMP("GT", "AC");
	REV_COMP("TA", "TA");
	REV_COMP("TC", "GA");
	REV_COMP("TG", "CA");
	REV_COMP("TT", "AA");


}

void testCompressUncompress() {
	TwoBitEncoding bin[1024];
	char fasta[1024];

	BaseLocationVectorType markup;
	markup = TwoBitSequence::compressSequence(fasta1, bin);
	BOOST_CHECK_EQUAL(markup.size(), 0ul);

	TwoBitSequence::uncompressSequence(bin, std::strlen(fasta1), fasta);
	BOOST_CHECK_EQUAL(fasta1, fasta);

	markup = TwoBitSequence::compressSequence(fasta2, bin);
	BOOST_CHECK_EQUAL(markup.size(), 0ul);

	TwoBitSequence::uncompressSequence(bin, std::strlen(fasta2), fasta);
	BOOST_CHECK_EQUAL(fasta2, fasta);

	markup = TwoBitSequence::compressSequence(fasta3, bin);
	BOOST_CHECK_EQUAL(markup.size(), 0ul);

	TwoBitSequence::uncompressSequence(bin, std::strlen(fasta3), fasta);
	BOOST_CHECK_EQUAL(fasta3, fasta);

	markup = TwoBitSequence::compressSequence(fasta4, bin);
	BOOST_CHECK_EQUAL(markup.size(), 0ul);

	TwoBitSequence::uncompressSequence(bin, std::strlen(fasta4), fasta);
	BOOST_CHECK_EQUAL(fasta4, fasta);

}

void testMarkup() {
	TwoBitEncoding bin[1024];
	std::string fasta(1024, '\0');

	BaseLocationVectorType markup;
	markup = TwoBitSequence::compressSequence(n1, bin);
	BOOST_CHECK_EQUAL(markup.size(), 1ul);
	BOOST_CHECK_EQUAL(markup[0].first, 'N');
	BOOST_CHECK_EQUAL(markup[0].second, 0ul);
	TwoBitSequence::uncompressSequence(bin, std::strlen(n1), fasta);
	TwoBitSequence::applyMarkup(fasta, markup);
	BOOST_CHECK_EQUAL(n1, fasta);

	markup = TwoBitSequence::compressSequence(n2, bin);
	BOOST_CHECK_EQUAL(markup.size(), 2ul);
	BOOST_CHECK_EQUAL(markup[0].first, 'N');
	BOOST_CHECK_EQUAL(markup[0].second, 0ul);
	BOOST_CHECK_EQUAL(markup[1].first, 'N');
	BOOST_CHECK_EQUAL(markup[1].second, 5ul);
	TwoBitSequence::uncompressSequence(bin, std::strlen(n2), fasta);
	TwoBitSequence::applyMarkup(fasta, markup);
	BOOST_CHECK_EQUAL(n2, fasta);

	markup = TwoBitSequence::compressSequence(n3, bin);
	BOOST_CHECK_EQUAL(markup.size(), 3ul);
	BOOST_CHECK_EQUAL(markup[0].first, 'N');
	BOOST_CHECK_EQUAL(markup[0].second, 0ul);
	BOOST_CHECK_EQUAL(markup[1].first, 'N');
	BOOST_CHECK_EQUAL(markup[1].second, 5ul);
	BOOST_CHECK_EQUAL(markup[2].first, 'N');
	BOOST_CHECK_EQUAL(markup[2].second, 10ul);
	TwoBitSequence::uncompressSequence(bin, std::strlen(n3), fasta);
	TwoBitSequence::applyMarkup(fasta, markup);
	BOOST_CHECK_EQUAL(n3, fasta);

	markup = TwoBitSequence::compressSequence(n4, bin);
	BOOST_CHECK_EQUAL(markup.size(), 4ul);
	BOOST_CHECK_EQUAL(markup[0].first, 'N');
	BOOST_CHECK_EQUAL(markup[0].second, 0ul);
	BOOST_CHECK_EQUAL(markup[1].first, 'N');
	BOOST_CHECK_EQUAL(markup[1].second, 5ul);
	BOOST_CHECK_EQUAL(markup[2].first, 'N');
	BOOST_CHECK_EQUAL(markup[2].second, 10ul);
	BOOST_CHECK_EQUAL(markup[3].first, 'N');
	BOOST_CHECK_EQUAL(markup[3].second, 15ul);
	TwoBitSequence::uncompressSequence(bin, std::strlen(n4), fasta);
	TwoBitSequence::applyMarkup(fasta, markup);
	BOOST_CHECK_EQUAL(n4, fasta);

	markup = TwoBitSequence::compressSequence(n5, bin);
	BOOST_CHECK_EQUAL(markup.size(), 5ul);
	BOOST_CHECK_EQUAL(markup[0].first, 'N');
	BOOST_CHECK_EQUAL(markup[0].second, 0ul);
	BOOST_CHECK_EQUAL(markup[1].first, 'N');
	BOOST_CHECK_EQUAL(markup[1].second, 5ul);
	BOOST_CHECK_EQUAL(markup[2].first, 'N');
	BOOST_CHECK_EQUAL(markup[2].second, 10ul);
	BOOST_CHECK_EQUAL(markup[3].first, 'N');
	BOOST_CHECK_EQUAL(markup[3].second, 15ul);
	BOOST_CHECK_EQUAL(markup[4].first, 'N');
	BOOST_CHECK_EQUAL(markup[4].second, 23ul);
	TwoBitSequence::uncompressSequence(bin, std::strlen(n5), fasta);
	TwoBitSequence::applyMarkup(fasta, markup);
	BOOST_CHECK_EQUAL(n5, fasta);

}

#define TEST_GC(inFasta,count)\
		sequenceLength = std::strlen(inFasta);\
		TwoBitSequence::compressSequence(inFasta, in);\
		BOOST_CHECK_EQUAL((long)TwoBitSequence::getGC(in,sequenceLength),(long)count);


void testGC()
{
	TEST_GC("A", 0);
	TEST_GC("C", 1);
	TEST_GC("G", 1);
	TEST_GC("T", 0);
	TEST_GC("AC", 1);
	TEST_GC("CG", 2);
	TEST_GC("GT", 1);
	TEST_GC("TA", 0);
	TEST_GC("ACG", 2);
	TEST_GC("CGT", 2);
	TEST_GC("GTA", 1);
	TEST_GC("TAC", 1);
	TEST_GC("ACGT", 2);
	TEST_GC("CGTA", 2);
	TEST_GC("GTAC", 2);
	TEST_GC("TACG", 2);
	TEST_GC("ACGTA", 2);
	TEST_GC("CGTAC", 3);
	TEST_GC("GTACG", 3);
	TEST_GC("TACGT", 2);

}

BOOST_AUTO_TEST_CASE( TwoBitSequenceTest )
{
	testLeftShift();
	testCompressUncompress();
	testMarkup();
	testReverseComplement();
	testPermuteBase();
	testGC();
}
