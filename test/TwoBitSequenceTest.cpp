// $Header: /repository/PI_annex/robsandbox/KoMer/test/TwoBitSequenceTest.cpp,v 1.1 2009-10-22 07:04:03 regan Exp $
//

#include "TwoBitSequence.h"
#define BOOST_TEST_MODULE TwoBitSequenceTest
#include <boost/test/unit_test.hpp>

const char fasta1[] = "ACGTCGTAGTACTACG";
const char rev1[]   = "CGTAGTACTACGACGT";
const char fasta2[] = "ACGTCGTAGTACTACGA";
const char rev2[]   = "TCGTAGTACTACGACGT";
const char fasta3[] = "CACGTCGTAGTACTACGA";
const char rev3[]   = "TCGTAGTACTACGACGTG";
const char fasta4[] = "GCACGTCGTAGTACTACGA";
const char rev4[]   = "TCGTAGTACTACGACGTGC";
const char n1[]     = "NCGTAGTACTACGACGTGCC";
const char n2[]     = "NCGTANTACTACGACGTGCCA";
const char n3[]     = "NCGTANTACTNCGACGTGCCAG";
const char n4[]     = "NCGTANTACTNCGACNTGCCAGT";
const char n5[]     = "NCGTANTACTNCGACNTGCCAGTN";

void testCompressUncompress()
{
  TwoBitEncoding bin[1024];
  char fasta[1024];
  
  BaseLocationVectorType markup;
  markup = TwoBitSequence::compressSequence(fasta1, bin);
  BOOST_CHECK_EQUAL( markup.size(), 0 );
  
  TwoBitSequence::uncompressSequence(bin, std::strlen(fasta1), fasta);
  BOOST_CHECK_EQUAL( fasta1, fasta );

  markup = TwoBitSequence::compressSequence(fasta2, bin);
  BOOST_CHECK_EQUAL( markup.size(), 0 );
  
  TwoBitSequence::uncompressSequence(bin, std::strlen(fasta2), fasta);
  BOOST_CHECK_EQUAL( fasta2, fasta );
  
  markup = TwoBitSequence::compressSequence(fasta3, bin);
  BOOST_CHECK_EQUAL( markup.size(), 0 );
  
  TwoBitSequence::uncompressSequence(bin, std::strlen(fasta3), fasta);
  BOOST_CHECK_EQUAL( fasta3, fasta );

  markup = TwoBitSequence::compressSequence(fasta4, bin);
  BOOST_CHECK_EQUAL( markup.size(), 0 );
  
  TwoBitSequence::uncompressSequence(bin, std::strlen(fasta4), fasta);
  BOOST_CHECK_EQUAL( fasta4, fasta );
  
}

void testMarkup()
{
  TwoBitEncoding bin[1024];
  char fasta[1024];
  
  BaseLocationVectorType markup;
  markup = TwoBitSequence::compressSequence(n1, bin);
  BOOST_CHECK_EQUAL( markup.size(), 1 );
  BOOST_CHECK_EQUAL( markup[0].first, 'N');
  BOOST_CHECK_EQUAL( markup[0].second, 0);
  TwoBitSequence::uncompressSequence(bin, std::strlen(n1), fasta);
  TwoBitSequence::applyMarkup(fasta, markup);
  BOOST_CHECK_EQUAL( n1, fasta );
  
  markup = TwoBitSequence::compressSequence(n2, bin);
  BOOST_CHECK_EQUAL( markup.size(), 2 );
  BOOST_CHECK_EQUAL( markup[0].first, 'N');
  BOOST_CHECK_EQUAL( markup[0].second, 0);
  BOOST_CHECK_EQUAL( markup[1].first, 'N');
  BOOST_CHECK_EQUAL( markup[1].second, 5);
  TwoBitSequence::uncompressSequence(bin, std::strlen(n2), fasta);
  TwoBitSequence::applyMarkup(fasta, markup);
  BOOST_CHECK_EQUAL( n2, fasta );
  
  markup = TwoBitSequence::compressSequence(n3, bin);
  BOOST_CHECK_EQUAL( markup.size(), 3 );
  BOOST_CHECK_EQUAL( markup[0].first, 'N');
  BOOST_CHECK_EQUAL( markup[0].second, 0);
  BOOST_CHECK_EQUAL( markup[1].first, 'N');
  BOOST_CHECK_EQUAL( markup[1].second, 5);
  BOOST_CHECK_EQUAL( markup[2].first, 'N');
  BOOST_CHECK_EQUAL( markup[2].second, 10);
  TwoBitSequence::uncompressSequence(bin, std::strlen(n3), fasta);
  TwoBitSequence::applyMarkup(fasta, markup);
  BOOST_CHECK_EQUAL( n3, fasta );
  
  markup = TwoBitSequence::compressSequence(n4, bin);
  BOOST_CHECK_EQUAL( markup.size(), 4 );
  BOOST_CHECK_EQUAL( markup[0].first, 'N');
  BOOST_CHECK_EQUAL( markup[0].second, 0);
  BOOST_CHECK_EQUAL( markup[1].first, 'N');
  BOOST_CHECK_EQUAL( markup[1].second, 5);
  BOOST_CHECK_EQUAL( markup[2].first, 'N');
  BOOST_CHECK_EQUAL( markup[2].second, 10);
  BOOST_CHECK_EQUAL( markup[3].first, 'N');
  BOOST_CHECK_EQUAL( markup[3].second, 15);
  TwoBitSequence::uncompressSequence(bin, std::strlen(n4), fasta);
  TwoBitSequence::applyMarkup(fasta, markup);
  BOOST_CHECK_EQUAL( n4, fasta );
  
  markup = TwoBitSequence::compressSequence(n5, bin);
  BOOST_CHECK_EQUAL( markup.size(), 5 );
  BOOST_CHECK_EQUAL( markup[0].first, 'N');
  BOOST_CHECK_EQUAL( markup[0].second, 0);
  BOOST_CHECK_EQUAL( markup[1].first, 'N');
  BOOST_CHECK_EQUAL( markup[1].second, 5);
  BOOST_CHECK_EQUAL( markup[2].first, 'N');
  BOOST_CHECK_EQUAL( markup[2].second, 10);
  BOOST_CHECK_EQUAL( markup[3].first, 'N');
  BOOST_CHECK_EQUAL( markup[3].second, 15);
  BOOST_CHECK_EQUAL( markup[4].first, 'N');
  BOOST_CHECK_EQUAL( markup[4].second, 23);
  TwoBitSequence::uncompressSequence(bin, std::strlen(n5), fasta);
  TwoBitSequence::applyMarkup(fasta, markup);
  BOOST_CHECK_EQUAL( n5, fasta );
  
}

BOOST_AUTO_TEST_CASE( TwoBitSequenceTest )
{
	testCompressUncompress();
	testMarkup();
}

//
// $Log: TwoBitSequenceTest.cpp,v $
// Revision 1.1  2009-10-22 07:04:03  regan
// added a few unit tests
// minor refactor
//
//
