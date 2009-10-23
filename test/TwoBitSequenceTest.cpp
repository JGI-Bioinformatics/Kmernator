// $Header: /repository/PI_annex/robsandbox/KoMer/test/TwoBitSequenceTest.cpp,v 1.4 2009-10-23 00:13:56 cfurman Exp $
//

#include "TwoBitSequence.h"
#define BOOST_TEST_MODULE TwoBitSequenceTest
#include <boost/test/unit_test.hpp>

const char fastaA[]  = "AAAA";
const char fastaC[]  = "CCCC";
const char fastaG[]  = "GGGG";
const char fastaT[]  = "TTTT";


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


TwoBitEncoding in[1024],out[1024],test[1024];
char fasta[1024];
SequenceLengthType twoBitLength,sequenceLength;

#define BIT_SHIFT(basesIn, targetBases, baseShift)\
    TwoBitSequence::compressSequence(basesIn, in);\
    out[0] = TwoBitSequence::bitShiftTable[((unsigned long)*((unsigned short*)in))*3ul+baseShift-1];\
    TwoBitSequence::uncompressSequence(out,4,fasta);\
    BOOST_CHECK_EQUAL(targetBases,fasta);

void testBitShift()
{
   BIT_SHIFT("AAAAAAAA",fastaA,1);
   BIT_SHIFT("AAAAAAAA",fastaA,2);
   BIT_SHIFT("AAAAAAAA",fastaA,3);

   BIT_SHIFT("CCCCCCCC",fastaC,1);
   BIT_SHIFT("CCCCCCCC",fastaC,2);
   BIT_SHIFT("CCCCCCCC",fastaC,3);
      
   BIT_SHIFT("GGGGGGGG",fastaG,1);
   BIT_SHIFT("GGGGGGGG",fastaG,2);
   BIT_SHIFT("GGGGGGGG",fastaG,3);

   BIT_SHIFT("TTTTTTTT",fastaT,1);
   BIT_SHIFT("TTTTTTTT",fastaT,2);
   BIT_SHIFT("TTTTTTTT",fastaT,3);

   BIT_SHIFT("ACGTACGT","CGTA",1);
   BIT_SHIFT("ACGTACGT","GTAC",2);
   BIT_SHIFT("ACGTACGT","TACG",3);
   
}

//static void shiftLeft(const TwoBitEncoding *in, TwoBitEncoding *out, SequenceLengthType twoBitLength, unsigned char shiftAmountInBases);

#define LEFT_SHIFT(inFasta,targetFasta,shiftAmount)\
  sequenceLength = std::strlen(inFasta);\
  twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(sequenceLength);\
  TwoBitSequence::compressSequence(inFasta, in);\
  TwoBitSequence::shiftLeft(in,out,twoBitLength,shiftAmount);\
  TwoBitSequence::uncompressSequence(out, sequenceLength, fasta);\
  BOOST_CHECK_EQUAL(targetFasta,fasta);
  




void testLeftShift()
{
   LEFT_SHIFT("ACGT","ACGT",0);
   LEFT_SHIFT("ACGT", "CGTA",1);
   LEFT_SHIFT("ACGT",  "GTAA",2);
   LEFT_SHIFT("ACGT",   "TAAA",3);

   LEFT_SHIFT("AAAAAAAA","AAAAAAAA",0);
   LEFT_SHIFT("AAAAAAAA", "AAAAAAAA",1);
   LEFT_SHIFT("AAAAAAAA",  "AAAAAAAA",2);
   LEFT_SHIFT("AAAAAAAA",   "AAAAAAAA",3);

   LEFT_SHIFT("CCCCCCCC","CCCCCCCA",1);
   LEFT_SHIFT("CCCCCCCC", "CCCCCCAA",2);
   LEFT_SHIFT("CCCCCCCC",  "CCCCCAAA",3);
      
   LEFT_SHIFT("ACGTACGT","CGTACGTA",1);
   LEFT_SHIFT("ACGTACGT", "GTACGTAA",2);
   LEFT_SHIFT("ACGTACGT",  "TACGTAAA",3);

   LEFT_SHIFT("TACGTACGT","ACGTACGTA",1);
   LEFT_SHIFT("TACGTACGT", "CGTACGTAA",2);
   LEFT_SHIFT("TACGTACGT",  "GTACGTAAA",3);
   
   LEFT_SHIFT("GTACGTACGT","TACGTACGTA",1);
   LEFT_SHIFT("GTACGTACGT", "ACGTACGTAA",2);
   LEFT_SHIFT("GTACGTACGT",  "CGTACGTAAA",3);

   LEFT_SHIFT("CGTACGTACGT","GTACGTACGTA",1);
   LEFT_SHIFT("CGTACGTACGT", "TACGTACGTAA",2);
   LEFT_SHIFT("CGTACGTACGT",  "ACGTACGTAAA",3);

}

#define REV_COMP(fwd,rev) \
  TwoBitSequence::compressSequence(fwd, in);\
  TwoBitSequence::compressSequence(rev, test); \
  sequenceLength = std::strlen(fwd);\
  twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(sequenceLength);\
  TwoBitSequence::reverseComplement(in,out,sequenceLength);\
  BOOST_CHECK_EQUAL( memcmp(out,test,twoBitLength), 0);\
  TwoBitSequence::uncompressSequence(out, sequenceLength, fasta);\
  BOOST_CHECK_EQUAL(rev,fasta);



void testReverseComplement()
{

  REV_COMP(fastaA,fastaT);
  REV_COMP(fastaC,fastaG);
  REV_COMP(fastaG,fastaC);
  REV_COMP(fastaT,fastaA);
  
  REV_COMP(fasta1,rev1);
  REV_COMP(fasta2,rev2);
  REV_COMP(fasta3,rev3);
  REV_COMP(fasta4,rev4);
 
}


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
    testBitShift();
    testLeftShift();
	testCompressUncompress();
    testMarkup();
    testReverseComplement();
}

//
// $Log: TwoBitSequenceTest.cpp,v $
// Revision 1.4  2009-10-23 00:13:56  cfurman
// reverse complement now works
//
// Revision 1.3  2009-10-22 21:46:47  regan
// fixed ushort to ulong conversion problems
//
// Revision 1.2  2009-10-22 20:49:18  cfurman
// tests added
//
// Revision 1.1  2009-10-22 07:04:03  regan
// added a few unit tests
// minor refactor
//
//
