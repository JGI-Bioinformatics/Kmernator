// $Header: /repository/PI_annex/robsandbox/KoMer/test/KmerTest.cpp,v 1.1 2009-10-23 07:06:57 regan Exp $
//
 

#include "TwoBitSequence.h"
#include "Kmer.h"

#define BOOST_TEST_MODULE KmerSetTest
#include <boost/test/unit_test.hpp>
#include <cstring>

using namespace std;

TwoBitEncoding twoBit1[1024],twoBit2[1024],twoBit3[1024];
Kmer    *kmer1     = (Kmer*)&twoBit1;
Kmer    *kmer2     = (Kmer*)&twoBit2;
Kmer    *kmer3     = (Kmer*)&twoBit3;

#define SET_KMERS(fasta1, fasta2, fasta3) \
    TwoBitSequence::compressSequence(fasta1,   twoBit1);\
    TwoBitSequence::compressSequence(fasta2,   twoBit2);\
    TwoBitSequence::compressSequence(fasta3,   twoBit3);\
    KmerSizer::set(std::strlen(fasta1));\
    BOOST_CHECK_EQUAL( fasta1, kmer1->toFasta());\
    BOOST_CHECK_EQUAL( fasta2, kmer2->toFasta());

void testKmerCompare()
{
  SET_KMERS("AAAA", "AAAA", "");
  BOOST_CHECK( *kmer1 == *kmer2 );
  
  SET_KMERS("AAAAA", "AAAAA", "");
  BOOST_CHECK( *kmer1 == *kmer2 );
  
  SET_KMERS("AAAAAA", "AAAAAA", "");
  BOOST_CHECK( *kmer1 == *kmer2 );
  
  SET_KMERS("AAAAAAA", "AAAAAAA", "");
  BOOST_CHECK( *kmer1 == *kmer2 );

  SET_KMERS("AAAAAAAA", "AAAAAAAA", "");
  BOOST_CHECK( *kmer1 == *kmer2 );

  SET_KMERS("ACGT", "ACGT", "");
  BOOST_CHECK( *kmer1 == *kmer2 );
  
  SET_KMERS("ACGTC", "ACGTC", "");
  BOOST_CHECK( *kmer1 == *kmer2 );
  
  SET_KMERS("ACGTCG", "ACGTCG", "");
  BOOST_CHECK( *kmer1 == *kmer2 );
  
  SET_KMERS("ACGTCGT", "ACGTCGT", "");
  BOOST_CHECK( *kmer1 == *kmer2 );

  SET_KMERS("ACGTCGTC", "ACGTCGTC", "");
  BOOST_CHECK( *kmer1 == *kmer2 );

  SET_KMERS("ACGT", "CGTA", "");
  BOOST_CHECK( *kmer1 != *kmer2 );
  
  SET_KMERS("ACGTC", "CGTCA", "");
  BOOST_CHECK( *kmer1 != *kmer2 );
  
  SET_KMERS("ACGTCG", "CGTCGA", "");
  BOOST_CHECK( *kmer1 != *kmer2 );
  
  SET_KMERS("ACGTCGT", "CGTCGTA", "");
  BOOST_CHECK( *kmer1 != *kmer2 );

  SET_KMERS("ACGTCGTC", "CGTCGTCA", "");
  BOOST_CHECK( *kmer1 != *kmer2 );

}

KmerPtr kptr1(kmer1);
KmerPtr kptr2(kmer2);
KmerPtr kptr3(kmer3);

void testKmerPtr()
{
  //verify initial conditions
  BOOST_CHECK_EQUAL( kptr1.get(), kmer1);
  BOOST_CHECK_EQUAL( kptr2.get(), kmer2);
  BOOST_CHECK_EQUAL( kptr3.get(), kmer3);
  
  KmerPtr a = kptr1;
  KmerPtr b = kptr2;
  KmerPtr c = kptr3;
  
  BOOST_CHECK( kptr1 == a);
  BOOST_CHECK( kptr2 == b);
  BOOST_CHECK( kptr3 == c);
  BOOST_CHECK( a != c);
  
  SET_KMERS("ACGTCGTA", "TACGACGT", "AAAACCCCGGGGTTTTACGTCGTAGTACTACG");
  BOOST_CHECK_EQUAL( "ACGTCGTA", a++->toFasta());
  BOOST_CHECK_EQUAL( "ACGTCGTA", (--a)->toFasta());
  
  BOOST_CHECK_EQUAL( "TACGACGT", b++->toFasta());
  BOOST_CHECK_EQUAL( "TACGACGT", (--b)->toFasta());
  
  BOOST_CHECK_EQUAL( "AAAACCCC", c++->toFasta());
  BOOST_CHECK_EQUAL( "GGGGTTTT", c++->toFasta());
  BOOST_CHECK_EQUAL( "ACGTCGTA", c++->toFasta());
  BOOST_CHECK_EQUAL( "GTACTACG", c++->toFasta());

  BOOST_CHECK_EQUAL( "GTACTACG", (--c)->toFasta());
  BOOST_CHECK_EQUAL( "ACGTCGTA", (--c)->toFasta());
  BOOST_CHECK_EQUAL( "GGGGTTTT", (--c)->toFasta());
  BOOST_CHECK_EQUAL( "AAAACCCC", (--c)->toFasta());

  BOOST_CHECK( kptr3 == c);
  
  BOOST_CHECK_EQUAL( "AAAACCCC", c[0].toFasta());
  BOOST_CHECK_EQUAL( "GGGGTTTT", c[1].toFasta());
  BOOST_CHECK_EQUAL( "ACGTCGTA", c[2].toFasta());
  BOOST_CHECK_EQUAL( "GTACTACG", c[3].toFasta());

  BOOST_CHECK( kptr3 == c);
  
  // and original have not changed...
  a++; b++; c++;
  BOOST_CHECK_EQUAL( kptr1.get(), kmer1);
  BOOST_CHECK_EQUAL( kptr2.get(), kmer2);
  BOOST_CHECK_EQUAL( kptr3.get(), kmer3);  
}

BOOST_AUTO_TEST_CASE( KmerSetTest )
{
  testKmerCompare();
  testKmerPtr();
}

//
// $Log: KmerTest.cpp,v $
// Revision 1.1  2009-10-23 07:06:57  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
//
