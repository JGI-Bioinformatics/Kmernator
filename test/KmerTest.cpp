// $Header: /repository/PI_annex/robsandbox/KoMer/test/KmerTest.cpp,v 1.3 2009-10-23 20:32:52 cfurman Exp $
//
 

#include "TwoBitSequence.h"
#include "Kmer.h"

#define BOOST_TEST_MODULE KmerSetTest
#include <boost/test/unit_test.hpp>
#include <cstring>
#include <iostream>
#include <cstdlib>

using namespace std;

TwoBitEncoding twoBit1[1024],twoBit2[1024],twoBit3[1024];
Kmer    *kmer1     = (Kmer*)&twoBit1;
Kmer    *kmer2     = (Kmer*)&twoBit2;
Kmer    *kmer3     = (Kmer*)&twoBit3;
Kmer    &kmer1ref  = *kmer1;
Kmer    &kmer2ref  = *kmer2;
Kmer    &kmer3ref  = *kmer3;

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
SequenceLengthType kmerBytesJump = 0;
#define SS(str, hops, len) str.substr(hops*kmerBytesJump*4, len)

void testKmerPtr(SequenceLengthType size)
{
  unsigned char bitshift = size % 4;
  kmerBytesJump = (size+3) /4;
  cout << "Executing testKmerPtr(" <<  size << "); kmerBytesJump == " << kmerBytesJump << endl;
  //verify initial conditions
  BOOST_CHECK_EQUAL( kptr1.get(), kmer1);
  BOOST_CHECK_EQUAL( kptr2.get(), kmer2);
  BOOST_CHECK_EQUAL( kptr3.get(), kmer3);
 
  BOOST_CHECK( kptr1 == *&kptr1);
  BOOST_CHECK( kptr2 == *&kptr2);
  BOOST_CHECK( kptr3 == *&kptr3); 
  
  BOOST_CHECK( *kptr1 == kmer1ref);
  BOOST_CHECK( *kptr2 == kmer2ref);
  BOOST_CHECK( *kptr3 == kmer3ref); 

  KmerPtr a = kptr1;
  KmerPtr b = kptr2;
  KmerPtr c = kptr3;  

  BOOST_CHECK( kptr1 == a);
  BOOST_CHECK( kptr2 == b);
  BOOST_CHECK( kptr3 == c);
  BOOST_CHECK( a != c);

  Kmer *a_ptr = kmer1;
  Kmer *b_ptr = (Kmer *)b.get();
  Kmer *c_ptr = (Kmer *)((void *)&(*c));
  BOOST_CHECK( a_ptr == kmer1 );
  BOOST_CHECK( b_ptr == kmer2 );
  BOOST_CHECK( c_ptr == kmer3 );

  Kmer &a_ref = (Kmer &)*kmer1;
  Kmer &b_ref = *b_ptr;
  Kmer &c_ref = *c;  
  // test Kmer == Kmer
  BOOST_CHECK( a_ref == *a_ptr );
  BOOST_CHECK( b_ref == *b_ptr );
  BOOST_CHECK( c_ref == *c_ptr );
  // test KmerPtr == KmerPtr
  BOOST_CHECK( &a_ref == a_ptr );
  BOOST_CHECK( &b_ref == b_ptr );
  BOOST_CHECK( &c_ref == c_ptr );
  
  std::string A("ACGTCGTAACGTCGTA"), B("TACGACGTTACGACGT"), C("AAAACCCCGGGGTTTTACGTCGTAGTACTACGAAAACCCCGGGGTTTTACGTCGTAGTACTACG");
  SET_KMERS(A.c_str(), B.c_str(), C.c_str());
  KmerSizer::set(size);
  
  // check KmerPtr
 //   BOOST_CHECK_EQUAL( SS(A,0,size), (*a++).toFasta());
  BOOST_CHECK_EQUAL( SS(A,0,size), a++->toFasta());
  BOOST_CHECK_EQUAL( SS(A,0,size), (--a)->toFasta());
  
  BOOST_CHECK_EQUAL( SS(B,0,size), b++->toFasta());
  BOOST_CHECK_EQUAL( SS(B,0,size), (--b)->toFasta());
  
  BOOST_CHECK_EQUAL( SS(C,0,size), c++->toFasta());
  BOOST_CHECK_EQUAL( SS(C,1,size), c++->toFasta());
  BOOST_CHECK_EQUAL( SS(C,2,size), c++->toFasta());
  BOOST_CHECK_EQUAL( SS(C,3,size), c++->toFasta());

  BOOST_CHECK_EQUAL( SS(C,3,size), (--c)->toFasta());
  BOOST_CHECK_EQUAL( SS(C,2,size), (--c)->toFasta());
  BOOST_CHECK_EQUAL( SS(C,1,size), (--c)->toFasta());
  BOOST_CHECK_EQUAL( SS(C,0,size), (--c)->toFasta());

  BOOST_CHECK( kptr3 == c);
  
  BOOST_CHECK_EQUAL( SS(C,0,size), c[0].toFasta());
  BOOST_CHECK_EQUAL( SS(C,1,size), c[1].toFasta());
  BOOST_CHECK_EQUAL( SS(C,2,size), c[2].toFasta());
  BOOST_CHECK_EQUAL( SS(C,3,size), c[3].toFasta());

  BOOST_CHECK( kptr3 == c);
  
  

  // check Kmer * (I do not know if this will ever work...)
  BOOST_CHECK_EQUAL( SS(A,0,size), a_ptr++->toFasta());
  BOOST_CHECK_EQUAL( SS(A,0,size), (--a_ptr)->toFasta());
  
  BOOST_CHECK_EQUAL( SS(B,0,size), b_ptr++->toFasta());
  BOOST_CHECK_EQUAL( SS(B,0,size), (--b_ptr)->toFasta());
  if(0)  
  {
    BOOST_CHECK_EQUAL( SS(C,0,size), c_ptr++->toFasta());
    BOOST_CHECK_EQUAL( SS(C,1,size), c_ptr++->toFasta());
    BOOST_CHECK_EQUAL( SS(C,2,size), c_ptr++->toFasta());
    BOOST_CHECK_EQUAL( SS(C,3,size), c_ptr++->toFasta());

    BOOST_CHECK_EQUAL( SS(C,3,size), (--c_ptr)->toFasta());
    BOOST_CHECK_EQUAL( SS(C,2,size), (--c_ptr)->toFasta());
    BOOST_CHECK_EQUAL( SS(C,1,size), (--c_ptr)->toFasta());
    BOOST_CHECK_EQUAL( SS(C,0,size), (--c_ptr)->toFasta());

    BOOST_CHECK( kptr3 == c_ptr);  
  }
  
  // and original have not changed...
  a++; b++; c++;
  BOOST_CHECK_EQUAL( kptr1.get(), kmer1);
  BOOST_CHECK_EQUAL( kptr2.get(), kmer2);
  BOOST_CHECK_EQUAL( kptr3.get(), kmer3);  
}

BOOST_AUTO_TEST_CASE( KmerSetTest )
{
  testKmerCompare();
  testKmerPtr(1);
  testKmerPtr(2);
  testKmerPtr(3);
  testKmerPtr(4);
  testKmerPtr(5);
  testKmerPtr(6);
  testKmerPtr(7);
  testKmerPtr(8);
  testKmerPtr(9);
}

//
// $Log: KmerTest.cpp,v $
// Revision 1.3  2009-10-23 20:32:52  cfurman
// more kmer changes
//
// Revision 1.2  2009-10-23 17:22:45  regan
// added more tests
//
// Revision 1.1  2009-10-23 07:06:57  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
//
