// $Header: /repository/PI_annex/robsandbox/KoMer/test/KmerTest.cpp,v 1.10 2009-10-26 23:04:35 regan Exp $
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
// do not want to support these....
/* Kmer    *kmer1     = (Kmer*)&twoBit1;
Kmer    *kmer2     = (Kmer*)&twoBit2;
Kmer    *kmer3     = (Kmer*)&twoBit3;
Kmer    &kmer1ref  = *kmer1;
Kmer    &kmer2ref  = *kmer2;
Kmer    &kmer3ref  = *kmer3;
 */
 
KmerPtr kmer1(&twoBit1), kmer2(&twoBit2), kmer3(&twoBit3);


#define SET_KMERS(fasta1, fasta2, fasta3) \
    TwoBitSequence::compressSequence(fasta1,   twoBit1);\
    TwoBitSequence::compressSequence(fasta2,   twoBit2);\
    TwoBitSequence::compressSequence(fasta3,   twoBit3);\
    KmerSizer::set(std::strlen(fasta1));\
    BOOST_CHECK_EQUAL( fasta1, kmer1.toFasta());\
    BOOST_CHECK_EQUAL( fasta2, kmer2.toFasta());

void testKmerCompare()
{
  SET_KMERS("AAAA", "AAAA", "");
  BOOST_CHECK( kmer1.equals(kmer2) );
  
  SET_KMERS("AAAAA", "AAAAA", "");
  BOOST_CHECK( kmer1.equals(kmer2) );
  
  SET_KMERS("AAAAAA", "AAAAAA", "");
  BOOST_CHECK( kmer1.equals(kmer2) );
  
  SET_KMERS("AAAAAAA", "AAAAAAA", "");
  BOOST_CHECK( kmer1.equals(kmer2) );

  SET_KMERS("AAAAAAAA", "AAAAAAAA", "");
  BOOST_CHECK( kmer1.equals(kmer2) );

  SET_KMERS("ACGT", "ACGT", "");
  BOOST_CHECK( kmer1.equals(kmer2) );
  
  SET_KMERS("ACGTC", "ACGTC", "");
  BOOST_CHECK( kmer1.equals(kmer2) );
  
  SET_KMERS("ACGTCG", "ACGTCG", "");
  BOOST_CHECK( kmer1.equals(kmer2) );
  
  SET_KMERS("ACGTCGT", "ACGTCGT", "");
  BOOST_CHECK( kmer1.equals(kmer2) );

  SET_KMERS("ACGTCGTC", "ACGTCGTC", "");
  BOOST_CHECK( kmer1.equals(kmer2) );

  SET_KMERS("ACGT", "CGTA", "");
  BOOST_CHECK( !kmer1.equals(kmer2) );
  
  SET_KMERS("ACGTC", "CGTCA", "");
  BOOST_CHECK(  !kmer1.equals(kmer2) );
  
  SET_KMERS("ACGTCG", "CGTCGA", "");
  BOOST_CHECK(  !kmer1.equals(kmer2) );
  
  SET_KMERS("ACGTCGT", "CGTCGTA", "");
  BOOST_CHECK(  !kmer1.equals(kmer2) );

  SET_KMERS("ACGTCGTC", "CGTCGTCA", "");
  BOOST_CHECK(  !kmer1.equals(kmer2) );

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
  //verify initial conditions
  BOOST_CHECK( kptr1.get() == kmer1);
  BOOST_CHECK( kptr2.get() == kmer2);
  BOOST_CHECK( kptr3.get() == kmer3);
 
  BOOST_CHECK( kptr1 == *&kptr1);
  BOOST_CHECK( kptr2 == *&kptr2);
  BOOST_CHECK( kptr3 == *&kptr3); 
  
  // do not want to support these...
  //BOOST_CHECK( *kptr1 == kmer1ref);
  //BOOST_CHECK( *kptr2 == kmer2ref);
  //BOOST_CHECK( *kptr3 == kmer3ref); 

  KmerPtr a = kptr1;
  KmerPtr b = kptr2;
  KmerPtr c = kptr3;  

  BOOST_CHECK( kptr1 == a);
  BOOST_CHECK( kptr2 == b);
  BOOST_CHECK( kptr3 == c);
  BOOST_CHECK( a != c);

/*   // Do not want to support theses...
  Kmer *a_ptr = kmer1;
  Kmer *b_ptr = (Kmer *)b.get();
  Kmer *c_ptr = (Kmer *)c.get();//(Kmer *)((void *)&(*c));
  BOOST_CHECK( a_ptr == kmer1 );
  BOOST_CHECK( b_ptr == kmer2 );
  BOOST_CHECK( c_ptr == kmer3 );

  // Do not want to support theses...
  Kmer &a_ref = (Kmer &)*kmer1;
  Kmer &b_ref = *b_ptr;
  Kmer &c_ref = *c_ptr; //*c;  
  // test Kmer == Kmer
  BOOST_CHECK( a_ref == *a_ptr );
  BOOST_CHECK( b_ref == *b_ptr );
  BOOST_CHECK( c_ref == *c_ptr );
  // test KmerPtr == KmerPtr
  BOOST_CHECK( &a_ref == a_ptr );
  BOOST_CHECK( &b_ref == b_ptr );
  BOOST_CHECK( &c_ref == c_ptr );
 */  
  std::string A("ACGTCGTAACGTCGTA"), B("TACGACGTTACGACGT"), C("AAAACCCCGGGGTTTTACGTCGTAGTACTACGAAAACCCCGGGGTTTTACGTCGTAGTACTACG");
  SET_KMERS(A.c_str(), B.c_str(), C.c_str());
  KmerSizer::set(size);
  
  
  // check KmerPtr

  BOOST_CHECK( *a     == a[0] );
  BOOST_CHECK( *(a+0) == a[0] ); 
  BOOST_CHECK( *(a+1) == a[1] ); 
  
  BOOST_CHECK( *b     == b[0] );
  BOOST_CHECK( *(b+0) == b[0] ); 
  BOOST_CHECK( *(b+1) == b[1] ); 
  
  BOOST_CHECK( *c     == c[0] );
  BOOST_CHECK( *(c+0) == c[0] ); 
  BOOST_CHECK( *(c+1) == c[1] ); 
    
  KmerPtr copy = a;
  BOOST_CHECK( a == copy );
  BOOST_CHECK( copy == a );


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
  
  

/*   // Do not want o support these
  if(0)  
  {  
  	// check Kmer * (I do not know if this will ever work...)
    BOOST_CHECK_EQUAL( SS(A,0,size), a_ptr++->toFasta());
    BOOST_CHECK_EQUAL( SS(A,0,size), (--a_ptr)->toFasta());
  
    BOOST_CHECK_EQUAL( SS(B,0,size), b_ptr++->toFasta());
    BOOST_CHECK_EQUAL( SS(B,0,size), (--b_ptr)->toFasta());

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
 */  
  // and original have not changed...
  a++; b++; c++;
  BOOST_CHECK( kptr1.get() == kmer1);
  BOOST_CHECK( kptr2.get() == kmer2);
  BOOST_CHECK( kptr3.get() == kmer3);  
}

#define SS2(str, offset, len) str.substr(offset, len)

void testKmerArray(SequenceLengthType size)
{
  kmerBytesJump = (size+3) /4;
	
  std::string A("ACGTCGTAACGTCGTA"), B("TACGACGTTACGACGT"), C("AAAACCCCGGGGTTTTACGTCGTAGTACTACGAAAACCCCGGGGTTTTACGTCGTAGTACTACG");
  KmerArray<WeakKmerTag>::releasePools();
  KmerArray<SolidKmerTag>::releasePools();
  SET_KMERS(A.c_str(), B.c_str(), C.c_str());
  KmerSizer::set(size);
  
  struct myTag {};
  KmerArray<WeakKmerTag> kmersA(twoBit1, A.length());
  KmerArray<SolidKmerTag> kmersB(twoBit2, B.length());
  KmerArray<> kmersC(twoBit3, C.length());
  
  for (int i=0; i< A.length() - size +1; i++) {
    BOOST_CHECK_EQUAL( kmersA[i].toFasta(), SS2(A,i,size));
  }
  for (int i=0; i< B.length() - size +1; i++) {
    BOOST_CHECK_EQUAL( kmersB[i].toFasta(), SS2(B,i,size));
  }
  for (int i=0; i< C.length() - size +1; i++) {
    BOOST_CHECK_EQUAL( kmersC[i].toFasta(), SS2(C,i,size));
  }
  
}

void testKmerNewDelete()
{
#if 0
  // This should never compile!
  Kmer *one = new Kmer;
  delete one;
  
  Kmer *many = new Kmer[100];
  delete [] many;
  
  Kmer bad;
#endif
}
/* 
void testKmerInstance()
{
  KmerInstance one;
  KmerInstance many[100];
  KmerSizer::set(8);
  
  std::string A("ACGTCGTAACGTCGTA"), B("TACGACGTTACGACGT"), C("AAAACCCCGGGGTTTTACGTCGTAGTACTACGAAAACCCCGGGGTTTTACGTCGTAGTACTACG");
  SET_KMERS(A.c_str(), B.c_str(), C.c_str()); 
  KmerSizer::set(8);
  
  one = *kmer1;
  many[0] = *kptr1;//kmer1ref;
  many[1] = *kptr2;
  many[2] = *kmer3;

  BOOST_CHECK_EQUAL( "ACGTCGTA", one.toFasta() );
  BOOST_CHECK_EQUAL( "ACGTCGTA", many[0].toFasta() );
  BOOST_CHECK_EQUAL( "TACGACGT", many[1].toFasta() );
  BOOST_CHECK_EQUAL( "AAAACCCC", many[2].toFasta() );
}
 */

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

  testKmerArray(1);
  testKmerArray(2);
  testKmerArray(3);
  testKmerArray(4);
  testKmerArray(5);
  testKmerArray(6);
  testKmerArray(7);
  testKmerArray(8);
  testKmerArray(9);
  
  testKmerNewDelete();
  
  //testKmerInstance();
}

//
// $Log: KmerTest.cpp,v $
// Revision 1.10  2009-10-26 23:04:35  regan
// checkpoint make Kmer private inner class
//
// Revision 1.9  2009-10-26 17:42:26  regan
// templated KmerArray
//
// Revision 1.8  2009-10-24 00:37:50  regan
// fixed tests
//
// Revision 1.7  2009-10-24 00:32:47  regan
// added bugs
//
// Revision 1.6  2009-10-24 00:03:51  regan
// checkpoint
//
// Revision 1.5  2009-10-23 23:22:44  regan
// checkpoint
//
// Revision 1.4  2009-10-23 21:54:48  regan
// checkpoint
//
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
