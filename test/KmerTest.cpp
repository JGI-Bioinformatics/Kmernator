//
// Kmernator/test/KmerTest.cpp
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

#include "config.h"
#include "TwoBitSequence.h"
#include "Kmer.h"

#define BOOST_TEST_MODULE KmerSetTest
#include <boost/test/unit_test.hpp>
#include <cstring>
#include <iostream>
#include <cstdlib>

// Note for versbosity: export BOOST_TEST_LOG_LEVEL=message

using namespace std;

TwoBitEncoding twoBit1[1024], twoBit2[1024], twoBit3[1024];
// do not want to support these....
/* Kmer    *kmer1     = (Kmer*)&twoBit1;
 Kmer    *kmer2     = (Kmer*)&twoBit2;
 Kmer    *kmer3     = (Kmer*)&twoBit3;
 Kmer    &kmer1ref  = *kmer1;
 Kmer    &kmer2ref  = *kmer2;
 Kmer    &kmer3ref  = *kmer3;
 */

Kmer &kmer1 = *((Kmer *) (twoBit1));
Kmer &kmer2 = *((Kmer *) (twoBit2));
Kmer &kmer3 = *((Kmer *) (twoBit3));

#define SET_KMERS(fasta1, fasta2, fasta3) \
		TwoBitSequence::compressSequence(fasta1,   twoBit1);\
		TwoBitSequence::compressSequence(fasta2,   twoBit2);\
		TwoBitSequence::compressSequence(fasta3,   twoBit3);\
		KmerSizer::set(std::strlen(fasta1));\
		BOOST_CHECK_EQUAL( fasta1, kmer1.toFasta());\
		BOOST_CHECK_EQUAL( fasta2, kmer2.toFasta());

void testKmerCompare() {

	SET_KMERS("AAAA", "AAAA", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("AAAAA", "AAAAA", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("AAAAAA", "AAAAAA", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("AAAAAAA", "AAAAAAA", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("AAAAAAAA", "AAAAAAAA", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("ACGT", "ACGT", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("ACGTC", "ACGTC", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("ACGTCG", "ACGTCG", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("ACGTCGT", "ACGTCGT", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("ACGTCGTC", "ACGTCGTC", "");
	BOOST_CHECK(kmer1 == (kmer2));

	SET_KMERS("ACGT", "CGTA", "");
	BOOST_CHECK(kmer1 != (kmer2));

	SET_KMERS("ACGTC", "CGTCA", "");
	BOOST_CHECK(kmer1 != (kmer2));

	SET_KMERS("ACGTCG", "CGTCGA", "");
	BOOST_CHECK(kmer1 != (kmer2));

	SET_KMERS("ACGTCGT", "CGTCGTA", "");
	BOOST_CHECK(kmer1 != (kmer2));

	SET_KMERS("ACGTCGTC", "CGTCGTCA", "");
	BOOST_CHECK(kmer1 != (kmer2));

}
#if 0
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
	BOOST_CHECK( kptr1.getTwoBitSequence() == kmer1);
	BOOST_CHECK( kptr2.getTwoBitSequence() == kmer2);
	BOOST_CHECK( kptr3.getTwoBitSequence() == kmer3);

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

	// Do not want to support theses...
	/*
	 Kmer *a_ptr = kmer1;
	 Kmer *b_ptr = (Kmer *)b.getTwoBitSequence();
	 Kmer *c_ptr = (Kmer *)c.getTwoBitSequence();//(Kmer *)((void *)&(*c));
	 BOOST_CHECK( a_ptr == kmer1 );
	 BOOST_CHECK( b_ptr == kmer2 );
	 BOOST_CHECK( c_ptr == kmer3 );
	 */
	// Do not want to support theses...
	/*
	 Kmer &a_ref = (Kmer &)*kmer1;
	 Kmer &b_ref = *b_ptr;
	 Kmer &c_ref = *c_ptr; // *c;
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

	BOOST_CHECK( *a == a[0] );
	BOOST_CHECK( *(a+0) == a[0] );
	BOOST_CHECK( *(a+1) == a[1] );

	BOOST_CHECK( *b == b[0] );
	BOOST_CHECK( *(b+0) == b[0] );
	BOOST_CHECK( *(b+1) == b[1] );

	BOOST_CHECK( *c == c[0] );
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
	BOOST_CHECK( kptr1.getTwoBitSequence() == kmer1);
	BOOST_CHECK( kptr2.getTwoBitSequence() == kmer2);
	BOOST_CHECK( kptr3.getTwoBitSequence() == kmer3);
}
#endif

#define SS2(str, offset, len) str.substr(offset, len)

void testKmerArray(SequenceLengthType size) {
	//   kmerBytesJump = (size+3) /4;
	//
	//BOOST_MESSAGE( size );
	std::string A("ACGTCGTAACGTCGTA"), B("TACGACGTTACGACGT"), C(
			"AAAACCCCGGGGTTTTACGTCGTAGTACTACGAAAACCCCGGGGTTTTACGTCGTAGTACTACG");

	SET_KMERS(A.c_str(), B.c_str(), C.c_str());
	KmerSizer::set(size);

	struct myTag {
	};
	KmerArrayPair<char> kmersA(twoBit1, A.length());
	KmerArrayPair<char> kmersB(twoBit2, B.length());
	KmerArrayPair<char> kmersC(twoBit3, C.length());

	if (size == 1) {
		for (Kmer::IndexType i = 0; i < kmersC.size(); i++) {
			char *ch = (char*) kmersC[i].getTwoBitSequence();
			BOOST_CHECK_EQUAL(*ch & 0x3f, 0x00);
		}
	}
	if (size == 2) {
		for (Kmer::IndexType i = 0; i < kmersC.size(); i++) {
			char *ch = (char*) kmersC[i].getTwoBitSequence();
			BOOST_CHECK_EQUAL(*ch & 0x0f, 0x00);
		}
	}
	if (size == 3) {
		for (Kmer::IndexType i = 0; i < kmersC.size(); i++) {
			char *ch = (char*) kmersC[i].getTwoBitSequence();
			BOOST_CHECK_EQUAL(*ch & 0x03, 0x00);
		}
	}
	KmerArrayPair<char> *kmersD = new KmerArrayPair<char> (twoBit3, C.length());

	for (size_t i = 0; i < A.length() - size + 1; i++) {
		BOOST_CHECK_EQUAL(kmersA[i].toFasta(), SS2(A,i,size));
	}
	for (size_t i = 0; i < B.length() - size + 1; i++) {
		BOOST_CHECK_EQUAL(kmersB[i].toFasta(), SS2(B,i,size));
	}
	for (size_t i = 0; i < C.length() - size + 1; i++) {
		BOOST_CHECK_EQUAL(kmersC[i].toFasta(), SS2(C,i,size));
	}
	for (size_t i = 0; i < C.length() - size + 1; i++) {
		BOOST_CHECK_EQUAL((*kmersD)[i].toFasta(), SS2(C,i,size));
	}
	delete kmersD;

	KmerArrayPair<float> kmersFloat(twoBit3, C.length());
	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		float &valRef = kmersFloat.valueAt(i);
		valRef = i * 2.0;
	}
	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		float &valRef = kmersFloat.valueAt(i);
		float val = i * 2.0;
		BOOST_CHECK_EQUAL(val, valRef);
		BOOST_CHECK_EQUAL(kmersC[i].toFasta(), kmersFloat[i].toFasta());
	}

	// test KmerArray assignment

	void *mem1, *mem2;
	mem1 = kmersFloat[0].getTwoBitSequence();
	KmerArrayPair<float> copy = kmersFloat;
	BOOST_CHECK_EQUAL(copy.size(), kmersFloat.size());
	BOOST_CHECK_EQUAL(mem1, (void*) kmersFloat[0].getTwoBitSequence());
	mem2 = copy[0].getTwoBitSequence();
	BOOST_CHECK(mem1 != mem2);

	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		float &valRef = copy.valueAt(i);
		float &valRef2 = kmersFloat.valueAt(i);
		BOOST_CHECK_EQUAL(valRef2, valRef);
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i].toFasta());
	}

	// test resize and []
	unsigned long oldSize = copy.size();
	copy.resize(oldSize + 1);
	mem1 = copy[0].getTwoBitSequence();
	BOOST_CHECK(mem1 != mem2);

	BOOST_CHECK_EQUAL(copy.size(), kmersFloat.size() + 1);
	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		float &valRef = copy.valueAt(i);
		float &valRef2 = kmersFloat.valueAt(i);
		BOOST_CHECK_EQUAL(valRef2, valRef);
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i].toFasta());
	}

	copy[oldSize] = kmersFloat[0];
	BOOST_CHECK_EQUAL(kmersFloat[0].toFasta(), copy[oldSize].toFasta());
	copy.valueAt(oldSize) = oldSize * 2.0;
	BOOST_CHECK_EQUAL(oldSize * 2.0, copy.valueAt(oldSize));
	BOOST_CHECK(copy[oldSize].getTwoBitSequence() != kmersFloat[0].getTwoBitSequence());

	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		float &valRef = copy.valueAt(i);
		float &valRef2 = kmersFloat.valueAt(i);
		BOOST_CHECK_EQUAL(valRef2, valRef);
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i].toFasta());
	}
	BOOST_CHECK_EQUAL(kmersFloat[0].toFasta(), copy[oldSize].toFasta());
	BOOST_CHECK_EQUAL(oldSize * 2.0, copy.valueAt(oldSize));

	// now reduce
	copy.resize(oldSize);
	mem2 = copy[0].getTwoBitSequence();

	BOOST_CHECK_EQUAL(copy.size(), kmersFloat.size());
	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		float &valRef = copy.valueAt(i);
		float &valRef2 = kmersFloat.valueAt(i);
		BOOST_CHECK_EQUAL(valRef2, valRef);
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i].toFasta());
	}

	// test find
	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		Kmer::IndexType idx = copy.findIndex(kmersFloat[i]);
		BOOST_REQUIRE(idx != KmerArrayPair<float>::MAX_INDEX);
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[idx].toFasta());
		// duplicates could be present and first one will be found
		BOOST_CHECK(idx <= i);
		BOOST_CHECK(copy.valueAt(idx) <= kmersFloat.valueAt(i));
	}

	// test insert
	BOOST_CHECK_EQUAL(kmersFloat.size(), copy.size());
	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		//BOOST_MESSAGE( i );

		BOOST_CHECK_EQUAL(kmersFloat.size() + i, copy.size());
		copy.insertAt(i * 2, kmersFloat[i]);
		copy.valueAt(i * 2) = kmersFloat.valueAt(i) - 1.0;
		BOOST_CHECK_EQUAL(kmersFloat.size() + i + 1, copy.size());

		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i * 2].toFasta());
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i * 2 + 1].toFasta());
		BOOST_CHECK_EQUAL(kmersFloat.valueAt(i), copy.valueAt(i * 2 + 1));
		BOOST_CHECK_EQUAL(kmersFloat.valueAt(i) - 1.0, copy.valueAt(i * 2));
	}
	BOOST_CHECK_EQUAL(kmersFloat.size() * 2, copy.size());

	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		float &valRef = copy.valueAt(i * 2 + 1);
		float &valRef2 = kmersFloat.valueAt(i);
		BOOST_CHECK_EQUAL(valRef2, valRef);
		BOOST_CHECK_EQUAL(valRef2 - 1.0, copy.valueAt(i * 2));
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i * 2].toFasta());
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i * 2 + 1].toFasta());
	}

	// test remove
	BOOST_CHECK_EQUAL(kmersFloat.size() * 2, copy.size());
	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		BOOST_CHECK_EQUAL(kmersFloat.size() * 2 - i, copy.size());
		copy.remove(i);
		BOOST_CHECK_EQUAL(kmersFloat.size() * 2 - i - 1, copy.size());

		for (Kmer::IndexType k = 0; k < i + 1; k++) {
			BOOST_CHECK_EQUAL(kmersFloat[k].toFasta(), copy[k].toFasta());
			BOOST_CHECK_EQUAL(kmersFloat.valueAt(k), copy.valueAt(k));
		}
	}

	BOOST_CHECK_EQUAL(copy.size(), kmersFloat.size());
	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		float &valRef = copy.valueAt(i);
		float &valRef2 = kmersFloat.valueAt(i);
		BOOST_CHECK_EQUAL(valRef2, valRef);
		BOOST_CHECK_EQUAL(kmersFloat[i].toFasta(), copy[i].toFasta());
	}

	BOOST_CHECK_EQUAL(copy.size(), kmersFloat.size());
	Kmer::IndexType count = 0;
	for (KmerArrayPair<float>::Iterator it = kmersFloat.begin(); it
	!= kmersFloat.end(); it++) {
		BOOST_CHECK_EQUAL(kmersFloat.valueAt(count), it->value());
		BOOST_CHECK_EQUAL(kmersFloat[count].toFasta(), it->key().toFasta());
		count++;
	}
	BOOST_CHECK_EQUAL(copy.size(), kmersFloat.size());
	BOOST_CHECK_EQUAL(copy.size(), count);

	/* print out the permutations table
	 if (size == 4)  {
	 for(int i=0; i<=255; i++) {
	 TwoBitEncoding tb = i;
	 KmerPtr k2( &tb );
	 std::stringstream ss;
	 ss << k2->toFasta() << ": ";
	 for (int j=0; j<12; j++) {
	 KmerPtr k( &(TwoBitSequence::permutations[i*12+j]) );
	 ss << k->toFasta() << ", ";
	 }
	 BOOST_MESSAGE( ss.str() );
	 }
	 }  else
	 return;
	 */

	for (Kmer::IndexType i = 0; i < kmersFloat.size(); i++) {
		KmerArrayPair<float> permutations = KmerArrayPair<float>::permuteBases(
				kmersFloat[i]);
		BOOST_CHECK_EQUAL(KmerSizer::getSequenceLength() * 3,
				permutations.size());
		//BOOST_MESSAGE( "start" );
		//BOOST_MESSAGE( kmersFloat[i].toFasta() );
		for (Kmer::IndexType j = 0; j < permutations.size(); j++) {
			BOOST_CHECK_NE(kmersFloat[i].toFasta(), permutations[j].toFasta());
			//BOOST_MESSAGE( permutations[j].toFasta() );
			for (Kmer::IndexType k = 0; k < j; k++) {
				BOOST_CHECK_NE(permutations[j].toFasta(),
						permutations[k].toFasta());
			}
		}
	}

	//permuteBases(const Kmer &kmer, KmerArray &kmers, short editDistance, bool leastComplement = false)

	for (Kmer::IndexType i = 0; i < kmersFloat.size() && i < 5; i++) {
		KmerArrayPair<float> permutations;
		KmerArrayPair<float>::permuteBases(	kmersFloat[i], permutations, 2);

		//BOOST_MESSAGE( "start" );
		//BOOST_MESSAGE( kmersFloat[i].toFasta() );
		for (Kmer::IndexType j = 0; j < permutations.size(); j++) {
			BOOST_CHECK_NE(kmersFloat[i].toFasta(), permutations[j].toFasta());
			//BOOST_MESSAGE( permutations[j].toFasta() );
			for (Kmer::IndexType k = 0; k < j; k++) {
				BOOST_CHECK_NE(permutations[j].toFasta(),
						permutations[k].toFasta());
			}
		}
	}
	for (Kmer::IndexType i = 0; i < kmersFloat.size() && i < 5; i++) {
		KmerArrayPair<float> permutations;
		KmerArrayPair<float>::permuteBases(	kmersFloat[i], permutations, 3);

		//BOOST_MESSAGE( "start" );
		//BOOST_MESSAGE( kmersFloat[i].toFasta() );
		for (Kmer::IndexType j = 0; j < permutations.size(); j++) {
			BOOST_CHECK_NE(kmersFloat[i].toFasta(), permutations[j].toFasta());
			//BOOST_MESSAGE( permutations[j].toFasta() );
			for (Kmer::IndexType k = 0; k < j; k++) {
				BOOST_CHECK_NE(permutations[j].toFasta(),
						permutations[k].toFasta());
			}
		}
	}



}

template<typename Map>
class Tester {
public:

	Map &kmerF;
	Tester(Map map) :
		kmerF(map) {

	}
	template<typename U>
	void operator()(U &e) {
		BOOST_CHECK(kmerF.exists(e.key()));
		BOOST_CHECK_EQUAL(kmerF[e.key()], e.value());
	}
};

template<typename MapV, typename MapPairV>
void initTestKmerMap(MapV &kmerF, MapPairV &kmerP, const KmerArrayPair<char> &kmersC) {
	typedef typename MapPairV::mapped_type Pair;

	for (Kmer::IndexType i = 0; i < kmersC.size(); i++) {
		kmerF[kmersC[i]] = i * 2.0;
		kmerP[kmersC[i]] = Pair(i, i * 3.0);
		//BOOST_MESSAGE(kmersC[i].toFasta());
		BOOST_CHECK_EQUAL(i * 2.0, kmerF[kmersC[i]]);
		BOOST_CHECK_EQUAL(i, kmerP[kmersC[i]].first);
		BOOST_CHECK_EQUAL(i * 3.0, kmerP[kmersC[i]].second);

		//BOOST_MESSAGE ( kmerF.toString() );
	}

}
template<typename MapV, typename MapPairV>
void testStore(SequenceLengthType size) {

	std::string
	A("ACGTCGTAACGTCGTA"),
	B("TACGACGTTACGACGT"),
	C("AAAACCCCGGGGTTTTTACGTCGTAGTACTACGAAAACCCCGGGGTTTTACGTCGTAGTACTACG");
	SET_KMERS(A.c_str(), B.c_str(), C.c_str());
	KmerSizer::set(size);

	// test KmerMap construction, destruction
	// test insert, find, delete

	KmerArrayPair<char> kmersC(twoBit3, C.length());
	Kmer::IndexType s = kmersC.size();
	BOOST_CHECK_EQUAL(s, kmersC.size());


	typedef typename MapV::Iterator MapVIterator;
	MapV kmerF(4);
	MapPairV kmerP(8);


	initTestKmerMap(kmerF, kmerP, kmersC);

	Kmernator::MmapFile mmapF = kmerF.store();
	BOOST_CHECK(mmapF.is_open());
	BOOST_CHECK(mmapF.size() > 0);
	Kmernator::MmapFile mmapP = kmerP.store();
	BOOST_CHECK(mmapP.is_open());
	BOOST_CHECK(mmapP.size() > 0);

	{

		MapV kmf = MapV::restore(mmapF.data());
		MapPairV  kmp = MapPairV::restore(mmapP.data());

		MapV kcf(mmapF.data());
		MapPairV  kcp(mmapP.data());

		for (Kmer::IndexType i = 0; i < kmersC.size(); i++) {
			Kmer &kmer = kmersC[i];
			BOOST_CHECK_EQUAL( kmerF[kmer], kmf[kmer] );
			BOOST_CHECK_EQUAL( kmerP[kmer].first,  kmp[kmer].first );
			BOOST_CHECK_EQUAL( kmerP[kmer].second, kmp[kmer].second );
			BOOST_CHECK_EQUAL( kmerF[kmer], kcf[kmer] );
			BOOST_CHECK_EQUAL( kmerP[kmer].first,  kcp[kmer].first );
			BOOST_CHECK_EQUAL( kmerP[kmer].second, kcp[kmer].second );
		}
	}

}
template<typename MapV, typename MapPairV>
void testKmerMap(SequenceLengthType size) {
	std::string
	A("ACGTCGTAACGTCGTA"),
	B("TACGACGTTACGACGT"),
	C("AAAACCCCGGGGTTTTTACGTCGTAGTACTACGAAAACCCCGGGGTTTTACGTCGTAGTACTACG");
	SET_KMERS(A.c_str(), B.c_str(), C.c_str());
	KmerSizer::set(size);

	// test KmerMap construction, destruction
	// test insert, find, delete

	KmerArrayPair<char> kmersC(twoBit3, C.length());
	Kmer::IndexType s = kmersC.size();
	BOOST_CHECK_EQUAL(s, kmersC.size());

	typedef typename MapPairV::mapped_type Pair;
	typedef typename MapV::Iterator MapVIterator;
	MapV kmerF(4);
	MapPairV kmerP(8);

	initTestKmerMap(kmerF, kmerP, kmersC);

	kmerF.clear();
	kmerP.clear();
	//BOOST_MESSAGE ("Testing exists()");
	for (Kmer::IndexType i = 0; i < kmersC.size(); i++) {
		BOOST_CHECK(!kmerF.exists(kmersC[i]));
		BOOST_CHECK(!kmerF.exists(kmersC[i]));
		BOOST_CHECK(!kmerP.exists(kmersC[i]));
		BOOST_CHECK(!kmerP.exists(kmersC[i]));
	}

	for (Kmer::IndexType i = 0; i < kmersC.size(); i++) {
		if (!kmerF.exists(kmersC[i]))
			kmerF[kmersC[i]] = i * 2.0;
		if (!kmerP.exists(kmersC[i]))
			kmerP[kmersC[i]] = Pair(i, i * 3.0);

		BOOST_CHECK(i * 2.0 >= kmerF[kmersC[i]]);
		BOOST_CHECK(i >= kmerP[kmersC[i]].first);
		BOOST_CHECK(i * 3.0 >= kmerP[kmersC[i]].second);
	}

	std::for_each(kmerF.begin(), kmerF.end(), Tester< MapV >(kmerF));

	Kmer::IndexType count = 0;
	for (MapVIterator it = kmerF.begin(); it != kmerF.end(); it++) {
		BOOST_CHECK(kmerF.exists(it->key()));
		BOOST_CHECK_EQUAL(kmerF[it->key()], it->value());
		count++;
	}

	BOOST_CHECK_EQUAL(kmerF.size(), count);

	Kmer::IndexType countThread;

	countThread = 0;
#pragma omp parallel num_threads(1) reduction(+: countThread)
	for(MapVIterator it = kmerF.beginThreaded(); it != kmerF.end(); it++) {
		BOOST_CHECK(kmerF.exists(it->key()));
		BOOST_CHECK_EQUAL(kmerF[it->key()], it->value());
		countThread++;
	}
	BOOST_CHECK_EQUAL(kmerF.size(), countThread);

	countThread = 0;
#pragma omp parallel num_threads(2) reduction(+: countThread)
	for(MapVIterator it = kmerF.beginThreaded(); it != kmerF.end(); it++) {
		BOOST_CHECK(kmerF.exists(it->key()));
		BOOST_CHECK_EQUAL(kmerF[it->key()], it->value());
		countThread++;
	}
	BOOST_CHECK_EQUAL(kmerF.size(), countThread);

	countThread = 0;
#pragma omp parallel num_threads(3) reduction(+: countThread)
	for(MapVIterator it = kmerF.beginThreaded(); it != kmerF.end(); it++) {
		BOOST_CHECK(kmerF.exists(it->key()));
		BOOST_CHECK_EQUAL(kmerF[it->key()], it->value());
		countThread++;
	}
	BOOST_CHECK_EQUAL(kmerF.size(), countThread);

	countThread = 0;
#pragma omp parallel num_threads(4) reduction(+: countThread)
	for(MapVIterator it = kmerF.beginThreaded(); it != kmerF.end(); it++) {
		BOOST_CHECK(kmerF.exists(it->key()));
		BOOST_CHECK_EQUAL(kmerF[it->key()], it->value());
		countThread++;
	}
	BOOST_CHECK_EQUAL(kmerF.size(), countThread);

	countThread = 0;
#pragma omp parallel num_threads(7) reduction(+: countThread)
	for(MapVIterator it = kmerF.beginThreaded(); it != kmerF.end(); it++) {
		BOOST_CHECK(kmerF.exists(it->key()));
		BOOST_CHECK_EQUAL(kmerF[it->key()], it->value());
		countThread++;
	}
	BOOST_CHECK_EQUAL(kmerF.size(), countThread);

	countThread = 0;
#pragma omp parallel num_threads(31) reduction(+: countThread)
	for(MapVIterator it = kmerF.beginThreaded(); it != kmerF.end(); it++) {
		BOOST_CHECK(kmerF.exists(it->key()));
		BOOST_CHECK_EQUAL(kmerF[it->key()], it->value());
		countThread++;
	}
	BOOST_CHECK_EQUAL(kmerF.size(), countThread);
}

BOOST_AUTO_TEST_CASE( KmerSetTest )
{
	testKmerCompare();
	/*
	 testKmerPtr(1);
	 testKmerPtr(2);
	 testKmerPtr(3);
	 testKmerPtr(4);
	 testKmerPtr(5);
	 testKmerPtr(6);
	 testKmerPtr(7);
	 testKmerPtr(8);
	 testKmerPtr(9);
	 */
	testKmerArray(1);
	testKmerArray(2);
	testKmerArray(3);
	testKmerArray(4);
	testKmerArray(5);
	testKmerArray(6);
	testKmerArray(7);
	testKmerArray(8);
	testKmerArray(9);
	testKmerArray(10);
	testKmerArray(11);
	testKmerArray(12);

	{
	typedef KmerMapByKmerArrayPair<float> M1;
	typedef KmerMapByKmerArrayPair< std::pair<unsigned int, float> > MP;
	testKmerMap<M1,MP>(1);
	testKmerMap<M1,MP>(2);
	testKmerMap<M1,MP>(3);
	testKmerMap<M1,MP>(5);
	testKmerMap<M1,MP>(6);
	testKmerMap<M1,MP>(7);
	testKmerMap<M1,MP>(8);
	testKmerMap<M1,MP>(9);
	testKmerMap<M1,MP>(10);
	testKmerMap<M1,MP>(11);
	testKmerMap<M1,MP>(12);

	testStore<M1,MP>(1);
	testStore<M1,MP>(2);
	testStore<M1,MP>(3);
	testStore<M1,MP>(5);
	testStore<M1,MP>(6);
	testStore<M1,MP>(7);
	testStore<M1,MP>(8);
	testStore<M1,MP>(9);
	testStore<M1,MP>(10);
	testStore<M1,MP>(11);
	testStore<M1,MP>(12);
	}

	{
	typedef KmerMapBoost<float> M1;
	typedef KmerMapBoost< std::pair<unsigned int, float> > MP;
	testKmerMap<M1,MP>(1);
	testKmerMap<M1,MP>(2);
	testKmerMap<M1,MP>(3);
	testKmerMap<M1,MP>(5);
	testKmerMap<M1,MP>(6);
	testKmerMap<M1,MP>(7);
	testKmerMap<M1,MP>(8);
	testKmerMap<M1,MP>(9);
	testKmerMap<M1,MP>(10);
	testKmerMap<M1,MP>(11);
	testKmerMap<M1,MP>(12);
	}

	{
	typedef KmerMapGoogleSparse<float> M1;
	typedef KmerMapGoogleSparse< std::pair<unsigned int, float> > MP;
	testKmerMap<M1,MP>(1);
	testKmerMap<M1,MP>(2);
	testKmerMap<M1,MP>(3);
	testKmerMap<M1,MP>(5);
	testKmerMap<M1,MP>(6);
	testKmerMap<M1,MP>(7);
	testKmerMap<M1,MP>(8);
	testKmerMap<M1,MP>(9);
	testKmerMap<M1,MP>(10);
	testKmerMap<M1,MP>(11);
	testKmerMap<M1,MP>(12);
	}

}
