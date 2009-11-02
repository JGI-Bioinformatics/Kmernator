// $Header: /repository/PI_annex/robsandbox/KoMer/test/ktest2.cpp,v 1.11 2009-11-02 18:50:36 regan Exp $
//

#include <iostream>

#include <cstdlib>
#include <cstring>

#include <tr1/unordered_map>

#include "Utils.h"
#include "ReadSet.h"
#include "Kmer.h"

using namespace std;

int main(int argc, char *argv[]) {
    
    ReadSet store;

    KmerSizer::set(atoi(argv[1]));
    for (int i = 2 ; i< argc ; i++) {
      cerr << "reading " << argv[i] << endl;
      store.appendFastq(argv[i]);
      cerr << "loaded " << store.getSize() << " Reads, " << store.getBaseCount() << " Bases " << endl;
    }
    
    unsigned long numBuckets = estimateWeakKmerBucketSize( store, 64 );

    cerr << "targetting " << numBuckets << endl;
    
    KmerCountMap kmerCounts( numBuckets ); // weak (possible) kmers are much larger than solid estimate
    KmerSolidMap solidKmers( numBuckets / 32 );
    
    TwoBitEncoding _kmer1[KmerSizer::getTwoBitLength()];
    TwoBitEncoding _kmer2[KmerSizer::getTwoBitLength()];
    KmerPtr kmer1(&_kmer1), kmer2(&_kmer2);
    KmerPtr least;
    
    unsigned long minDepth = 10;
    
    
    for (int i=0 ; i < store.getSize(); i++)
    {
       KmerWeights kmers = buildWeightedKmers(store.getRead(i));
  
       for (int j=0; j < kmers.size(); j++)
       {
       	  
       	  *kmer2 = kmers[j];
       	  TwoBitSequence::reverseComplement((TwoBitEncoding*)kmer1.get(), (TwoBitEncoding*)kmer2.get(), KmerSizer::getSequenceLength());
       	  
       	  bool keepDirection = *kmer1 < *kmer2;
       	  least = kmer1;
       	  if (!keepDirection)
       	     least = kmer2;
       	  
       	  if ( solidKmers.exists( *least ) ) {
       	  	// track stats
       	  	solidKmers[ *least ].value.track( kmers.valueAt(j), keepDirection );
       	  	
       	  } else if (++kmerCounts[ *least ] > minDepth) {
          	// track stats and pop out of weak hash
          	solidKmers[ *least ].value.track( kmers.valueAt(j), keepDirection );
          	
          	kmerCounts.remove(*least);
          };
       }
       if (i % 1000000 == 0) {
       	 KmerSolidMap::Iterator it = solidKmers.begin();
         cerr << i << " reads, " << solidKmers.size() << " / " << kmerCounts.size() << " kmers so far " 
              << it.bucketIndex() << endl ;//<< " : " << it.bucket().toString() << endl;
         
       }
    }
     KmerCountMap::Iterator it = kmerCounts.begin();
    cerr << store.getSize() << " reads, " << kmerCounts.size() << " kmers so far " << it.bucketIndex() << endl;//<< " : " << it.bucket().toString() << endl;
    
}


//
// $Log: ktest2.cpp,v $
// Revision 1.11  2009-11-02 18:50:36  regan
// more changes
//
// Revision 1.10  2009-10-31 00:16:38  regan
// minor changes and optimizations
//
// Revision 1.9  2009-10-30 00:51:37  regan
// bug fix and working on executable
//
// Revision 1.8  2009-10-26 17:42:26  regan
// templated KmerArray
//
// Revision 1.7  2009-10-22 07:04:03  regan
// added a few unit tests
// minor refactor
//
// Revision 1.6  2009-10-22 01:39:46  cfurman
// bug fix in kmer.h
//
// Revision 1.5  2009-10-21 06:51:37  regan
// bug fixes
// build lookup tables for twobitsequence
//
// Revision 1.4  2009-10-21 00:02:02  cfurman
// working on kmers....
//
// Revision 1.3  2009-10-20 20:56:29  cfurman
// Got it to compile!
//
// Revision 1.2  2009-10-20 17:25:53  regan
// added CVS tags
//
//
