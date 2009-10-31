// $Header: /repository/PI_annex/robsandbox/KoMer/test/ktest2.cpp,v 1.10 2009-10-31 00:16:38 regan Exp $
//

#include <iostream>

#include <cstdlib>
#include <cstring>

#include <tr1/unordered_map>

#include "ReadSet.h"
#include "Kmer.h"

using namespace std;

unsigned long estimateBucketSize( ReadSet &store, unsigned long assumedReadCoverage = 25 ) {
	unsigned long baseCount = store.getBaseCount();
	unsigned long avgSequenceLength = baseCount / store.getSize();
	unsigned long kmersPerRead = (avgSequenceLength - KmerSizer::getSequenceLength() + 1);
	unsigned long rawKmers = kmersPerRead * store.getSize();
	unsigned long estimatedUniqueKmers =  store.getSize() / assumedReadCoverage  * kmersPerRead;
	unsigned long targetBuckets = estimatedUniqueKmers; // put avg 256 kmers per bucket
	unsigned long maxBuckets = 128*1024*1024;
	unsigned long minBuckets = 128;
	return targetBuckets > maxBuckets ? maxBuckets : (targetBuckets < minBuckets ? minBuckets : targetBuckets);
}
int main(int argc, char *argv[]) {
    
    ReadSet store;

    KmerSizer::set(atoi(argv[1]));
    for (int i = 2 ; i< argc ; i++) {
      cerr << "reading " << argv[i] << endl;
      store.appendFastq(argv[i]);
      cerr << "loaded " << store.getSize() << " Reads, " << store.getBaseCount() << " Bases " << endl;
    }
    
    unsigned long numBuckets = estimateBucketSize( store );

    cerr << "targetting " << numBuckets << endl;
    
    KmerCountMap kmerCounts( numBuckets ); // weak (possible) kmers are much larger than solid estimate
    
    for (int i=0 ; i < store.getSize(); i++)
    {
       KmerArray<WeakKmerTag> kmers(store.getRead(i).getTwoBitSequence(),store.getRead(i).getLength());
       for (int j=0; j < kmers.size(); j++)
       {
          kmerCounts[ kmers[j] ]++;
       }
       if (i % 1000000 == 0) {
       	 KmerCountMap::Iterator it = kmerCounts.begin();
         cerr << i << " reads, " << kmerCounts.size() << " kmers so far " << it.bucketIndex() << endl;//<< " : " << it.bucket().toString() << endl;
         
       }
    }
     KmerCountMap::Iterator it = kmerCounts.begin();
    cerr << store.getSize() << " reads, " << kmerCounts.size() << " kmers so far " << it.bucketIndex() << endl;//<< " : " << it.bucket().toString() << endl;
    
}


//
// $Log: ktest2.cpp,v $
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
