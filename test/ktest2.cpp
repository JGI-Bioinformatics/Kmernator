// $Header: /repository/PI_annex/robsandbox/KoMer/test/ktest2.cpp,v 1.14 2009-11-04 18:26:18 regan Exp $
//

#include <iostream>

#include <cstdlib>
#include <cstring>

#include <tr1/unordered_map>

#include "ReadSet.h"
#include "Kmer.h"
#include "Utils.h"
#include "MemoryUtils.h"

using namespace std;

int main(int argc, char *argv[]) {
    
    ReadSet store;

    cerr << MemoryUtils::getMemoryUsage() << endl;
    
    KmerSizer::set(atoi(argv[1]));
    for (int i = 2 ; i< argc ; i++) {
      cerr << "reading " << argv[i] << endl;
      store.appendFastq(argv[i]);
      cerr << "loaded " << store.getSize() << " Reads, " << store.getBaseCount() << " Bases " << endl;
      cerr << MemoryUtils::getMemoryUsage() << endl;
    }
    
    unsigned long numBuckets = estimateWeakKmerBucketSize( store, 64 );
    cerr << "targetting " << numBuckets << " buckets " << endl;
    
    KmerSpectrum spectrum(numBuckets);
    cerr << MemoryUtils::getMemoryUsage() << endl;

    TrackingData::minimumDepth = 10;
    TrackingData::minimumWeight = 0.25;
        
    buildKmerSpectrum( store, spectrum );
    
    cerr << MemoryUtils::getMemoryUsage() << endl;
    
}


//
// $Log: ktest2.cpp,v $
// Revision 1.14  2009-11-04 18:26:18  regan
// refactored
// added statistics calculations and histograms
//
// Revision 1.13  2009-11-03 17:15:43  regan
// minor refactor
//
// Revision 1.12  2009-11-02 21:19:28  regan
// fixed types and boundary tests
//
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
