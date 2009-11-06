// $Header: /repository/PI_annex/robsandbox/KoMer/test/ktest2.cpp,v 1.15 2009-11-06 04:10:23 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "ReadSet.h"
#include "Kmer.h"
#include "Utils.h"
#include "MemoryUtils.h"
#include "Options.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;


int main(int argc, char *argv[]) {
    
    if (!Options::parseOpts(argc,argv))
      throw std::invalid_argument("Please fix the command line arguments");
      
    ReadSet refReads, reads;

    cerr << MemoryUtils::getMemoryUsage() << endl;
    
    KmerSizer::set(Options::getKmerSize());
    Options::FileListType references = Options::getReferenceFiles();
    Options::FileListType inputs     = Options::getInputFiles();
    
    cerr << "Reading Reference Files" << endl;
    foreach( string referenceFile, references ) {
    	cerr << "Reading Reference: " << referenceFile << endl;
    	refReads.appendFastq( referenceFile.c_str() );
    	cerr << "loaded " << refReads.getSize() << " Reads, " << refReads.getBaseCount() << " Bases " << endl;
        cerr << MemoryUtils::getMemoryUsage() << endl; 
    }
    unsigned long numBuckets = estimateWeakKmerBucketSize( refReads, 256 );
    cerr << "targetting " << numBuckets << " buckets for reference " << endl;
    KmerSpectrum refSpectrum( numBuckets );
    refSpectrum.weak.clear();
    buildKmerSpectrum( refReads, refSpectrum, true );
    cerr << MemoryUtils::getMemoryUsage() << endl;
    
    cerr << "Reading Input Files" << endl;
    foreach( string inputFile, inputs ) {
      cerr << "reading " << inputFile << endl;
      reads.appendFastq(inputFile);
      cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount() << " Bases " << endl;
      cerr << MemoryUtils::getMemoryUsage() << endl;
    }
    
    numBuckets = estimateWeakKmerBucketSize( reads, 64 );
    cerr << "targetting " << numBuckets << " buckets for reads " << endl;
    
    KmerSpectrum spectrum(numBuckets);
    cerr << MemoryUtils::getMemoryUsage() << endl;

    TrackingData::minimumDepth = 10;
    TrackingData::minimumWeight = 0.25;
        
    buildKmerSpectrum( reads, spectrum );
    cerr << MemoryUtils::getMemoryUsage() << endl;
    spectrum.promote( Options::getSolidQuantile() );
    
    if (refReads.getSize() > 0) {
    	cerr << "Contrasted to reference kmer-spectrum:" << endl;
        cerr << spectrum.contrastSpectrums( refSpectrum ) << endl;
    }
}


//
// $Log: ktest2.cpp,v $
// Revision 1.15  2009-11-06 04:10:23  regan
// refactor of cmd line option handling
// added methods to evaluate spectrums
//
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
