// $Header: /repository/PI_annex/robsandbox/KoMer/test/ktest2.cpp,v 1.23 2009-12-23 07:16:50 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "ReadSet.h"
#include "Kmer.h"
#include "Utils.h"
#include "KmerSpectrum.h"
#include "MemoryUtils.h"
#include "Options.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;

typedef KmerSpectrum<TrackingData,TrackingDataWithAllReads> KS;

int main(int argc, char *argv[]) {
    
    if (!Options::parseOpts(argc,argv))
      throw std::invalid_argument("Please fix the command line arguments");
      
    ReadSet refReads, reads;

    cerr << MemoryUtils::getMemoryUsage() << endl;
    
    KmerSizer::set(Options::getKmerSize());
    Options::FileListType references = Options::getReferenceFiles();
    Options::FileListType inputs     = Options::getInputFiles();
 
    TrackingDataWithAllReads test;
    test.track(0.99, true, 1, 2);
    KmerArray<TrackingDataWithAllReads> test2;
    TEMP_KMER( blah );
    test2.insertAt(0, blah);
    test2.valueAt(0).track(0.98,true,2,3);
        
    
    cerr << "Reading Reference Files" << endl;
    refReads.appendAllFiles( references );
    cerr << "loaded " << refReads.getSize() << " Reads, " << refReads.getBaseCount() << " Bases " << endl;
    cerr << MemoryUtils::getMemoryUsage() << endl; 
    
    unsigned long numBuckets = estimateWeakKmerBucketSize( refReads, 256 );
    cerr << "targetting " << numBuckets << " buckets for reference " << endl;
    KS refSpectrum( numBuckets );
    refSpectrum.weak.clear();
    refSpectrum.buildKmerSpectrum( refReads, true );
    TrackingData::resetGlobalCounters();
    cerr << MemoryUtils::getMemoryUsage() << endl;
    
    cerr << "Reading Input Files" << endl;
    reads.appendAllFiles(inputs);
    cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount() << " Bases " << endl;
    cerr << MemoryUtils::getMemoryUsage() << endl;
    
    
    numBuckets = estimateWeakKmerBucketSize( reads, 64 );
    cerr << "targetting " << numBuckets << " buckets for reads " << endl;
    
    KS spectrum(numBuckets);
    cerr << MemoryUtils::getMemoryUsage() << endl;

    TrackingData::minimumDepth = 10;
    TrackingData::minimumWeight = 0.25;
        
    spectrum.buildKmerSpectrum( reads );
    cerr << MemoryUtils::getMemoryUsage() << endl;
    
    if (refReads.getSize() > 0) {
    	cerr << "Getting real error rate" << endl;
    	spectrum.getErrorRates(refSpectrum.solid);
    }
    
    unsigned long promoted = spectrum.autoPromote();//spectrum.promote( Options::getSolidQuantile() );
    
    cerr << "Promoted " << promoted << " kmers" << endl;
    
    if (refReads.getSize() > 0) {
    	cerr << "Contrasted to reference kmer-spectrum:" << endl;
        cerr << spectrum.contrastSpectrums( refSpectrum ) << endl;
    } else if (Options::getVerbosity() > 0){
    	cerr << "Dumping kmer spectrum" << endl;
    	cerr << "Solid:" << endl;
    	for( KS::SolidMapType::Iterator it = spectrum.solid.begin(); it != spectrum.solid.end(); it++) {
    	  if (it->value().getCount() > 15 && (it->value().getNormalizedDirectionBias() > 0.9 || it->value().getNormalizedDirectionBias() < 0.1) )
    		cerr << "\t" << spectrum.pretty( it->key(), it->value().toString() );
    	}
    	cerr << "Weak:" << endl;
    	for( KS::WeakMapType::Iterator it = spectrum.weak.begin(); it != spectrum.weak.end(); it++) {
    	  if (it->value().getCount() > 15 && (it->value().getNormalizedDirectionBias() > 0.9 || it->value().getNormalizedDirectionBias() < 0.1) )
    		cerr << "\t" << spectrum.pretty( it->key(), it->value().toString() );
    	}
    }
}


//
// $Log: ktest2.cpp,v $
// Revision 1.23  2009-12-23 07:16:50  regan
// fixed reading of fasta files
// parallelized reading of multiple files
//
// Revision 1.22  2009-11-28 01:00:10  regan
// fixed bugs and warnings
//
// Revision 1.21  2009-11-27 01:53:43  regan
// refactored and got first pass at error rate by position
//
// Revision 1.20  2009-11-26 09:03:34  regan
// refactored and stuff
//
// Revision 1.19  2009-11-22 08:16:43  regan
// some fixes some bugs... optimized vs debug vs deb4/5 give different results
//
// Revision 1.18  2009-11-21 15:58:31  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.17  2009-11-11 07:57:26  regan
// built framework for autoPromote (not working) - make_heap is broken
//
// Revision 1.16  2009-11-09 19:37:20  regan
// enhanced some debugging / analysis output
//
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
