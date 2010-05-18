// $Header: /repository/PI_annex/robsandbox/KoMer/test/HashTester.cpp,v 1.4 2010-05-18 20:50:21 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "Options.h"
#include "ReadSet.h"
#include "KmerSpectrum.h"
#include "Utils.h"

#include <boost/lexical_cast.hpp>

using namespace std;

typedef TrackingDataMinimal4f DataType;
typedef KmerSpectrum<DataType, DataType> KS;

int main(int argc, char *argv[]) {
	if (!Options::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	MemoryUtils::getMemoryUsage();
	cerr << MemoryUtils::getMemoryUsage() << endl;

	ReadSet reads;
	KmerSizer::set(Options::getKmerSize());

	Options::FileListType inputs = Options::getInputFiles();
	cerr << "Reading Input Files" << endl;
	reads.appendAllFiles(inputs);
	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
			<< " Bases " << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	KS spectrumSolid(0), spectrumNormal(0), spectrumParts(0);

	if (Options::getKmerSize() > 0) {

	  long numBuckets = 64*64;
	  cerr << "targeting " << numBuckets << " buckets for reads " << endl;


	  cerr << MemoryUtils::getMemoryUsage() << endl;

	  TrackingData::minimumWeight = 0;
	  TrackingData::resetGlobalCounters();

	  cerr << "building solid spectrum " << endl << MemoryUtils::getMemoryUsage() << endl;
	  spectrumSolid = KS(numBuckets);
	  spectrumSolid.buildKmerSpectrum(reads, true);
	  spectrumSolid.printHistograms();
	  cerr << MemoryUtils::getMemoryUsage() << endl;


	  for(Kmer::IndexType i = 0; i < spectrumSolid.solid.getNumBuckets(); i++) {
		  cerr << i << ": " << spectrumSolid.solid.getBucket(i).size() << endl;
	  }

	  TrackingData::resetGlobalCounters();
	  cerr << "building normal spectrum " << endl << MemoryUtils::getMemoryUsage() << endl;
	  spectrumNormal = KS(numBuckets);
	  spectrumNormal.buildKmerSpectrum(reads);
	  cerr << MemoryUtils::getMemoryUsage() << endl;
	  for(Kmer::IndexType i = 0; i < spectrumNormal.weak.getNumBuckets(); i++) {
		  cerr << i << ": " << spectrumNormal.weak.getBucket(i).size() << endl;
	  }

	  TrackingData::resetGlobalCounters();
	  cerr << "building normal spectrum in parts" << endl << MemoryUtils::getMemoryUsage() << endl;
	  spectrumParts = KS(numBuckets);
	  KoMer::MmapFileVector mmaps = spectrumParts.buildKmerSpectrumInParts(reads, Options::getBuildPartitions());
	  cerr << MemoryUtils::getMemoryUsage() << endl;
	  for(Kmer::IndexType i = 0; i < spectrumParts.weak.getNumBuckets(); i++) {
		  cerr << i << ": " << spectrumParts.weak.getBucket(i).size() << endl;
	  }
	  cerr << MemoryUtils::getMemoryUsage() << endl;

//
//	  int threads = 8;
//	  for(int thread = 0; thread < threads; thread++) {
//		  spectrum.solid.clear(false);
//		  cerr << MemoryUtils::getMemoryUsage() << endl;
//		  spectrum.buildKmerSpectrum(reads, true, thread, threads);
//		  for(Kmer::NumberType i = 0; i < spectrum.solid.getNumBuckets(); i++) {
//			  cerr << i << ": " << spectrum.solid.getBucket(i).size() << endl;
//		  }
//	  }
//
//	  for(Kmer::NumberType i = 0; i < spectrum.solid.getNumBuckets(); i++) {
//		  cerr << i << ": " << spectrum.solid.getBucket(i).size() << endl;
//	  }
	  if (KmerSizer::getSequenceLength() <= 8) {
		  for(KS::WeakIterator it = spectrumNormal.weak.begin(); it != spectrumNormal.weak.end(); it++) {
			  cerr << it->key().toFasta() << "\t" << it->value().getCount() << endl;
		  }
	  }
	}
}

// $Log: HashTester.cpp,v $
// Revision 1.4  2010-05-18 20:50:21  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.3.2.1  2010-05-07 22:59:41  regan
// refactored base type declarations
//
// Revision 1.3  2010-05-06 21:46:51  regan
// merged changes from PerformanceTuning-20100501
//
//
