//
// Kmernator/test/HashTester.cpp
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

class _HashTesterOptions : public OptionsBaseInterface {
public:
	void _resetDefaults() {
		KmerBaseOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		p.add("kmer-size", 1);
		p.add("input-file", -1);
		KmerBaseOptions::_setOptions(desc,p);
		KmerSpectrumOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= KmerBaseOptions::_parseOptions(vm);
		ret &= KmerSpectrumOptions::_parseOptions(vm);
		ret &= GeneralOptions::_parseOptions(vm);
		return ret;
	}
};
typedef OptionsBaseTemplate< _HashTesterOptions > HashTesterOptions;

int main(int argc, char *argv[]) {
	HashTesterOptions::parseOpts(argc, argv);

	MemoryUtils::getMemoryUsage();
	cerr << MemoryUtils::getMemoryUsage() << endl;

	ReadSet reads;

	OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
	cerr << "Reading Input Files" << endl;
	reads.appendAllFiles(inputs);
	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
					<< " Bases " << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	KS spectrumSolid(0), spectrumNormal(0), spectrumParts(0);

	if (KmerBaseOptions::getOptions().getKmerSize() > 0) {

		long numBuckets = 64*64;
		cerr << "targeting " << numBuckets << " buckets for reads " << endl;


		cerr << MemoryUtils::getMemoryUsage() << endl;

		TrackingData::setMinimumWeight(0);
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
		Kmernator::MmapFileVector mmaps = spectrumParts.buildKmerSpectrumInParts(reads, KmerSpectrumOptions::getOptions().getBuildPartitions());
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
