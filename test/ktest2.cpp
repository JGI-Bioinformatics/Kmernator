//
// Kmernator/test/ktest2.cpp
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
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerReadUtils.h"
#include "Utils.h"
#include "FilterKnownOddities.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "MemoryUtils.h"
#include "Options.h"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;

typedef KmerSpectrum<TrackingData, TrackingDataWithAllReads> KS;

class _Ktest2Options : public OptionsBaseInterface {
public:
	void _resetDefaults() {
		KmerBaseOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();
		FilterKnownOdditiesOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		KmerBaseOptions::_setOptions(desc,p);
		KmerSpectrumOptions::_setOptions(desc,p);
		FilterKnownOdditiesOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc, p);
		p.add("kmer-size", 1);
		p.add("input-file", -1);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= KmerBaseOptions::_parseOptions(vm);
		ret &= KmerSpectrumOptions::_parseOptions(vm);
		ret &= FilterKnownOdditiesOptions::_parseOptions(vm);
		ret &= GeneralOptions::_parseOptions(vm);
		return ret;
	}
};
typedef OptionsBaseTemplate< _Ktest2Options > Ktest2Options;

int main(int argc, char *argv[]) {

	Ktest2Options::parseOpts(argc, argv);

	ReadSet refReads, reads;

	MemoryUtils::getMemoryUsage();
	cerr << MemoryUtils::getMemoryUsage() << endl;

	OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();

	TrackingDataWithAllReads test;
	test.track(0.99, true, 1, 2);
	KmerArray<TrackingDataWithAllReads> test2;
	TEMP_KMER(blah);
	test2.insertAt(0, blah);
	test2.valueAt(0).track(0.98, true, 2, 3);

	cerr << MemoryUtils::getMemoryUsage() << endl;
	cerr << "loaded " << refReads.getSize() << " Reads, "
			<< refReads.getBaseCount() << " Bases " << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	unsigned long numBuckets = KS::estimateWeakKmerBucketSize(refReads, 256);
	cerr << "targeting " << numBuckets << " buckets for reference " << endl;
	KS refSpectrum(numBuckets);
	refSpectrum.buildKmerSpectrum(refReads, true);
	TrackingData::resetGlobalCounters();
	cerr << MemoryUtils::getMemoryUsage() << endl;

	cerr << "Reading Input Files" << endl;
	reads.appendAllFiles(inputs);
	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
					<< " Bases " << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	{
		cerr << "Preparing filter" << endl;
		FilterKnownOddities filter;
		cerr << MemoryUtils::getMemoryUsage() << endl;

		cerr << "Applying filter to Input Files" << endl;
		unsigned long filtered = filter.applyFilter(reads);
		cerr << "filter affected " << filtered << " Reads " << endl;
		cerr << MemoryUtils::getMemoryUsage() << endl;
		cerr << "Purging filter" << endl;
	}
	cerr << MemoryUtils::getMemoryUsage() << endl;

	numBuckets = KS::estimateWeakKmerBucketSize(reads, 64);
	cerr << "targeting " << numBuckets << " buckets for reads " << endl;

	KS spectrum(numBuckets);
	cerr << MemoryUtils::getMemoryUsage() << endl;


	cerr << "building spectrum" << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;
	spectrum.buildKmerSpectrum(reads);
	cerr << MemoryUtils::getMemoryUsage() << endl;

	if (refReads.getSize() > 0) {
		//cerr << "Getting real error rate" << endl;
		//spectrum.getErrorRates(refSpectrum.solid);
	}

	unsigned long promoted = 0; //spectrum.autoPromote();//spectrum.promote( Options::getOptions().getSolidQuantile() );

	cerr << "Promoted " << promoted << " kmers" << endl;

	if (refReads.getSize() > 0) {
		cerr << "Contrasted to reference kmer-spectrum:" << endl;
		cerr << spectrum.contrastSpectrums(cerr, refSpectrum) << endl;
	} else if (Options::getOptions().getVerbose() > 0) {
		cerr << "Dumping kmer spectrum" << endl;
		cerr << "Solid:" << endl;
		for (KS::SolidMapType::Iterator it = spectrum.solid.begin(); it
		!= spectrum.solid.end(); it++) {
			if (it->value().getCount() > 15
					&& (it->value().getNormalizedDirectionBias() > 0.9
							|| it->value().getNormalizedDirectionBias() < 0.1))
				cerr << "\t" << spectrum.pretty(it->key(),
						it->value().toString());
		}
		cerr << "Weak:" << endl;
		for (KS::WeakMapType::Iterator it = spectrum.weak.begin(); it
		!= spectrum.weak.end(); it++) {
			if (it->value().getCount() > 15
					&& (it->value().getNormalizedDirectionBias() > 0.9
							|| it->value().getNormalizedDirectionBias() < 0.1))
				cerr << "\t" << spectrum.pretty(it->key(),
						it->value().toString());
		}
	}
}
