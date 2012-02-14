//
// Kmernator/test/ktest2.cpp
//
// Author: Rob Egan, Craig Furman
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
		KmerOptions::_resetDefaults();
		FilterKnownOdditiesOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		KmerOptions::_setOptions(desc,p);
		FilterKnownOdditiesOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc, p);
		p.add("kmer-size", 1);
		p.add("input-file", -1);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= KmerOptions::_parseOptions(vm);
		ret &= FilterKnownOdditiesOptions::_parseOptions(vm);
		ret &= GeneralOptions::_parseOptions(vm);
		return ret;
	}
};
typedef OptionsBaseTemplate< _Ktest2Options > Ktest2Options;

int main(int argc, char *argv[]) {

	if (!Ktest2Options::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	ReadSet refReads, reads;

	MemoryUtils::getMemoryUsage();
	cerr << MemoryUtils::getMemoryUsage() << endl;

	OptionsBaseInterface::FileListType references = Options::getOptions().getReferenceFiles();
	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();

	TrackingDataWithAllReads test;
	test.track(0.99, true, 1, 2);
	KmerArray<TrackingDataWithAllReads> test2;
	TEMP_KMER(blah);
	test2.insertAt(0, blah);
	test2.valueAt(0).track(0.98, true, 2, 3);

	cerr << MemoryUtils::getMemoryUsage() << endl;
	cerr << "Reading Reference Files" << endl;
	refReads.appendAllFiles(references);
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

//
// $Log: ktest2.cpp,v $
// Revision 1.31  2010-05-01 21:57:51  regan
// merged head with serial threaded build partitioning
//
// Revision 1.30.8.1  2010-04-26 22:52:44  regan
// more testing
//
// Revision 1.30  2010-03-02 15:04:03  regan
// modified to use options
//
// Revision 1.29  2010-02-26 13:01:21  regan
// reformatted
//
// Revision 1.28  2010-01-16 01:07:40  regan
// refactored
//
// Revision 1.27  2010-01-13 23:49:11  regan
// refactored
//
// Revision 1.26  2010-01-13 07:20:10  regan
// refactored filter
// checkpoint on read picker
//
// Revision 1.25  2010-01-08 06:24:53  regan
// refactored some code
//
// Revision 1.24  2010-01-06 15:20:27  regan
// code to screen out primers
//
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
