//
// Kmernator/apps/CompareSpectrums.cpp
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
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "config.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "KmerReadUtils.h"
#include "Options.h"
#include "Utils.h"
#include "Log.h"

typedef TrackingDataMinimal4 DataType;
typedef KmerMap<DataType> KmerSolidMap;
typedef KmerSpectrum<DataType, DataType> KS;

using namespace std;


class _CS_Options : public OptionsBaseInterface {

	// cache of variables (for inline lookup and defaults)

public:
	virtual ~_CS_Options() {}
	static bool getCircularReference() { return  getVarMap()["circular-reference"].as<unsigned int>() != 0; }
	static bool getPerRead() { return getVarMap()["per-read"].as<unsigned int>() != 0; }

	void _resetDefaults() {
		KmerOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
	}

	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		// set options specific to this program
		p.add("kmer-size", 1);
		p.add("reference-file", 1);
		p.add("input-file", -1);

		po::options_description opts("CompareSpectrum options");

		opts.add_options()
		 ("circular-reference", po::value<unsigned int>()->default_value(0),
				 "reference file should be treated as circular")
		 ("per-read", po::value<unsigned int>()->default_value(0),
				 "if set, each read in readset1 will be compared to the entire readset2 separately")
	    ;
		desc.add(opts);
		KmerOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= KmerOptions::_parseOptions(vm);
		return ret;
	}

};

typedef OptionsBaseTemplate< _CS_Options > CS_Options;

typedef std::vector<unsigned long> NumbersVector;

// returns:
//  0: the common count of unique kmers
//  1: the cumulative count of the common kmers from map 1
//  2: the cumulative count of the common kmers from map 2
//  3: the total cumulative count of all kmers from map 1
//  4: the total cumulative count of all kmers from map 2
static NumbersVector countCommonKmers(KmerSolidMap &m1, KmerSolidMap &m2) {
	NumbersVector ret(5);

	for (KmerSolidMap::Iterator it(m1.begin()), itEnd(m1.end()); it != itEnd; it++) {
		ret[3] += it->value().getCount();
		KmerSolidMap::ElementType element = m2.getElementIfExists(it->key());
		if (element.isValid()) {
			ret[0]++;
			ret[1] += it->value().getCount();
			ret[2] += element.value().getCount();
		}
	}
	for (KmerSolidMap::Iterator it(m2.begin()), itEnd(m2.end()); it != itEnd; it++) {
		ret[4] += it->value().getCount();
	}
	return ret;
}

void evaluate(std::ostream &os, KS &ks1, KS &ks2, std::string label = "");
void evaluatePerRead(std::ostream &os, KS &ks1, KS &ks2, ReadSet &readSet1);

int main(int argc, char *argv[]) {

	if (!CS_Options::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	MemoryUtils::getMemoryUsage();
	ReadSet readSet1, readSet2;

	OfstreamMap om(Options::getOptions().getOutputFile(), "");
	std::ostream *outPtr = &std::cout;
	if (! Options::getOptions().getOutputFile().empty() ) {
		outPtr = &om.getOfstream("");
	}


	KmerSizer::set(KmerOptions::getOptions().getKmerSize());
	OptionsBaseInterface::FileListType fileList1 = Options::getOptions().getReferenceFiles();
	OptionsBaseInterface::FileListType fileList2 = Options::getOptions().getInputFiles();

	LOG_VERBOSE(1, "Reading 1st file set:");
	readSet1.appendAllFiles(fileList1);

	LOG_VERBOSE(1, " loaded " << readSet1.getSize() << " Reads, "
			<< readSet1.getBaseCount() << " Bases ");

	if (CS_Options::getOptions().getCircularReference())
		readSet1.circularize(KmerSizer::getSequenceLength());

	LOG_VERBOSE(1, "Reading 2nd file set:");
	readSet2.appendAllFiles(fileList2);
	LOG_VERBOSE(1, " loaded " << readSet2.getSize() << " Reads, "
			<< readSet2.getBaseCount() << " Bases ");

	long buckets = std::max(KS::estimateWeakKmerBucketSize(readSet1),
			KS::estimateWeakKmerBucketSize(readSet2)) * 64;

	LOG_DEBUG(1, "Estimated bucket size: " << buckets );

	KS ks1(buckets);
	KS ks2(buckets);
	ks1.setSolidOnly();
	ks2.setSolidOnly();

	LOG_VERBOSE(1, "Building map 2");
	ks2.buildKmerSpectrum(readSet2, true);

	*outPtr << endl;
	*outPtr << "Set 1\tSet 2\tCommon\t%Uniq1\t%Tot1\t%Uniq2\t%Tot2\n";


	if (CS_Options::getOptions().getPerRead()) {
		evaluatePerRead(*outPtr, ks1, ks2, readSet1);
	} else {
		LOG_VERBOSE(1, "Building map 1");
		ks1.buildKmerSpectrum(readSet1, true);
		evaluate(*outPtr, ks1, ks2);
	}
}

void evaluate(std::ostream &os, KS &ks1, KS &ks2, std::string label) {
	LOG_VERBOSE(1, "Counting common Kmers\n");
	KmerSolidMap &m1 = ks1.solid;
	KmerSolidMap &m2 = ks2.solid;

	NumbersVector common = countCommonKmers(m1, m2);

	// TODO fix check common to iterater through same-bucketed KmerMaps in sorted order (fast)
	// TODO if perRead1, iterate through common matches and output non-zero % matches

	os << m1.size() << '\t' << m2.size() << '\t' << common[0] << '\t'
			<< setprecision(4) << (common[0] * 100.0) / m1.size() << '\t'
			<< (common[1] * 100.0 / common[3]) << "\t" << (common[0]
			* 100.0) / m2.size() << "\t" << (common[2] * 100.0 / common[4])
			<< '\t' << label
			<< endl;

}
void evaluatePerRead(std::ostream &os, KS &ks1, KS &ks2, ReadSet &readSet1) {
	for(ReadSet::ReadSetSizeType readIdx = 0; readIdx < readSet1.getSize(); readIdx++ ) {

		TrackingData::resetGlobalCounters();
		ks1.reset(false);

		ReadSet a;
		Read &read = readSet1.getRead(readIdx);
		a.append(read);
		ks1.buildKmerSpectrum(a, true);

		evaluate(os, ks1,ks2,read.getName());
	}
}
//
// $Log: CompareSpectrums.cpp,v $
// Revision 1.11  2010-05-24 21:48:49  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.10.2.2  2010-05-19 23:41:02  regan
// re-added memory optimizations
//
// Revision 1.10.2.1  2010-05-19 23:38:15  regan
// bug and performance fixes
//
// Revision 1.10  2010-05-18 20:50:18  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.9.20.4  2010-05-13 20:29:13  regan
// minor refactor
//
// Revision 1.9.20.3  2010-05-12 22:44:36  regan
// reworked option handling
//
// Revision 1.9.20.2  2010-05-10 17:57:11  regan
// fixing types
//
// Revision 1.9.20.1  2010-05-07 22:59:29  regan
// refactored base type declarations
//
// Revision 1.9  2010-03-04 06:36:36  regan
// fixed compiler warnings for non-openmp compilers
//
// Revision 1.8  2010-02-26 13:01:15  regan
// reformatted
//
// Revision 1.7  2010-02-22 14:41:31  regan
// checkpoint
//
// Revision 1.6  2010-01-14 19:27:43  regan
// bugfixes
//
// Revision 1.5  2010-01-14 00:46:51  regan
// refactor and other changes
//
// Revision 1.4  2010-01-08 18:34:51  regan
// refactored a bit
// enabled openmp parallelizations
//
// Revision 1.3  2010-01-05 06:44:37  regan
// fixed warnings
//
// Revision 1.2  2009-12-24 00:54:08  regan
// fixed reading of fasta files
// parallelized reading of multiple files
//
// Revision 1.1  2009-11-27 23:28:07  cfurman
// CompareSpectrum application added
//

