//
// Kmernator/apps/CompareSpectrums.cpp
//
// Author: Rob Egan
//
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/

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
typedef KmerSpectrum<KmerSolidMap, KmerSolidMap> KS;

using namespace std;


class _CS_Options : public OptionsBaseInterface {

	// cache of variables (for inline lookup and defaults)

public:
	_CS_Options() : referenceFiles(), circularReference(false), perRead(false) {}
	virtual ~_CS_Options() {}

	void _resetDefaults() {
		KmerBaseOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
	}


	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		// set options specific to this program
		p.add("kmer-size", 1);
		p.add("reference-file", 1);
		p.add("input-file", -1);

		po::options_description opts("CompareSpectrum options");

		opts.add_options()
					 ("reference-file", po::value<FileListType>(), "set reference file(s)")

					 ("circular-reference", po::value<bool>()->default_value(circularReference), "if set, reference file should be treated as circular")

					 ("per-read", po::value<bool>()->default_value(perRead), "if set, each read in readset1 will be compared to the entire readset2 separately");


		desc.add(opts);
		KmerBaseOptions::_setOptions(desc,p);
		KmerSpectrumOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		setOpt2("reference-file", referenceFiles);
		setOpt("circular-reference", circularReference);
		setOpt("per-read", perRead);

		ret &= GeneralOptions::_parseOptions(vm);
		ret &= KmerBaseOptions::_parseOptions(vm);
		ret &= KmerSpectrumOptions::_parseOptions(vm);
		return ret;
	}
	FileListType &getReferenceFiles()
	{
		return referenceFiles;
	}
	bool &getCircularReference() {
		return  circularReference;
	}
	bool &getPerRead() {
		return perRead;
	}

private:
	FileListType referenceFiles;
	bool circularReference, perRead;

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

	CS_Options::parseOpts(argc, argv);

	MemoryUtils::getMemoryUsage();
	ReadSet readSet1, readSet2;

	OfstreamMap om(Options::getOptions().getOutputFile(), "");
	std::ostream *outPtr = &std::cout;
	if (! Options::getOptions().getOutputFile().empty() ) {
		outPtr = &om.getOfstream("");
	}


	KmerSizer::set(KmerBaseOptions::getOptions().getKmerSize());
	OptionsBaseInterface::FileListType &fileList1 = CS_Options::getOptions().getReferenceFiles();
	OptionsBaseInterface::FileListType &fileList2 = Options::getOptions().getInputFiles();

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

	long estimatedRawKmers = std::max(KS::estimateRawKmers(readSet1),
			KS::estimateRawKmers(readSet2));

	KS ks1(estimatedRawKmers);
	KS ks2(estimatedRawKmers);
	ks1.setSolidOnly();
	ks2.setSolidOnly();

	LOG_VERBOSE(1, "Building map 2");
	ks2.buildKmerSpectrum(readSet2, true);
	ks2.optimize();

	*outPtr << endl;
	*outPtr << "Set 1\tSet 2\tCommon\t%Uniq1\t%Tot1\t%Uniq2\t%Tot2\n";


	if (CS_Options::getOptions().getPerRead()) {
		evaluatePerRead(*outPtr, ks1, ks2, readSet1);
	} else {
		LOG_VERBOSE(1, "Building map 1");
		ks1.buildKmerSpectrum(readSet1, true);
		ks1.optimize();
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
		ks1.optimize();

		evaluate(os, ks1,ks2,read.getName());
	}
}
