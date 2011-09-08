//
// Kmernator/apps/RandomlySample.cpp
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
#include <cmath>

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "Options.h"
#include "Utils.h"
#include "Log.h"

using namespace std;

class _RSOptions : public OptionsBaseInterface {
public:
	static int getByPair() {
		return getVarMap()["by-pair"].as<int> ();
	}
	static int getNumSamples() {
		return getVarMap()["num-samples"].as<int> ();
	}
	static int getMinBytesPerRecord() {
		return getVarMap()["min-bytes-per-record"].as<int> ();
	}
	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		GeneralOptions::getOptions().getVerbose() = 0;
		GeneralOptions::getOptions().getMmapInput() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		p.add("input-file", -1);
		po::options_description opts("Randomly Sample Options");
		opts.add_options()
				("by-pair", po::value<int>()->default_value(1), "If set, pairs are sampled, if not set, reads are sampled")
				("num-samples",  po::value<int>()->default_value(1000), "The number of samples to output")
				("min-bytes-per-record", po::value<int>()->default_value(900), "The minimum number of bytes between two records (should be >2x greatest record size)"
				);
		desc.add(opts);
		GeneralOptions::_setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {

		bool ret = GeneralOptions::_parseOptions(vm);

		if (Options::getOptions().getInputFiles().empty() || Options::getOptions().getInputFiles().size() > 1) {
			ret = false;
			LOG_ERROR(1, "Please specify at a single input file");
		}
		return ret;
	}
};
typedef OptionsBaseTemplate< _RSOptions > RSOptions;


int main(int argc, char *argv[]) {
	if (!RSOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	std::srand(static_cast<unsigned>(std::time(0)));
	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();
	std::string file = Options::getOptions().getInputFiles()[0];
	ReadSet reads;
	LOG_VERBOSE(1, "Selecting Input File positions");
	ReadFileReader rfr(file, "");

	unsigned long fileSize = rfr.getFileSize();
	LOG_DEBUG(1, "FileSize of " << file << " is " << fileSize);
	unsigned long minBytes = RSOptions::getOptions().getMinBytesPerRecord();
	unsigned long numSamples = RSOptions::getOptions().getNumSamples();
	if (fileSize * 0.333 < minBytes * numSamples) {
		minBytes = fileSize * 0.333 / numSamples;
		LOG_DEBUG(1, "Overriding minBytes: " << minBytes);
	}
	std::vector<unsigned long> positions;
	positions.reserve(numSamples);
	unsigned long attempts = 0;
	unsigned long maxAttempts =  (log(numSamples)/log(2)+3)*2;
	while (positions.size() < numSamples && attempts++ < maxAttempts) {
		long newSamples = numSamples - positions.size();
		for(long i = 0; i < newSamples; i++)
			positions.push_back( LongRand::rand() % fileSize);
		std::sort(positions.begin(), positions.end());
		long lastPos = positions[0];
		std::vector<long> deleteThese;
		for(long i = 1 ; i < (long) positions.size(); i++) {
			if (positions[i] < lastPos + minBytes || positions[i] >= fileSize - minBytes*50) {
				deleteThese.push_back(i);
			} else {
				lastPos = positions[i];
			}
		}
		for(long i = deleteThese.size() - 1 ; i >= 0; i--) {
			std::swap(positions[deleteThese[i]], positions.back());
			LOG_DEBUG(2, "attempt " << attempts << " size " << positions.size() << " removing " << positions.back() << " from " << deleteThese[i]);
			positions.pop_back();
		}
		LOG_DEBUG(1, "Selected " <<  positions.size() << " after removing " << deleteThese.size() << " attempt " << attempts << " of " << maxAttempts);
	}

	if (positions.size() < numSamples) {
		LOG_WARN(1, "Could not find " << numSamples << ", attempting only " << positions.size());
	}

	if (Log::isDebug(1)) {
		std::stringstream ss;
		for(long i = 0 ; i < (long) positions.size(); i++)
			ss << "\t" << positions[i];
		std::string s = ss.str();
		LOG_DEBUG(2, "Picked positions(" << positions.size() << "):" << s);
	}

	bool byPair = (RSOptions::getOptions().getByPair() == 1);

	if (rfr.seekToNextRecord(0, true)) {
		LOG_DEBUG(1, "Reading first two records to determine inherent pairing");
		std::string name1, name2, bases1, bases2, quals1, quals2;
		rfr.nextRead(name1, bases1, quals1);
		rfr.nextRead(name2, bases2, quals2);

		bool isFilePaired= ReadSet::isPair(name1,name2);
		LOG_DEBUG(1, "Reading first two records to determine inherent pairing: " << isFilePaired << " " << name1 << " " << name2);
		byPair &= isFilePaired;
	}	

	LOG_DEBUG(1, "detecting by pair: " << byPair);
	OfstreamMap *ofm = NULL;
	if (!Options::getOptions().getOutputFile().empty()) {
		ofm = new OfstreamMap(Options::getOptions().getOutputFile(), "");
	}
	ostream &output = (ofm == NULL ? std::cout : ofm->getOfstream(""));

	unsigned long count = 0;
	unsigned long lastPos = fileSize;
	for(long i = 0; i < (long) positions.size(); i++) {
		if (rfr.eof())
			break;
		rfr.seekToNextRecord(positions[i], byPair);
		unsigned long myPos = rfr.getPos();
		while (lastPos < fileSize && lastPos >= myPos) {
			LOG_DEBUG(2, "Re-seeking from " << myPos << " because lastPos is larger " << lastPos);
			rfr.seekToNextRecord(myPos + 1, byPair);
			if (rfr.eof())
				break;
			myPos = rfr.getPos();
		}
		if (rfr.eof())
			break;
		lastPos = myPos;
		std::string name, bases, quals;
		rfr.nextRead(name, bases, quals);
		Read read(name, bases, quals);
		read.write(output);
		if (byPair) {
			rfr.nextRead(name, bases, quals);
			Read read2(name, bases, quals);
			read2.write(output);
		}
		count++;
	}
	if (ofm != NULL)
		delete ofm;

	if (count < numSamples) {
		LOG_WARN(1, "Could not select all samples. " << count << " selected.");
	}

}

// $Log: FixPair.cpp,v $
// Revision 1.4  2010-05-18 20:50:18  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.3.2.1  2010-05-07 22:59:29  regan
// refactored base type declarations
//
// Revision 1.3  2010-05-06 21:46:57  regan
// merged changes from PerformanceTuning-20100501
//
//
