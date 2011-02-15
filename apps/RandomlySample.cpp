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

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "Options.h"
#include "Log.h"

using namespace std;

class RSOptions : public Options {
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
	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("input-file", -1);
		getDesc().add_options()("help", "produce help message")
				("by-pair", po::value<int>()->default_value(1), "If set, pairs are sampled, if not set, reads are sampled")
				("num-samples",  po::value<int>()->default_value(1000), "The number of samples to output")
				("min-bytes-per-record", po::value<int>()->default_value(2000), "The minimum number of bytes between two records (should be >2x greatest record size)"
				);


		bool ret = Options::parseOpts(argc, argv);
		if (getInputFiles().empty() || getInputFiles().size() > 1) {
			ret = false;
			LOG_ERROR(1, "Please specify at a single input file");
		}
		return ret;
	}
};

unsigned long longRand() {
	return (((unsigned long)(std::rand() & 0xFF)) << 56) |
		   (((unsigned long)(std::rand() & 0xFF)) << 48) |
		   (((unsigned long)(std::rand() & 0xFF)) << 40) |
		   (((unsigned long)(std::rand() & 0xFF)) << 32) |
		   (((unsigned long)(std::rand() & 0xFF)) << 24) |
		   (((unsigned long)(std::rand() & 0xFF)) << 16) |
		   (((unsigned long)(std::rand() & 0xFF)) << 8) |
		   (((unsigned long)(std::rand() & 0xFF)) );
}

int main(int argc, char *argv[]) {
	Options::getVerbosity() = 0;
	Options::getMmapInput() = 0;
	if (!RSOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	std::srand(static_cast<unsigned>(std::time(0)));
	Options::FileListType inputs = Options::getInputFiles();
	std::string file = Options::getInputFiles()[0];
	ReadSet reads;
	LOG_VERBOSE(1, "Selecting Input File positions");
	ReadFileReader rfr(file, "");

	unsigned long fileSize = rfr.getFileSize();
	LOG_DEBUG(1, "FileSize of " << file << " is " << fileSize);
	unsigned long minBytes = RSOptions::getMinBytesPerRecord();
	unsigned long numSamples = RSOptions::getNumSamples();

	if ( (numSamples * 2) > fileSize / minBytes ) {
		minBytes = 1.1 * fileSize / (numSamples * 2);
		LOG_DEBUG(1, "File seems small, resetting minBytesPerRecord to " << minBytes);
	}
	fileSize -= minBytes * 2;
	std::vector<unsigned long> positions;
	positions.reserve(numSamples);
	unsigned long attempts = 0;
	while (positions.size() < numSamples && attempts++ < numSamples) {
		long newSamples = numSamples - positions.size();
		for(long i = 0; i < newSamples; i++)
			positions.push_back( longRand() % fileSize);
		std::sort(positions.begin(), positions.end());
		long lastPos = positions[0];
		std::vector<long> deleteThese;
		for(long i = 1 ; i < (long) positions.size(); i++) {
			if (positions[i] < lastPos + minBytes) {
				deleteThese.push_back(i);
			} else {
				lastPos = positions[i];
			}
		}
		for(long i = deleteThese.size() - 1 ; i >= 0; i--) {
			std::swap(positions[deleteThese[i]], positions.back());
			LOG_DEBUG(1, "attempt " << attempts << " size " << positions.size() << " removing " << positions.back() << " from " << deleteThese[i]);
			positions.pop_back();
		}
	}

	if (positions.size() < numSamples) {
		LOG_WARN(1, "Could not find " << numSamples << ", attempting only " << positions.size());
	}

	if (Log::isDebug(1)) {
		std::stringstream ss;
		for(long i = 0 ; i < (long) positions.size(); i++)
			ss << "\t" << positions[i];
		std::string s = ss.str();
		LOG_DEBUG(1, "Picked positions(" << positions.size() << "):" << s);
	}

	bool byPair = (RSOptions::getByPair() == 1);

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

	unsigned long numRecords = 0;
	unsigned long lastPos = fileSize;
	for(long i = 0; i < (long) positions.size(); i++) {
		if (!rfr.seekToNextRecord(positions[i], byPair)) {
			break;
		}
		if (numRecords > 0 && lastPos >= rfr.getPos()) {
			if (!rfr.seekToNextRecord(std::max(rfr.getPos(), lastPos+1), byPair)) {
				break;
			}
		}

		lastPos = rfr.getPos();
		std::string name, bases, quals;
		rfr.nextRead(name, bases, quals);
		Read read(name, bases, quals);
		read.write(std::cout);
		if (byPair) {
			rfr.nextRead(name, bases, quals);
			Read read2(name, bases, quals);
			read2.write(std::cout);
		}

		numRecords++;
	}

	if (numRecords < numSamples) {
		LOG_WARN(1, "Only " << numRecords << " " << (byPair?"pairs":"reads") << " were selected, perhaps the input file was too small or irregular?");
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
