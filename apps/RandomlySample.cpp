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
	static int getMaxPercentForFseek() {
		return getVarMap()["max-percent-for-fseek"].as<int> ();
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
								("min-bytes-per-record", po::value<int>()->default_value(2000), "The minimum number of bytes between two records (should be >2x greatest record size)")
								("max-percent-for-fseek", po::value<int>()->default_value(15), "The estimated maximum % of reads to select by fseek instead of blocked reads");
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

typedef std::vector<unsigned long> Positions;
Positions selectRandomDense(unsigned long numSamples, unsigned long limit) {

	if (numSamples > limit) {
		LOG_WARN(1, "Requesting more samples than available.  requested samples: " << numSamples << " available: " << limit);
		Positions positions;
		positions.reserve(limit);
		for(unsigned long i = 0; i < limit; i++)
			positions.push_back(i);
		return positions;
	}
	std::set<unsigned long> positions;
	if (numSamples * 2 > limit) {
		// select positions to remove
		for(unsigned long i = 0; i < limit ; i++) {
			positions.insert(i);
		}
		while (positions.size() > numSamples) {
			positions.erase( (unsigned long) (LongRand::rand() % limit) );
		}

	} else {
		// select positions to add

		while (positions.size() < numSamples) {
			positions.insert( (unsigned long) (LongRand::rand() % limit) );
		}
	}
	Positions list(positions.begin(), positions.end());
	return list;
}

Positions selectRandom(unsigned long numSamples, unsigned long limit, unsigned long minSpacing, unsigned long maxEdgeSpacing) {
	if (minSpacing == 1 && numSamples > limit / 2) {
		return selectRandomDense(numSamples, limit);
	}
	LOG_DEBUG(1, "selectRandom(" << numSamples << ", " << limit << ", " << minSpacing << ", " << maxEdgeSpacing << ")");
	Positions positions;
	positions.reserve(numSamples);
	unsigned long attempts = 0;
	unsigned long maxAttempts =  (log(numSamples)/log(2)+3)*2;
	while (positions.size() < numSamples && attempts++ < maxAttempts) {
		long newSamples = numSamples - positions.size();
		for(long i = 0; i < newSamples; i++)
			positions.push_back( LongRand::rand() % limit );
		std::sort(positions.begin(), positions.end());
		unsigned long lastPos = positions[0];
		std::vector<long> deleteThese;
		for(long i = 1 ; i < (long) positions.size(); i++) {
			if (positions[i] < lastPos + minSpacing || positions[i] >= limit - (maxEdgeSpacing)) {
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
	std::sort(positions.begin(), positions.end());

	if (positions.size() < numSamples) {
		LOG_WARN(1, "Could not find " << numSamples << ", attempting only " << positions.size());
	}
	return positions;
}

long pickByBlock(ReadFileReader &rfr, long numSamples) {
	bool byPair = (RSOptions::getOptions().getByPair() == 1);
	long numBlocks = std::min((long) 100, numSamples / 5);
	long numPicksPerBlock = (numSamples / numBlocks);
	long count = 0;

	LOG_DEBUG(1, "pickByBlock: " << numSamples << " blocks: " << numBlocks << " picksPerBlock: " << numPicksPerBlock);

	OfstreamMap *ofm = NULL;
	if (!Options::getOptions().getOutputFile().empty()) {
		ofm = new OfstreamMap(Options::getOptions().getOutputFile(), "");
	}
	ostream &output = (ofm == NULL ? std::cout : ofm->getOfstream(""));

	for(long block = 0; block < numBlocks; block++) {
		if (numSamples <= 0)
			break;

		ReadSet reads;
		LOG_DEBUG(1, "Reading block " << block << " of " << numBlocks);
		reads.appendFasta(rfr, block, numBlocks);
		if (byPair)
			reads.identifyPairs();
		long numRead = byPair ? reads.getPairSize() : reads.getSize();
		numPicksPerBlock = std::min(numPicksPerBlock+1, numSamples);
		long numPicks = std::min(numRead, numPicksPerBlock);
		LOG_DEBUG(1, "Read " << numRead << (byPair?" pairs" : " reads") << " selecting " << numPicks << " remaining: " << numSamples);
		Positions positions = selectRandom(numPicks, numRead, 1, 0);
		// write
		LOG_DEBUG(1, "Writing " << positions.size());
		for(Positions::iterator it = positions.begin(); it != positions.end(); it++) {
			if (byPair)
				reads.write(output, reads.getPair(*it));
			else
				reads.write(output, *it);
		}
		count += positions.size();
		numSamples -= positions.size();
	}
	if (ofm != NULL)
		delete ofm;

	return count;
}
long pickBySeeks(ReadFileReader &rfr, unsigned long numSamples, unsigned long minBytes) {

	LOG_DEBUG(1, "pickBySeeks: " << numSamples << ", " << minBytes);
	unsigned long fileSize = rfr.getFileSize();
	Positions positions = selectRandom(numSamples, fileSize, minBytes, minBytes * 10);

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

	long count = 0;
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

	return count;
}

int main(int argc, char *argv[]) {
	RSOptions::parseOpts(argc, argv);

	std::srand(static_cast<unsigned>(std::time(0)));
	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();
	std::string file = inputs[0];
	LOG_VERBOSE(1, "Selecting Input File positions");
	ReadFileReader rfr(file, "");

	unsigned long fileSize = rfr.getFileSize();
	LOG_DEBUG(1, "FileSize of " << file << " is " << fileSize);
	unsigned long minBytes = RSOptions::getOptions().getMinBytesPerRecord();
	unsigned long numSamples = RSOptions::getOptions().getNumSamples();

	// if more than 15% of the file is estimated to be requested, read it sequentially
	// and pick reads randomly
	unsigned long count = 0;
	if (numSamples * minBytes * 100 > fileSize * RSOptions::getOptions().getMaxPercentForFseek()) {
		count = pickByBlock(rfr, numSamples);
	} else {
		count = pickBySeeks(rfr, numSamples, minBytes);
	}
	if (count < numSamples) {
		LOG_WARN(1, "Could not select all samples. " << count << " selected.");
	}
	LOG_DEBUG(1, "wrote " << count);


}
