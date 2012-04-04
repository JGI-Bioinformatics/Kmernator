//
// Kmernator/apps/FixPair.cpp
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

class _FixPair : public OptionsBaseInterface {
public:
	static int getSplitSizeMegaBase() {
		return getVarMap()["split-size-mbase"].as<int> ();
	}
	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		p.add("input-file", -1);
		GeneralOptions::_setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = GeneralOptions::_parseOptions(vm);
		return ret;
	}
};
typedef OptionsBaseTemplate< _FixPair > FixPair;

int main(int argc, char *argv[]) {
	FixPair::parseOpts(argc, argv);

	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();
	std::string outputFilename = Options::getOptions().getOutputFile();
	if (outputFilename.empty())
		LOG_THROW("Invalid: Please specify an --ouput-file");
	ReadSet reads;
	LOG_VERBOSE(1, "Reading Input Files");
	reads.appendAllFiles(inputs);
	LOG_VERBOSE(1,"loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
			<< " Bases ");
	reads.identifyPairs();

	OfstreamMap ofmap = OfstreamMap(outputFilename, ".fastq");

	string filekey;
	for(ReadSet::ReadSetSizeType pairIdx = 0 ; pairIdx < reads.getPairSize(); pairIdx++) {
		const ReadSet::Pair &pair = reads.getPair(pairIdx);

		std::string read1Label = "";
		if (reads.isValidRead(pair.read1)) {
			const Read &read = reads.getRead(pair.read1);
			if (! read.isMmaped() )
				LOG_THROW("Invalid read: " << read.getName());
			Kmernator::RecordPtr record = read.getRecord();
			SequenceRecordParser::nextLine(read1Label, record);
			size_t pos = read1Label.find_first_of(" \t");
			if (pos == std::string::npos)
				read1Label.erase();
			else
				read1Label.erase(0, pos +1);
		}
		std::string read2Label = "";
		if (reads.isValidRead(pair.read2)) {
			const Read &read = reads.getRead(pair.read2);
			if (! read.isMmaped() )
				LOG_THROW("Invalid: read: " << read.getName());
			Kmernator::RecordPtr record = read.getRecord();
			SequenceRecordParser::nextLine(read2Label, record);
			size_t pos = read2Label.find_first_of(" \t");
			if (pos == std::string::npos)
				read2Label.erase();
			else
				read2Label.erase( 0, pos +1);
		}

		reads.write(ofmap, pair, 0, MAX_SEQUENCE_LENGTH, read1Label, 0, MAX_SEQUENCE_LENGTH, read2Label, FormatOutput::FastqUnmasked(), true);
	}

}
