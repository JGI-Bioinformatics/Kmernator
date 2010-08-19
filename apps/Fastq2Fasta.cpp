//
// Kmernator/apps/Fastq2Fasta.cpp
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
#include "Options.h"

using namespace std;

class Fastq2FastaOptions : Options {
public:
	static int getSplitSizeMegaBase() {
		return getVarMap()["split-size-mbase"].as<int> ();
	}
	static int getSplitPairs() {
		return getVarMap()["split-pairs"].as<int>();
	}
	static bool parseOpts(int argc, char *argv[]) {
		// override the default output format!
		Options::getFormatOutput() = 3;

		// set options specific to this program
		getPosDesc().add("input-file", -1);

		getDesc().add_options()

		("split-pairs", po::value<int>()->default_value(0),
				"if set, pairs will be directed into separate files")

		("split-size-mbase", po::value<int>()->default_value(0),
				"maximum size of output fastas.  requires --output-file");

		bool ret = Options::parseOpts(argc, argv);

		return ret;
	}
};

int main(int argc, char *argv[]) {
	if (!Fastq2FastaOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	Options::FileListType inputs = Options::getInputFiles();
	long splitSizeBase = Fastq2FastaOptions::getSplitSizeMegaBase() * 1000000;

	ReadSet reads;
	cerr << "Reading Input Files" << endl;
	reads.appendAllFiles(inputs);

	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
		<< " Bases " << endl;

	reads.identifyPairs();

	long currentBase = 0;
	OfstreamMap ofmap;
	string outputFilename = Options::getOutputFile();
	bool hasOfMap = false;
	ostream *out = &cout;

	int partitionNum = 1;
	if (!outputFilename.empty()) {
		ofmap = OfstreamMap(outputFilename);
		hasOfMap = true;
	} else {
		splitSizeBase = 0; // do not support splitting when no output is specified
	}

	bool splitPairs = Fastq2FastaOptions::getSplitPairs() != 0;
	string filekey;
	for(ReadSet::ReadSetSizeType pairIdx = 0 ; pairIdx < reads.getPairSize(); pairIdx++) {
		ReadSet::Pair pair = reads.getPair(pairIdx);

		ReadSet::ReadSetSizeType lesserIdx  = std::min(pair.read1, pair.read2);

		if (hasOfMap) {
		    filekey = reads.getReadFileNamePrefix(lesserIdx);
		} else {
			filekey.clear();
		}

		if (splitSizeBase > 0) {
			SequenceLengthType len = reads.getRead(lesserIdx).getLength();
			currentBase += len;
			if (currentBase > splitSizeBase) {
			  // new output handle
			   partitionNum++;
			   currentBase = len;
			}
			filekey += "-" + boost::lexical_cast<string>( partitionNum );
		}


		if (reads.isValidRead(pair.read1) && reads.isValidRead(pair.read2)) {

			const Read read = reads.getRead(pair.read1);
			if (hasOfMap) {
				if (splitPairs) {
					filekey += "-1";
				}
				out = &( ofmap.getOfstream(filekey) );
			}

			reads.getRead(pair.read1).write(*out);
			if (splitPairs) {
				filekey[filekey.length()-1] = '2';
				out = &( ofmap.getOfstream(filekey) );
			}
			reads.getRead(pair.read2).write(*out);

		} else {
			if (hasOfMap) {
				out = &( ofmap.getOfstream(filekey) );
			}
			reads.getRead(lesserIdx).write(*out);
		}

	}

}

// $Log: Fastq2Fasta.cpp,v $
// Revision 1.7  2010-05-24 21:48:50  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.6.2.1  2010-05-19 00:20:49  regan
// refactored fomat output options
// added options to fastq2fasta
//
// Revision 1.6  2010-05-18 20:50:18  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.5.2.1  2010-05-07 22:59:29  regan
// refactored base type declarations
//
// Revision 1.5  2010-05-06 21:46:57  regan
// merged changes from PerformanceTuning-20100501
//
//
