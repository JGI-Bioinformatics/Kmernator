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
#include "Log.h"

using namespace std;

class _Fastq2FastaOptions : public OptionsBaseInterface {
public:
	static int getSplitSizeMegaBase() {
		return getVarMap()["split-size-mbase"].as<int> ();
	}
	static int getSplitPairs() {
		return getVarMap()["split-pairs"].as<int>();
	}
	void _resetDefaults() {
		GeneralOptions::getOptions()._resetDefaults();
		GeneralOptions::getOptions().getFormatOutput() = 3;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		// override the default output format!

		// set options specific to this program
		p.add("input-file", -1);

		po::options_description opts("Fastq to Fasta Options");
		opts.add_options()

				("split-pairs", po::value<int>()->default_value(0), "if set, pairs will be directed into separate files")

				("split-size-mbase", po::value<int>()->default_value(0), "maximum size of output fastas.  requires --output-file");

		desc.add(opts);

		GeneralOptions::getOptions()._setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {
		return GeneralOptions::getOptions()._parseOptions(vm);
	}
};
typedef OptionsBaseTemplate< _Fastq2FastaOptions > Fastq2FastaOptions;

int main(int argc, char *argv[]) {
	Fastq2FastaOptions::parseOpts(argc, argv);

	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();
	long splitSizeBase = Fastq2FastaOptions::getOptions().getSplitSizeMegaBase() * 1000000;

	ReadSet reads;
	LOG_VERBOSE(1, "Reading Input Files" );
	reads.appendAllFiles(inputs);

	LOG_VERBOSE(1, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
			<< " Bases ");

	reads.identifyPairs();

	long currentBase = 0;
	OfstreamMap ofmap;
	string outputFilename = Options::getOptions().getOutputFile();
	bool hasOfMap = false;
	ostream *out = &cout;

	int partitionNum = 1;
	if (!outputFilename.empty()) {
		ofmap = OfstreamMap(outputFilename);
		hasOfMap = true;
	} else {
		splitSizeBase = 0; // do not support splitting when no output is specified
	}

	bool splitPairs = Fastq2FastaOptions::getOptions().getSplitPairs() != 0;
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
