//
// Kmernator/apps/Fastq2Fasta.cpp
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
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "ReadSet.h"
#include "Options.h"
#include "Log.h"

using namespace std;

class _Fastq2FastaOptions : public OptionsBaseInterface {
public:
	_Fastq2FastaOptions(): splitPairs(false), splitSizeMbase(0) {}
	virtual ~_Fastq2FastaOptions() {}

	int &getSplitSizeMegaBase() {
		return splitSizeMbase;
	}
	bool &getSplitPairs() {
		return splitPairs;
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

				("split-pairs", po::value<int>()->default_value(splitPairs), "if set, pairs will be directed into separate files")

				("split-size-mbase", po::value<int>()->default_value(splitSizeMbase), "maximum size of output fastas.  requires --output-file");

		desc.add(opts);

		GeneralOptions::getOptions()._setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {
		setOpt("split-pairs", splitPairs);
		setOpt("split-size-mbase", splitSizeMbase);

		return GeneralOptions::getOptions()._parseOptions(vm);
	}
private:
	bool splitPairs;
	int splitSizeMbase;
};
typedef OptionsBaseTemplate< _Fastq2FastaOptions > Fastq2FastaOptions;

int main(int argc, char *argv[]) {

	if (!Fastq2FastaOptions::parseOpts(argc, argv)) exit(1);

	Cleanup::prepare();

	OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
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
