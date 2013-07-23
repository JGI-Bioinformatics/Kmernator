//
// Kmernator/apps/FixPair.cpp
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

	if (!FixPair::parseOpts(argc, argv)) exit(1);

	Cleanup::prepare();

	OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
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
