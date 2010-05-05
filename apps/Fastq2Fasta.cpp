// $Header: /repository/PI_annex/robsandbox/KoMer/apps/Fastq2Fasta.cpp,v 1.4 2010-05-05 06:28:40 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "ReadSet.h"

using namespace std;

class Fastq2FastaOptions : Options {
public:
	static int getSplitSizeMegaBase() {
		return getVarMap()["split-size-mbase"].as<int> ();
	}
	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("input-file", -1);

		getDesc().add_options()

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

	long currentBase = 0;
	OfstreamMap ofmap;
	string outputFilename = Options::getOutputFile();
	bool hasOfMap = false;
	ostream *out = &cout;

	int partitionNum = 1;
	if (!outputFilename.empty()) {
		ofmap = OfstreamMap(outputFilename, ".fasta");
		hasOfMap = true;
	} else {
		splitSizeBase = 0; // do not support splitting when no output is specified
	}

	string filekey;
	for(ReadSet::ReadSetSizeType idx = 0 ; idx < reads.getSize(); idx++) {
		const Read &read = reads.getRead(idx);

		if (hasOfMap)
		    filekey = reads.getReadFileNamePrefix(idx);

		if (splitSizeBase > 0) {
			currentBase += read.getLength();
			if (currentBase > splitSizeBase) {
			  // new output handle
			   partitionNum++;
			   currentBase = read.getLength();
			}
			filekey += "-" + boost::lexical_cast<string>( partitionNum );
		}
		if (hasOfMap)
			 out = &( ofmap.getOfstream(filekey) );

		read.write(*out, Sequence::MAX_SEQUENCE_LENGTH, "", 3);
	}

}
