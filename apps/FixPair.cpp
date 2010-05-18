// $Header: /repository/PI_annex/robsandbox/KoMer/apps/FixPair.cpp,v 1.4 2010-05-18 20:50:18 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "Options.h"

using namespace std;

class FixPair : Options {
public:
	static int getSplitSizeMegaBase() {
		return getVarMap()["split-size-mbase"].as<int> ();
	}
	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("input-file", -1);

		bool ret = Options::parseOpts(argc, argv);

		return ret;
	}
};


int main(int argc, char *argv[]) {
	if (!FixPair::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	Options::FileListType inputs = Options::getInputFiles();
	std::string outputFilename = Options::getOutputFile();
	if (outputFilename.empty())
		throw std::invalid_argument("Please specify an --ouput-file");
	ReadSet reads;
	cerr << "Reading Input Files" << endl;
	reads.appendAllFiles(inputs);
	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
		<< " Bases " << endl;
	reads.identifyPairs();

	OfstreamMap ofmap = OfstreamMap(outputFilename, ".fastq");

	string filekey;
	for(ReadSet::ReadSetSizeType pairIdx = 0 ; pairIdx < reads.getPairSize(); pairIdx++) {
		const ReadSet::Pair &pair = reads.getPair(pairIdx);

		std::string read1Label = "";
		if (reads.isValidRead(pair.read1)) {
			const Read &read = reads.getRead(pair.read1);
			if (! read.isMmaped() )
				throw std::invalid_argument(read.getName());
			KoMer::RecordPtr record = read.getRecord();
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
				throw std::invalid_argument(read.getName());
			KoMer::RecordPtr record = read.getRecord();
			SequenceRecordParser::nextLine(read2Label, record);
			size_t pos = read2Label.find_first_of(" \t");
			if (pos == std::string::npos)
				read2Label.erase();
			else
			    read2Label.erase( 0, pos +1);
		}

		reads.write(ofmap, pair, MAX_SEQUENCE_LENGTH, read1Label, MAX_SEQUENCE_LENGTH, read2Label, 2, true);
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
