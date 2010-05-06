// $Header: /repository/PI_annex/robsandbox/KoMer/apps/Fastq2FastaQual.cpp,v 1.2 2010-05-06 21:46:57 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "ReadSet.h"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc <= 1)
	    throw std::invalid_argument("Please specify one or more input files.\nFasta will be output to STDOUT and quals will be output to STDERR");

	ReadSet reads;
	for(int i = 1 ; i < argc; i++) {
	   reads.appendAnyFile(std::string(argv[i]));
	}

	for(ReadSet::ReadSetSizeType idx = 0 ; idx < reads.getSize(); idx++) {
		const Read read = reads.getRead(idx);

		cout << ">" << read.getName() << endl << read.getFasta() << endl;
		cerr << ">" << read.getName() << endl << read.getFormattedQuals() << endl;
	}

}

// $Log: Fastq2FastaQual.cpp,v $
// Revision 1.2  2010-05-06 21:46:57  regan
// merged changes from PerformanceTuning-20100501
//
//
