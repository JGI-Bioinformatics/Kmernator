// $Header: /repository/PI_annex/robsandbox/KoMer/apps/Fastq2FastaQual.cpp,v 1.1 2010-04-14 15:58:34 regan Exp $
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
