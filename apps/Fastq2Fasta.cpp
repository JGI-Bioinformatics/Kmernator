// $Header: /repository/PI_annex/robsandbox/KoMer/apps/Fastq2Fasta.cpp,v 1.1 2010-03-16 18:58:10 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "ReadSet.h"

using namespace std;

int main(int argc, char *argv[]) {
	if (argc <= 1)
	    throw std::invalid_argument("Please specify one or more input files");

	ReadSet reads;
	for(int i = 1 ; i < argc; i++) {
	   cerr << "Reading Input Files" << endl;
	   reads.appendAnyFile(std::string(argv[i]));
	}
	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
		<< " Bases " << endl;

	for(ReadSet::ReadSetSizeType idx = 0 ; idx < reads.getSize(); idx++) {
		const Read read = reads.getRead(idx);

		cout << ">" << read.getName() << endl << read.getFasta() << endl;
	}

}
