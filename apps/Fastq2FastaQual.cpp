//
// Kmernator/apps/Fastq2FastaQual.cpp
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
