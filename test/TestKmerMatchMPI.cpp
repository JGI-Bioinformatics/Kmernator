//
// Kmernator/test/KmerMatchTest.cpp
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

#include "DistributedFunctions.h"
#include "KmerMatch.h"

// Note for versbosity: export BOOST_TEST_LOG_LEVEL=message

using namespace std;

void testKmerMatch() {


}

class _KmerMatchTestOptions : public _MatcherInterfaceOptions, public _KmerMatchOptions, public _MPIOptions {
public:
	void _resetDefaults() {}
	void _setOptions(po::options_description &desc,
				po::positional_options_description &p) {
		GeneralOptions::getOptions()._setOptions(desc,p);
		_MatcherInterfaceOptions::_setOptions(desc,p);
		_KmerMatchOptions::_setOptions(desc,p);
		_MPIOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::getOptions()._parseOptions(vm);
		ret &= _MatcherInterfaceOptions::_parseOptions(vm);
		ret &= _KmerMatchOptions::_parseOptions(vm);
		ret &= _MPIOptions::_parseOptions(vm);
		return ret;
	}
};

typedef OptionsBaseTemplate< _KmerMatchTestOptions > KmerMatchTestOptions;

int main(int argc, char **argv)
{

	mpi::communicator world = initializeWorldAndOptions< KmerMatchTestOptions >(argc, argv);


	MPI_Finalize();

}

