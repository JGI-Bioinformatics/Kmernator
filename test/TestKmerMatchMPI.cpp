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

class _KmerMatchTestOptions : public OptionsBaseInterface {
public:
	void _resetDefaults() {}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {
		GeneralOptions::getOptions()._setOptions(desc,p);
		MatcherInterfaceOptions::_setOptions(desc,p);
		KmerMatchOptions::_setOptions(desc,p);
		MPIOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::getOptions()._parseOptions(vm);
		ret &= MatcherInterfaceOptions::_parseOptions(vm);
		ret &= KmerMatchOptions::_parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);
		return ret;
	}
};

typedef OptionsBaseTemplate< _KmerMatchTestOptions > KmerMatchTestOptions;

bool containsRead(const ReadSet &rs, const Read &read) {
	bool hasRead = false;
	if (read.getFirstMarkupLength() != read.getLength())
		return true; // can not expect a trimmed read to return as a match!
	std::string name = read.getName();
	for(ReadSet::ReadSetSizeType i = 0; i < rs.getSize(); i++) {
		if (name.compare(rs.getRead(i).getName()) == 0)
			hasRead = true;
	}
	return hasRead;
}

bool testMatchesSelf(mpi::communicator &world, ReadSet &q, ReadSet &t, bool includePair) {
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "running testMatchesSelf on " << q.getSize() << " " << t.getSize());
	bool passed = true;
	KmerMatch matcher(world, t, includePair);
	//std::cout << matcher.getKmerSpectrum().solid.toString() << std::endl;
	MatcherInterface::MatchReadResults results = matcher.match(q);
	passed &= results.size() == q.getSize();
	for(ReadSet::ReadSetSizeType idx = 0; idx < q.getSize(); idx++) {
		ReadSet &rs = results[idx];
		Read &r = q.getRead(idx);
		//std::cout << passed << "\t" << rs.getSize() << "\t" << r.getName();
		if (!containsRead( rs, r )) {
			//std::cout << "\tmissing self match!";
			passed = false;
		}
		//std::cout << std::endl;
	}
	if (!passed)
		LOG_WARN(1, "Failed testMatchesSelf!");
	return passed;
}
int main(int argc, char **argv)
{

	mpi::communicator world = initializeWorldAndOptions< KmerMatchTestOptions >(argc, argv);

	ReadSet reads, reads2, greads, greads2;
	reads.appendAnyFile("10.fastq");
	reads.identifyPairs();
	greads.appendAnyFile("10.fastq", "", world.rank(), world.size());
	greads.identifyPairs();
	setGlobalReadSetOffsets(world, greads);

	reads2.appendAnyFile("1000.fastq");
	reads2.identifyPairs();
	greads2.appendAnyFile("1000.fastq", "", world.rank(), world.size());
	greads2.identifyPairs();
	setGlobalReadSetOffsets(world, greads2);

	KmerSizer::set(21);

	bool passed = true;

	passed &= testMatchesSelf(world, greads, greads, false);
	//passed &= testMatchesSelf(world, greads, greads, true);
	//passed &= testMatchesSelf(world, greads2, greads2, false);
	passed &= testMatchesSelf(world, greads2, greads2, true);

	MPI_Finalize();
	return passed ? 0 : -1;
}

