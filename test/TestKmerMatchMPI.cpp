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
#include "KmerAlign.h"

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
		KmerAlignOptions::_setOptions(desc,p);
		MPIOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::getOptions()._parseOptions(vm);
		ret &= MatcherInterfaceOptions::_parseOptions(vm);
		ret &= KmerMatchOptions::_parseOptions(vm);
		ret &= KmerAlignOptions::_parseOptions(vm);
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

bool testMatchesSelf(mpi::communicator &world, ReadSet &q, ReadSet &t) {
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "running testMatchesSelf on " << q.getSize() << " " << t.getSize());
	bool passed = true;
	bool isPair = MatcherInterfaceOptions::getOptions().getIncludeMate();
	KmerMatch matcher(world, t);
	//std::cout << matcher.getKmerSpectrum().solid.toString() << std::endl;
	MatcherInterface::MatchReadResults results = matcher.match(q);
	passed &= results.size() == q.getSize();
	if (!passed) {
		LOG_WARN(1, "result size not the same as query size!" << results.size() << " vs " << q.getSize());
	}
	for(ReadSet::ReadSetSizeType idx = 0; idx < q.getSize(); idx++) {
		ReadSet &rs = results[idx];
		Read &r = q.getRead(idx);
		//std::cout << passed << "\t" << rs.getSize() << "\t" << r.getName();
		if (!containsRead( rs, r )) {
			//std::cout << "\tmissing self match!";
			stringstream ss;
			for(ReadSet::ReadSetSizeType j = 0; j < rs.getSize(); j++)
				ss << rs.getRead(j).getName() << ", ";
			std::string s =ss.str();
			LOG_WARN(1, "readset " << idx << " did not contain read: " << r.toString() << " results: " << rs.getSize() << " : " << s);
			passed = false;
		}
		//std::cout << std::endl;

		KmerAlign testAlign(r);
		bool pairMatched = false;
		int pairNum = 0;
		std::string commonName;
		Alignment lastAln;
		for(ReadSet::ReadSetSizeType j = 0 ; j <rs.getSize(); j++) {
			Read query = rs.getRead(j);
			std::string thisCommonName = SequenceRecordParser::commonName(query.getName());
			if (isPair && (commonName.empty() || commonName.compare(thisCommonName) != 0)) {
				if ((!commonName.empty()) && (!pairMatched)) {
					LOG_WARN(1, "neither pair matched or only one pair! " << lastAln.toString() << " of " << r.toString() << " to " << commonName << " " << rs.getRead(j-1).toString() << " now on " << thisCommonName);
				}
				commonName = thisCommonName;
				pairNum = 0;
				pairMatched = false;
			}
			Alignment aln = testAlign.getAlignment(query);
			if ( aln.getOverlap() >= KmerSizer::getSequenceLength() ) {
				pairMatched = true;
			}
			if (isPair) {
				if (pairNum != 0 && !pairMatched) {
					LOG_WARN(1, "neither pair matched! " << aln.toString() << " of " << r.toString() << " to " << query.toString() << " or " << rs.getRead(j-1).toString() << " with " << lastAln.toString());
					passed = false;
				}
			} else {
				if (!pairMatched) {
					LOG_WARN(1, "read did not match! " << aln.toString() << " of " << r.toString() << " to " << query.toString());
					passed = false;
				}
				pairMatched = false;
			}
			lastAln = aln;
			pairNum++;
		}
	}
	if (!passed)
		LOG_WARN(1, "Failed testMatchesSelf!");
	LOG_VERBOSE(1, "Done with testMatchesSelf()");
	return passed;
}
int main(int argc, char **argv)
{

	mpi::communicator world = initializeWorldAndOptions< KmerMatchTestOptions >(argc, argv);

	TrackingData::setMinimumWeight(0.0);
	KmerOptions::getOptions().getMinKmerQuality() = 0;
	GeneralOptions::getOptions().getMinQuality() = 2;
	KmerSizer::set(21);
	MatcherInterfaceOptions::getOptions().getMaxReadMatches() = 10000;
	MatcherInterfaceOptions::getOptions().getMinOverlap() = KmerSizer::getSequenceLength();

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


	bool passed = true;

	MatcherInterfaceOptions::getOptions().setIncludeMate(true);
	passed &= testMatchesSelf(world, greads, greads);
	passed &= testMatchesSelf(world, greads2, greads2);

	MatcherInterfaceOptions::getOptions().setIncludeMate(false);
	LOG_VERBOSE(1, "Ignoring pairs for matching");
	passed &= testMatchesSelf(world, greads, greads);
	passed &= testMatchesSelf(world, greads2, greads2);

	MPI_Finalize();
	return passed ? 0 : -1;
}

