//
// Kmernator/test/KmerMatchTest.cpp
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


#include "DistributedFunctions.h"
#include "KmerMatch.h"
#include "KmerAlign.h"

// Note for versbosity: export BOOST_TEST_LOG_LEVEL=message

using namespace std;

void testKmerMatch() {


}

class _KmerMatchTestOptions : public OptionsBaseInterface {
public:
	void _resetDefaults() {
		KmerBaseOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
		MatcherInterfaceOptions::_resetDefaults();

		KmerBaseOptions::getOptions().getKmerSize() = 21;
		KmerSpectrumOptions::getOptions().getMinKmerQuality() = 0;
		KmerSpectrumOptions::getOptions().getMinDepth() = 1;
		GeneralOptions::getOptions().getMinQuality() = 2;
		MatcherInterfaceOptions::getOptions().getMaxReadMatches() = 10000;
		MatcherInterfaceOptions::getOptions().getMinOverlap() = KmerSizer::getSequenceLength();
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {
		GeneralOptions::getOptions()._setOptions(desc,p);
		MatcherInterfaceOptions::_setOptions(desc,p);
		KmerBaseOptions::_setOptions(desc,p);
		KmerSpectrumOptions::_setOptions(desc,p);
		KmerMatchOptions::_setOptions(desc,p);
		KmerAlignOptions::_setOptions(desc,p);
		MPIOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::getOptions()._parseOptions(vm);
		ret &= MatcherInterfaceOptions::_parseOptions(vm);
		ret &= KmerBaseOptions::_parseOptions(vm);
		ret &= KmerSpectrumOptions::_parseOptions(vm);
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
	LOG_VERBOSE_GATHER(1, "Done with testMatchesSelf()");
	return passed;
}
int main(int argc, char **argv)
{

	ScopedMPIComm< KmerMatchTestOptions > world(argc, argv);

	TrackingData::setMinimumWeight(0.0);

	ReadSet reads, reads2, greads, greads2;
	reads.appendAnyFile("10.fastq");
	reads.identifyPairs();
	greads.appendAnyFile("10.fastq", "", world.rank(), world.size());
	greads.identifyPairs();
	setGlobalReadSetConstants(world, greads);
	reads2.appendAnyFile("1000.fastq");
	reads2.identifyPairs();
	greads2.appendAnyFile("1000.fastq", "", world.rank(), world.size());
	greads2.identifyPairs();
	setGlobalReadSetConstants(world, greads2);


	bool passed = true;

	MatcherInterfaceOptions::getOptions().setIncludeMate(true);
	passed &= testMatchesSelf(world, greads, greads);
	passed &= testMatchesSelf(world, greads2, greads2);

	MatcherInterfaceOptions::getOptions().setIncludeMate(false);
	LOG_VERBOSE(1, "Ignoring pairs for matching");
	passed &= testMatchesSelf(world, greads, greads);
	passed &= testMatchesSelf(world, greads2, greads2);

	return passed ? 0 : -1;
}

