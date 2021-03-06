/*
 * Meraculous.h
 *
 *  Created on: Sep 16, 2011
 *      Author: regan
 */
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

#ifndef MERACULOUS_H_
#define MERACULOUS_H_

#include "KmerTrackingData.h"
#include "DistributedFunctions.h"

class _MeraculousOptions : public OptionsBaseInterface {
public:
	// use to set/overrided any defaults on options that are stored persistently
	void _resetDefaults() {
		TrackingData::setMinimumWeight(0.0);
	}
	// use to set the description of all options
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {

		p.add("input-file", -1);
		po::options_description opts("Meraculous Options");
		//opts.add_options();

		desc.add(opts);
	}
	// use to post-process options, returning true if everything is okay
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		return ret;
	}
};
typedef OptionsBaseTemplate< _MeraculousOptions > MeraculousOptions;
typedef ExtensionTrackingData DataType;
typedef ExtensionTrackingDataSingleton SingletonDataType;



//typedef BucketExposedMap<KmerInstance, DataType, boost::unordered_map<KmerInstance, DataType, KmerHasher>, KmerHasher > MapType;
//typedef KmerMap< DataType > DefaultMapType;
//typedef KmerMap< SingletonDataType > SDefaultMapType;
//typedef KmerMapBoost< DataType > DefaultMapType;
//typedef KmerMapBoost< SingletonDataType > SDefaultMapType;
typedef KmerMapByKmerArrayPair< DataType > DefaultMapType;
typedef KmerMapByKmerArrayPair< SingletonDataType > SDefaultMapType;

typedef DefaultMapType MapType;
typedef SDefaultMapType SMapType;


typedef DistributedKmerSpectrum<MapType, MapType, SMapType> _MeraculousDistributedKmerSpectrum;
class MeraculousDistributedKmerSpectrum : public _MeraculousDistributedKmerSpectrum {
public:
	typedef _MeraculousDistributedKmerSpectrum::WeakMapType WeakMapType;
	typedef WeakMapType::Iterator KmerIter;
	MeraculousDistributedKmerSpectrum(mpi::communicator &_world, unsigned long buckets = 0, bool separateSingletons = true) :
		_MeraculousDistributedKmerSpectrum(_world, buckets, separateSingletons) {

	}
	~MeraculousDistributedKmerSpectrum() {}

	void dumpCounts(std::string filename, int minDepth) {
		DistributedOfstreamMap dom(getWorld(), filename, "");
		ostream &os = dom.getOfstream("");
		for(KmerIter it = weak.begin(); it != weak.end(); it++) {
			if ((int) it.value().getCount() < minDepth)
				continue;
			long count = it.value().getDirectionBias();
			long revCount = it.value().getCount() - count;
			os << it.key().toFasta() << "\t" << count+revCount << std::endl;
			TEMP_KMER(rev);
			it.key().buildReverseComplement(rev);
			os << rev.toFasta() << "\t" << revCount+count << std::endl;
		}
	}
	void dumpGraphs(std::string filename, int minDepth) {
		DistributedOfstreamMap dom(getWorld(), filename, "");
		ostream &os = dom.getOfstream("");
		for(KmerIter it = weak.begin(); it != weak.end(); it++) {
			if ((int) it.value().getCount() < minDepth)
				continue;
			ExtensionTracking extensions = it.value().getExtensionTracking();
			os << it.key().toFasta() << "\t" << extensions.toTextValues() << std::endl;
			TEMP_KMER(rev);
			it.key().buildReverseComplement(rev);
			os << rev.toFasta() << "\t" << extensions.getReverseComplement().toTextValues() << std::endl;
		}
	}

};

#endif /* MERACULOUS_H_ */
