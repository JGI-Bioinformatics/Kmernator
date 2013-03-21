/*
 * KmerMatch.h
 *
 *  Created on: Nov 1, 2011
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

#ifndef KMERMATCH_H_
#define KMERMATCH_H_

#include "Kmer.h"
#include "KmerTrackingData.h"
#include "KmerSpectrum.h"
#include "DistributedFunctions.h"
#include "ReadSet.h"
#include "Options.h"
#include "MatcherInterface.h"
#include "KmerAlign.h"

class _KmerMatchOptions  : public OptionsBaseInterface {
public:
	_KmerMatchOptions() : maxPositionsFromEdge(500) {}
	~_KmerMatchOptions() {}
	int &getMatchMaxPositionsFromEdge() {
		return maxPositionsFromEdge;
	}
	void _resetDefaults() {
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("Kmer-Match Options");
		opts.add_options()

				("match-max-positions-from-edge", po::value<int>()->default_value(maxPositionsFromEdge), "if >0 then match only reads with max-postitions-from-edge bases of either end")

				;
		desc.add(opts);

	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		setOpt("match-max-positions-from-edge", maxPositionsFromEdge);
		return ret;
	}
protected:
	int maxPositionsFromEdge;
};
typedef OptionsBaseTemplate< _KmerMatchOptions > KmerMatchOptions;


class KmerMatch : public MatcherInterface {
public:
	typedef MatcherInterface::MatchResults MatchResults;

	typedef DistributedKmerSpectrum< KmerMap< TrackingDataWithAllReads >, KmerMap< TrackingDataWithAllReads >, KmerMap< TrackingDataSingletonWithReadPosition > > KS;

// TODO make work with GSH
//	typedef KmerMapGoogleSparse< TrackingDataWithAllReads > MapType;
//	typedef DistributedKmerSpectrum< MapType, MapType, KmerMapGoogleSparse< TrackingDataSingletonWithReadPosition > > KS;

	KmerMatch(mpi::communicator &world, const ReadSet &target)
	: MatcherInterface(world, target), _spectrum(world, KS::estimateRawKmers(world, target)) {
		assert(target.isGlobal());
		_spectrum._buildKmerSpectrumMPI(target, false);
		_spectrum.purgeMinDepth(KmerSpectrumOptions::getOptions().getMinDepth());
		_spectrum.optimize(true);
	}
	virtual ~KmerMatch() {}
	MatchResults matchLocalImpl(std::string queryFile) {
		ReadSetStream query(queryFile);;
		return _matchLocal(query);
	}

	void _addResults(MatchHitSet &mhs, TrackingData::ReadPositionWeightVector &rpwv) {
		for(TrackingData::ReadPositionWeightVector::iterator it = rpwv.begin(); it != rpwv.end(); it++) {
			ReadSet::ReadSetSizeType globalReadIdx = it->readId;
			LOG_DEBUG(5, "KmerMatch::matchLocal: globalTarget " << globalReadIdx << "@" << it->position);
			mhs.insert( globalReadIdx );
		}
	}
	// works on the full copy of the query set
	MatchResults _matchLocal(ReadSetStream &query) {
		MatchResults matchResults;
		int maxPositionsFromEdge = KmerMatchOptions::getOptions().getMatchMaxPositionsFromEdge();
		int maxKmersFromEdge = maxPositionsFromEdge - KmerSizer::getSequenceLength() + 1;

		long totalMatches = 0;
		long maxMatch = 0;
		ReadSet::ReadSetSizeType contigIdx = 0;
		while(query.hasNext()) {
			Read read = query.getRead();
			matchResults.push_back(MatchHitSet());

			KmerWeights kmers(read.getTwoBitSequence(), read.getLength(), true);
			unsigned int lowerMaxKmer = maxPositionsFromEdge > 0 ? maxKmersFromEdge : kmers.size();
			unsigned int upperMinKmer = maxPositionsFromEdge > 0 ? kmers.size() - maxKmersFromEdge : 0;
			LOG_DEBUG(4, "KmerMatch::matchLocal: " << contigIdx << " " << read.toString() << " kmers: " << kmers.size());

			for(unsigned int j = 0; j < kmers.size(); j++) {
				if (j > lowerMaxKmer && j < upperMinKmer ) {
					LOG_DEBUG(5, "Skipping match to middle of " << contigIdx << " len:" << read.getLength() << " kmer:" << j);
					continue;
				}

				KS::WeakElementType element = _spectrum.getIfExistsWeak( kmers[j] );
				if (element.isValid()) {
					TrackingData::ReadPositionWeightVector rpwv = element.value().getEachInstance();
					_addResults(matchResults[contigIdx], rpwv);
				} else if (_spectrum.hasSingletons) {
					KS::SingletonElementType element = _spectrum.getIfExistsSingleton( kmers[j] );
					if (element.isValid()) {
						TrackingData::ReadPositionWeightVector rpwv = element.value().getEachInstance();
						_addResults(matchResults[contigIdx], rpwv);
					}
				}
			}
			totalMatches += matchResults[contigIdx].size();
			if (maxMatch < (long) matchResults[contigIdx].size())
				maxMatch = matchResults[contigIdx].size();
			LOG_DEBUG_OPTIONAL(2, (int) (contigIdx % (100*getWorld().size())) == (100*getWorld().rank()), "KmerMatch::_matchLocal(): processed " << contigIdx << " with " << matchResults[contigIdx].size() << " total: " << totalMatches);
			contigIdx++;
		}
		LOG_DEBUG(1, "KmerMatch::_matchLocal(): processed " << contigIdx << " total: " << totalMatches << " max: " << maxMatch<< ". " << MemoryUtils::getMemoryUsage());

		return matchResults;
	}
	const KS getKmerSpectrum() const { return _spectrum; }
protected:
	KS _spectrum;
};
#endif /* KMERMATCH_H_ */
