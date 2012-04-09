/*
 * KmerMatch.h
 *
 *  Created on: Nov 1, 2011
 *      Author: regan
 */

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

class _KmerMatchOptions  : public _KmerBaseOptions {
public:
	_KmerMatchOptions() : _KmerBaseOptions(), maxPositionsFromEdge(0) {}
	~_KmerMatchOptions() {}
	int &getMatchMaxPositionsFromEdge() {
		return maxPositionsFromEdge;
	}
	void _resetDefaults() {
		_KmerBaseOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		_KmerBaseOptions::_setOptions(desc,p);
		po::options_description opts("Kmer-Match Options");
		opts.add_options()

				("match-max-positions-from-edge", po::value<int>()->default_value(maxPositionsFromEdge), "if >0 then match only reads with max-postitions-from-edge bases of either end")

				;
		desc.add(opts);

	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = _KmerBaseOptions::_parseOptions(vm);
		setOpt<int>("match-max-positions-from-edge", maxPositionsFromEdge);
		return ret;
	}
protected:
	int maxPositionsFromEdge;
};
typedef OptionsBaseTemplate< _KmerMatchOptions > KmerMatchOptions;


class KmerMatch : public MatcherInterface {
public:
	typedef MatcherInterface::MatchResults MatchResults;
	typedef DistributedKmerSpectrum< TrackingDataWithAllReads, TrackingDataWithAllReads, TrackingDataSingletonWithReadPosition > KS;

	KmerMatch(mpi::communicator &world, const ReadSet &target)
	: MatcherInterface(world, target), _spectrum(world, KS::estimateWeakKmerBucketSize(target)) {
		assert(target.isGlobal());
		_spectrum._buildKmerSpectrumMPI(target, true);
		_spectrum.optimize();
	}
	MatchResults matchLocalImpl(std::string queryFile) {
		ReadSetStream query(queryFile);;
		return _matchLocal(query);
	}

	// works on the full copy of the query set
	MatchResults _matchLocal(ReadSetStream &query) {
		MatchResults matchResults;
		int maxPositionsFromEdge = KmerMatchOptions::getOptions().getMatchMaxPositionsFromEdge();
		int maxKmersFromEdge = maxPositionsFromEdge - KmerSizer::getSequenceLength() + 1;

		long totalMatches = 0;
		ReadSet::ReadSetSizeType contigIdx = 0;
		while(query.readNext()) {
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

				KS::SolidElementType element = _spectrum.getIfExistsSolid( kmers[j] );
				if (element.isValid()) {
					TrackingData::ReadPositionWeightVector rpwv = element.value().getEachInstance();
					for(TrackingData::ReadPositionWeightVector::iterator it = rpwv.begin(); it != rpwv.end(); it++) {
						ReadSet::ReadSetSizeType globalReadIdx = it->readId;
						LOG_DEBUG(5, "KmerMatch::matchLocal: localRead " << contigIdx << "@" << j << " globalTarget " << globalReadIdx << "@" << it->position << " " << kmers[j].toFasta());
						matchResults[contigIdx].insert( globalReadIdx );
					}
				}
			}
			totalMatches += matchResults[contigIdx].size();
			LOG_DEBUG_OPTIONAL(1, contigIdx % (100*getWorld().size()) == (100*getWorld().rank()), "KmerMatch::_matchLocal(): processed " << contigIdx << " with " << matchResults[contigIdx-1].size() << " total: " << totalMatches);
			contigIdx++;
		}

		return matchResults;
	}
	const KS getKmerSpectrum() const { return _spectrum; }
protected:
	KS _spectrum;
};
#endif /* KMERMATCH_H_ */
