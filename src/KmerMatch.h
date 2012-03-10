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

class _KmerMatchOptions  : public _KmerBaseOptions {
public:
	_KmerMatchOptions() : _KmerBaseOptions(), maxPositionsFromEdge(0),includeMate(1) {}
	~_KmerMatchOptions() {}
	int &getMaxPositionsFromEdge() {
		return maxPositionsFromEdge;
	}
	bool getIncludeMate() {
		return includeMate == 1;
	}
	void _resetDefaults() {
		_KmerBaseOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		_KmerBaseOptions::_setOptions(desc,p);
		po::options_description opts("Kmer-Match Options");
		opts.add_options()

				("max-positions-from-edge", po::value<int>()->default_value(maxPositionsFromEdge), "if >0 then match only reads with max-postitions-from-edge bases of either end")

				("include-mate", po::value<int>()->default_value(includeMate), "1 - include mates, 0 - do not")
				;
		desc.add(opts);

	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = _KmerBaseOptions::_parseOptions(vm);
		setOpt<int>("max-positions-from-edge", maxPositionsFromEdge);
		setOpt<int>("include-mate", includeMate);

		return ret;
	}
protected:
	int maxPositionsFromEdge, includeMate;
};
typedef OptionsBaseTemplate< _KmerMatchOptions > KmerMatchOptions;


class KmerMatch : public MatcherInterface {
public:
	typedef MatcherInterface::MatchResults MatchResults;
	typedef DistributedKmerSpectrum< TrackingDataWithAllReads, TrackingDataWithAllReads, TrackingDataSingletonWithReadPosition > KS;

	KmerMatch(mpi::communicator &world, const ReadSet &target, bool returnPairedMatches = KmerMatchOptions::getOptions().getIncludeMate())
	: MatcherInterface(world, target, returnPairedMatches), _spectrum(world, KS::estimateWeakKmerBucketSize(target)) {
		assert(target.isGlobal());
		_spectrum._buildKmerSpectrumMPI(target, true);
		_spectrum.optimize();
	}
	MatchResults matchLocal(std::string queryFile) {
		ReadSet query;
		query.appendAnyFile(queryFile);
		return this->_matchLocal(query);
	}
	virtual MatchResults matchLocal(const ReadSet &query) {
		return MatcherInterface::matchLocal(query);
	}

	// works on the full copy of the query set
	MatchResults _matchLocal(const ReadSet &query) {
		assert(!query.isGlobal());
		MatchResults matchResults;
		matchResults.resize(query.getSize());
		bool includeMates = getTarget().hasPairs() && isReturnPairedMatches();
		int maxPositionsFromEdge = KmerMatchOptions::getOptions().getMaxPositionsFromEdge();
		int maxKmersFromEdge = maxPositionsFromEdge - KmerSizer::getSequenceLength() + 1;

		for(unsigned int i = 0; i < query.getSize(); i++) {
			const Read &read = query.getRead(i);
			KmerWeights kmers(read.getTwoBitSequence(), read.getLength(), true);
			unsigned int lowerMaxKmer = maxPositionsFromEdge > 0 ? maxKmersFromEdge : kmers.size();
			unsigned int upperMinKmer = maxPositionsFromEdge > 0 ? kmers.size() - maxKmersFromEdge : 0;

			for(unsigned int j = 0; j < kmers.size(); j++) {
				if (j > lowerMaxKmer && j < upperMinKmer ) {
					LOG_DEBUG(5, "Skipping match to middle of " << i << " len:" << read.getLength() << " kmer:" << j);
					continue;
				}
				KS::SolidElementType element = _spectrum.getIfExistsSolid( kmers[j] );
				if (element.isValid()) {
					TrackingData::ReadPositionWeightVector rpwv = element.value().getEachInstance();
					for(TrackingData::ReadPositionWeightVector::iterator it = rpwv.begin(); it != rpwv.end(); it++) {
						ReadSet::ReadSetSizeType globalReadIdx = it->readId;
						matchResults[i].insert( globalReadIdx );
						LOG_DEBUG(5, "KmerMatch::matchLocal: localRead " << i << "@" << j << " globalTarget " << globalReadIdx << "@" << it->position << " " << kmers[j].toFasta());
						if (includeMates) {
							matchResults[i].insert( getTarget().getGlobalPairIdx( globalReadIdx ) );
						}
					}
				}
			}
		}

		return matchResults;
	}
	const KS getKmerSpectrum() const { return _spectrum; }
protected:
	KS _spectrum;
};
#endif /* KMERMATCH_H_ */
