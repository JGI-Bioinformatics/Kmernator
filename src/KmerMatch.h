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
};
typedef OptionsBaseTemplate< _KmerMatchOptions > KmerMatchOptions;


class KmerMatch : public MatcherInterface {
public:
	typedef MatcherInterface::MatchResults MatchResults;
	typedef DistributedKmerSpectrum< TrackingDataWithAllReads, TrackingDataWithAllReads, TrackingDataSingletonWithReadPosition > KS;

	KmerMatch(mpi::communicator &world, const ReadSet &target, bool returnPairedMatches = true)
	: MatcherInterface(world, target, returnPairedMatches), _spectrum(world) {
		_spectrum._buildKmerSpectrumMPI(target, true);
	}
	MatchResults matchLocal(std::string queryFile) {
		ReadSet query;
		query.appendAnyFile(queryFile);
		return this->matchLocal(query);
	}
	virtual MatchResults matchLocal(const ReadSet &query) {
		// TODO
		MatchResults matchResults;
		matchResults.resize(query.getSize());
		for(unsigned int i = 0; i < query.getSize(); i++) {
			const Read &read = query.getRead(i);
			KmerWeights kmers(read.getTwoBitSequence(), read.getLength(), true);

			for(unsigned int j = 0; j < kmers.size(); j++) {
				KS::SolidElementType element = _spectrum.getIfExistsSolid( kmers[j] );
				if (element.isValid()) {
					TrackingData::ReadPositionWeightVector rpwv = element.value().getEachInstance();
					for(TrackingData::ReadPositionWeightVector::iterator it = rpwv.begin(); it != rpwv.end(); it++) {
						// TODO add positional logic to include or not
						matchResults[i].insert( getTarget().getGlobalReadIdx( it->readId ) );
					}
				}
			}
		}

		return matchResults;
	}
protected:
	KS _spectrum;
};
#endif /* KMERMATCH_H_ */
