/*
 * Meraculous.h
 *
 *  Created on: Sep 16, 2011
 *      Author: regan
 */

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

typedef DistributedKmerSpectrum<ExtensionTrackingData, ExtensionTrackingData, ExtensionTrackingDataSingleton> _MeraculousDistributedKmerSpectrum;
class MeraculousDistributedKmerSpectrum : public _MeraculousDistributedKmerSpectrum {
public:
	typedef _MeraculousDistributedKmerSpectrum::WeakMapType WeakMapType;
	typedef WeakMapType::Iterator KmerIter;
	MeraculousDistributedKmerSpectrum(mpi::communicator &_world, unsigned long buckets = 0, bool separateSingletons = true) :
		_MeraculousDistributedKmerSpectrum(_world, buckets, separateSingletons) {

	}
	~MeraculousDistributedKmerSpectrum() {}

	void dumpCounts(std::string filename) {
		DistributedOfstreamMap dom(getWorld(), filename, "");
		ostream &os = dom.getOfstream("");
		for(KmerIter it = weak.begin(); it != weak.end(); it++) {
			long count = it.value().getDirectionBias();
			long revCount = it.value().getCount() - count;
			os << it.key().toFasta() << "\t" << count+revCount << std::endl;
			TEMP_KMER(rev);
			it.key().buildReverseComplement(rev);
			os << rev.toFasta() << "\t" << revCount+count << std::endl;
		}
	}
	void dumpGraphs(std::string filename) {
		DistributedOfstreamMap dom(getWorld(), filename, "");
		ostream &os = dom.getOfstream("");
		for(KmerIter it = weak.begin(); it != weak.end(); it++) {
			ExtensionTracking extensions = it.value().getExtensionTracking();
			os << it.key().toFasta() << "\t" << extensions.toTextValues() << std::endl;
			TEMP_KMER(rev);
			it.key().buildReverseComplement(rev);
			os << rev.toFasta() << "\t" << extensions.getReverseComplement().toTextValues() << std::endl;
		}
	}

};

#endif /* MERACULOUS_H_ */
