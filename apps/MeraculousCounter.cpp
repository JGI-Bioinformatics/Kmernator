/*
 * MeraculousCounter-P.cpp
 *
 *  Created on: Sep 20, 2011
 *      Author: regan
 */

#include "config.h"
#include "Options.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "KmerTrackingData.h"
#include "Utils.h"
#include "Log.h"
#include "DistributedFunctions.h"
#include "Meraculous.h"

class _MeraculousCounterOptions : public OptionsBaseInterface {
public:
	// use to set/overrided any defaults on options that are stored persistently
	void _resetDefaults() {
		MeraculousOptions::_resetDefaults();
		MPIOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
		KmerOptions::_resetDefaults();

		GeneralOptions::getOptions().getMmapInput() = 0;
		GeneralOptions::getOptions().getVerbose() = 2;
		GeneralOptions::getOptions().getMinQuality() = 2;

		KmerOptions::getOptions().getMinKmerQuality() = 0;
		KmerOptions::getOptions().getSaveKmerMmap() = 0;
	}
	// use to set the description of all options
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {

		MeraculousOptions::_setOptions(desc, p);
		MPIOptions::_setOptions(desc,p);
		KmerOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc, p);

	}
	// use to post-process options, returning true if everything is okay
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= MeraculousOptions::_parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);
		ret &= KmerOptions::_parseOptions(vm);
		if (KmerOptions::getOptions().getKmerSize() == 0) {
			setOptionsErrorMsg("The Kmer size can not be 0");
		}
		if (GeneralOptions::getOptions().getInputFiles().empty()) {
			setOptionsErrorMsg("You must specify at least one input file");
		}
		return ret;
	}
};
typedef OptionsBaseTemplate< _MeraculousCounterOptions > MeraculousCounterOptions;
typedef MeraculousDistributedKmerSpectrum KS;

int main(int argc, char *argv[]) {

	mpi::communicator world = initializeWorldAndOptions< MeraculousCounterOptions > (argc, argv);

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	ReadSet reads;

	try {
		OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files");

		// TODO save memory! read file and build spectrum in 100MB chunks
		reads.appendAllFiles(inputs, world.rank(), world.size());

		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		setGlobalReadSetOffsets(world, reads);
		long numBuckets = 0;
		numBuckets = KS::estimateWeakKmerBucketSize(reads);

		numBuckets = all_reduce(world, numBuckets, mpi::maximum<int>());
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "targeting " << numBuckets << " buckets for reads");

		KS spectrum(world, numBuckets);

		spectrum.buildKmerSpectrum(reads);
		if (Log::isVerbose(1)) {
			std::string hist = spectrum.getHistogram(false);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Kmer Histogram\n" << hist);
		}
		std::string outputFilenameBase = outputFilename + ".mercount.m" + boost::lexical_cast<std::string>(KmerSizer::getSequenceLength());
		spectrum.dumpCounts(outputFilenameBase);
		outputFilenameBase = outputFilename + ".mergraph.m" + boost::lexical_cast<std::string>(KmerSizer::getSequenceLength()) + ".D" + boost::lexical_cast<std::string>(KmerOptions::getOptions().getMinDepth());
		spectrum.dumpGraphs(outputFilenameBase);
	} catch (...) {
		LOG_ERROR(1, "caught an error!" << StackTrace::getStackTrace());
	}
	world.barrier();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");

	MPI_Finalize();

	return 0;

}

