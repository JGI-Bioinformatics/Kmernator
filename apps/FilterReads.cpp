// $Header: /repository/PI_annex/robsandbox/KoMer/apps/FilterReads.cpp,v 1.19 2010-05-06 16:43:59 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "Options.h"
#include "ReadSet.h"
#include "FilterKnownOddities.h"
#include "KmerSpectrum.h"
#include "ReadSelector.h"
#include "Utils.h"

#include <boost/lexical_cast.hpp>

using namespace std;

typedef TrackingDataMinimal4f DataType;
typedef KmerSpectrum<DataType, DataType> KS;
typedef ReadSelector<DataType> RS;

// TODO add outputformat of fasta
class FilterReadsOptions : Options {
public:
	static int getMaxKmerDepth() {
		return getVarMap()["max-kmer-depth"].as<int> ();
	}
	static int getPartitionByDepth() {
		return getVarMap()["partition-by-depth"].as<int> ();
	}
	static bool getBothPairs() {
		return getVarMap()["min-passing-in-pair"].as<int>() == 2;
	}
	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("kmer-size", 1);
		getPosDesc().add("input-file", -1);

		getDesc().add_options()

		("max-kmer-depth", po::value<int>()->default_value(-1),
				"maximum number of times a kmer will be represented among the selected reads (mutually exclusive with partition-by-depth)")

		("partition-by-depth", po::value<int>()->default_value(-1),
				"partition filtered reads by powers-of-two coverage depth (mutually exclusive with max-kmer-depth)")

		("min-passing-in-pair", po::value<int>()->default_value(1),
				"1 or 2 reads in a pair must pass filters");

		bool ret = Options::parseOpts(argc, argv);

		if (ret) {
			// verify mutually exclusive options are not set
			if ( (getMaxKmerDepth() > 0 && getPartitionByDepth() >  0) )
			{
				throw std::invalid_argument("You can not specify both max-kmer-depth and partition-by-depth");
			}
		}
		return ret;
	}
};

int main(int argc, char *argv[]) {
	if (!FilterReadsOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	MemoryUtils::getMemoryUsage();

	ReadSet reads;
	KmerSizer::set(Options::getKmerSize());

	Options::FileListType inputs = Options::getInputFiles();
	cerr << "Reading Input Files" << endl;
	reads.appendAllFiles(inputs);
	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
			<< " Bases " << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	cerr << "Identifying Pairs: ";
	long numPairs = reads.identifyPairs();
	cerr << "Pairs + single = " << numPairs << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	if (Options::getSkipArtifactFilter() == 0) {

	  cerr << "Preparing artifact filter: ";
      FilterKnownOddities filter;
      cerr << MemoryUtils::getMemoryUsage() << endl;

	  cerr << "Applying sequence artifact filter to Input Files" << endl;
	  unsigned long filtered = filter.applyFilter(reads);
	  cerr << "filter affected (trimmed/removed) " << filtered << " Reads " << endl;;
	  cerr << MemoryUtils::getMemoryUsage() << endl;

	  cerr << "Applying DuplicateFragmentPair Filter to Input Files" << endl;
	  unsigned long duplicateFragments = filter.filterDuplicateFragmentPairs(reads);
	  cerr << "filter affected  (removed) " << duplicateFragments << endl;
	  cerr << MemoryUtils::getMemoryUsage() << endl;
	}

	KS spectrum(0);

	KoMer::MmapFileVector spectrumMmaps;

	if (Options::getKmerSize() > 0) {

	  long numBuckets = KS::estimateWeakKmerBucketSize(reads, 64);
	  cerr << "targeting " << numBuckets << " buckets for reads " << endl;

	  spectrum = KS(numBuckets);
	  cerr << MemoryUtils::getMemoryUsage() << endl;

	  TrackingData::minimumWeight = Options::getMinKmerQuality();

	  spectrumMmaps = spectrum.buildKmerSpectrumInParts(reads, Options::getBuildPartitions());
	  cerr << MemoryUtils::getMemoryUsage() << endl;

	  if (Options::getMinDepth() > 1) {
        cerr << "Clearing singletons from memory" << endl;
        spectrum.singleton.clear();
	    cerr << MemoryUtils::getMemoryUsage() << endl;
	  }
	}

	cerr << "Trimming reads: ";
	RS selector(reads, spectrum.weak, Options::getMinDepth());
	cerr << MemoryUtils::getMemoryUsage() << endl;
	cerr << "Picking reads: " << endl;

	long oldPicked = 0;
	long picked = 0;

	int maximumKmerDepth = FilterReadsOptions::getMaxKmerDepth();

	std::string outputFilename = Options::getOutputFile();
	OfstreamMap ofmap(outputFilename, ".fastq");

	if (maximumKmerDepth > 0) {
		for (int depth = 1; depth < maximumKmerDepth; depth++) {
			cerr << "Picking depth " << depth << " layer of reads" << endl;
			if (reads.hasPairs())
				picked += selector.pickBestCoveringSubsetPairs(depth,
						Options::getMinDepth(), Options::getMinReadLength(), FilterReadsOptions::getBothPairs());
			else
				picked += selector.pickBestCoveringSubsetReads(depth,
						Options::getMinDepth(), Options::getMinReadLength());
            cerr << MemoryUtils::getMemoryUsage() << endl;
		}

		if (picked > 0 && !outputFilename.empty()) {
			cerr << "Writing " << picked << " reads to output file(s)" << endl;
			selector.writePicks(ofmap, oldPicked);
		}
        cerr << MemoryUtils::getMemoryUsage() << endl;


	} else {

		int maxDepth = FilterReadsOptions::getPartitionByDepth();
		if (maxDepth < 0) {
			maxDepth = 1;
		}

		for (unsigned int depth = maxDepth; depth >= 1; depth /= 2) {

			string ofname = outputFilename + '-';
			if (maxDepth > 1) {
				ofname += "-MinDepth" + boost::lexical_cast< string >( depth );
			}
			ofmap = OfstreamMap(ofname, ".fastq");
			float minDepth = std::max(Options::getMinDepth(), depth);
			if (Options::getKmerSize() == 0) {
				minDepth = 0;
				depth = 0;
			}
			cerr << "Selecting reads over depth: " << depth << " (" << minDepth << ") " << endl;

			if (reads.hasPairs()) {
				picked = selector.pickAllPassingPairs(minDepth,
						Options::getMinReadLength(),
						FilterReadsOptions::getBothPairs());
			} else {
				picked = selector.pickAllPassingReads(minDepth,
						Options::getMinReadLength());
			}
			cerr << "At or above coverage: " << depth << " Picked " << picked
					<< " / " << reads.getSize() << " reads" << endl;
            cerr << MemoryUtils::getMemoryUsage() << endl;

			if (picked > 0 && !outputFilename.empty()) {
				cerr << "Writing " << picked << " reads  to output files" << endl;
				selector.writePicks(ofmap, oldPicked);
			}
			oldPicked += picked;

			if (Options::getMinDepth() > depth) {
				break;
			}
		}
	}
}
