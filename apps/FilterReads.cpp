// $Header: /repository/PI_annex/robsandbox/KoMer/apps/FilterReads.cpp,v 1.6 2010-02-26 13:01:15 regan Exp $
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

using namespace std;

typedef TrackingDataMinimal8 DataType;
typedef KmerSpectrum<DataType, DataType> KS;
typedef ReadSelector<DataType> RS;

// TODO add outputformat of fasta
// TODO min-read-length 'FULL'
class FilterReadsOptions: Options {
public:
	static int getMaxKmerDepth() {
		return getVarMap()["max-kmer-depth"].as<unsigned int> ();
	}
	static int getPartitionByDepth() {
		return getVarMap()["partition-by-depth"].as<unsigned int> ();
	}
	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("kmer-size", 1);
		getPosDesc().add("input-file", -1);

		getDesc().add_options()(
				"max-kmer-depth",
				po::value<unsigned int>()->default_value(0),
				"maximum number of times a kmer will be represented among the selected reads (mutually exclusive with partition-by-depth)")(
				"partition-by-depth",
				po::value<unsigned int>()->default_value(1),
				"partition filtered reads by powers-of-two coverage depth (mutually exclusive with max-kmer-depth)");

		return Options::parseOpts(argc, argv);
	}
};

int main(int argc, char *argv[]) {
	if (!FilterReadsOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	cerr << MemoryUtils::getMemoryUsage() << endl;

	ReadSet reads;
	FilterKnownOddities filter;
	KmerSizer::set(Options::getKmerSize());

	Options::FileListType inputs = Options::getInputFiles();
	cerr << "Reading Input Files" << endl;
	reads.appendAllFiles(inputs);
	cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
			<< " Bases " << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	cerr << "Identifying Pairs: ";
	long numPairs = reads.identifyPairs();
	cerr << numPairs << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	cerr << "Applying filter to Input Files" << endl;
	unsigned long filtered = filter.applyFilter(reads);
	cerr << "filter affected " << filtered << " Reads " << endl;
	cerr << MemoryUtils::getMemoryUsage() << endl;

	long numBuckets = KS::estimateWeakKmerBucketSize(reads, 64);
	cerr << "targetting " << numBuckets << " buckets for reads " << endl;

	KS spectrum(numBuckets);
	cerr << MemoryUtils::getMemoryUsage() << endl;

	TrackingData::minimumWeight = 0.25;

	spectrum.buildKmerSpectrum(reads);
	cerr << MemoryUtils::getMemoryUsage() << endl;

	cerr << "Picking reads: " << endl;
	RS selector(reads, spectrum.weak, Options::getMinDepth());

	long oldPicked = 0;
	long picked = 0;

	int maximumKmerDepth = FilterReadsOptions::getMaxKmerDepth();

	if (maximumKmerDepth > 0) {
		for (int depth = 1; depth < maximumKmerDepth; depth++) {
			cerr << "Picking depth " << depth << " layer of reads" << endl;
			if (reads.hasPairs())
				picked += selector.pickBestCoveringSubsetPairs(depth,
						Options::getMinDepth(), Options::getMinReadLength());
			else
				picked += selector.pickBestCoveringSubsetReads(depth,
						Options::getMinDepth(), Options::getMinReadLength());
		}
		std::string outputFile = Options::getOutputFile();
		if (picked > 0 && !outputFile.empty()) {
			cerr << "Writing reads to output file: " << outputFile << endl;
			ofstream out;
			out.open(outputFile.c_str());
			selector.writePicks(out, oldPicked);
			out.close();
		}
	} else {

		unsigned int maxDepth = FilterReadsOptions::getPartitionByDepth();

		for (unsigned int depth = maxDepth; depth >= 1; depth /= 2) {
			if (reads.hasPairs())
				picked = selector.pickAllPassingPairs(depth,
						Options::getMinReadLength());
			else
				picked = selector.pickAllPassingReads(depth,
						Options::getMinReadLength());
			cerr << "At or above coverage: " << depth << " Picked " << picked
					<< " / " << reads.getSize() << " reads" << endl;

			std::string outputFile = Options::getOutputFile();
			if (picked > 0 && !outputFile.empty()) {
				outputFile += "-" + boost::lexical_cast<std::string>(depth);
				cerr << "Writing reads to output file: " << outputFile << endl;
				ofstream out;
				out.open(outputFile.c_str());
				selector.writePicks(out, oldPicked);
				out.close();
			}
			oldPicked += picked;
		}
	}
}
