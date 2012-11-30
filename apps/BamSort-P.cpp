/*
 * BamSort-P.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: regan
 */

#include <string>

#include "mpi.h"
#include "SamUtils.h"
#include "Log.h"
#include "Options.h"
#include "MPIUtils.h"

class _BamSortOptions : public OptionsBaseInterface {
public:
	_BamSortOptions() : unmappedReads(), unmappedReadPairs(), numPartitions(0) {
	}
	virtual ~_BamSortOptions() {}
	std::string &getUnmappedReads() {
		return unmappedReads;
	}
	std::string &getUnmappedReadPairs() {
		return unmappedReadPairs;
	}
	std::string &getOutputBam() {
		return outputBam;
	}
	FileListType &getInputBams() {
		return inputBams;
	}
	int &getNumPartitions() {
		return numPartitions;
	}
	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		GeneralOptions::getOptions().getDebug() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		p.add("output-bam", 1);
		p.add("input-bams", -1);

		po::options_description opts("BamSort-P <options> output.bam input.bam [...]\n\nOptions");
		opts.add_options()
				("output-bam", po::value<std::string>())
				("input-bams", po::value<FileListType>())
				("unmapped-read-pairs", po::value<std::string>()->default_value(unmappedReadPairs), "gzipped file to place unmapped read Pairs Fastqs (can be same as --unmapped-reads)")
				("unmapped-reads", po::value<std::string>()->default_value(unmappedReads), "gzipped file to place unmapped reads Fastqs (can be same as --unmapped-read-pairs)")
				("num-partitions", po::value<int>()->default_value(numPartitions), "The number of alignment-index partitions to merge. Input bams expected to come ordered grouped by in batches of num-partitions where each group has the exact same read counts in the exact same order");
		desc.add(opts);

		GeneralOptions::_setOptions(desc, p);
		// Other *::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		ret |= GeneralOptions::_parseOptions(vm);

		setOpt("unmapped-read-pairs", unmappedReadPairs);
		setOpt("unmapped-reads", unmappedReads);
		setOpt("output-bam", outputBam);
		setOpt2("input-bams", inputBams);
		setOpt("num-partitions", numPartitions);

		if (outputBam.empty())
			setOptionsErrorMsg("You must specify at least the outputBam");

		if (inputBams.empty())
			setOptionsErrorMsg("You must specify at least one input bam");

		// Other ret &= *::_parseOptions(vm);
		return ret;
	}
private:
	std::string outputBam;
	FileListType inputBams;
	std::string unmappedReads;
	std::string unmappedReadPairs;
	int numPartitions;
};
typedef OptionsBaseTemplate< _BamSortOptions > BamSortOptions;

int main(int argc, char **argv)
{
	ScopedMPIComm< BamSortOptions > world(argc, argv);

	int rank = world.rank(), size = world.size();

	BamVector reads;
	BamHeaderPtr header;

	std::string outputBam = BamSortOptions::getOptions().getOutputBam();
	OptionsBaseInterface::FileListType inputBams = BamSortOptions::getOptions().getInputBams();
	int partitions = BamSortOptions::getOptions().getNumPartitions();

	if (partitions > 1) {
		if (size != inputBams.size()) {
			if (rank == 0)
				std::cerr << "The number of files is a mismatch to the job size.  When merging partitions, this is necessary." << BamSortOptions::getDesc() << std::endl;
			exit(1);
		}
		if (size % partitions != 0) {
			if (rank == 0)
				std::cerr << "The partitions " << partitions << " is not a factor of the size " << size << BamSortOptions::getDesc() << std::endl;
			exit(1);
		}
	}

	unlink(outputBam.c_str());

	if (partitions > 1) {
		std::string myInputFile = inputBams[rank];
		int color = rank / partitions;
		mpi::communicator partitionWorld = world.split(color);
		LOG_VERBOSE(1, "Input: " << myInputFile << " split color: " << color << " rank: " << partitionWorld.rank() << " of " << partitionWorld.size());
		SamUtils::MPIMergeSam mergeSam(partitionWorld, myInputFile, reads);

		if (color == 0)
			mergeSam.outputMergedHeader(outputBam);

		world.barrier();

	} else {
		header = BamStreamUtils::readBamFile(world, inputBams, reads);
	}

	bool needsCollapse = false;
	std::string unmappedReadPairFile = BamSortOptions::getOptions().getUnmappedReadPairs();
	if (!unmappedReadPairFile.empty()) {
		BamVector unmappedReads;
		SamUtils::splitUnmapped(reads, unmappedReads, true);
		if (!unmappedReads.empty())
			needsCollapse = true;
		if (unmappedReadPairFile.compare("/dev/null") == 0) {
			BamManager::destroyOrRecycleBamVector(unmappedReads);
		} else {
			SamUtils::writeFastqGz(world, unmappedReads, unmappedReadPairFile, true);
		}
		assert(unmappedReads.empty());
	}

	std::string unmappedReadsFile  = BamSortOptions::getOptions().getUnmappedReads();
	if (!unmappedReadsFile.empty()) {
		BamVector unmappedReads;
		SamUtils::splitUnmapped(reads, unmappedReads, false);
		if (!unmappedReads.empty())
			needsCollapse = true;
		if (unmappedReadsFile.compare("/dev/null") == 0) {
			BamManager::destroyOrRecycleBamVector(unmappedReads);
		} else {
			SamUtils::writeFastqGz(world, unmappedReads, unmappedReadsFile, true);
		}
		assert(unmappedReads.empty());
	}

	if (needsCollapse)
		SamUtils::collapseVector(reads);

	{
		BamStreamUtils::distributeReadsFinal(world, reads);

		LOG_VERBOSE(1, "Sorting myreads: " << reads.size());

		SamUtils::MPISortBam sortem(world, reads, outputBam, header.get());
	}

	header.reset();

	LOG_VERBOSE(1, "Finished");

	return 0;
}
