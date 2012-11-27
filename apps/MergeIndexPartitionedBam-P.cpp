/*
 * MergeIndexPartitionedBam-P.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: regan
 */

#undef _USE_OPENMP

#include <string>
#include "mpi.h"
#include "SamUtils.h"
#include "Log.h"
#include "Options.h"
#include "MPIUtils.h"

class _MergeIndexPartitionedBamOptions : public OptionsBaseInterface {
public:
	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		GeneralOptions::getOptions().getDebug() = 0;
		GeneralOptions::getOptions().getMaxThreads() = 1;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		p.add("output-bam", 1);
		p.add("num-partitions", 1);
		p.add("input-bams", -1);

		po::options_description opts("MergeIndexPartitionedBam-P <options> output.bam num_partitions input.bam [...]\n\nOptions");
		opts.add_options()
				("output-bam", po::value<std::string>())
				("input-bams", po::value<FileListType>())
				("num-partitions", po::value<int>()->default_value(0), "The number of index partitions");
		desc.add(opts);

		GeneralOptions::_setOptions(desc, p);
		// Other *::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		ret |= GeneralOptions::_parseOptions(vm);

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
	std::string outputBam;
	FileListType inputBams;
	int numPartitions;
};
typedef OptionsBaseTemplate< _MergeIndexPartitionedBamOptions > MergeIndexPartitionedBamOptions;
int main(int argc, char **argv)
{


	ScopedMPIComm< MergeIndexPartitionedBamOptions > world(argc, argv);

	int rank = world.rank(), size = world.size();

	std::string ourOutputBam = MergeIndexPartitionedBamOptions::getOptions().outputBam;
	int partitions = MergeIndexPartitionedBamOptions::getOptions().numPartitions;
	std::string myInputFile = MergeIndexPartitionedBamOptions::getOptions().inputBams[rank];

	if (size != (int) MergeIndexPartitionedBamOptions::getOptions().inputBams.size()) {
		if (rank == 0)
			std::cerr << " Your run size must be equal to the number of inputs\n " << MergeIndexPartitionedBamOptions::getDesc() << std::endl;
		exit(1);
	}

	if (partitions == 0 || size % partitions != 0) {
		if (rank == 0)
			std::cerr << "The partitions " << partitions << " is not a factor of the size " << size << MergeIndexPartitionedBamOptions::getDesc() << std::endl;
		exit(1);
	}
	BamVector myReads;
	{
		int color = rank / partitions;
		mpi::communicator partitionWorld = world.split(color);
		LOG_VERBOSE(1, "Input: " << myInputFile << " split color: " << color << " rank: " << partitionWorld.rank() << " of " << partitionWorld.size());
		SamUtils::MPIMergeSam mergeSam(partitionWorld, myInputFile, myReads);
		BamStreamUtils::distributeReadsFinal(world, myReads);
		unlink(ourOutputBam.c_str());
		mergeSam.outputMergedHeader(ourOutputBam);
		SamUtils::MPISortBam sortBam(world, myReads, ourOutputBam, NULL);
	}

	LOG_VERBOSE(1, "Finished");

	return 0;
}
