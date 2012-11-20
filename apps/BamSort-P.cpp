/*
 * BamSort-P.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: regan
 */

#undef _USE_OPENMP

#include <string>

#include "mpi.h"
#include "SamUtils.h"
#include "Log.h"
#include "Options.h"
#include "MPIUtils.h"

class _BamSortOptions : public OptionsBaseInterface {
public:
	_BamSortOptions() : unmappedReads(), unmappedReadPairs() {
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
	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		GeneralOptions::getOptions().getDebug() = 0;
		GeneralOptions::getOptions().getMaxThreads() = 1;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		p.add("output-bam", 1);
		p.add("input-bams", -1);

		po::options_description opts("BamSort-P <options> output.bam input.bam [...]\n\nOptions");
		opts.add_options()
				("output-bam", po::value<std::string>())
				("input-bams", po::value<FileListType>())
				("unmapped-read-pairs", po::value<std::string>()->default_value(unmappedReadPairs), "gzipped file to place unmapped read Pairs Fastqs (can be same as --unmapped-reads)")
				("unmapped-reads", po::value<std::string>()->default_value(unmappedReads), "gzipped file to place unmapped reads Fastqs (can be same as --unmapped-read-pairs)");
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
};
typedef OptionsBaseTemplate< _BamSortOptions > BamSortOptions;

int main(int argc, char **argv)
{
	ScopedMPIComm< BamSortOptions > world(argc, argv);

	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	BamVector reads;
	std::string outputBam = BamSortOptions::getOptions().getOutputBam();
	OptionsBaseInterface::FileListType inputBams = BamSortOptions::getOptions().getInputBams();

	BamHeaderPtr header = BamStreamUtils::readBamFile(world, inputBams, reads);

	bool needsCollapse = false;
	if (!BamSortOptions::getOptions().getUnmappedReadPairs().empty()) {
		BamVector unmappedReads;
		SamUtils::splitUnmapped(reads, unmappedReads, true);
		if (!unmappedReads.empty())
			needsCollapse = true;
		SamUtils::writeFastqGz(world, unmappedReads, BamSortOptions::getOptions().getUnmappedReadPairs(), true);
		assert(unmappedReads.empty());
	}
	if (!BamSortOptions::getOptions().getUnmappedReads().empty()) {
		BamVector unmappedReads;
		SamUtils::splitUnmapped(reads, unmappedReads, false);
		if (!unmappedReads.empty())
			needsCollapse = true;
		needsCollapse = true;
		SamUtils::writeFastqGz(world, unmappedReads, BamSortOptions::getOptions().getUnmappedReads(), true);
		assert(unmappedReads.empty());

	}

	if (needsCollapse)
		SamUtils::collapseVector(reads);

	unlink(outputBam.c_str());
	LOG_VERBOSE(1, "Sorting myreads: " << reads.size());

	{
		SamUtils::MPISortBam sortem(world, reads, outputBam, header.get());
	}

	header.reset();

	LOG_VERBOSE(1, "Finished");

	return 0;
}
