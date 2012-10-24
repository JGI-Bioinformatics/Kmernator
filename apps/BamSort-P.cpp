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
	_BamSortOptions() : unmappedReads(), unmappedReadPairs() {
	}
	virtual ~_BamSortOptions() {}
	void _resetDefaults() {
		// Other *::_resetDefaults();
	}
	std::string &getUnmappedReads() {
		return unmappedReads;
	}
	std::string &getUnmappedReadPairs() {
		return unmappedReadPairs;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("BamSort-P Options");
		opts.add_options()
				("unmapped-read-pairs", po::value<std::string>()->default_value(unmappedReadPairs), "gzipped file to place unmapped read Pairs Fastqs (can be same as --unmapped-reads)")
				("unmapped-reads", po::value<std::string>()->default_value(unmappedReads), "gzipped file to place unmapped reads Fastqs (can be same as --unmapped-read-pairs)");
		desc.add(opts);
		// Other *::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		setOpt("unmapped-read-pairs", unmappedReadPairs);
		setOpt("unmapped-reads", unmappedReads);
		// Other ret &= *::_parseOptions(vm);
		return ret;
	}
private:
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
	if (argc <= 2) {
		if (rank == 0)
			std::cerr << "Usage: mpirun BamSort-P out.bam [in.[sb]am ...]\n\n";
		MPI_Finalize();
		exit(1);
	}

	GeneralOptions::getOptions().getDebug() = 0;
	BamVector reads;
	std::string outputBam(argv[1]);
	samfile_t *fh = NULL;
	for(int i = 2; i < argc; i++) {
		if ((i-2) % size == rank) {
			std::string inputFile(argv[i]);
			LOG_VERBOSE_OPTIONAL(1, true, "Reading " << inputFile);
			if (fh != NULL)
				BamStreamUtils::closeSamOrBam(fh);
			fh = BamStreamUtils::openSamOrBam(inputFile);
			long s = reads.size();
			BamStreamUtils::readBamFile(fh, reads);
			LOG_VERBOSE_OPTIONAL(1, true, "Loaded " << (reads.size() - s) << " new bam records");
		}
	}

	unlink(outputBam.c_str());
	LOG_VERBOSE(1, "Sorting myreads: " << reads.size());

	{
		SamUtils::MPISortBam sortem(world, reads, outputBam, (fh != NULL) ? fh->header : NULL);
	}

	if (fh != NULL)
		BamStreamUtils::closeSamOrBam(fh);

	LOG_VERBOSE(1, "Finished");

	return 0;
}
