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

int main(int argc, char **argv)
{
	if (argc <= 2) {
		std::cerr << "Usage: mpirun BamSort-P out.bam [in.[sb]am ...]";
		exit(1);
	}
	if (MPI_SUCCESS != MPI_Init(&argc, &argv))
		LOG_THROW("MPI_Init() failed: ");

	mpi::communicator world;
	Logger::setWorld(&world);

	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	GeneralOptions::getOptions().getDebug() = 1;
	BamVector reads;
	std::string outputBam(argv[1]);
	samfile_t *fh = NULL;
	for(int i = 2; i < argc; i++) {
		if ((i-2) % size == rank) {
			std::string inputFile(argv[i]);
			LOG_VERBOSE_OPTIONAL(1, true, "Reading " << inputFile);
			if (fh != NULL)
				samclose(fh);
			fh = BamStreamUtils::openSamOrBam(inputFile);
			long s = reads.size();
			BamStreamUtils::readBamFile(fh, reads);
			LOG_VERBOSE_OPTIONAL(1, true, "Loaded " << (reads.size() - s) << " new bam records");
		}
	}
	LOG_VERBOSE(1, "Sorting myreads: " << reads.size());

	{
		SamUtils::MPISortBam sortem(MPI_COMM_WORLD, reads);
	}

	LOG_VERBOSE(1, "Writing our sorted reads: " << reads.size());
	SamUtils::writeBamVector(MPI_COMM_WORLD, outputBam, reads, (fh != NULL) ? fh->header : NULL, true);

	if (fh!=NULL)
		samclose(fh);

	if (MPI_SUCCESS != MPI_Finalize())
		LOG_THROW("MPI_Finalize() failed: ");
	return 0;
}
