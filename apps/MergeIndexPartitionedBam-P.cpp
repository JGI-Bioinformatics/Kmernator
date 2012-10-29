/*
 * MergeIndexPartitionedBam-P.cpp
 *
 *  Created on: Oct 17, 2012
 *      Author: regan
 */

#include <string>
#include "mpi.h"
#include "SamUtils.h"
#include "Log.h"
#include "Options.h"

int main(int argc, char **argv)
{

	if (MPI_SUCCESS != MPI_Init(&argc, &argv))
		LOG_THROW("MPI_Init() failed: ");

	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (argc != 2 + size) {
		if (rank == 0)
			std::cerr << "Usage: mpirun MergeIndexPartitionedBam-P out.bam [in1.[sb]am in2.[sb]am ...]\n\tYour run size must be equal to the number of inputs\n\n";
		MPI_Finalize();
		exit(1);
	}
	mpi::communicator world;
	Logger::setWorld(&world);
	GeneralOptions::getOptions().getDebug() = 2;

	std::string ourOutputBam(argv[1]);
	std::string myInputFile(argv[2+rank]);
	BamVector myReads;
	{
		SamUtils::MPIMergeSam mergeSam(MPI_COMM_WORLD, myInputFile, myReads);
		unlink(ourOutputBam.c_str());
		mergeSam.outputMergedHeader(ourOutputBam);
		SamUtils::MPISortBam sortBam(MPI_COMM_WORLD, myReads, ourOutputBam, NULL);
	}

	LOG_VERBOSE(1, "Finished");

	Logger::setWorld(NULL);
	if (MPI_SUCCESS != MPI_Finalize())
		LOG_THROW("MPI_Finalize() failed: ");
	return 0;
}
