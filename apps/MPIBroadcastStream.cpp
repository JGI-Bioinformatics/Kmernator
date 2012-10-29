/*
 * MPIBroadcastStream.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: regan
 */

#include <string>
#include "mpi.h"

#include <boost/lexical_cast.hpp>
#include <boost/iostreams/copy.hpp>

#include "Log.h"
#include "MPIUtils.h"
#include "Options.h"
#include "BroadcastOstream.h"


int main(int argc, char **argv)
{
	GeneralOptions::getOptions().getMaxThreads() = 1;
	GeneralOptions::getOptions().getDebug() = 0;
	ScopedMPIComm<> world(argc, argv);

	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (argc != 3) {
		if (rank == 0)
			std::cerr << "Usage: mpirun MPIBroadcastStream input.file output.file[.rank]\n\tThis will create n copies of the input.file\n\n";
		MPI_Finalize();
		exit(1);
	}

	std::ifstream ifs(argv[1]);
	std::string output(argv[2]);
	output += "." + boost::lexical_cast<std::string>(rank);
	LOG_VERBOSE(1, "Copying " << argv[1] << " to " << output);
	{
		std::ofstream os(output.c_str());
		int root = 0;
		BroadcastOstream bcastos(world, root, os);

		if (rank == root)
			boost::iostreams::copy(ifs, bcastos);
	}
	LOG_VERBOSE_OPTIONAL(1, true, "Done");


	return 0;
}
