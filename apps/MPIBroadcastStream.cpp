/*
 * MPIBroadcastStream.cpp
 *
 *  Created on: Oct 22, 2012
 *      Author: regan
 */
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/

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
