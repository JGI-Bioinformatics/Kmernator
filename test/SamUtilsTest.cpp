//
// Kmernator/test/SamUtilsTest.cpp
//
// Author: Rob Egan
//
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

#include <fstream>
#include "mpi.h"
#include "SamUtils.h"
#include "Options.h"

int main(int argc, char **argv)
{
	bool passed = true;
	if (MPI_SUCCESS != MPI_Init(&argc, &argv))
		LOG_THROW("MPI_Init() failed: ");

	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	mpi::communicator world;
	Logger::setWorld(&world);
	GeneralOptions::getOptions().getDebug() = 2;

	if (argc != 4)
		LOG_THROW("Usage: SamUtilsTest my.bam testout.bam testout-sort.bam");

	samfile_t *fh = NULL;
	// BamStreamUtils::openSamOrBam(f)
	BamVector reads;
	BamHeaderPtr header = BamStreamUtils::readBamFile(world, argv[1], reads);
	if (0) {
		fh = BamStreamUtils::openSamOrBam(argv[1]);

		// readBamFile(f, BamVector &)
		BamStreamUtils::readBamFile(fh, reads);
		LOG_VERBOSE(1, "reads.size(): " << reads.size());
		// HACK to support without bamseek
		int blocksize = (reads.size() + size - 1) / size;
		BamVector keep;
		keep.reserve(blocksize);
		for(int i = blocksize * rank; i < blocksize * (rank+1); i++) {
			if (i == (int) reads.size())
				break;
			keep.push_back(reads[i]);
			reads[i] = NULL;
		}
		BamManager::destroyBamVector(reads);
		reads.swap(keep);
	}
	LOG_VERBOSE(1, "reads.size(): " << reads.size());

	LOG_VERBOSE(1, "Copying bam to " << argv[2]);
	unlink(argv[2]);

	SamUtils::writeBamVector(MPI_COMM_WORLD, argv[2], reads, header.get(), false);

	unlink(argv[3]);
	{
		SamUtils::MPISortBam sortem(MPI_COMM_WORLD, reads);
	}
	SamUtils::writeBamVector(MPI_COMM_WORLD, argv[3], reads, header.get(), false);

	BamManager::destroyBamVector(reads);
	if (fh != NULL)
		samclose(fh);


	if (MPI_SUCCESS != MPI_Finalize())
		LOG_THROW("MPI_Finalize() failed: ");
	return passed ? 0 : -1;
}
