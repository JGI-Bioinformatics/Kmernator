//
// Kmernator/test/SamUtilsTest.cpp
//
// Author: Rob Egan
//
// Copyright 2012 The Regents of the University of California.
// All rights reserved.
//
// The United States Government has rights in this work pursuant
// to contracts DE-AC03-76SF00098, W-7405-ENG-36 and/or
// W-7405-ENG-48 between the United States Department of Energy
// and the University of California.
//
// Redistribution and use in source and binary forms are permitted
// provided that: (1) source distributions retain this entire
// copyright notice and comment, and (2) distributions including
// binaries display the following acknowledgement:  "This product
// includes software developed by the University of California,
// JGI-PSF and its contributors" in the documentation or other
// materials provided with the distribution and in all advertising
// materials mentioning features or use of this software.  Neither the
// name of the University nor the names of its contributors may be
// used to endorse or promote products derived from this software
// without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE.
//

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
	GeneralOptions::getOptions().getDebug() = 1;

	if (argc != 4)
		LOG_THROW("Usage: SamUtilsTest my.bam testout.bam testout-sort.bam");

	// BamStreamUtils::openSamOrBam(f)
	samfile_t *fh = BamStreamUtils::openSamOrBam(argv[1]);

	// readBamFile(f, BamVector &)
	BamVector reads;
	BamStreamUtils::readBamFile(fh, reads);
	LOG_VERBOSE(1, "reads.size(): " << reads.size());

	// HACK until bamseek is available
	int blocksize = (reads.size() + size - 1) / size;
	{
		BamVector keep;
		keep.reserve(blocksize);
		for(int i = blocksize * rank; i < blocksize * (rank+1); i++) {
			if (i == (int) reads.size())
				break;
			keep.push_back(reads[i]);
			reads[i] = NULL;
		}
		SamUtils::destroyBamVector(reads);
		reads.swap(keep);
	}
	LOG_VERBOSE(1, "reads.size(): " << reads.size());

	LOG_VERBOSE(1, "Copying bam to " << argv[2]);
	unlink(argv[2]);

	SamUtils::writeBamVector(MPI_COMM_WORLD, argv[2], reads, fh->header, false);

	unlink(argv[3]);
	{
		SamUtils::MPISortBam sortem(MPI_COMM_WORLD, reads);
	}
	SamUtils::writeBamVector(MPI_COMM_WORLD, argv[3], reads, fh->header, false);

	SamUtils::destroyBamVector(reads);
	samclose(fh);




	if (MPI_SUCCESS != MPI_Finalize())
		LOG_THROW("MPI_Finalize() failed: ");
	return passed ? 0 : -1;
}
