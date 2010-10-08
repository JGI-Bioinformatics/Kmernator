//
// Kmernator/src/DistributedFunctions.h
//
// Author: Rob Egan, Craig Furman
//
// Copyright 2010 The Regents of the University of California.
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

#ifndef DISTRIBUTED_FUNCTIONS_H_
#define DISTRIBUTED_FUNCTIONS_H_

#include "config.h"
#include "KmerSpectrum.h"

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

void reduceOMPThreads(mpi::communicator &world) {
	int numThreads = omp_get_max_threads();
	numThreads = all_reduce(world, numThreads, mpi::minimum<int>());
	omp_set_num_threads(numThreads);
	if (world.rank() == 0)
		LOG_DEBUG_MT(1, world.rank() << ": set OpenMP threads to " << numThreads);
}

template<typename So, typename We, typename Si = TrackingDataSingleton>
class DistributedKmerSpectrum : public KmerSpectrum<So, We, Si>
{
	typedef KmerSpectrum<So, We, Si> KS;
public:
	DistributedKmerSpectrum(unsigned long buckets = 0, bool separateSingletons = true)
	: KS(buckets, separateSingletons) {}
	~DistributedKmerSpectrum() {}
	DistributedKmerSpectrum &operator=(const KS &other) {
		*((KS*) this) = other;
		return *this;
	}
};

#endif /* DISTRIBUTED_FUNCTIONS_H_ */
