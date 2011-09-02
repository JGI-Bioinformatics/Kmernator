/*
 * MPIBase.h
 *
 *  Created on: Feb 2, 2011
 *      Author: regan
 */

#ifndef MPIBASE_H_
#define MPIBASE_H_

#include "config.h"

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

#include "boost/optional.hpp"

typedef boost::optional< mpi::request > MPIOptionalRequest;
typedef boost::optional< mpi::status > MPIOptionalStatus;

#include <boost/thread/thread.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#define WAIT_MS 1
#define WAIT_AND_WARN( iterations, warningMessage ) \
	{ if ((iterations % (30000/WAIT_MS)) == 0) LOG_WARN(1, warningMessage << " waiting in loop: " << iterations);  \
	  boost::this_thread::sleep( boost::posix_time::milliseconds(WAIT_MS) ); \
	}


#endif /* MPIBASE_H_ */
