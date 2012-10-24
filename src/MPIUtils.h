/*
 * MPI.h
 *
 *  Created on: Oct 24, 2012
 *      Author: regan
 */

#ifndef MPI_UTILS_H_
#define MPI_UTILS_H_

#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>

#include "boost/thread.hpp"

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

#include "mpi.h"
#include "Options.h"
#include "Log.h"


// collective
void reduceOMPThreads(mpi::communicator &world) {
	Options::getOptions().validateOMPThreads();
#ifdef _USE_OPENMP
	int numThreads = omp_get_max_threads();
	numThreads = all_reduce(world, numThreads, mpi::minimum<int>());
	omp_set_num_threads(numThreads);
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "set OpenMP threads to " << numThreads);
#endif
}

// collective
void validateMPIWorld(mpi::communicator &world) {
	int provided;
	MPI_Query_thread(&provided);
#ifdef _USE_OPENMP
	if (provided != MPI_THREAD_FUNNELED && omp_get_max_threads() > 1) {
		if (world.rank() == 0)
			LOG_WARN(1, "Your version of MPI does not support MPI_THREAD_FUNNELED (" << provided << "), reducing OpenMP threads to 1")
			omp_set_num_threads(1);
	}
#endif
	reduceOMPThreads(world);
}

template< typename OptionsTempl = NullOptions >
class ScopedMPIComm {
public:
	typedef boost::shared_ptr< mpi::communicator > Comm;

	ScopedMPIComm(int argc, char *argv[]) {
		_world = initializeWorldAndOptions(argc, argv);
	}
	~ScopedMPIComm() {
		LOG_DEBUG_OPTIONAL(1, Logger::isMaster(), "Finishing");
		Logger::setWorld(NULL);
		_world.reset();
		MPI_Finalize();
	}
	operator mpi::communicator() const {
		return *_world;
	}
	operator mpi::communicator&() {
		return *_world;
	}
	operator const mpi::communicator&() const {
		return *_world;
	}

	operator MPI_Comm() const {
		return (MPI_Comm) *_world;
	}
	int rank() const {
		return _world->rank();
	}
	int size() const {
		return _world->size();
	}
	void barrier() const {
		_world->barrier();
	}
	void abort(int i) const {
		_world->abort(i);
	}
	mpi::communicator split(int color) {
		return _world->split(color);
	}
	mpi::communicator split(int color, int key) {
		return _world->split(color, key);
	}

protected:
	Comm initializeWorldAndOptions(int argc, char *argv[]) {
		int threadProvided;
		int threadRequest = omp_get_max_threads() == 1 ? MPI_THREAD_SINGLE : MPI_THREAD_FUNNELED;
		MPI_Init_thread(&argc, &argv, threadRequest, &threadProvided);
		mpi::environment env(argc, argv);
		Comm comm(new mpi::communicator(MPI_COMM_WORLD, mpi::comm_duplicate));
		MPI_Comm_set_errhandler( *comm, MPI::ERRORS_THROW_EXCEPTIONS );

		try {
			Logger::setWorld(comm.get());

			OptionsTempl::parseOpts(argc, argv);

			if (GeneralOptions::getOptions().getGatheredLogs())
				Logger::setWorld(comm.get(), Options::getOptions().getDebug() >= 2);

			validateMPIWorld(*comm);

		} catch (...) {
			std::cerr << "Could not initializeWorldAndOptions!" << std::endl;
			MPI_Finalize();
			exit(1);
		}
		comm->barrier();
		return comm;
	}

private:
	Comm _world;

};


class DistributedDirectoryManagement {
public:

	// collective
	static std::string _getRankSubDir(mpi::communicator &world, std::string prefix) {
		int subRank = world.rank() / 256;
		std::stringstream ss;
		ss << prefix << "/rank-subdirs-" << world.size() << "-0x" << std::hex << subRank;
		std::string subDir = ss.str();
		return subDir;
	}
	static std::string _makeRankSubDir(mpi::communicator &world, std::string prefix) {
		mkdir(prefix.c_str(), 0777); // all ranks must mkdir if writing to local disks
		int subRank = world.rank() / 256;
		std::string subDir = _getRankSubDir(world, prefix);
		LOG_DEBUG_OPTIONAL(1, world.rank() == subRank * 256, "Making rank-subdirs: " << subDir);
		mkdir(subDir.c_str(), 0777); // all ranks must mkdir if writing to local disks
		return subDir;
	}
	static void _rmRankSubDir(mpi::communicator &world, std::string prefix) {
		std::string subDir = _getRankSubDir(world, prefix);
		rmdir(subDir.c_str());  // all ranks must rmdir if writing to local disks
	}

	static std::string getRankSubDir(mpi::communicator &world, std::string prefix) {
		std::string subSubDir = _getRankSubDir(world, prefix);
		std::string subDir = subSubDir + "/" + boost::lexical_cast<std::string>(world.rank()) + "of" + boost::lexical_cast<std::string>(world.size()) + "/";
		return subDir;
	}
	static std::string makeRankSubDir(mpi::communicator &world, std::string prefix) {
		std::string subSubDir = _makeRankSubDir(world, prefix);
		std::string subDir = getRankSubDir(world, prefix);
		LOG_DEBUG(2, "getRankSubDir(" << prefix << "): " << subDir);
		mkdir(subDir.c_str(), 0777);
		return subDir;
	}

	static void rmRankSubDir(mpi::communicator &world, std::string prefix) {
		std::string subSubDir = _getRankSubDir(world,prefix);
		std::string subDir = getRankSubDir(world, prefix);
		rmdir(subDir.c_str());
		rmdir(subSubDir.c_str());
	}

};

// collective
void niceBarrier(mpi::communicator &world, int waitMs = 1) {
	mpi::communicator tmpWorld(world, mpi::comm_duplicate);
	int rank = tmpWorld.rank();
	int size = tmpWorld.size();
	char buf, buf2;
	buf  = '\0';
	buf2 = '\0';
	mpi::request rreq, rreq2;
	mpi::request sreq, sreq2;

	LOG_DEBUG_OPTIONAL(2, true, "Entering niceBarrier");
	if (rank == 0)
		sreq = tmpWorld.isend((rank+1) % size, 0, buf);

	rreq = tmpWorld.irecv((rank+size-1) % size, 0, buf2);
	while (! rreq.test() ) {
		boost::this_thread::sleep( boost::posix_time::milliseconds(waitMs) );
	}

	if (rank != 0)
		sreq = tmpWorld.isend((rank+1) % size, 0, buf);

	while (! sreq.test() ) {
		boost::this_thread::sleep( boost::posix_time::milliseconds(waitMs) );
	}

	if (rank == 0)
		sreq2 = tmpWorld.isend((rank+1) % size, 1, buf);

	rreq2 = tmpWorld.irecv((rank+size-1) % size, 1, buf2);
	while (! rreq2.test() ) {
		boost::this_thread::sleep( boost::posix_time::milliseconds(waitMs) );
	}

	if (rank != 0)
		sreq2 = tmpWorld.isend((rank+1) % size, 1, buf);

	while (! sreq2.test() ) {
		boost::this_thread::sleep( boost::posix_time::milliseconds(waitMs) );
	}

	LOG_DEBUG_OPTIONAL(2, true, "Exiting niceBarrier");

}


#endif /* MPI_UTILS_H_ */
