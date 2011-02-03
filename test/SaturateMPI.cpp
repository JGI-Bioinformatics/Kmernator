/*
 * SaturateMPI.cpp
 *
 *  Created on: Nov 12, 2010
 *      Author: regan
 */


#include <stdio.h>
#include <omp.h>
#include "config.h"
#include "Log.h"
#include "Options.h"
#include "MPIBase.h"

#include <boost/date_time/posix_time/posix_time.hpp>

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

void all2all(mpi::communicator &world, int numMessages, int messageSize)
{
	std::vector< std::vector<char> > inbuf(world.size(), std::vector<char>(messageSize, ' ')),
			outbuf(world.size(), std::vector<char>(messageSize, ' '));
	std::vector< MPIOptionalRequest > recvRequests(world.size());
	std::vector< MPIOptionalRequest > sendRequests(world.size());
	std::vector< long > recvCounts(world.size(),0);
	std::vector< long > sendCounts(world.size(),0);

	world.barrier();
	boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::ptime lastActivity = start, lastStallMessage = start;
	// initialize and send/receive first message
	for(int i=1; i < world.size(); i++) {
		int dest = (i + world.rank()) % world.size();
		for(int j = 0; j < messageSize ; j++)
			outbuf[dest][j] = (char) (j%26)+65;
		recvRequests[dest] = world.irecv(dest, 0, &inbuf[dest][0], messageSize);
		sendRequests[dest] = world.isend(dest, 0, &outbuf[dest][0], messageSize);
	}

	LOG_DEBUG(2, "Starting loop");
	long inactiveIterations = 0;
	long activeIterations = 0;
	bool hasActivity = false;
	lastActivity = boost::posix_time::microsec_clock::local_time();
	bool done = false;
	while (! done) {
		if (hasActivity) {
			inactiveIterations = 0;
			hasActivity = false;
			if (++activeIterations > 10000) {
				activeIterations = 0;
				lastActivity = boost::posix_time::microsec_clock::local_time();
			}
		}
		int countSentDone = 1;
		int countRecvDone = 1;
		for(int i=1; i < world.size(); i++) {
			int dest = (i + world.rank()) % world.size();

			if (!!recvRequests[dest]) {
				MPIOptionalStatus status = recvRequests[dest].get().test();
				if (!! status) {
					if (++recvCounts[dest] < numMessages) {
						hasActivity = true;
						recvRequests[dest] = world.irecv(dest, 0, &inbuf[dest][0], messageSize);
					} else {
						recvRequests[dest] = MPIOptionalRequest();
						countRecvDone++;
					}
				}
			} else {
				countRecvDone++;
			}
			if (!!sendRequests[dest]) {
				MPIOptionalStatus status = sendRequests[dest].get().test();
				if (!! status) {
					if (++sendCounts[dest] < numMessages) {
						hasActivity = true;
						sendRequests[dest] = world.isend(dest, 0, &outbuf[dest][0], messageSize);
					} else {
						sendRequests[dest] = MPIOptionalRequest();
						countSentDone++;
					}
				}
			} else {
				countSentDone++;
			}
		}
		if (! hasActivity) {
			if (++inactiveIterations > 10000) {
				inactiveIterations = 0;
				activeIterations += 10000; // reset on the next activity
				//boost::this_thread::sleep( boost::posix_time::milliseconds(1) );
				boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
				boost::posix_time::time_duration elapsed = now - lastActivity;
				boost::posix_time::time_duration elapsedMessage = now - lastStallMessage;
				if (elapsed.total_milliseconds() > 1000 && elapsedMessage.total_milliseconds() > 2500) {
					std::stringstream ss;
					ss << "Sending still pending: ";
					for(int i=1; i < world.size(); i++) {
						int dest = (i + world.rank()) % world.size();
						if (sendCounts[dest] < numMessages)
							ss << dest << ", ";
					}
					ss << "Receiving still pending: ";
					for(int i=1; i < world.size(); i++) {
						int dest = (i + world.rank()) % world.size();
						if (recvCounts[dest] < numMessages)
							ss << dest << ", ";
					}
					std::string s = ss.str();

					LOG_WARN(1, "Detected stalling for: " << elapsed.total_milliseconds() << " ms. " << s);
					lastStallMessage = boost::posix_time::microsec_clock::local_time();
				}
			}
		}
		LOG_DEBUG(2, "Recv: " << countRecvDone << " Sent: " << countSentDone);
		done = (countSentDone == world.size() && countRecvDone == world.size());
	}

	world.barrier();
	boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration elapsed = end - start;

	double megaBytes = (double) messageSize * (double) numMessages * double (world.size() - 1) / 1024.0 / 1024.0;
	double rate = megaBytes / (double) elapsed.total_milliseconds() * 1000.0;
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, numMessages << " message rounds, " << messageSize << " bytes: "
			<< numMessages * (world.size()-1) << " total messages, " << megaBytes << " MB, " << rate << " MB/s, " << rate * 8 << " Mbit/s, " << elapsed.total_milliseconds() << " ms");

}

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	mpi::environment env(argc, argv);
	mpi::communicator world;

	int maxCycles = 100;
	if (argc > 1)
		maxCycles = atoi(argv[1]);
	int maxMessageSize = 1024*1024*64; // 65535 kb
	if (argc > 2)
		maxMessageSize = atoi(argv[2]) * 1024;
	int maxNumMessages = 1024*1024*256; // 256 million+ messages
	if (argc > 3)
		maxNumMessages = atoi(argv[3]);
	int maxDataSize = 1024*1024*256; // 256 MB
	if (argc > 4)
		maxDataSize = atoi(argv[4])*1024;

	int minDataSize = 1024*1024*16;  // 16 MB

	//Options::getDebug() = 2;
	Logger::setWorld(&world);
	{
		char hostname[128];
		gethostname(hostname, 128);
		LOG_VERBOSE(1, "Starting on " << hostname << " with Thread Support: " << provided);
	}
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Starting with " << maxCycles << " cycles");


	for(int cycle = 0 ; cycle < maxCycles ; cycle++) {
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Starting cycle " << cycle);
		for(int num = 1; num < maxNumMessages ; num *= 2) {
			for(int size = 512; size < maxMessageSize ; size *= 2) {
				if (num * size * (world.size()-1) < minDataSize || num * size * (world.size()-1) > maxDataSize)
					continue;

				for (int i = 0 ; i < 3 ; i++)
					all2all(world, num, size);
			}
		}
	}

	MPI_Finalize();
	return 0;
}

