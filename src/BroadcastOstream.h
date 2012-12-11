/*
 * BroadcastOstream.h
 *
 *  Created on: Sep 21, 2012
 *      Author: regan
 *
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

#ifndef BROADCASTOSTREAM_H_
#define BROADCASTOSTREAM_H_

#include <algorithm>

#include <boost/iostreams/filter/symmetric.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

#include "MemoryBufferStream.h"

template< int BUFFER_SIZE = 131072 >
class BroadcastOstreamDetail {
public:
	const static int padding = sizeof(int) * 2; // (int)offset, (int)instance, (char[]) streambuffer
	const static int MAX_SIZE = BUFFER_SIZE - padding - 128;

	BroadcastOstreamDetail(const MPI_Comm &_comm, int _broadcastRank, std::ostream &_destination) :
		broadcastRank(_broadcastRank), destination(_destination) {
		myComm = mpi::communicator(_comm, mpi::comm_duplicate);
		assert(MPI::Is_thread_main());
		assert(!omp_in_parallel());
		MPI_Comm_rank(myComm, &myRank);
		init();
	}
	~BroadcastOstreamDetail() {
		assert(MPI::Is_thread_main());
		assert(!omp_in_parallel());
		_close();
		free(buf1); buf1 = NULL;
		free(buf2); buf2 = NULL;
	}

	bool isActive() {
		return activeIn != NULL;
	}
	std::streamsize write(const char* s, std::streamsize n) {
		std::streamsize wrote = 0, remaining = n;
		while (remaining > 0) {
			wrote = writeSome(s, remaining);
			s+=wrote;
			remaining -= wrote;
		}
		return n;
	}

	std::streamsize writeSome(const char* s, std::streamsize n) {
		assert(MPI::Is_thread_main());
		assert(myRank == broadcastRank);
		assert(omp_get_thread_num() == 0);
		assert(isActive());
		int &offset = _getOffset(activeIn);
		std::streamsize wrote = std::min((std::streamsize) (MAX_SIZE-offset), n);
		memcpy(activeIn + offset + padding, s, wrote);
		offset += wrote;
		if (offset == MAX_SIZE)
			_flush();
		return wrote;
	}

protected:
	int &_getOffset(char *buf) {
		return *((int*)buf);
	}
	int &_getInstance(char *buf) {
		return *(((int*)buf) + 1);
	}
	void _swapBuffers() {
		if (activeIn == buf1) {
			activeIn = buf2;
		} else {
			activeIn = buf1;
		}
	}
	void init() {
		buf1 = (char*) calloc(MAX_SIZE + padding, 1);
		buf2 = (char*) calloc(MAX_SIZE + padding, 1);
		activeIn = buf1;
		activeOut = NULL;
		_getOffset(buf1) = 0;
		_getOffset(buf2) = 0;
	}
	void _close() {
		LOG_DEBUG_OPTIONAL(1, true, "Entered BroadcastOstream::_close(): " << this << " isActive: " << isActive());

		if (buf1 == NULL)
			return;

		if (activeIn != NULL) {
			int &offset = _getOffset(activeIn);
			if (myRank == broadcastRank) {
				if (offset == MAX_SIZE) {
					_flush(); _flush();
				} else {
					_flush();
				}
			} else {
				while (isActive())
					_flush();
			}
		}

		activeIn = NULL;
	}
	void _flush() {
		LOG_DEBUG_OPTIONAL(2, true, "BroadcastOstream::Entered flush(): " << this);
		if (omp_in_parallel()) {
			_flushThreaded();
		} else
			_flushSerial();
	}
	void _outputBuffer(char *buf) {
		int &myoffset = _getOffset(buf);
		LOG_DEBUG_OPTIONAL(2, true, "BroadcastOstream::_outputBuffer(): writing output buffer " << myoffset << " bytes: " << this);
		if (myoffset > 0) {
			destination.write(buf+padding, myoffset);
		}
		if (myRank != broadcastRank) {
			if (myoffset < MAX_SIZE) {
				LOG_DEBUG_OPTIONAL(2, true, "BroadcastOstream detected End of transmission");
				activeIn = NULL;
			}
		}
		myoffset = 0;
	}

	void _flushSerial() {
		assert(MPI::Is_thread_main());
		assert(omp_get_thread_num() == 0);
		assert(activeOut == NULL);
		assert(isActive());
		LOG_DEBUG_OPTIONAL(2, true, "BroadcastOstream::_flushSerial() Entered: " << ((myRank == broadcastRank) ? _getOffset(activeIn) : 0) );

		activeOut = activeIn;
		_swapBuffers();

		MPI_Bcast(activeOut, MAX_SIZE+padding, MPI_BYTE, broadcastRank, myComm);

		_outputBuffer(activeOut);
		activeOut = NULL;

	}
	void _flushThreaded() {
		assert(omp_in_parallel());
		if (MPI::Is_thread_main()) { // (omp_get_thread_num() == 0) {
			{
				boost::mutex::scoped_lock mylock(mutex);
				while (activeOut != NULL) {
					LOG_DEBUG_OPTIONAL(1, true, "BroadcastOstream::_flushThreaded(): Waiting for flushReady.");
					flushReady.wait(mylock);
				}
				_flushSerial();
			}
			broadcastReady.notify_one();
		} else {
			{
				boost::mutex::scoped_lock mylock(mutex);
				while (activeOut == NULL) {
					LOG_DEBUG_OPTIONAL(1, true, "BroadcastOstream::_flushThreaded(): Waiting for broadcastReady.");
					broadcastReady.wait(mylock);
				}
				_outputBuffer(activeOut);
				activeOut = NULL;
			}
			flushReady.notify_one();
		}
	}


private:
	mpi::communicator myComm;
	int broadcastRank;
	std::ostream &destination;
	int myRank;
	int myInstance;
	char *buf1, *buf2, *activeIn, *activeOut;
	boost::mutex mutex;
	boost::condition_variable flushReady, broadcastReady;

};

typedef stream_impl_template< BroadcastOstreamDetail< > > BroadcastOstream;

#endif /* BROADCASTOSTREAM_H_ */
