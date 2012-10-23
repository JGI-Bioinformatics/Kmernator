/*
 * BroadcastOstream.h
 *
 *  Created on: Sep 21, 2012
 *      Author: regan
 */

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

template< int BUFFER_SIZE = 524288 >
class BroadcastOstreamDetail {
public:
	typedef std::map< int, BroadcastOstreamDetail *> InstanceMap;
	const static int padding = sizeof(int) * 2; // (int)offset, (int)instance, (char[]) streambuffer
	const static int MAX_SIZE = BUFFER_SIZE - padding - 128;
	inline static boost::mutex &getMutex() {
		static boost::mutex mutex;
		return mutex;
	}
	inline static boost::condition_variable &getFlushReady() {
		static boost::condition_variable flushReady;
		return flushReady;
	}
	inline static boost::condition_variable &getBroadcastReady() {
		static boost::condition_variable broadcastReady;
		return broadcastReady;
	}

	BroadcastOstreamDetail(MPI_Comm &_comm, int _broadcastRank, std::ostream &_destination) :
		broadcastRank(_broadcastRank), destination(_destination) {
		assert(!omp_in_parallel());
		if (MPI_Comm_dup(_comm, &comm) != MPI_SUCCESS)
			LOG_THROW("Could not duplicate MPI Communicator");
		MPI_Comm_rank(comm, &myRank);
		init();
	}
	~BroadcastOstreamDetail() {
		assert(!omp_in_parallel());
		_close();
		free(buf1); buf1 = NULL;
		free(buf2); buf2 = NULL;
		if (MPI_Comm_free(&comm) != MPI_SUCCESS)
			LOG_THROW("Could not free MPI communicator in ~BroadcastOstream");
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

/*
	std::ostream &operator<<(std::string str) {
		assert(myRank == broadcastRank);
		assert(isActive());
		int len = str.length();
		assert(len < MAX_SIZE);
		int &offset = _getOffset(activeIn);
		if (len + offset >= MAX_SIZE) {
			_flush();
		}
		memcpy(activeIn+offset+padding, str.c_str(), len);
		offset += len;
		destination << str;
		return destination;
	}
*/

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
		assert(omp_get_thread_num() == 0);
		assert(activeOut == NULL);
		assert(isActive());
		LOG_DEBUG_OPTIONAL(2, true, "BroadcastOstream::_flushSerial() Entered: " << ((myRank == broadcastRank) ? _getOffset(activeIn) : 0) );

		activeOut = activeIn;
		_swapBuffers();

		MPI_Bcast(activeOut, MAX_SIZE+padding, MPI_BYTE, broadcastRank, comm);

		_outputBuffer(activeOut);
		activeOut = NULL;

	}
	void _flushThreaded() {
		assert(omp_in_parallel());
		if (omp_get_thread_num() == 0) {
			{
				boost::mutex::scoped_lock mylock(getMutex());
				while (activeOut != NULL) {
					LOG_DEBUG_OPTIONAL(1, true, "BroadcastOstream::_flushThreaded(): Waiting for flushReady.");
					getFlushReady().wait(mylock);
				}
				_flushSerial();
			}
			getBroadcastReady().notify_one();
		} else {
			{
				boost::mutex::scoped_lock mylock(getMutex());
				while (activeOut == NULL) {
					LOG_DEBUG_OPTIONAL(1, true, "BroadcastOstream::_flushThreaded(): Waiting for broadcastReady.");
					getBroadcastReady().wait(mylock);
				}
				_outputBuffer(activeOut);
				activeOut = NULL;
			}
			getFlushReady().notify_one();
		}
	}


private:
	MPI_Comm comm;
	int broadcastRank;
	std::ostream &destination;
	int myRank;
	int myInstance;
	char *buf1, *buf2, *activeIn, *activeOut;

};

typedef stream_impl_template< BroadcastOstreamDetail< > > BroadcastOstream;

#endif /* BROADCASTOSTREAM_H_ */
