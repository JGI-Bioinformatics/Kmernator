//
// Kmernator/src/MPIBuffer.h
//
// Author: Rob Egan
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

#ifndef MPIBUFFER_H_
#define MPIBUFFER_H_

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

#ifdef _USE_OPENMP
#define X_OPENMP_CRITICAL_MPI
#endif

#include "boost/optional.hpp"
#include <boost/thread/thread.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <vector>

#define MPI_BUFFER_DEFAULT_SIZE (16 * 1024)

template <typename C, int BufferSize = MPI_BUFFER_DEFAULT_SIZE>
class MPIMessageBufferBase {
public:
	static const int MESSAGE_BUFFER_SIZE = BufferSize;
	typedef C MessageClass;
	typedef boost::optional< mpi::request > OptionalRequest;
	typedef boost::optional< mpi::status > OptionalStatus;

protected:
	mpi::communicator _world;
	int _messageSize;
	char *_message;
	long _count;

public:
	MPIMessageBufferBase(mpi::communicator &world, int messageSize)
	  : _world(world), _messageSize(messageSize), _count(0) {
		_message = new char[ getMessageSize() ];
		assert(getMessageSize() >= (int) sizeof(MessageClass));
	}
	~MPIMessageBufferBase() {
		delete [] _message;
	}
	MessageClass *getTmpMessage() {
		return (MessageClass*) _message;
	}
	inline mpi::communicator &getWorld() {
		return _world;
	}
	inline int getMessageSize() {
		return _messageSize;
	}
	inline long getCount() {
		return _count;
	}
	inline void newMessage() {
		#pragma omp atomic
		_count++;
	}
};

template <typename C, int BufferSize = MPI_BUFFER_DEFAULT_SIZE>
class MPIRecvMessageBuffer: public MPIMessageBufferBase<C, BufferSize> {
public:
	typedef MPIMessageBufferBase<C, BufferSize> BufferBase;
	typedef typename BufferBase::OptionalRequest OptionalRequest;
	typedef typename BufferBase::OptionalStatus OptionalStatus;
	typedef C MessageClass;
protected:
	char *_recvBuffers;
	OptionalRequest *_requests;
	int _numCheckpoints;
	int _tag;
	int _source;
	int _size;
	std::vector<int> _requestAttempts;

public:
	MPIRecvMessageBuffer(mpi::communicator &world, int messageSize, int tag = mpi::any_tag) :
		BufferBase(world, messageSize), _numCheckpoints(0), _tag(tag), _source(mpi::any_source), _size(0) {
		_recvBuffers = new char[ BufferBase::MESSAGE_BUFFER_SIZE * this->getWorld().size() ];
		_requests = new OptionalRequest[ this->getWorld().size() ];
		_requestAttempts.resize(this->getWorld().size(), 0);
		for(int destRank = 0; destRank < this->getWorld().size(); destRank++) {
			_requests[destRank] = irecv(destRank);
		}
	}
	~MPIRecvMessageBuffer() {
		cancelAllRequests();
		delete [] _recvBuffers;
		delete [] _requests;
	}

	inline int getSource() {
		return _source;
	}
	inline int getTag() {
		return _tag;
	}
	inline int getSize() {
		return _size;
	}
	void cancelAllRequests() {
		for(int source = 0 ; source < this->getWorld().size() ; source++) {
			if (!!_requests[source]) {
				_requests[source].get().cancel();
				_requests[source] = OptionalRequest();
			}
		}
	}
	bool receiveIncomingMessage(int rankSource) {
		bool returnValue = false;
		OptionalStatus optionalStatus;
		OptionalRequest &optionalRequest = _requests[rankSource];

		if (!optionalRequest) {
			LOG_WARN(1, "Detected non-pending request for tag: " << _tag);
			optionalRequest = irecv(rankSource);
		}
		if (!!optionalRequest) {
			bool retry = ++_requestAttempts[rankSource] > 1000;

			optionalStatus = optionalRequest.get().test();
			if (!optionalStatus && retry) {
				LOG_WARN(1, "Canceling pending request that looks to be stuck tag: " << _tag);
				optionalRequest.get().cancel();
				optionalStatus = optionalRequest.get().test();
				if (!optionalStatus) {
					optionalRequest = irecv(rankSource);
					optionalStatus = optionalRequest.get().test();
				}
			}

			if (!!optionalStatus) {
				mpi::status status = optionalStatus.get();
				process(status, rankSource);
				optionalRequest = irecv(rankSource);
				returnValue = true;
			}
		}
		return returnValue;
	}
	int receiveAllIncomingMessages() {
		int messages = 0;

		LOG_DEBUG(3, _tag << ": receiving all messages for tag. cp: " << getNumCheckpoints() << " msgCount: " << this->getCount());
		bool mayHavePending = true;
		while (mayHavePending) {
			mayHavePending = false;
			for(int rankSource = 0 ; rankSource < this->getWorld().size(); rankSource++) {
				if (receiveIncomingMessage(rankSource)) {
					messages++;
					mayHavePending = true;
				}
			}
		}
		if (messages > 0)
			LOG_DEBUG(3, _tag << ": processed messages: " << messages);

		return messages;
	}
	void finalize(int checkpointFactor = 1) {
		LOG_DEBUG(2, _tag << ": Entering recv checkpoint: " << getNumCheckpoints());
		while ( ! reachedCheckpoint(checkpointFactor) ) {
			receiveAllIncomingMessages();
			boost::this_thread::sleep( boost::posix_time::milliseconds(2) );
		}
		_numCheckpoints = 0;
	}

	void checkpoint() {
		#pragma omp atomic
		_numCheckpoints++;
		LOG_DEBUG(2, _tag << ": checkpoint received:" << _numCheckpoints);
	}
	inline int getNumCheckpoints() const {
		return _numCheckpoints;
	}
	bool reachedCheckpoint(int checkpointFactor = 1) {
		return _numCheckpoints == this->getWorld().size() * checkpointFactor;
	}
private:
	OptionalRequest irecv(int sourceRank) {
		OptionalRequest oreq;
		LOG_DEBUG(4, "Starting irecv for " << sourceRank << "," << _tag);
#ifdef OPENMP_CRITICAL_MPI
		#pragma omp critical (MPI_buffer)
#endif
		oreq = this->getWorld().irecv(sourceRank, _tag, _recvBuffers + sourceRank*BufferBase::MESSAGE_BUFFER_SIZE, BufferBase::MESSAGE_BUFFER_SIZE);
		_requestAttempts[sourceRank] = 0;
		assert(!!oreq);
		return oreq;
	}
	int processMessages(int sourceRank, int size) {
		char *msg, *end, *start = _recvBuffers + sourceRank*BufferBase::MESSAGE_BUFFER_SIZE;
		msg = start;
		end = (start+size);
		int count = 0;
		while (msg != end) {
			((MessageClass*) msg)->process(this);
			msg += this->getMessageSize();
		}
		return count;
	}
	bool process(mpi::status &status, int rankSource) {

		bool wasMessage = false;
		int source = status.source();
		int tag =  status.tag();
		int size = status.count<char>().get();
		if (status.cancelled()) {
			LOG_WARN(1, _tag << ": request was successfully canceled from " << source << " size " << size);
		} else {

			assert( rankSource == source );
			assert( _tag == tag );
			_source = source;
			_size = size;

			this->newMessage();
			LOG_DEBUG(3, _tag << ": received message " << this->getCount() << " from " << source << "," << tag << " size " << size << " probe attempts: " << _requestAttempts[rankSource]);

			if (size == 0) {
				checkpoint();
				LOG_DEBUG(2, _tag << ": got checkpoint from " << source);
			} else {
				processMessages(source, size);
			}
			wasMessage = true;
		}
		return wasMessage;
	}

};


template <typename C, int BufferSize = MPI_BUFFER_DEFAULT_SIZE>
class MPISendMessageBuffer: public MPIMessageBufferBase<C, BufferSize> {
public:
	typedef MPIMessageBufferBase<C, BufferSize> BufferBase;
	typedef MPIRecvMessageBuffer<C, BufferSize> RecvBuffer;
	typedef typename BufferBase::OptionalRequest OptionalRequest;
	typedef typename BufferBase::OptionalStatus OptionalStatus;
	typedef C MessageClass;
	typedef std::vector< RecvBuffer* > RecvBufferCallbackVector;
protected:
	char *_sendBuffers;
	int  *_offsets;
	RecvBufferCallbackVector _receivingBufferCallbacks;

public:
	MPISendMessageBuffer(mpi::communicator &world, int messageSize) :
		BufferBase(world, messageSize) {
		_sendBuffers = new char[ BufferBase::MESSAGE_BUFFER_SIZE * world.size() ];
		_offsets = new int[ this->getWorld().size() ];
		for(int destRank = 0; destRank < this->getWorld().size(); destRank++) {
			_offsets[destRank] = 0;
		}
	}
	~MPISendMessageBuffer() {
		delete [] _sendBuffers;
		delete [] _offsets;
	}

	// receive buffers to flush before and/or during send
	void addCallback( RecvBuffer& receiveBuffer ) {
		_receivingBufferCallbacks.push_back( &receiveBuffer );
		receiveAnyIncoming();
		assert(receiveBuffer.getCount() == 0);
	}

	MessageClass *bufferMessage(int rankDest, int tagDest) {
		bool wasSent;
		return bufferMessage(rankDest, tagDest, wasSent);
	}
	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest, bool &wasSent) {
		char *buffStart = _sendBuffers + BufferBase::MESSAGE_BUFFER_SIZE * rankDest;
		int &offset = _offsets[rankDest];
		if (offset + this->getMessageSize() > BufferBase::MESSAGE_BUFFER_SIZE) {
			flushMessageBuffer(rankDest, tagDest, buffStart, offset);
			wasSent = true;
		} else {
			wasSent = false;
		}
		MessageClass *buf = (MessageClass *) (buffStart+offset);
		offset += this->getMessageSize();

		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg) {
		char *buf = (char*) bufferMessage(rankDest, tagDest);
		memcpy(buf, (char*) msg, this->getMessageSize());
	}

protected:
	void flushMessageBuffer(int rankDest, int tagDest, char *buffer, int &offset, bool sendZeroMessage = false) {
		receiveAnyIncoming();
		if (offset > 0 || sendZeroMessage) {
			mpi::request request;
			LOG_DEBUG(3, "sending message to " << rankDest << ", " << tagDest << " size " << offset);
#ifdef OPENMP_CRITICAL_MPI
#pragma omp critical (MPI_buffer)
#endif
			request = this->getWorld().isend(rankDest, tagDest, buffer, offset);
			this->newMessage();
			int count = 0;
			OptionalStatus optStatus = request.test();
			while ( ! optStatus ) {
				if (++count > 1000) {
					request.cancel();
					if ( ! request.test() ) {
						request = this->getWorld().isend(rankDest, tagDest, buffer, offset);
						count = 0;
						LOG_WARN(1, "Canceled and retried pending message to " << rankDest << ", " << tagDest << " size " << offset << " msgCount " << this->getCount());
					}
				}
				LOG_DEBUG(3, "waiting for send to finish to " << rankDest << ", " << tagDest << " size " << offset << " msgCount " << this->getCount() << " attempt count " << count);
				receiveAnyIncoming();
				boost::this_thread::sleep( boost::posix_time::milliseconds(2) );
				optStatus = request.test();
			}
			mpi::status status = optStatus.get();
			// hmmm looks like error is sometimes populated with junk...
			//if (status.error() > 0)
			//	LOG_WARN(1, "sending message returned an error: " << status.error());
			LOG_DEBUG(3, "finished sending message to " << rankDest << ", " << tagDest << " size " << offset << " msgCount " << this->getCount());
		}
		offset = 0;
	}
public:
	void flushMessageBuffer(int rankDest, int tagDest, bool sendZeroMessage = false) {
		int &offset = _offsets[rankDest];
		char *buffer = _sendBuffers + BufferBase::MESSAGE_BUFFER_SIZE * rankDest;
		flushMessageBuffer(rankDest, tagDest, buffer, offset, sendZeroMessage);
	}
	void flushAllMessageBuffers(int tagDest, bool sendZeroMessage = false) {
		for(int rankDest = 0 ; rankDest < this->getWorld().size(); rankDest++)
			flushMessageBuffer(rankDest, tagDest, sendZeroMessage);
	}
	int receiveAnyIncoming() {
		int count = 0;
		for(unsigned int i = 0; i < _receivingBufferCallbacks.size(); i++)
			count += _receivingBufferCallbacks[i]->receiveAllIncomingMessages();
		return count;
	}
	void finalize(int tagDest) {
		LOG_DEBUG(3, "entering finalize() for tag " << tagDest);
		flushAllMessageBuffers(tagDest);
		// send zero message buffer as checkpoint signal to stop
		LOG_DEBUG(3, "sending checkpoint tag " << tagDest);
		flushAllMessageBuffers(tagDest, true);
		LOG_DEBUG(2, "sent checkpoints tag " << tagDest);
	}
};

#if 0
template <typename C, int BufferSize = 16 * 1024>
class MPIMessageBuffersOMP {
public:
	static const int MESSAGE_BUFFER_SIZE = BufferSize;
	typedef C MessageClass;
	typedef boost::optional< mpi::request > OptionalRequest;
	typedef boost::optional< mpi::status > OptionalStatus;

protected:
	mpi::communicator &_world;
	int _messageSize;
	int _numCheckpoints;
	int _numThreads;
	char *_message;
	char **_recvBuffers;
	char **_sendBuffers;
	int  **_offsets;
	OptionalRequest **_requests;

public:
	MPIMessageBuffersOMP(mpi::communicator &world, int messageSize) : _world(world), _messageSize(messageSize), _numCheckpoints(0) {
		_numThreads = omp_get_max_threads();
		_message = new char * [ _numThreads ];
		_recvBuffers = new char * [ _numThreads ];
		_sendBuffers = new char * [ _numThreads ];
		_offsets = new int * [ _numThreads ];
		_requests = new OptionalRequest * [ _numThreads ];

		#pragma omp parallel num_threads(_numThreads)
		{
			int threadId = omp_get_thread_num();
			_message[threadId] = new char[ _messageSize ];
			_recvBuffers[threadId] = new char[ MESSAGE_BUFFER_SIZE * world.size() ];
			_sendBuffers[threadId] = new char[ MESSAGE_BUFFER_SIZE * world.size() ];
			_offsets[threadId] = new int[ _world.size() ];

			_requests[threadId] = new OptionalRequest[ _world.size() ];
			for(int destRank = 0; destRank < _world.size(); destRank++) {
				_offsets[threadId][destRank] = 0;
				_requests[threadId][destRank] = _world.irecv(destRank, threadId, _recvBuffers[threadId] + destRank*MESSAGE_BUFFER_SIZE, MESSAGE_BUFFER_SIZE);
			}
		}
	}
	~MPIMessageBuffersOMP() {
		#pragma omp parallel num_threads(_numThreads)
		{
			int threadId = omp_get_thread_num();
			delete [] _message[threadId];
			delete [] _recvBuffers[threadId];
			delete [] _sendBuffers[threadId];
			delete [] _offsets[threadId];
			delete [] _requests[threadId];
		}
		delete [] _message;
		delete [] _recvBuffers;
		delete [] _sendBuffers;
		delete [] _offsets;
		delete [] _requests;
	}
	inline int getMessageSize() {
		return _messageSize;
	}
	MessageClass *getTmpMessage() {
		return getTmpMessage(omp_get_thread_num());
	}
	MessageClass *getTmpMessage(int threadId) {
		return _message[threadId];
	}
	int processMessages(char *start, int size) {
		char *msg, *end;
		msg = start;
		end = (start+size);
		int count = 0;
		while (msg != end) {
			_processMessage((MessageClass*) msg);
			msg += getMessageSize();
		}
		return count;
	}
	virtual void _processMessage(MessageClass *msg) = 0;

	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest) {
		int threadId = omp_get_thread_num();
		char *buffStart = _sendBuffers[threadId] + MESSAGE_BUFFER_SIZE * rankDest;
		int &offset = _offsets[threadId][rankDest];
		if (offset + getMessageSize() > MESSAGE_BUFFER_SIZE)
			flushMessageBuffer(rankDest, tagDest, buffStart, offset);
		MessageClass *buf = (MessageClass *) (buffStart+offset);
		offset += getMessageSize();
		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg) {
		char *buf = (char*) bufferMessage(rankDest, tagDest);
		memcpy(buf, (char*) msg, getMessageSize());
	}

protected:
	void flushMessageBuffer(int rankDest, int tagDest, char *buffer, int &offset, bool checkReceive = true, bool sendZeroMessage = false) {
		int threadId = omp_get_thread_num();
		if (checkReceive)
			receiveAllIncomingMessages(threadId);
		if (offset > 0 || sendZeroMessage) {
			_world.send(rankDest, tagDest, buffer, offset);
		}
		offset = 0;
	}
public:
	void flushMessageBuffer(int rankDest, int tagDest, bool checkReceive = true, bool sendZeroMessage = false) {
		int threadId = omp_get_thread_num();
		int &offset = _offsets[threadId][rankDest];
		char *buffer = _sendBuffers[threadId] + MESSAGE_BUFFER_SIZE * rankDest;
		flushMessageBuffer(rankDest, tagDest, buffer, offset, checkReceive, sendZeroMessage);
	}
	void flushAllMessageBuffers(int tagDest, bool checkReceive = true, bool sendZeroMessage = false) {
		int threadId = omp_get_thread_num();
		if (checkReceive)
			receiveAllIncomingMessages(threadId);
		for(int i = 0 ; i < _world.size(); i++)
			flushMessageBuffer(i, tagDest, checkReceive, sendZeroMessage);
	}
	bool receiveIncomingMessage(int rankSource, int msgTag) {
		int threadId = omp_get_thread_num();
		assert(threadId == msgTag);
		OptionalStatus optionalStatus;
		OptionalRequest &optionalRequest = _requests[threadId][rankSource];
		if (!!optionalRequest) {
			optionalStatus = optionalRequest.get().test();
			if (!!optionalStatus) {
				int source = optionalStatus.get().source();
				int tag = optionalStatus.get().tag();
				int size = optionalStatus.get().count<char>().get();

				assert( rankSource == source );
				assert( threadId == tag );

				LOG_DEBUG(4, "receiving message from " << source << " size " << size << " t=" << threadId);

				char *buffer = _recvBuffers[ threadId ] + ( source * MESSAGE_BUFFER_SIZE );
				if (size == 0) {
					checkpoint();
					LOG_DEBUG(2, "got checkpoint from " << source);
					optionalRequest = OptionalRequest();
				} else {
					processMessages(buffer, size);
					optionalRequest = _world.irecv(source, threadId, buffer, MESSAGE_BUFFER_SIZE);
				}
				return true;
			}
		}
		return false;
	}
	int receiveAllIncomingMessages(int msgTag) {
		int threadId = omp_get_thread_num();
		assert(threadId == msgTag);
		int messages = 0;

		LOG_DEBUG(3, "receiving all messages for thread " << msgTag);
		bool mayHavePending = true;
		while (mayHavePending) {
			LOG_DEBUG(4, "calling iprobe");

			mayHavePending = false;
			for(int rankSource = 0 ; rankSource < _world.size(); rankSource++) {
				if (receiveIncomingMessage(rankSource, msgTag)) {
					messages++;
					mayHavePending = true;
				}
			}

		}
		LOG_DEBUG(3, "processed messages: " << messages);
		return messages;
	}
	void checkpoint() {
		#pragma omp atomic
		_numCheckpoints++;
		LOG_DEBUG(2, "checkpoint received:" << _numCheckpoints);
	}
	inline int getNumCheckpoints() const {
		return _numCheckpoints;
	}
};

#endif

#endif /* MPIBUFFER_H_ */
