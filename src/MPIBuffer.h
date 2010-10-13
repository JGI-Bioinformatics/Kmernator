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

#include "boost/optional.hpp"
#include <boost/thread/thread.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

template <typename C, int BufferSize = 16 * 1024>
class MPIMessageBuffers {
public:
	static const int MESSAGE_BUFFER_SIZE = BufferSize;
	typedef C MessageClass;
	typedef boost::optional< mpi::request > OptionalRequest;
	typedef boost::optional< mpi::status > OptionalStatus;

protected:
	mpi::communicator &_world;
	int _tag;
	int _messageSize;
	int _numCheckpoints;
	char *_message;
	char *_recvBuffers;
	char *_sendBuffers;
	int  *_offsets;
	OptionalRequest *_requests;

public:
	MPIMessageBuffers(mpi::communicator &world, int messageSize, int tag = 0) : _world(world), _messageSize(messageSize), _tag(tag), _numCheckpoints(0) {
		_message = new char[ _messageSize ];
		_recvBuffers = new char[ MESSAGE_BUFFER_SIZE * world.size() ];
		_sendBuffers = new char[ MESSAGE_BUFFER_SIZE * world.size() ];
		_offsets = new int[ _world.size() ];
		_requests = new OptionalRequest[ _world.size() ];
		for(int destRank = 0; destRank < _world.size(); destRank++) {
			_offsets[destRank] = 0;
			_requests[destRank] = _world.irecv(destRank, _tag, _recvBuffers + destRank*MESSAGE_BUFFER_SIZE, MESSAGE_BUFFER_SIZE);
		}
	}
	~MPIMessageBuffers() {
		delete [] _message;
		delete [] _recvBuffers;
		delete [] _sendBuffers;
		delete [] _offsets;
		delete [] _requests;
	}
	MessageClass *getTmpMessage() {
		return (MessageClass*) _message;
	}
	int processMessages(char *start, int size) {
		char *msg, *end;
		msg = start;
		end = (start+size);
		int count = 0;
		while (msg != end) {
			_processMessage((MessageClass*) msg);
			msg += _messageSize;
		}
		return count;
	}
	virtual void _processMessage(MessageClass *msg) = 0;

	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest) {
		char *buffStart = _sendBuffers + MESSAGE_BUFFER_SIZE * rankDest;
		int &offset = _offsets[rankDest];
		if (offset + _messageSize > MESSAGE_BUFFER_SIZE)
			flushMessageBuffer(rankDest, tagDest, buffStart, offset);
		MessageClass *buf = (MessageClass *) (buffStart+offset);
		offset += _messageSize;
		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg) {
		char *buf = (char*) bufferMessage(rankDest, tagDest);
		memcpy(buf, (char*) msg, _messageSize);
	}

protected:
	void flushMessageBuffer(int rankDest, int tagDest, char *buffer, int &offset, bool checkReceive = true, bool sendZeroMessage = false) {
		if (checkReceive)
			receiveAllIncomingMessages();
		if (offset > 0 || sendZeroMessage) {
			_world.send(rankDest, tagDest, buffer, offset);
		}
		offset = 0;
	}
public:
	void flushMessageBuffer(int rankDest, int tagDest, bool checkReceive = true, bool sendZeroMessage = false) {
		int &offset = _offsets[rankDest];
		char *buffer = _sendBuffers + MESSAGE_BUFFER_SIZE * rankDest;
		flushMessageBuffer(rankDest, tagDest, buffer, offset, checkReceive, sendZeroMessage);
	}
	void flushAllMessageBuffers(int tagDest, bool checkReceive = true, bool sendZeroMessage = false) {
		if (checkReceive)
			receiveAllIncomingMessages();
		for(int i = 0 ; i < _world.size(); i++)
			flushMessageBuffer(i, tagDest, checkReceive, sendZeroMessage);
	}
	bool receiveIncomingMessage(int rankSource) {
		OptionalStatus optionalStatus;
		OptionalRequest &optionalRequest = _requests[rankSource];
		if (!!optionalRequest) {
			optionalStatus = optionalRequest.get().test();
			if (!!optionalStatus) {
				int source = optionalStatus.get().source();
				int tag = optionalStatus.get().tag();
				int size = optionalStatus.get().count<char>().get();

				assert( rankSource == source );
				assert( _tag == tag );

				LOG_DEBUG_MT(4, _world.rank() << ": " << _tag << ": receiving message from " << source << " size " << size);

				char *buffer = _recvBuffers + ( source * MESSAGE_BUFFER_SIZE );
				if (size == 0) {
					checkpoint();
					LOG_DEBUG_MT(2,_world.rank() << ": " << _tag << ": got checkpoint from " << source);
					optionalRequest = OptionalRequest();
				} else {
					processMessages(buffer, size);
					optionalRequest = _world.irecv(source, _tag, buffer, MESSAGE_BUFFER_SIZE);
				}
				return true;
			}
		}
		return false;
	}
	int receiveAllIncomingMessages() {
		int messages = 0;

		LOG_DEBUG_MT(3, _world.rank() << ": " << _tag << ": receiving all messages for thread ");
		bool mayHavePending = true;
		while (mayHavePending) {
			LOG_DEBUG_MT(4, _world.rank() << ": " << _tag << ": calling iprobe");

			mayHavePending = false;
			for(int rankSource = 0 ; rankSource < _world.size(); rankSource++) {
				if (receiveIncomingMessage(rankSource)) {
					messages++;
					mayHavePending = true;
				}
			}

		}
		LOG_DEBUG_MT(3, _world.rank() << ": " << _tag << ": processed messages: " << messages);
		return messages;
	}
	void checkpoint() {
		#pragma omp atomic
		_numCheckpoints++;
		LOG_DEBUG_MT(2, _world.rank() << ": " << _tag << ": checkpoint received:" << _numCheckpoints);
	}
	inline int getNumCheckpoints() const {
		return _numCheckpoints;
	}
	bool reachedCheckpoint() const {
		return _numCheckpoints == _world.rank();
	}
	void finalize() {
		flushAllMessageBuffers(true);
		// send zero message buffer as checkpoint signal to stop
		LOG_DEBUG_MT(2, _world.rank() << ": " << _tag << ": sending stop message");
		flushAllMessageBuffers(true, true);
		checkpoint();
		LOG_DEBUG_MT(2, _world.rank() << ": " << _tag << ": sent checkpoints: " << getNumCheckpoints());
		while ( ! reachedCheckpoint() ) {
			receiveAllIncomingMessages();
			boost::this_thread::sleep( boost::posix_time::milliseconds(2) );
		}

	}
};



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
			msg += _messageSize;
		}
		return count;
	}
	virtual void _processMessage(MessageClass *msg) = 0;

	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest) {
		int threadId = omp_get_thread_num();
		char *buffStart = _sendBuffers[threadId] + MESSAGE_BUFFER_SIZE * rankDest;
		int &offset = _offsets[threadId][rankDest];
		if (offset + _messageSize > MESSAGE_BUFFER_SIZE)
			flushMessageBuffer(rankDest, tagDest, buffStart, offset);
		MessageClass *buf = (MessageClass *) (buffStart+offset);
		offset += _messageSize;
		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg) {
		char *buf = (char*) bufferMessage(rankDest, tagDest);
		memcpy(buf, (char*) msg, _messageSize);
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

				LOG_DEBUG_MT(4, _world.rank() << ": " << threadId << ": receiving message from " << source << " size " << size << " t=" << threadId);

				char *buffer = _recvBuffers[ threadId ] + ( source * MESSAGE_BUFFER_SIZE );
				if (size == 0) {
					checkpoint();
					LOG_DEBUG_MT(2,_world.rank() << ": " << threadId << ": got checkpoint from " << source);
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

		LOG_DEBUG_MT(3, _world.rank() << ": " << threadId << ": receiving all messages for thread " << msgTag);
		bool mayHavePending = true;
		while (mayHavePending) {
			LOG_DEBUG_MT(4, _world.rank() << ": " << threadId << ": calling iprobe");

			mayHavePending = false;
			for(int rankSource = 0 ; rankSource < _world.size(); rankSource++) {
				if (receiveIncomingMessage(rankSource, msgTag)) {
					messages++;
					mayHavePending = true;
				}
			}

		}
		LOG_DEBUG_MT(3, _world.rank() << ": " << threadId << ": processed messages: " << messages);
		return messages;
	}
	void checkpoint() {
		#pragma omp atomic
		_numCheckpoints++;
		LOG_DEBUG_MT(2, _world.rank() << ": " << omp_get_thread_num() << ": checkpoint received:" << _numCheckpoints);
	}
	inline int getNumCheckpoints() const {
		return _numCheckpoints;
	}
};


#endif /* MPIBUFFER_H_ */
