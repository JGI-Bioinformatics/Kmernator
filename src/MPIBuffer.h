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

#define _RETRY_MESSAGES false
#define _RETRY_THRESHOLD 10000

#include "boost/optional.hpp"
#include <boost/thread/thread.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <vector>

#define MPI_BUFFER_DEFAULT_SIZE (256 * 1024)
#define WAIT_MS 10

class _MPIMessageBufferBase
{
public:
	typedef boost::optional< mpi::request > OptionalRequest;
	typedef boost::optional< mpi::status > OptionalStatus;
	typedef std::pair< _MPIMessageBufferBase*, int> CallbackBase;
	typedef std::vector< CallbackBase > CallbackVector;

protected:
	CallbackVector _flushAllCallbacks;
	CallbackVector _receiveAllCallbacks;
	long _deliveries;
	long _numMessages;

public:
	_MPIMessageBufferBase() : _deliveries(0), _numMessages(0) {}
	virtual ~_MPIMessageBufferBase() {}

	// receive buffers to flush before and/or during send
	void addReceiveAllCallback( _MPIMessageBufferBase *receiveAllBuffer ) {
		_receiveAllCallbacks.push_back( CallbackBase(receiveAllBuffer, -1) );
		receiveAll();
		assert(receiveAllBuffer->getNumDeliveries() == 0);
	}
	void addFlushAllCallback( _MPIMessageBufferBase *flushAllBuffer, int tagDest ) {
		_flushAllCallbacks.push_back( CallbackBase(flushAllBuffer, tagDest) );
		flushAll( );
		assert(flushAllBuffer->getNumDeliveries() == 0);
	}
	virtual long receiveAllIncomingMessages() { return 0; }
	virtual long flushAllMessageBuffers(int tagDest) { return 0; }

	long receiveAll() {
		long count = 0;
		for(unsigned int i = 0; i < _receiveAllCallbacks.size(); i++)
			count += _receiveAllCallbacks[i].first->receiveAllIncomingMessages();
		LOG_DEBUG(4, "receiveAll() with " << count);
		return count;
	}
	long flushAll() {
		long count = 0;
		for(unsigned int i = 0; i < _flushAllCallbacks.size(); i++)
			count += _flushAllCallbacks[i].first->flushAllMessageBuffers(_flushAllCallbacks[i].second);
		LOG_DEBUG(4, "flushAll() with count " << count);
		return count;
	}
	inline long getNumDeliveries() {
		return _deliveries;
	}
	inline void newMessageDelivery() {
		#pragma omp atomic
		_deliveries++;
	}
	inline long getNumMessages() {
		return _numMessages;
	}
	inline void newMessage() {
		#pragma omp atomic
		_numMessages++;
	}

};

template <typename C, int BufferSize = MPI_BUFFER_DEFAULT_SIZE>
class MPIMessageBufferBase : public _MPIMessageBufferBase {
public:
	static const int MESSAGE_BUFFER_SIZE = BufferSize;
	typedef C MessageClass;
	typedef _MPIMessageBufferBase::OptionalRequest OptionalRequest;
	typedef _MPIMessageBufferBase::OptionalStatus  OptionalStatus;
	typedef char * Buffer;
	typedef std::vector< Buffer > FreeBufferCache;

protected:
	mpi::communicator _world;
	int _messageSize;
	Buffer _message;
	FreeBufferCache freeBuffers;

public:
	MPIMessageBufferBase(mpi::communicator &world, int messageSize)
	: _world(world), _messageSize(messageSize) {
		_message = new char[ getMessageSize() ];
		assert(getMessageSize() >= (int) sizeof(MessageClass));
	}
	~MPIMessageBufferBase() {
		delete [] _message;
		for(FreeBufferCache::iterator it = freeBuffers.begin(); it != freeBuffers.end(); it++)
			delete [] *it;
		freeBuffers.clear();
	}
	inline mpi::communicator &getWorld() {
		return _world;
	}
	inline int getMessageSize() {
		return _messageSize;
	}
	MessageClass *getTmpMessage() {
		return (MessageClass*) this->_message;
	}
	Buffer getNewBuffer() {
		if (freeBuffers.empty()) {
			return new char[ MESSAGE_BUFFER_SIZE ];
		} else {
			Buffer buf;
			#pragma omp critical (mpi_buffer_base_free_buffers)
			{
				buf = freeBuffers.back();
				freeBuffers.pop_back();
			}
			return buf;
		}	
	}
	void returnBuffer(Buffer buf) {
		if (freeBuffers.size() >= (size_t) (_world.size() * 3) ) {
			delete [] buf;
		} else {
			#pragma omp critical (mpi_buffer_base_free_buffers)
			freeBuffers.push_back(buf);
		}
	}
};

template <typename C, int BufferSize = MPI_BUFFER_DEFAULT_SIZE>
class MPIRecvMessageBuffer: public MPIMessageBufferBase<C, BufferSize> {
public:
	typedef MPIMessageBufferBase<C, BufferSize> BufferBase;
	typedef typename BufferBase::OptionalRequest OptionalRequest;
	typedef typename BufferBase::OptionalStatus OptionalStatus;
	typedef typename BufferBase::Buffer Buffer;
	class BufferReceived {
	public:
		Buffer buffer;
		int size;
		int source;
		int tag;
		BufferReceived(Buffer b, int s, int src, int t) : buffer(b), size(s), source(src), tag(t) {}
	};
	typedef std::list<BufferReceived> ProcessBufferQueue;
	typedef C MessageClass;

protected:
	Buffer *_recvBuffers;
	OptionalRequest *_requests;
	int _numCheckpoints;
	int _tag;
	int _recvSource;
	int _recvTag;
	int _recvSize;
	std::vector<int> _requestAttempts;
	ProcessBufferQueue _processBufferQueue;
	bool _isProcessing;

public:
	MPIRecvMessageBuffer(mpi::communicator &world, int messageSize, int tag = mpi::any_tag) :
		BufferBase(world, messageSize), _numCheckpoints(0), _tag(tag), _recvSource(mpi::any_source), _recvTag(mpi::any_tag), _recvSize(0), _isProcessing(false) {
		_recvBuffers = new Buffer[ this->getWorld().size() ];
		_requests = new OptionalRequest[ this->getWorld().size() ];
		_requestAttempts.resize(this->getWorld().size(), 0);
		for(int destRank = 0; destRank < this->getWorld().size(); destRank++) {
			_recvBuffers[destRank] = NULL;
			_requests[destRank] = irecv(destRank);
		}
	}
	~MPIRecvMessageBuffer() {
		assert(_processBufferQueue.empty());
		cancelAllRequests();
		for(int i = 0 ; i < this->getWorld().size(); i++) {
			Buffer &buf = _recvBuffers[i];
			if (buf != NULL)
				delete [] buf;
			buf = NULL;
		}
		delete [] _recvBuffers;
		delete [] _requests;
	}
	Buffer &getBuffer(int rank) {
		Buffer &buf = _recvBuffers[rank];
		if (buf == NULL) 
			buf = this->getNewBuffer();
		return buf;
	}

	inline int getRecvSource() {
		return _recvSource;
	}
	inline int getRecvTag() {
		return _recvTag;
	}
	inline int getRecvSize() {
		return _recvSize;
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
			++_requestAttempts[rankSource];
			bool retry = _RETRY_MESSAGES && _requestAttempts[rankSource] > _RETRY_THRESHOLD;

			optionalStatus = optionalRequest.get().test();
			if (retry && !optionalStatus) {
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
				queueMessage(status, rankSource);
				optionalRequest = irecv(rankSource);
				returnValue = true;
			}
		}
		processQueue();
		return returnValue;
	}
	long receiveAllIncomingMessages() {
		long messages = this->receiveAll();

		LOG_DEBUG(4, _tag << ": receiving all messages for tag. cp: " << getNumCheckpoints() << " msgCount: " << this->getNumDeliveries());
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
			LOG_DEBUG(4, _tag << ": processed messages: " << messages);

		return messages;
	}
private:
	long _finalize(int checkpointFactor) {
		long messages;
		do {
			messages = 0;
			messages += processQueue();
			messages += this->receiveAllIncomingMessages();
			messages += this->flushAll();
			if (reachedCheckpoint(checkpointFactor) && messages != 0) {
				LOG_DEBUG_OPTIONAL(2, true, "Recv " << _tag << ": Achieved checkpoint but more messages are pending");
			}
			if (messages == 0)
				boost::this_thread::sleep( boost::posix_time::milliseconds(WAIT_MS) );
		} while (!reachedCheckpoint(checkpointFactor));
		return messages;
	}
public:
	void finalize(int checkpointFactor = 1) {
		LOG_DEBUG_OPTIONAL(1, true, "Recv " << _tag << ": Entering finalize checkpoint: " << getNumCheckpoints() << " out of " << (checkpointFactor*this->getWorld().size()));

		long globalMessages = 1;
		while(globalMessages != 0) {
			long messages = _finalize(checkpointFactor);
			LOG_DEBUG_OPTIONAL(2, true, "waiting for globalMessage to be 0: " << messages);
			all_reduce(this->getWorld(), messages, globalMessages, mpi::maximum<long>());
		}

		LOG_DEBUG_OPTIONAL(2, true, "Recv " << _tag << ": Finished finalize checkpoint: " << getNumCheckpoints());
		_numCheckpoints = 0;
	}

	void checkpoint() {
		#pragma omp atomic
		_numCheckpoints++;
		LOG_DEBUG(3, _tag << ": checkpoint received:" << _numCheckpoints);
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
		LOG_DEBUG(5, "Starting irecv for " << sourceRank << "," << _tag);
#ifdef OPENMP_CRITICAL_MPI
#pragma omp critical (MPI_buffer)
#endif
		oreq = this->getWorld().irecv(sourceRank, _tag, getBuffer(sourceRank), BufferBase::MESSAGE_BUFFER_SIZE);
		assert(!!oreq);
		_requestAttempts[sourceRank] = 0;
		return oreq;
	}
	int processMessages(BufferReceived bufferReceived) {
		Buffer msg, end, start = bufferReceived.buffer;
		this->_recvSize = bufferReceived.size;
		this->_recvSource = bufferReceived.source;
		this->_recvTag = bufferReceived.tag;

		msg = start;
		end = (start+bufferReceived.size);
		int count = 0;
		while (msg != end) {
			((MessageClass*) msg)->process(this);
			msg += this->getMessageSize();
			this->newMessage();
			count++;
		}
		return count;
	}
	bool queueMessage(mpi::status &status, int rankSource) {

		bool wasMessage = false;
		int source = status.source();
		int tag =  status.tag();
		int size = status.count<char>().get();
		if (status.cancelled()) {
			LOG_WARN(1, _tag << ": request was successfully canceled from " << source << " size " << size);
		} else {

			assert( rankSource == source );
			assert( _tag == tag );

			Buffer &bufLoc = getBuffer(source);
			Buffer buf = bufLoc;
			bufLoc = NULL;
			BufferReceived bufferReceived( buf, size, source, tag );

			#pragma omp critical (mpi_recv_buffer_process_queue)
			_processBufferQueue.push_back( bufferReceived );

			this->newMessageDelivery();
			LOG_DEBUG(4, _tag << ": received delivery " << this->getNumDeliveries() << " from " << source << "," << tag << " size " << size << " probe attempts: " << _requestAttempts[rankSource]);

			wasMessage = true;
		}
		return wasMessage;

	}
	long processQueue() {
		long processCount = 0;
		// code to allow recursive message generation
		if (!_isProcessing) {
			_isProcessing = true;

			while( !_processBufferQueue.empty() ) {
				BufferReceived bufferReceived = _processBufferQueue.front();

				#pragma omp critical (mpi_recv_buffer_process_queue)
				_processBufferQueue.pop_front();

				if (bufferReceived.size == 0) {
					checkpoint();
					LOG_DEBUG(3, _tag << ": got checkpoint from " << bufferReceived.source << "/" << bufferReceived.tag << ": " << _numCheckpoints);
				} else {
					processMessages(bufferReceived);
				}
				this->returnBuffer(bufferReceived.buffer);
				processCount++;
			}

			_isProcessing = false;
		}
		return processCount;
	}

};


template <typename C, int BufferSize = MPI_BUFFER_DEFAULT_SIZE>
class MPISendMessageBuffer: public MPIMessageBufferBase<C, BufferSize> {
public:
	typedef MPIMessageBufferBase<C, BufferSize> BufferBase;
	typedef MPIRecvMessageBuffer<C, BufferSize> RecvBuffer;
	typedef typename BufferBase::OptionalRequest OptionalRequest;
	typedef typename BufferBase::OptionalStatus OptionalStatus;
	typedef typename BufferBase::Buffer Buffer;
	class SentBuffer {
	public:
		mpi::request request;
		Buffer buffer;
		int destRank;
		int destTag;
		int size;
		long deliveryNum;
		long pollCount;

		SentBuffer(Buffer &_buffer, int _rank, int _tag, int _size, int _id) : request(), buffer(_buffer), destRank(_rank), destTag(_tag), size(_size), deliveryNum(_id), pollCount(0) {}
		SentBuffer() : request(), buffer(NULL), destRank(mpi::any_source), destTag(mpi::any_tag), size(0), deliveryNum(0), pollCount(0) {}
		SentBuffer(const SentBuffer &copy) {
			*this = copy;
		}
		SentBuffer &operator=(const SentBuffer &copy) {
			request = copy.request;
			buffer = copy.buffer;
			destRank = copy.destRank;
			destTag = copy.destTag;
			size = copy.size;
			deliveryNum = copy.deliveryNum;
			pollCount = copy.pollCount;
			return *this;
		}
		void reset() {
			*this = SentBuffer();
		}
		// predicate test for removal
		bool operator() (const SentBuffer &sent) {
			return (sent.buffer == NULL);
		}
	};
	typedef std::list< SentBuffer > SentBuffers;
	typedef typename SentBuffers::iterator SentBuffersIterator;
	typedef C MessageClass;
	typedef std::vector< _MPIMessageBufferBase* > RecvBufferCallbackVector;

protected:
	Buffer *_sendBuffers;
	int  *_offsets;
	SentBuffers _sentBuffers;

public:
	MPISendMessageBuffer(mpi::communicator &world, int messageSize) :
		BufferBase(world, messageSize) {
		_sendBuffers = new Buffer[ world.size() ];
		_offsets = new int[ this->getWorld().size() ];
		for(int destRank = 0; destRank < this->getWorld().size(); destRank++) {
			_sendBuffers[destRank] = NULL;
			_offsets[destRank] = 0;
		}
	}
	~MPISendMessageBuffer() {
		for(int i = 0 ; i < this->getWorld().size(); i++) {
			Buffer &buf = _sendBuffers[i];
			if (buf != NULL)
				delete [] buf;
			buf = NULL;
		}
		delete [] _sendBuffers;
		delete [] _offsets;
	}

	Buffer &getBuffer(int rank) {
		Buffer &buf = _sendBuffers[rank];
		if (buf == NULL) {
			buf = this->getNewBuffer();
			_offsets[rank] = 0;
		}
		return buf;
	}

	MessageClass *bufferMessage(int rankDest, int tagDest) {
		bool wasSent;
		long messages = 0;
		return bufferMessage(rankDest, tagDest, wasSent, messages);
	}
	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest, bool &wasSent, long &messages) {
		int &offset = _offsets[rankDest];
		wasSent = false;
		while (offset + this->getMessageSize() > BufferBase::MESSAGE_BUFFER_SIZE) {
			messages += flushMessageBuffer(rankDest, tagDest);
			wasSent = true;
		}
		Buffer &buffStart = getBuffer(rankDest);
		assert(buffStart != NULL);
		MessageClass *buf = (MessageClass *) (buffStart+offset);
		offset += this->getMessageSize();
		this->newMessage();

		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg) {
		char *buf = (char *) bufferMessage(rankDest, tagDest);
		memcpy(buf, (char *) msg, this->getMessageSize());
	}

	long receiveAllIncomingMessages() {
		return this->receiveAll();
	}

public:
	long flushMessageBuffer(int rankDest, int tagDest, bool sendZeroMessage = false) {
		long messages = 0;
		long newMessages = 0;

		int &offsetLocation = _offsets[rankDest];
		bool waitForSend = sendZeroMessage == true;

		messages += checkSentBuffers( waitForSend );
		while (sendZeroMessage && (offsetLocation > 0 || newMessages != 0) ) {
			// flush all messages until there is nothing in the buffer
			newMessages = 0;
			boost::this_thread::sleep( boost::posix_time::milliseconds(WAIT_MS) );
			newMessages = flushMessageBuffer(rankDest, tagDest, false);
			newMessages += checkSentBuffers( waitForSend );
			messages += newMessages;
		}

		if (offsetLocation > 0 || sendZeroMessage) {
			Buffer &bufferLocation = getBuffer(rankDest);

			Buffer buffer = bufferLocation;
			bufferLocation = NULL;
			int offset = offsetLocation;
			offsetLocation = 0;

			this->newMessageDelivery();
			SentBuffer sent(buffer, rankDest, tagDest, offset, this->getNumDeliveries());
			LOG_DEBUG(3, "sending message to " << rankDest << ", " << tagDest << " size " << offset);
#ifdef OPENMP_CRITICAL_MPI
#pragma omp critical (MPI_buffer)
#endif
			sent.request = this->getWorld().isend(rankDest, tagDest, buffer, offset);

			messages++;

			recordSentBuffer(sent);

		}
		messages += checkSentBuffers( waitForSend );
		return messages;
	}

	void recordSentBuffer(SentBuffer &sent) {
		#pragma omp critical (mpi_send_sent_buffers)
		_sentBuffers.push_back(sent);
	}
	long checkSentBuffers(bool wait = false) {
		long messages = 0;
		long iterations = 0;
		while (wait || iterations++ == 0) {
			messages += this->receiveAllIncomingMessages();
			for(SentBuffersIterator it = _sentBuffers.begin(); it != _sentBuffers.end(); it++) {
				SentBuffer &sent = *it;

				OptionalStatus optStatus = sent.request.test();

				if ( _RETRY_MESSAGES && (!optStatus) && ++sent.pollCount > _RETRY_THRESHOLD) {
					sent.request.cancel();
					if ( ! sent.request.test() ) {
						sent.request = this->getWorld().isend(sent.destRank, sent.destTag, sent.buffer, sent.size);
						sent.pollCount = 0;
						LOG_WARN(1, "Canceled and retried pending message to " << sent.destRank << ", " << sent.destTag << " size " << sent.size << " deliveryCount " << sent.deliveryNum);
					}
					optStatus = sent.request.test();
				}

				if (!!optStatus) {

					mpi::status status = optStatus.get();
					// hmmm looks like error is sometimes populated with junk...
					//if (status.error() > 0)
					//	LOG_WARN(1, "sending message returned an error: " << status.error());
					LOG_DEBUG(4, "finished sending message to " << sent.destRank << ", " << sent.destTag << " size " << sent.size << " deliveryCount " << sent.deliveryNum);

					this->returnBuffer(sent.buffer);
					sent.reset();

				}
			}
			#pragma omp critical (mpi_send_sent_buffers)
			_sentBuffers.remove_if(SentBuffer());

			if (wait) {
				boost::this_thread::sleep( boost::posix_time::milliseconds(WAIT_MS) );
				wait = !_sentBuffers.empty();
			}
		}
		return messages;
	}
	long flushAllMessageBuffers(int tagDest) {
		return flushAllMessageBuffers(tagDest, false);
	}
	long flushAllMessageBuffers(int tagDest, bool sendZeroMessage) {
		long messages = 0;
		for(int rankDest = 0 ; rankDest < this->getWorld().size(); rankDest++)
			messages += flushMessageBuffer(rankDest, tagDest, sendZeroMessage);
		messages += this->flushAll();
		return messages;
	}
	void flushAllMessagesUntilEmpty(int tagDest) {
		long messages = flushAllMessageBuffers(tagDest);
		while (messages != 0) {
			boost::this_thread::sleep( boost::posix_time::milliseconds(WAIT_MS) );
			messages = flushAllMessageBuffers(tagDest);
			messages += checkSentBuffers(true);
		}

	}
	void finalize(int tagDest) {
		LOG_DEBUG_OPTIONAL(1, true, "Send " << tagDest << ": entering finalize()");

		// first clear the buffer and in-flight messages
		flushAllMessagesUntilEmpty(tagDest);
		receiveAllIncomingMessages();

		LOG_DEBUG_OPTIONAL(3, true,"Send " << tagDest << ": entering finalize stage2()");
		// send zero-message as checkpoint signal to stop
		flushAllMessageBuffers(tagDest, true);
		boost::this_thread::sleep( boost::posix_time::milliseconds(2) );

		LOG_DEBUG_OPTIONAL(3, true,"Send " << tagDest << ": entering finalize stage3()");
		// now continue to flush until there is nothing left in the buffer;
		flushAllMessagesUntilEmpty(tagDest);

		LOG_DEBUG_OPTIONAL(2, true,"Send " << tagDest << ": finished finalize()");

	}
	long processPending() {
		return checkSentBuffers(false);
	}
};

#endif /* MPIBUFFER_H_ */
