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

#include "config.h"
#include "MPIBase.h"

#ifdef _USE_OPENMP
#define X_OPENMP_CRITICAL_MPI
#endif

#define _RETRY_MESSAGES false
#define _RETRY_THRESHOLD 10000

#include <vector>

#define MPI_BUFFER_DEFAULT_SIZE (512 * 1024)
#define WAIT_MS 1
#define WAIT_AND_WARN( iterations, warningMessage ) \
	if ((iterations % (60000/WAIT_MS)) == 0) LOG_WARN(1, warningMessage << " waiting in loop: " << iterations);  \
	boost::this_thread::sleep( boost::posix_time::milliseconds(WAIT_MS) );

class MPIMessageBufferBase
{
public:
	typedef std::pair< MPIMessageBufferBase*, int> CallbackBase;
	typedef std::vector< CallbackBase > CallbackVector;

protected:
	CallbackVector _flushAllCallbacks;
	CallbackVector _receiveAllCallbacks;
	CallbackVector _sendReceiveCallbacks;
	long _deliveries;
	long _numMessages;

public:
	MPIMessageBufferBase() : _deliveries(0), _numMessages(0) {}
	virtual ~MPIMessageBufferBase() {}
	static bool isThreadSafe() {
		int threadlevel;
		MPI_Query_thread(&threadlevel);
		return threadlevel >= MPI_THREAD_FUNNELED;
	}

	// receive buffers to flush before and/or during send
	void addReceiveAllCallback( MPIMessageBufferBase *receiveAllBuffer ) {
		_receiveAllCallbacks.push_back( CallbackBase(receiveAllBuffer, -1) );
		receiveAll();
		assert(receiveAllBuffer->getNumDeliveries() == 0);
	}
	void addFlushAllCallback( MPIMessageBufferBase *flushAllBuffer, int tagDest ) {
		_flushAllCallbacks.push_back( CallbackBase(flushAllBuffer, tagDest) );
		flushAll( );
		assert(flushAllBuffer->getNumDeliveries() == 0);
	}
	void addSendReceiveCallback( MPIMessageBufferBase *sendReceiveBuffer ) {
		_sendReceiveCallbacks.push_back( CallbackBase(sendReceiveBuffer, -1) );
	}
	virtual long receiveAllIncomingMessages(bool untilFlushed = true) { return 0; }
	virtual long flushAllMessageBuffers(int tagDest) { return 0; }
	virtual long sendReceive(bool isFinalized) { return 0; }

	long receiveAll(bool untilFlushed = true) {
		long count = 0;
		for(unsigned int i = 0; i < _receiveAllCallbacks.size(); i++)
			count += _receiveAllCallbacks[i].first->receiveAllIncomingMessages(untilFlushed);
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
	long sendReceiveAll(bool isFinalized) {
		long count = 0;
		for(unsigned int i = 0; i < _sendReceiveCallbacks.size(); i++)
			count += _sendReceiveCallbacks[i].first->sendReceive(isFinalized);
		LOG_DEBUG(4, "sendReceiveAll() with count " << count);
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

template <typename C, typename B = MPIMessageBufferBase>
class DummyProcessor {
public:
	int process(C *msg, B *buffer) { return 0; }
};
template <typename C, typename CProcessor = DummyProcessor<C>, int BufferSize = MPI_BUFFER_DEFAULT_SIZE >
class MPIMessageBuffer : public MPIMessageBufferBase {
public:
	static const int MESSAGE_BUFFER_SIZE = BufferSize;
	static const int BUFFER_QUEUE_SOFT_LIMIT = 5;

	typedef C MessageClass;
	typedef CProcessor MessageClassProcessor;
	typedef char * Buffer;
	typedef std::vector< Buffer > FreeBufferCache;
	typedef std::vector< FreeBufferCache > ThreadedFreeBufferCache;
	class MessagePackage {
	public:
		Buffer buffer;
		int size;
		int source;
		int tag;
		MessagePackage(Buffer b, int s, int src, int t) : buffer(b), size(s), source(src), tag(t) {}
	};
	typedef std::list<MessagePackage> MessagePackageQueue;

protected:
	mpi::communicator _world;
	int _messageSize;
	MessageClassProcessor _processor;
	int _softMaxBufferSize;
	int _numCheckpoints;
	Buffer _message;
	ThreadedFreeBufferCache _freeBuffers;
	int _bufferSize;
	int _numThreads;
	int _worldSize;
	int _numWorldThreads;
	int _recvSource;
	int _recvTag;
	int _recvSize;

public:
	MPIMessageBuffer(mpi::communicator &world, int messageSize, MessageClassProcessor processor = MessageClassProcessor())
	: _world(world), _messageSize(messageSize), _processor(processor), _softMaxBufferSize(MESSAGE_BUFFER_SIZE), _numCheckpoints(0), _bufferSize(MESSAGE_BUFFER_SIZE), _recvSource(mpi::any_source), _recvTag(mpi::any_tag), _recvSize(0)  {
		_message = new char[ MESSAGE_BUFFER_SIZE ];
		assert(getMessageSize() >= (int) sizeof(MessageClass));
		if (omp_in_parallel())
			_numThreads = omp_get_num_threads();
		else
			_numThreads = omp_get_max_threads();
		_worldSize = _world.size();
		_numWorldThreads = _worldSize * _numThreads;
		_freeBuffers.resize(_numThreads);
		for(int threadId = 0; threadId < _numThreads; threadId++)
			_freeBuffers[threadId].reserve(BUFFER_QUEUE_SOFT_LIMIT * _numWorldThreads + 1);
		setSoftMaxBufferSize( Options::getBatchSize() * messageSize );
	}
	~MPIMessageBuffer() {
		delete [] _message;
		for(ThreadedFreeBufferCache::iterator it = _freeBuffers.begin(); it != _freeBuffers.end(); it++)
			for(FreeBufferCache::iterator it2 = it->begin(); it2 != it->end() ; it2++)
				delete [] *it2;
		_freeBuffers.clear();
	}
	inline mpi::communicator &getWorld() {
		return _world;
	}
	inline int getNumThreads() {
		return _numThreads;
	}
	inline int getNumWorldThreads() {
		return _numWorldThreads;
	}
	inline int getWorldSize () {
		return _worldSize;
	}
	inline int getMessageSize() {
		return _messageSize;
	}
	inline int &getRecvSource() {
		return _recvSource;
	}
	inline int &getRecvTag() {
		return _recvTag;
	}
	inline int &getRecvSize() {
		return _recvSize;
	}
	inline MessageClassProcessor &getProcessor() {
		return _processor;
	}
	void setMessageProcessor(MessageClassProcessor processor) {
		_processor = processor;
	}
	inline int getMaxBufferSize() const {
		return MESSAGE_BUFFER_SIZE;
	}
	inline int getSoftMaxBufferSize() const {
		return _softMaxBufferSize;
	}
	void setSoftMaxBufferSize(int size) {
		_softMaxBufferSize = std::max(getMessageSize(), std::min(size, getMaxBufferSize()));
	}
	MessageClass *getTmpMessage() {
		return (MessageClass*) this->_message;
	}
	void setBufferSize(int bufferSize) {
		_bufferSize = bufferSize;
	}
	Buffer getNewBuffer(int numLinearBuffers = 1, int padding = 0) {
		int threadId = omp_get_thread_num();
		if (_freeBuffers[threadId].empty()) {
			return new char[ _bufferSize ];
		} else {
			Buffer buf;
			buf = _freeBuffers[threadId].back();
			_freeBuffers[threadId].pop_back();
			return buf;
		}	
	}
	void returnBuffer(Buffer buf) {
		int threadId = omp_get_thread_num();
		if (_freeBuffers[threadId].size() >= (size_t) (BUFFER_QUEUE_SOFT_LIMIT * _numWorldThreads) ) {
			delete [] buf;
		} else {
			_freeBuffers[threadId].push_back(buf);
		}
	}

	long processMessagePackage(MessagePackage MessagePackage) {
		Buffer msg, start = MessagePackage.buffer;
		this->_recvSize = MessagePackage.size;
		this->_recvSource = MessagePackage.source;
		this->_recvTag = MessagePackage.tag;

		long count = 0;
		int offset = 0;
		while (offset < MessagePackage.size) {
			msg = start + offset;
			int trailingBytes = _processor.process((MessageClass*) msg, this);
			offset += this->getMessageSize() + trailingBytes;
			this->newMessage();
			count++;
		}
		return count;
	}

	inline int getNumCheckpoints() {
		return _numCheckpoints;
	}
	bool reachedCheckpoint(int checkpointFactor = 1) {
		return getNumCheckpoints() == this->getWorld().size() * checkpointFactor;
	}
	void resetCheckpoints() {
		_numCheckpoints = 0;
	}
	void checkpoint() {
		#pragma omp atomic
		_numCheckpoints++;
		LOG_DEBUG_OPTIONAL(3, true, "checkpoint received:" << _numCheckpoints);
	}

};

template <typename C, typename CProcessor, int BufferSize = MPI_BUFFER_DEFAULT_SIZE >
class MPIAllToAllMessageBuffer : public MPIMessageBuffer<C, CProcessor, BufferSize> {
public:
	typedef MPIMessageBuffer<C, CProcessor, BufferSize> BufferBase;
	typedef typename BufferBase::Buffer Buffer;
	typedef typename BufferBase::MessagePackage MessagePackage;
	typedef typename BufferBase::MessagePackageQueue MessagePackageQueue;
	typedef C MessageClass;
	typedef CProcessor MessageClassProcessor;

	class MessageHeader {
	public:
		int offset, threadSource, tag, dummy;
		MessageHeader() : offset(0), threadSource(0), tag(0) {}
		MessageHeader(const MessageHeader &copy) : offset(copy.offset) , threadSource(copy.threadSource), tag(copy.tag) {}
		MessageHeader &operator=(const MessageHeader &other) {
			if (this == &other)
				return *this;
			offset = other.offset ; threadSource = other.threadSource; tag = other.tag;
			return *this;
		}
		void reset(int _threadSource = 0, int _tag = 0) {
			offset = 0;
			threadSource = _threadSource;
			tag = _tag;
		}

	};
	class BuildBuffer
	{
	public:
		Buffer buffer;
		MessageHeader header;
		BuildBuffer() : buffer(NULL), header(MessageHeader()) {}
		~BuildBuffer() {
			reset();
		}
		void reset() {
			header.reset();
			if (buffer != NULL)
				delete [] buffer;
		}
	};
	static const int headerSize = sizeof(MessageHeader);

private:
	Buffer out,in;
	std::vector< std::vector< std::vector<BuildBuffer> > > buildsTWT;
	int threadsSending;
	int buildSize;
	int inoutDataSize;

public:
	MPIAllToAllMessageBuffer(mpi::communicator &world, int messageSize, MessageClassProcessor processor = MessageClassProcessor())
	: BufferBase(world, messageSize, processor), threadsSending(0) {
		assert(!omp_in_parallel());
		int worldSize = this->getWorldSize();
		int numThreads = this->getNumThreads();

		inoutDataSize = sizeof(int) + numThreads * numThreads * (headerSize + BufferSize);
		int inoutSize = worldSize * sizeof(int) * 2 + worldSize * inoutDataSize;
		buildsTWT.resize(numThreads);
		for(int threadId = 0 ; threadId < numThreads; threadId++) {
			buildsTWT[threadId].resize(worldSize);
		}
		out =   new char[ inoutSize ];
		in =    new char[ inoutSize ];
		for(int rankDest = 0 ; rankDest < worldSize ; rankDest++) {
			for(int threadId = 0 ; threadId < numThreads; threadId++) {
				buildsTWT[threadId][rankDest].resize(numThreads, BuildBuffer());
				for(int threadDest = 0; threadDest < numThreads; threadDest++) {
					BuildBuffer &bb = buildsTWT[threadId][rankDest][threadDest];
					bb.header.reset(threadId,threadDest);
					bb.buffer = this->getNewBuffer();
				}
			}
			// set in/out
			getInOutSize(rankDest, out) = 0;
			getInOutSize(rankDest, in) = inoutDataSize;
			// in/out displacements do not change
			getInOutOffset(rankDest, out) = rankDest * inoutDataSize;
			getInOutOffset(rankDest, in)  = rankDest * inoutDataSize;
		}
		buildSize = 0;
	}
	~MPIAllToAllMessageBuffer() {
		assert(!omp_in_parallel());
		delete [] in;
		delete [] out;
	}
	inline int getJump(int rank, int thread, int threadId = 0) {
		return threadId * this->getNumWorldThreads() + rank * this->getNumThreads() + thread;
	}

	inline int &getInOutSize(int rankDest, Buffer buf) {
		int *o = (int*) buf;
		return *(o + rankDest);
	}
	inline int &getInOutOffset(int rankDest, Buffer buf) {
		int *o = (int*) buf;
		return *(o + this->getWorldSize() + rankDest);
	}
	inline int &getInOutDataSize(int rankDest, Buffer buf) {
		int *o = (int*) buf;
		buf = (Buffer) (o + this->getWorldSize() * 2);
		return *((int*) (buf + rankDest * inoutDataSize));
	}
	inline Buffer getInOutBuffer(int rankDest, Buffer buf) {
		return (Buffer) (&getInOutDataSize(rankDest, buf)+1);
	}
	long sendReceive(bool isFinalized) {
		assert(omp_get_max_threads() == 1 || omp_in_parallel());
		int threadId = omp_get_thread_num();
		int numThreads = this->getNumThreads();
		int worldSize = this->getWorldSize();
		int numWorldThreads = this->getNumWorldThreads();

		Buffer outBuffers[numThreads * numWorldThreads];
		long numReceived = 0;

#pragma omp atomic
		threadsSending++;

#pragma omp critical
		{
			// allocate the out message headers
			for(int rankDest = 0 ; rankDest < worldSize ; rankDest++) {
				for(int threadDest = 0; threadDest < numThreads; threadDest++) {
					BuildBuffer &bb = buildsTWT[threadId][rankDest][threadDest];
					int &size = bb.header.offset;

					int jump = getJump(rankDest, threadDest, threadId);
					int &outSize = getInOutSize(rankDest, out);

					outBuffers[jump] = getInOutBuffer(rankDest, out) + outSize;
					buildSize += size;
					outSize += size + headerSize;
				}
			}
		}
		// copy any headers and messages to out
		for(int rankDest = 0 ; rankDest < worldSize ; rankDest++) {
			for(int threadDest = 0; threadDest < numThreads; threadDest++) {
				int jump = getJump(rankDest, threadDest, threadId);
				BuildBuffer &bb = buildsTWT[threadId][rankDest][threadDest];

				MessageHeader *outHeader = (MessageHeader *) outBuffers[jump];
				*outHeader = bb.header;
				if (bb.header.offset > 0)
					memcpy(outHeader + 1, bb.buffer, bb.header.offset);

				// reset the counter
				bb.header.offset = 0;
			}
		}

#pragma omp barrier

#pragma omp master
		{
			if (isFinalized && buildSize == 0) {
				for(int destRank = 0; destRank < worldSize ; destRank++) {
					getInOutDataSize(destRank, out) = 0;
					getInOutSize(destRank,out) = sizeof(int);
				}
			} else {
				// set transmitted sizes
				for(int destRank = 0; destRank < worldSize ; destRank++) {
					getInOutDataSize(destRank, out) = getInOutSize(destRank,out);
					getInOutSize(destRank,out) += sizeof(int);
				}
			}
			// mpi_alltoall
			MPI_Alltoallv(out+worldSize*2*sizeof(int), &getInOutSize(0,out), &getInOutOffset(0, out), MPI_BYTE, in+worldSize*2*sizeof(int), &getInOutSize(0,in), &getInOutOffset(0,in), MPI_BYTE, this->getWorld());

			// reset out sizes
			for(int destRank = 0; destRank < worldSize ; destRank++)
				getInOutSize(destRank,out) = 0;
			buildSize = 0;
			this->resetCheckpoints();
			threadsSending = 0;
		}
#pragma omp barrier

		for(int sourceRank = 0 ; sourceRank < worldSize ; sourceRank++) {
			int &size = getInOutDataSize(sourceRank, in);
			assert(size >= 0);
			if (size == 0) {
				this->checkpoint();
				continue;
			}
			Buffer buf = getInOutBuffer(sourceRank, in);
			Buffer begin = buf;
			Buffer end = buf + size;			while (begin != end) {
				assert(begin < end);
				MessageHeader *header = (MessageHeader*) begin;
				if (threadId == header->tag) {
					if (header->offset > 0)
						numReceived += processMessages(sourceRank, header);
				}
				begin = ((Buffer)(header)) + headerSize + header->offset;
			}
		}
#pragma omp barrier

		numReceived += this->sendReceiveAll(isFinalized);

		return numReceived;
	}

	int processMessages(int sourceRank, MessageHeader *header) {
		Buffer buf = (Buffer) (header+1);
		MessagePackage msgPkg(buf, header->offset, sourceRank, header->tag);
		return this->processMessagePackage(msgPkg);
	}

	void finalize() {
		LOG_DEBUG(2, "Entering finalize()");
		while (!this->reachedCheckpoint(this->getNumThreads())) {
			sendReceive(true);
		}
	}
	bool isReadyToSend(int offset, int trailingBytes) {
		return offset > 0 && (offset + trailingBytes + this->getMessageSize() >= this->getSoftMaxBufferSize() || threadsSending > this->getNumThreads() - 1);
	}

	MessageClass *bufferMessage(int rankDest, int tagDest, int trailingBytes = 0) {
		bool wasSent;
		long messages = 0;
		return bufferMessage(rankDest, tagDest, wasSent, messages, trailingBytes);
	}
	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest, bool &wasSent, long &messages, int trailingBytes = 0) {
		int threadId = omp_get_thread_num();
		BuildBuffer &bb = buildsTWT[threadId][rankDest][tagDest];
		int &offset = bb.header.offset;
		if (isReadyToSend(offset, trailingBytes))
			messages += sendReceive(false);

		assert(bb.buffer != NULL);
		MessageClass *buf = (MessageClass *) (bb.buffer+offset);
		offset += this->getMessageSize() + trailingBytes;
		this->newMessage();

		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg, int trailingBytes = 0) {
		char *buf = (char *) bufferMessage(rankDest, tagDest, trailingBytes);
		memcpy(buf, (char *) msg, this->getMessageSize() + trailingBytes);
	}


};

template <typename C, typename CProcessor, int BufferSize = MPI_BUFFER_DEFAULT_SIZE>
class MPIRecvMessageBuffer: public MPIMessageBuffer<C, CProcessor, BufferSize> {
public:
	typedef MPIMessageBuffer<C, CProcessor, BufferSize> BufferBase;

	typedef typename BufferBase::Buffer Buffer;
	typedef typename BufferBase::MessagePackage MessagePackage;
	typedef typename BufferBase::MessagePackageQueue MessagePackageQueue;
	typedef C MessageClass;
	typedef CProcessor MessageClassProcessor;

protected:
	Buffer *_recvBuffers;
	MPIOptionalRequest *_requests;
	int _tag;
	std::vector<int> _requestAttempts;
	MessagePackageQueue _MessagePackageQueue;
	bool _isProcessing;

public:
	MPIRecvMessageBuffer(mpi::communicator &world, int messageSize, int tag = mpi::any_tag, MessageClassProcessor processor = MessageClassProcessor())
	: BufferBase(world, messageSize, processor), _tag(tag), _isProcessing(false) {
		_recvBuffers = new Buffer[ this->getWorld().size() ];
		_requests = new MPIOptionalRequest[ this->getWorld().size() ];
		_requestAttempts.resize(this->getWorld().size(), 0);
		for(int destRank = 0; destRank < this->getWorld().size(); destRank++) {
			_recvBuffers[destRank] = NULL;
			_requests[destRank] = irecv(destRank);
		}
	}
	~MPIRecvMessageBuffer() {
		assert(_MessagePackageQueue.empty());
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
	void cancelAllRequests() {
		for(int source = 0 ; source < this->getWorld().size() ; source++) {
			if (!!_requests[source]) {
				_requests[source].get().cancel();
				_requests[source] = MPIOptionalRequest();
			}
		}
	}
	bool receiveIncomingMessage(int rankSource) {
		bool returnValue = false;
		MPIOptionalStatus MPIOptionalStatus;
		MPIOptionalRequest &MPIOptionalRequest = _requests[rankSource];

		if (!MPIOptionalRequest) {
			LOG_WARN(1, "Detected non-pending request for tag: " << _tag);
			MPIOptionalRequest = irecv(rankSource);
		}
		if (!!MPIOptionalRequest) {
			++_requestAttempts[rankSource];
			bool retry = _RETRY_MESSAGES && _requestAttempts[rankSource] > _RETRY_THRESHOLD;
#ifdef OPENMP_CRITICAL_MPI
#pragma omp critical (MPIBUFFER_RECV_REQUEST_TEST)
#endif
			{
			MPIOptionalStatus = MPIOptionalRequest.get().test();
			if (retry && !MPIOptionalStatus) {
				LOG_WARN(1, "Canceling pending request that looks to be stuck tag: " << _tag);
				MPIOptionalRequest.get().cancel();
				MPIOptionalStatus = MPIOptionalRequest.get().test();
				if (!MPIOptionalStatus) {
					MPIOptionalRequest = irecv(rankSource);
					MPIOptionalStatus = MPIOptionalRequest.get().test();
				}
			}
			}
			if (!!MPIOptionalStatus) {
				mpi::status status = MPIOptionalStatus.get();
				queueMessage(status, rankSource);
				MPIOptionalRequest = irecv(rankSource);
				returnValue = true;
			}
		}
		processQueue();
		return returnValue;
	}
	long receiveAllIncomingMessages(bool untilFlushed = true) {
		long messages = this->receiveAll(untilFlushed);

		LOG_DEBUG(4, _tag << ": receiving all messages for tag. cp: " << this->getNumCheckpoints() << " msgCount: " << this->getNumDeliveries());
		bool mayHavePending = true;
		while (mayHavePending) {
			mayHavePending = false;
			for(int rankSource = 0 ; rankSource < this->getWorld().size(); rankSource++) {
				if (receiveIncomingMessage(rankSource)) {
					messages++;
					mayHavePending = true;
				}
			}
			mayHavePending &= untilFlushed;
		}
		if (messages > 0)
			LOG_DEBUG(4, _tag << ": processed messages: " << messages);

		return messages;
	}
private:
	Buffer &getBuffer(int rank) {
		Buffer &buf = _recvBuffers[rank];
		if (buf == NULL)
			buf = this->getNewBuffer();
		return buf;
	}

	long _finalize(int checkpointFactor) {
		long iterations = 0;
		long messages;
		do {
			messages = 0;
			messages += processQueue();
			messages += this->receiveAllIncomingMessages();
			messages += this->flushAll();
			if (this->reachedCheckpoint(checkpointFactor) && messages != 0) {
				LOG_DEBUG_OPTIONAL(3, true, "Recv " << _tag << ": Achieved checkpoint but more messages are pending");
			}
			if (messages == 0)
				WAIT_AND_WARN(++iterations, "_finalize(" << checkpointFactor << ") with messages: " << messages << " checkpoint: " << this->getNumCheckpoints());
		} while (!this->reachedCheckpoint(checkpointFactor));
		return messages;
	}
public:
	void finalize(int checkpointFactor = 1) {
		LOG_DEBUG_OPTIONAL(2, true, "Recv " << _tag << ": Entering finalize checkpoint: " << this->getNumCheckpoints() << " out of " << (checkpointFactor*this->getWorld().size()));

		long messages = 0;
		do {
			messages = _finalize(checkpointFactor);
			LOG_DEBUG_OPTIONAL(3, true, "waiting for message to be 0: " << messages);
		} while (messages != 0);

		LOG_DEBUG_OPTIONAL(3, true, "Recv " << _tag << ": Finished finalize checkpoint: " << this->getNumCheckpoints());
		this->resetCheckpoints();
	}


private:
	MPIOptionalRequest irecv(int sourceRank) {
		MPIOptionalRequest oreq;
		LOG_DEBUG(5, "Starting irecv for " << sourceRank << "," << _tag);
#ifdef OPENMP_CRITICAL_MPI
#pragma omp critical (MPI_buffer_irecv)
#endif
		oreq = this->getWorld().irecv(sourceRank, _tag, getBuffer(sourceRank), BufferBase::MESSAGE_BUFFER_SIZE);
		assert(!!oreq);
		_requestAttempts[sourceRank] = 0;
		return oreq;
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
			MessagePackage MessagePackage( buf, size, source, tag );

			_MessagePackageQueue.push_back( MessagePackage );

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

			while( !_MessagePackageQueue.empty() ) {
				MessagePackage MessagePackage = _MessagePackageQueue.front();

				_MessagePackageQueue.pop_front();

				if (MessagePackage.size == 0) {
					this->checkpoint();
					LOG_DEBUG(3, _tag << ": got checkpoint from " << MessagePackage.source << "/" << MessagePackage.tag << ": " << this->getNumCheckpoints());
				} else {
					this->processMessagePackage(MessagePackage);
				}
				this->returnBuffer(MessagePackage.buffer);
				processCount++;
			}

			_isProcessing = false;
		}
		return processCount;
	}


};


template <typename C, typename CProcessor, int BufferSize = MPI_BUFFER_DEFAULT_SIZE>
class MPISendMessageBuffer: public MPIMessageBuffer<C, CProcessor, BufferSize> {
public:
	typedef MPIMessageBuffer<C, CProcessor, BufferSize> BufferBase;
	typedef MPIRecvMessageBuffer<C, CProcessor, BufferSize> RecvBuffer;
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
	typedef CProcessor MessageClassProcessor;
	typedef std::vector< MPIMessageBufferBase* > RecvBufferCallbackVector;

protected:
	Buffer *_sendBuffers;
	int  *_offsets;
	SentBuffers _sentBuffers;

public:
	MPISendMessageBuffer(mpi::communicator &world, int messageSize, MessageClassProcessor processor = MessageClassProcessor()) :
		BufferBase(world, messageSize, processor) {
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

	MessageClass *bufferMessage(int rankDest, int tagDest, int trailingBytes = 0) {
		bool wasSent;
		long messages = 0;
		return bufferMessage(rankDest, tagDest, wasSent, messages, trailingBytes);
	}
	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest, bool &wasSent, long &messages, int trailingBytes = 0) {
		int &offset = _offsets[rankDest];
		wasSent = flushIfFull(rankDest, tagDest, messages, trailingBytes);
		Buffer &buffStart = getBuffer(rankDest);
		assert(buffStart != NULL);
		MessageClass *buf = (MessageClass *) (buffStart+offset);
		offset += this->getMessageSize() + trailingBytes;
		this->newMessage();

		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg, int trailingBytes = 0) {
		char *buf = (char *) bufferMessage(rankDest, tagDest, trailingBytes);
		memcpy(buf, (char *) msg, this->getMessageSize() + trailingBytes);
	}

	long receiveAllIncomingMessages(bool untilFlushed = true) {
		return this->receiveAll(untilFlushed);
	}

	long flushIfFull(int rankDest, int tagDest, int trailingBytes = 0) {
		long messages = 0;
		flushIfFull(rankDest, tagDest, messages, trailingBytes);
		return messages;
	}
	bool flushIfFull(int rankDest, int tagDest, long &messages, int trailingBytes) {
		bool wasSent = false;
		int &offset = _offsets[rankDest];
		while (offset != 0 && offset + this->getMessageSize() + trailingBytes >= this->getSoftMaxBufferSize()) {
			messages += flushMessageBuffer(rankDest, tagDest);
			wasSent = true;
		}
		return wasSent;
	}
	long flushMessageBuffer(int rankDest, int tagDest, bool sendZeroMessage = false) {
		long messages = 0;
		long newMessages = 0;

		int &offsetLocation = _offsets[rankDest];
		bool waitForSend = sendZeroMessage == true;

		long iterations = 0;
		messages += checkSentBuffers( waitForSend );
		while (sendZeroMessage && (offsetLocation > 0 || newMessages != 0) ) {
			// flush all messages until there is nothing in the buffer
			newMessages = 0;
			WAIT_AND_WARN(++iterations, "flushMessageBuffer(" << rankDest << ", " << tagDest << ", " << sendZeroMessage << ")");
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
#pragma omp critical (MPI_buffer_send)
#endif
			sent.request = this->getWorld().isend(rankDest, tagDest, buffer, offset);

			messages++;

			recordSentBuffer(sent);

		}

		// do not let the sent queue get too large
		iterations = 0;
		newMessages = 0;
		while (newMessages == 0 && _sentBuffers.size() >= (size_t) BufferBase::BUFFER_QUEUE_SOFT_LIMIT * this->getWorld().size()) {
			if (newMessages != 0)
				WAIT_AND_WARN(++iterations, "flushMessageBuffer(" << rankDest << ", " << tagDest << ", " << sendZeroMessage << ") in checkSentBuffers loop");
			newMessages = checkSentBuffers( true );
			messages += newMessages;
		}

		return messages;
	}

	void recordSentBuffer(SentBuffer &sent) {
		_sentBuffers.push_back(sent);
	}
	void checkSent(SentBuffer &sent) {

		MPIOptionalStatus optStatus;
#ifdef OPENMP_CRITICAL_MPI
#pragma omp critical (MPI_BUFFER_SEND_REQUEST_TEST)
#endif
		{
		optStatus = sent.request.test();

		if ( _RETRY_MESSAGES && (!optStatus) && ++sent.pollCount > _RETRY_THRESHOLD) {
			sent.request.cancel();
			if ( ! sent.request.test() ) {
				sent.request = this->getWorld().isend(sent.destRank, sent.destTag, sent.buffer, sent.size);
				sent.pollCount = 0;
				LOG_WARN(1, "Canceled and retried pending message to " << sent.destRank << ", " << sent.destTag << " size " << sent.size << " deliveryCount " << sent.deliveryNum);
			}
			optStatus = sent.request.test();
		}
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
	long checkSentBuffers(bool wait = false) {
		long messages = 0;
		long iterations = 0;
		while (wait || iterations++ == 0) {
			messages += this->receiveAllIncomingMessages(wait);
			for(SentBuffersIterator it = _sentBuffers.begin(); it != _sentBuffers.end(); it++) {
				SentBuffer &sent = *it;
				checkSent(sent);
			}
			_sentBuffers.remove_if(SentBuffer());

			if (wait) {
				WAIT_AND_WARN(iterations+1, "checkSentBuffers()" );
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
		long messages = 0;
		long iterations = 0;
		while (iterations == 0 || messages != 0) {
			if (iterations++ > 0 )
				WAIT_AND_WARN(iterations, "flushAllMessageBuffersUntilEmpty(" << tagDest << ")");
			messages = flushAllMessageBuffers(tagDest);
			messages += checkSentBuffers(true);
		}

	}
	void finalize(int tagDest) {
		LOG_DEBUG_OPTIONAL(2, true, "Send " << tagDest << ": entering finalize()");

		// first clear the buffer and in-flight messages
		flushAllMessagesUntilEmpty(tagDest);
		receiveAllIncomingMessages();

		LOG_DEBUG_OPTIONAL(3, true,"Send " << tagDest << ": entering finalize stage2()");
		// send zero-message as checkpoint signal to stop
		flushAllMessageBuffers(tagDest, true);

		LOG_DEBUG_OPTIONAL(3, true,"Send " << tagDest << ": entering finalize stage3()");
		// now continue to flush until there is nothing left in the buffer;
		flushAllMessagesUntilEmpty(tagDest);

		LOG_DEBUG_OPTIONAL(3, true,"Send " << tagDest << ": finished finalize()");

	}
	long processPending() {
		return checkSentBuffers(false);
	}
};

#endif /* MPIBUFFER_H_ */
