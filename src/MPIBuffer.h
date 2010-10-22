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
    virtual int receiveAllIncomingMessages() { return 0; }
    virtual void flushAllMessageBuffers(int tagDest) { }

    int receiveAll() {
    	int count = 0;
    	for(unsigned int i = 0; i < _receiveAllCallbacks.size(); i++)
    		count += _receiveAllCallbacks[i].first->receiveAllIncomingMessages();
    	return count;
    }
    void flushAll() {
    	for(unsigned int i = 0; i < _flushAllCallbacks.size(); i++)
    		_flushAllCallbacks[i].first->flushAllMessageBuffers(_flushAllCallbacks[i].second);
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

protected:
        mpi::communicator _world;
        int _messageSize;
        char *_message;

public:
        MPIMessageBufferBase(mpi::communicator &world, int messageSize)
        : _world(world), _messageSize(messageSize) {
        	_message = new char[ getMessageSize() ];
        	assert(getMessageSize() >= (int) sizeof(MessageClass));
        }
        ~MPIMessageBufferBase() {
        	delete [] _message;
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
			++_requestAttempts[rankSource];
			bool retry =  _RETRY_MESSAGES && _requestAttempts[rankSource] > _RETRY_THRESHOLD;

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
				process(status, rankSource);
				optionalRequest = irecv(rankSource);
				returnValue = true;
			}
		}
		return returnValue;
	}
	int receiveAllIncomingMessages() {
		int messages = this->receiveAll();

		LOG_DEBUG(3, _tag << ": receiving all messages for tag. cp: " << getNumCheckpoints() << " msgCount: " << this->getNumDeliveries());
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
		LOG_DEBUG(1, _tag << ": Entering finalize checkpoint: " << getNumCheckpoints());
		while ( ! reachedCheckpoint(checkpointFactor) ) {
			this->receiveAllIncomingMessages();
			this->flushAll();
			boost::this_thread::sleep( boost::posix_time::milliseconds(2) );
		}
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
		LOG_DEBUG(4, "Starting irecv for " << sourceRank << "," << _tag);
#ifdef OPENMP_CRITICAL_MPI
		#pragma omp critical (MPI_buffer)
#endif
		oreq = this->getWorld().irecv(sourceRank, _tag, _recvBuffers + sourceRank*BufferBase::MESSAGE_BUFFER_SIZE, BufferBase::MESSAGE_BUFFER_SIZE);
		assert(!!oreq);
		_requestAttempts[sourceRank] = 0;
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
			this->newMessage();
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

			this->newMessageDelivery();
			LOG_DEBUG(3, _tag << ": received delivery " << this->getNumDeliveries() << " from " << source << "," << tag << " size " << size << " probe attempts: " << _requestAttempts[rankSource]);

			if (size == 0) {
				checkpoint();
				LOG_DEBUG(2, _tag << ": got checkpoint from " << source << ": " << _numCheckpoints);
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
	typedef std::vector< _MPIMessageBufferBase* > RecvBufferCallbackVector;

protected:
	char *_sendBuffers;
	int  *_offsets;

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
		this->newMessage();

		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg) {
		char *buf = (char*) bufferMessage(rankDest, tagDest);
		memcpy(buf, (char*) msg, this->getMessageSize());
	}

	int receiveAllIncomingMessages() {
		return this->receiveAll();
	}

protected:
	void flushMessageBuffer(int rankDest, int tagDest, char *buffer, int &offset, bool sendZeroMessage = false) {
		this->receiveAllIncomingMessages();
		if (offset > 0 || sendZeroMessage) {
			mpi::request request;
			LOG_DEBUG(3, "sending message to " << rankDest << ", " << tagDest << " size " << offset);
#ifdef OPENMP_CRITICAL_MPI
#pragma omp critical (MPI_buffer)
#endif
			request = this->getWorld().isend(rankDest, tagDest, buffer, offset);
			this->newMessageDelivery();
			int count = 0;
			OptionalStatus optStatus = request.test();
			while ( ! optStatus ) {
				if (_RETRY_MESSAGES && ++count > _RETRY_THRESHOLD) {
					request.cancel();
					if ( ! request.test() ) {
						request = this->getWorld().isend(rankDest, tagDest, buffer, offset);
						count = 0;
						LOG_WARN(1, "Canceled and retried pending message to " << rankDest << ", " << tagDest << " size " << offset << " deliveryCount " << this->getNumDeliveries());
					}
				}
				LOG_DEBUG(3, "waiting for send to finish to " << rankDest << ", " << tagDest << " size " << offset << " deliveryCount " << this->getNumDeliveries() << " attempt count " << count);
				this->receiveAllIncomingMessages();
				boost::this_thread::sleep( boost::posix_time::milliseconds(2) );
				optStatus = request.test();
			}
			mpi::status status = optStatus.get();
			// hmmm looks like error is sometimes populated with junk...
			//if (status.error() > 0)
			//	LOG_WARN(1, "sending message returned an error: " << status.error());
			LOG_DEBUG(3, "finished sending message to " << rankDest << ", " << tagDest << " size " << offset << " deliveryCount " << this->getNumDeliveries());
		}
		offset = 0;
	}
public:
	void flushMessageBuffer(int rankDest, int tagDest, bool sendZeroMessage = false) {
		int &offset = _offsets[rankDest];
		char *buffer = _sendBuffers + BufferBase::MESSAGE_BUFFER_SIZE * rankDest;
		flushMessageBuffer(rankDest, tagDest, buffer, offset, sendZeroMessage);
	}
	void flushAllMessageBuffers(int tagDest) {
		flushAllMessageBuffers(tagDest, false);
	}
	void flushAllMessageBuffers(int tagDest, bool sendZeroMessage) {
		for(int rankDest = 0 ; rankDest < this->getWorld().size(); rankDest++)
			flushMessageBuffer(rankDest, tagDest, sendZeroMessage);
		this->flushAll();
	}
	void finalize(int tagDest) {
		LOG_DEBUG(3, "entering finalize() for tag " << tagDest);
		flushAllMessageBuffers(tagDest);
		// send zero message buffer as checkpoint signal to stop
		LOG_DEBUG(3, "sending checkpoint tag " << tagDest);
		flushAllMessageBuffers(tagDest, true);
		LOG_DEBUG(3, "sent checkpoints tag " << tagDest);
	}
};

#endif /* MPIBUFFER_H_ */
