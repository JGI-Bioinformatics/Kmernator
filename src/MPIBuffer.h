//
// Kmernator/src/MPIBuffer.h
//
// Author: Rob Egan
//
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

#ifndef MPIBUFFER_H_
#define MPIBUFFER_H_

#include "config.h"
#include "Options.h"
#include "MPIBase.h"

#ifdef _USE_OPENMP
#define X_OPENMP_CRITICAL_MPI
#endif

#define _RETRY_MESSAGES false
#define _RETRY_THRESHOLD 10000

#include <vector>

// use about 32MB of memory total to batch & queue up messages between communications
// this is split up across world * thread * thread arrays
#define MPI_BUFFER_DEFAULT_SIZE (32 * 1024 * 1024)
#define MPI_MIN_TRANSMIT_DEFAULT_SIZE 2048

class _MPIOptions : public OptionsBaseInterface {
public:
	_MPIOptions() : mpiBufferSize(MPI_BUFFER_DEFAULT_SIZE), mpiMinTransmitSize(MPI_MIN_TRANSMIT_DEFAULT_SIZE) {}
	virtual ~_MPIOptions() {}
	int &getTotalBufferSize() {
		return mpiBufferSize;
	}
	int &getMinTransmitSize() {
		return mpiMinTransmitSize;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("MPI Options");
		opts.add_options()

					("mpi-buffer-size", po::value<int>()->default_value(mpiBufferSize),
							"total amount of RAM to devote to MPI message batching buffers in bytes")
					("mpi-min-transmit-size", po::value<int>()->default_value(mpiMinTransmitSize), "the minimum inter rank-thread buffer size")
							;
		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		setOpt("mpi-buffer-size", mpiBufferSize);
		setOpt("mpi-min-transmit-size", mpiMinTransmitSize);
		return true;
	}
protected:
	int mpiBufferSize, mpiMinTransmitSize;
};
typedef OptionsBaseTemplate< _MPIOptions > MPIOptions;

class MPIMessageBufferBase {
public:
	typedef std::pair<MPIMessageBufferBase*, int> CallbackBase;
	typedef std::vector<CallbackBase> CallbackVector;

protected:
	CallbackVector _flushAllCallbacks;
	CallbackVector _receiveAllCallbacks;
	CallbackVector _sendReceiveCallbacks;
	long _deliveries;
	long _numMessages;
	long _sendBytes, _recvBytes;
	long _syncPoints;
	long _syncAttempts;
	double _transitTime;
	double _threadWaitingTime, _masterWaitingTime;
	double _threadProcessTime, _masterProcessTime;
	double _startTime;

public:
	MPIMessageBufferBase() :
		_deliveries(0), _numMessages(0), _sendBytes(0), _recvBytes(0), _syncPoints(0), _syncAttempts(0), _transitTime(0), _threadWaitingTime(0), _masterWaitingTime(0), _threadProcessTime(0), _masterProcessTime(0) {
		_startTime = MPI_Wtime();
	}
	virtual ~MPIMessageBufferBase() {
		LOG_VERBOSE_GATHER(2, "~MPIMessageBufferBase(): " << _deliveries
				<< " deliveries, " << _numMessages << " messages, "
				<< _sendBytes << " b sent, " << _recvBytes << " b recv, "
				<< _syncPoints << " / " << _syncAttempts << " syncs in " << (MPI_Wtime() - _startTime)
				<< " seconds (" << _transitTime << " transit, " << _threadProcessTime << " threadProcess, " << _masterProcessTime << " masterProcess, " << _threadWaitingTime << " threadWait, " << _masterWaitingTime << " masterWait)");
	}
	static bool isThreadSafe() {
		int threadlevel;
		MPI_Query_thread(&threadlevel);
		return threadlevel >= MPI_THREAD_FUNNELED;
	}

	// receive buffers to flush before and/or during send
	void addReceiveAllCallback(MPIMessageBufferBase *receiveAllBuffer) {
		_receiveAllCallbacks.push_back(CallbackBase(receiveAllBuffer, -1));
		receiveAll();
		assert(receiveAllBuffer->getNumDeliveries() == 0);
	}
	void addFlushAllCallback(MPIMessageBufferBase *flushAllBuffer, int tagDest) {
		_flushAllCallbacks.push_back(CallbackBase(flushAllBuffer, tagDest));
		flushAll();
		assert(flushAllBuffer->getNumDeliveries() == 0);
	}
	void addSendReceiveCallback(MPIMessageBufferBase *sendReceiveBuffer) {
		_sendReceiveCallbacks.push_back(CallbackBase(sendReceiveBuffer, -1));
	}
	virtual long receiveAllIncomingMessages(bool untilFlushed = true) {
		return 0;
	}
	virtual long flushAllMessageBuffers(int tagDest) {
		return 0;
	}
	virtual long sendReceive(bool isFinalized) {
		return 0;
	}

	long receiveAll(bool untilFlushed = true) {
		long count = 0;
		for (unsigned int i = 0; i < _receiveAllCallbacks.size(); i++)
			count += _receiveAllCallbacks[i].first->receiveAllIncomingMessages(
					untilFlushed);
		LOG_DEBUG(4, "receiveAll() with " << count);
		return count;
	}
	long flushAll() {
		long count = 0;
		for (unsigned int i = 0; i < _flushAllCallbacks.size(); i++)
			count += _flushAllCallbacks[i].first->flushAllMessageBuffers(
					_flushAllCallbacks[i].second);
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
	inline void syncPoint() {
#pragma omp atomic
		_syncPoints++;
	}
	inline void syncAttempt() {
#pragma omp atomic
		_syncAttempts++;
	}
	inline void transit(double transit) {
#pragma omp atomic
		_transitTime += transit;
	}
	inline void threadWait(double wait) {
		if (omp_get_thread_num() == 0) {
#pragma omp atomic
			_masterWaitingTime += wait;
		} else {
#pragma omp atomic
			_threadWaitingTime += wait;
		}
	}
	inline void threadProcess(double duration) {
		if (omp_get_thread_num() == 0) {
#pragma omp atomic
			_masterProcessTime += duration;
		} else {
#pragma omp atomic
			_threadProcessTime += duration;
		}

	}

};

class MessagePackage {
public:
	typedef char * Buffer;
	Buffer buffer;
	int size;
	int source;
	int tag;
	MPIMessageBufferBase *bufferCallback;
	MessagePackage(Buffer b, int s, int src, int t, MPIMessageBufferBase *m) :
		buffer(b), size(s), source(src), tag(t), bufferCallback(m) {
	}
};

template<typename C>
class DummyProcessor {
public:
	int process(C *msg, MessagePackage &buffer) {
		return 0;
	}
};

template<typename C, typename CProcessor = DummyProcessor<C> >
class MPIMessageBuffer: public MPIMessageBufferBase {
public:
	static const int BUFFER_QUEUE_SOFT_LIMIT = 5;
	static const int BUFFER_INSTANCES = 5;

	typedef C MessageClass;
	typedef CProcessor MessageClassProcessor;
	typedef typename MessagePackage::Buffer Buffer;
	typedef std::vector<Buffer> FreeBufferCache;
	typedef std::vector<FreeBufferCache> ThreadedFreeBufferCache;
	typedef std::list<MessagePackage> MessagePackageQueue;

protected:
	mpi::communicator _world;
	int _bufferSize;
	int _messageSize;
	MessageClassProcessor _processor;
	float _softRatio;
	int _softMaxBufferSize;
	int _softNumThreads;
	int _numCheckpoints;
	ThreadedFreeBufferCache _freeBuffers;
	int _numThreads;
	int _worldSize;
	int _numWorldThreads;

public:
	MPIMessageBuffer(mpi::communicator &world, int messageSize,
			MessageClassProcessor processor = MessageClassProcessor(),
			int totalBufferSize = MPIOptions::getOptions().getTotalBufferSize(),
			float softRatio = 0.90) :
				_world(world, mpi::comm_duplicate), _bufferSize(totalBufferSize / _world.size() / BUFFER_INSTANCES),
				_messageSize(messageSize), _processor(processor), _softRatio(
						softRatio), _softNumThreads(0), _numCheckpoints(0) {
		assert(getMessageSize() >= (int) sizeof(MessageClass));
		assert(_softRatio >= 0.0 && _softRatio <= 1.0);
		_worldSize = _world.size();
		if (omp_in_parallel())
			_numThreads = omp_get_num_threads();
		else
			_numThreads = omp_get_max_threads();
		_softNumThreads = std::max(_numThreads * 3 / 4, 1);
		_bufferSize /= (_numThreads*_numThreads);
		_bufferSize = std::max(std::max(messageSize,MPIOptions::getOptions().getMinTransmitSize()/_numThreads), _bufferSize);
		setSoftMaxBufferSize();

		_numWorldThreads = _worldSize * _numThreads;
		_freeBuffers.resize(_numThreads);
		for (int threadId = 0; threadId < _numThreads; threadId++)
			_freeBuffers[threadId].reserve(BUFFER_QUEUE_SOFT_LIMIT
					* _numWorldThreads + 1);
		LOG_DEBUG_OPTIONAL(1, Logger::isMaster(), "MPIMessageBuffer(): " << _numThreads << " threads, " << _bufferSize << " size " << getMessageSize() << " message size");
	}
	~MPIMessageBuffer() {
		for (ThreadedFreeBufferCache::iterator it = _freeBuffers.begin(); it != _freeBuffers.end(); it++)
			for (FreeBufferCache::iterator it2 = it->begin(); it2 != it->end(); it2++)
				delete [] *it2;
		_freeBuffers.clear();
	}
	inline mpi::communicator &getWorld() {
		return _world;
	}
	inline int getNumThreads() const {
		return _numThreads;
	}
	inline int getSoftNumThreads() const {
		return _softNumThreads;
	}
	inline int getNumWorldThreads() const {
		return _numWorldThreads;
	}
	inline int getWorldSize() const {
		return _worldSize;
	}
	inline int getMessageSize() const {
		return _messageSize;
	}
	inline MessageClassProcessor &getProcessor() {
		return _processor;
	}
	void setMessageProcessor(MessageClassProcessor processor) {
		_processor = processor;
	}
	inline int getSoftMaxBufferSize() const {
		return _softMaxBufferSize;
	}
	void setSoftMaxBufferSize() {
		_softMaxBufferSize = std::max(getMessageSize(), (int) (_bufferSize
				* _softRatio));
		assert(_softMaxBufferSize >= 0 && _softMaxBufferSize <= _bufferSize);
	}
	inline int getBufferSize() const {
		return _bufferSize;
	}
	void setBufferSize(int bufferSize) {
		_bufferSize = bufferSize;
		setSoftMaxBufferSize();
	}
	Buffer getNewBuffer() {
		int threadId = omp_get_thread_num();
		Buffer buf;
		if (_freeBuffers[threadId].empty()) {
			buf = new char[_bufferSize];
			LOG_DEBUG(1, "getNewBuffer(): " << (void*) buf << " of " << _bufferSize << " ending: " << (void*) (buf + _bufferSize));
		} else {
			buf = _freeBuffers[threadId].back();
			_freeBuffers[threadId].pop_back();
		}
		return buf;

	}
	void returnBuffer(Buffer buf) {
		int threadId = omp_get_thread_num();
		if (_freeBuffers[threadId].size() >= (size_t) (BUFFER_QUEUE_SOFT_LIMIT
				* _numWorldThreads)) {
			delete [] buf;
			buf = NULL;
		} else {
			_freeBuffers[threadId].push_back(buf);
		}
	}

	long processMessagePackage(MessagePackage &messagePackage) {
		double duration = MPI_Wtime();
		Buffer msg, start = messagePackage.buffer;

		long count = 0;
		int offset = 0;
		while (offset < messagePackage.size) {
			msg = start + offset;
			LOG_DEBUG(4, "processMessagePackage(): (" << messagePackage.source
					<< ", " << messagePackage.tag << "): " << (void*) msg
					<< " offset: " << offset);
			int trailingBytes = _processor.process((MessageClass*) msg,
					messagePackage);
			offset += this->getMessageSize() + trailingBytes;
			this->newMessage();
			count++;
		}
		this->threadProcess(MPI_Wtime() - duration);
		return count;
	}

	inline int getNumCheckpoints() const {
		return _numCheckpoints;
	}
	bool reachedCheckpoint(int checkpointFactor = 1) const {
		return getNumCheckpoints() == _worldSize * checkpointFactor;
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

template<typename C, typename CProcessor >
class MPIAllToAllMessageBuffer: public MPIMessageBuffer<C, CProcessor > {
public:
	typedef MPIMessageBuffer<C, CProcessor > BufferBase;
	typedef typename BufferBase::Buffer Buffer;
	typedef typename BufferBase::MessagePackageQueue MessagePackageQueue;
	typedef C MessageClass;
	typedef CProcessor MessageClassProcessor;
	static const int NUM_BUFFERS = 4;

	class MessageHeader {
	public:
		MessageHeader() :
			offset(0), threadSource(0), tag(0) {
			setDummy();
		}
		MessageHeader(const MessageHeader &copy) :
			offset(copy.offset), threadSource(copy.threadSource),
			tag(copy.tag), dummy(copy.dummy) {
		}
		MessageHeader &operator=(const MessageHeader &other) {
			if (this == &other)
				return *this;
			offset = other.offset;
			threadSource = other.threadSource;
			tag = other.tag;
			dummy = other.dummy;
			return *this;
		}
		inline void reset(int _threadSource = 0, int _tag = 0) {
			LOG_DEBUG(5, "MessageHeader::reset(" << _threadSource << ", "
					<< _tag << "):" << (void*) this);
			offset = 0;
			threadSource = _threadSource;
			tag = _tag;
			setDummy();
		}
		inline void resetOffset() {
			LOG_DEBUG(5, "MessageHeader::resetOffset():" << (void*) this
					<< " from: " << offset);
			offset = 0;
			setDummy();
		}
		int append(int dataSize) {
			LOG_DEBUG(5, "MessageHeader::append(" << dataSize << "):"
					<< (void*) this << " " << offset);
			offset += dataSize;
			setDummy();
			return offset;
		}
		inline bool validate() const {
			return dummy == _getDummy();
		}
		inline int getOffset() const {
			return offset;
		}
		inline int getThreadSource() const {
			return threadSource;
		}
		inline int getTag() const {
			return tag;
		}
		std::string toString() const {
			std::stringstream ss;
			ss << "MessageHeader(" << (void*) this << "): offset: " << offset
					<< " threadSource: " << threadSource << " tag: " << tag
					<< " dummy: " << dummy << " valid: " << validate();
			return ss.str();
		}
	private:
		int offset, threadSource, tag, dummy;
		void setDummy() {
			dummy = _getDummy();
		}
		int _getDummy() const {
			return (offset + threadSource + tag);
		}

	};
	class BuildBuffer {
	public:
		Buffer buffer;
		MessageHeader header;
		BuildBuffer() :
			buffer(NULL), header(MessageHeader()) {
		}
		~BuildBuffer() {
			reset();
		}
		void reset() {
			header.reset();
			if (buffer != NULL) {
				delete [] buffer;
				buffer = NULL;
			}
		}
	};
	class TransmitBuffer {
	public:
		enum StateType { EMPTY_OUT, BUILDING_OUT, READY_OUT, EMPTY_IN, BUILDING_IN, READY_IN, DRAINING_IN, UNUSED };
		TransmitBuffer(int _numThreads, int _worldSize, int _numTags, int bufferSize) :
			numThreads(_numThreads), worldSize(_worldSize), numTags(_numTags), buildSize(0), finalCount(0) {
			dataSize = sizeof(int) + numThreads * (numThreads * numTags)
							* (sizeof(MessageHeader) + bufferSize);
			totalSize = getHeaderSize() + worldSize * dataSize;
			jumps = new int[numThreads * worldSize * numThreads * numTags];
			xmit = new char[totalSize];
			// in/out displacements do not change
			for (int rankDest = 0; rankDest < worldSize; rankDest++) {
				getOffset(rankDest) = rankDest * dataSize;
			}
			threadStates.resize(numThreads);
			setAllStates(UNUSED);
		}
		~TransmitBuffer() {
			delete [] jumps;
			delete [] xmit;
		}
		inline int getJump(int rankDest, int threadDest, int threadId = 0) {
			return threadId * worldSize * numThreads * numTags + rankDest
					* numThreads * numTags + threadDest;
		}
		inline Buffer getJumpBuffer(int rankDest, int threadDest, int threadId) {
			return getBuffer(rankDest) + jumps[getJump(rankDest,threadDest,threadId)];
		}
		inline int &getSize(int rankDest) {
			// Transmit size for each rank (includes extra int of datasize)
			int *o = (int*) xmit;
			return *(o + rankDest);
		}
		inline int &getOffset(int rankDest) {
			// displacement for each rank
			int *o = (int*) xmit;
			return *(o + worldSize + rankDest);
		}
		inline int &getDataSize(int rankDest) {
			// total size of the data sent
			int *o = (int*) xmit;
			Buffer buf = (Buffer) (o + worldSize * 2);
			return *((int*) (buf + rankDest * dataSize));
		}
		inline Buffer getBuffer(int rankDest) {
			// the data (immediately after datasize int)
			return (Buffer) (&getDataSize(rankDest) + 1);
		}
		inline int getHeaderSize() {
			return worldSize * sizeof(int) * 2;
		}
		inline int getBuildSize() const {
			return buildSize;
		}
		void prepOut() {
			assert( areAllInState(UNUSED) );
			for (int rankDest = 0; rankDest < worldSize; rankDest++) {
				getSize(rankDest) = 0;
				getDataSize(rankDest) = 0;
			}
			buildSize = 0;
			setAllStates(EMPTY_OUT);
			finalCount = 0;
		}
		void prepIn() {
			assert( areAllInState(UNUSED) );
			for (int rankDest = 0; rankDest < worldSize; rankDest++) {
				getSize(rankDest) = dataSize;
			}
			buildSize = 0;
			setAllStates(EMPTY_IN);
			finalCount = 0;
		}
		void reset() {
			assert( areAllInState(UNUSED) || areAllInState(EMPTY_IN) || areAllInState(READY_OUT));
			buildSize = 0;
			setAllStates(UNUSED);
			finalCount = 0;
		}
		static void AllToAll(TransmitBuffer &out, TransmitBuffer &in, mpi::communicator &world) {
			// mpi_alltoall
			assert(out.areAllInState(READY_OUT));
			assert(in.areAllInState(EMPTY_IN));
			in.setAllStates(BUILDING_IN);
			MPI_Alltoallv(out.xmit + out.getHeaderSize(), &out.getSize(0),
					&out.getOffset(0), MPI_BYTE, in.xmit
					+ in.getHeaderSize(), &in.getSize(0),
					&in.getOffset(0), MPI_BYTE, world);

			out.setAllStates(UNUSED);
			in.setAllStates(READY_IN);
		}
		// must be in critical section!
		void setSize(int rankDest, int threadDest, int threadId, BuildBuffer &bb) {

			int jump = getJump(rankDest, threadDest, threadId);
			int &outSize = getSize(rankDest);

			jumps[jump] = outSize;
			int size = bb.header.getOffset();
			buildSize += size;
			outSize += size + sizeof(MessageHeader);
			LOG_DEBUG(4, "sendReceive(): sending (" << rankDest << ", "
					<< threadDest << "): " << size
					<< " bytes, outSize: " << outSize << " buildSize: "
					<< buildSize);
		}
		void setTransmitSizes(bool isFinalized) {
			long bytesSend = 0, bytesRecv = 0;
			setTransmitSizes(isFinalized, bytesSend, bytesRecv);
		}
		void setTransmitSizes(bool isFinalized, long &bytesSend, long &bytesRecv) {
			if (isFinalized && buildSize == 0) {
				for (int destRank = 0; destRank < worldSize; destRank++) {
					getDataSize(destRank) = 0; // send no MessageHeaders
					// add the (int) datasize to the total length sent so something is always sent
					bytesSend += getSize(destRank) = sizeof(int);
				}
			} else {
				// set transmitted sizes
				for (int destRank = 0; destRank < worldSize; destRank++) {
					int dataSize = getSize(destRank);
					getDataSize(destRank) = dataSize; // includes MessageHeaders
					// add the (int) datasize to the total length sent
					bytesSend += getSize(destRank) += sizeof(int);
				}
			}
		}
		bool inState(StateType _state, int threadNum) const {
			return threadStates[threadNum] == _state;
		}
		bool inState(StateType _state) const {
			return inState(_state, omp_get_thread_num());
		}
		bool areAllInState(StateType _state) const {
			for(int i = 0; i < numThreads; i++)
				if ( !inState(_state, i) ) {
					LOG_DEBUG(5, "areAllInState() Failed: " << toString() << " i:" << i << " ts:" << threadStates[i] << " _state: " << _state);
					return false;
				}
			return true;
		}
		bool areAllReadyOut() const {
			return areAllInState(READY_OUT);
		}
		bool isEmptyOut() const {
			return inState(EMPTY_OUT);
		}
		bool isReadyOut() const {
			return inState(READY_OUT);
		}
		bool isEmptyIn() const {
			return inState(EMPTY_IN);
		}
		bool isBuildingIn() const {
			return inState(BUILDING_IN);
		}
		bool isReadyIn() const {
			return inState(READY_IN);
		}
		bool isUnused() const {
			return inState(UNUSED);
		}
		void setFinal() {
#pragma omp atomic
			finalCount++;
		}
		bool areAllFinal() {
			return finalCount == numThreads;
		}
		void setThreadBuildingOut() {
			assert(inState(EMPTY_OUT));
			setThreadState(BUILDING_OUT);
		}
		void setThreadReadyOut() {
			assert(inState(BUILDING_OUT));
			setThreadState(READY_OUT);
		}
		void setThreadDrainingIn() {
			assert(inState(READY_IN));
			setThreadState(DRAINING_IN);
		}
		void setThreadFinishedIn() {
			assert(inState(DRAINING_IN));
			setThreadState(UNUSED);
		}
		void clearFinal() {
			finalCount = 0;
		}
		std::string toString() const {
			std::stringstream ss;
			ss << "TransmitBuffer(" << (void*) this << "): {buildSize: " << buildSize << ", ";
			ss << "finalCount: " << finalCount << ", ";
			for(int i = 0; i < numThreads; i++)
				ss << i << ": " << threadStates[i] << ", ";
			ss << "}";
			return ss.str();
		}
	private:
		int numThreads, worldSize, numTags, buildSize, dataSize, totalSize, finalCount;
		Buffer xmit;
		int *jumps;
		std::vector< StateType > threadStates;
		void setThreadState(StateType _newState) {
			threadStates[omp_get_thread_num()] = _newState;
		}
		void setAllStates(StateType _newState) {
			for(int i = 0; i < numThreads; i++)
				threadStates[i] = _newState;
		}
	};

private:
	TransmitBuffer* buffers[NUM_BUFFERS];
	int currentBuffer;
	std::vector<std::vector<std::vector<BuildBuffer> > > buildsTWT;
	int numTags;
	int threadsSending;

public:

	MPIAllToAllMessageBuffer(mpi::communicator &world, int messageSize,
			MessageClassProcessor processor = MessageClassProcessor(),
			int _numTags = 1, int totalBufferSize = MPIOptions::getOptions().getTotalBufferSize(), double softRatio = 0.90) :
				BufferBase(world, messageSize, processor, totalBufferSize, softRatio),
				numTags(_numTags), threadsSending(0) {
		assert(!omp_in_parallel());
		assert(omp_get_thread_num() == 0);
		assert(numTags > 0);
		int worldSize = this->getWorldSize();
		int numThreads = this->getNumThreads();

		buildsTWT.resize(numThreads);
		for (int threadId = 0; threadId < numThreads; threadId++) {
			buildsTWT[threadId].resize(worldSize);
		}
		for(int i = 0; i < NUM_BUFFERS; i++) {
			buffers[i] = new TransmitBuffer(numThreads, worldSize, numTags, this->getBufferSize());
		}
		currentBuffer = 0;

		for (int rankDest = 0; rankDest < worldSize; rankDest++) {
			for (int threadId = 0; threadId < numThreads; threadId++) {
				buildsTWT[threadId][rankDest].resize(numThreads * numTags, BuildBuffer());
				for (int threadDest = 0; threadDest < (numThreads * numTags); threadDest++) {
					BuildBuffer &bb = buildsTWT[threadId][rankDest][threadDest];
					bb.header.reset(threadId, threadDest);
					bb.buffer = this->getNewBuffer();
				}
			}
		}
	}
	~MPIAllToAllMessageBuffer() {
		assert(!omp_in_parallel());
		for(int i = 0; i < NUM_BUFFERS; i++) {
			if (buffers[i] != NULL) {
				delete buffers[i];
				buffers[i] = NULL;
			}
		}
	}

private:
	long sendReceive(bool isFinalized) {
		assert(omp_get_max_threads() == 1 || omp_in_parallel());

		if (isFinalized)
			this->syncAttempt();
		long iterations = 0;
		double waitTime;

		int thisBuffer = currentBuffer;

		TransmitBuffer &last  = *buffers[(thisBuffer+0) % NUM_BUFFERS];
		TransmitBuffer &in    = *buffers[(thisBuffer+1) % NUM_BUFFERS];
		TransmitBuffer &out   = *buffers[(thisBuffer+2) % NUM_BUFFERS];
		//TransmitBuffer &last2 = *buffers[(thisBuffer+3) % NUM_BUFFERS];

		LOG_DEBUG(3, "sendReceive(" << isFinalized << ") "
				<< this->getNumMessages() << " messages "
				<< this->getNumDeliveries() << " deliveries "
				<< thisBuffer << " bufferCount "
				<< threadsSending << " threadsSending");
		// reset checkpoints before processing 'last'
		if (omp_get_thread_num() == 0) {
			assert(out.areAllInState( TransmitBuffer::UNUSED ));
			out.prepOut();
			assert(in.areAllInState( TransmitBuffer::UNUSED ));;
			in.prepIn();
		}

#pragma omp atomic
		threadsSending++;

		// ensure in, out & last buffer are all prepared
		// and all threads have reset checkpoints
		LOG_DEBUG(3, "sendReceive() starting barrier1, buffer: " << thisBuffer);
		waitTime = MPI_Wtime();
#pragma omp barrier
		this->threadWait( MPI_Wtime() - waitTime );
		// make sure all threads get here before incrementing the buffer or resetting checkpoints

		if (omp_get_thread_num() == 0) {
			this->resetCheckpoints();
#pragma omp atomic
			currentBuffer++;
		}
		if (isFinalized)
			out.setFinal();

		iterations = 1;
		while (! out.isEmptyOut() )
			WAIT_AND_WARN(iterations, "sendReceive() waiting for out buffer to be prepared!!! (" << TransmitBuffer::EMPTY_OUT << "): buffer: " << thisBuffer << " " << out.toString()); // spin lock

		assert(out.isEmptyOut());
		assert(in.isEmptyIn());

#pragma omp atomic
		threadsSending--;

		waitTime = MPI_Wtime();
#pragma omp barrier
		this->threadWait( MPI_Wtime() - waitTime );
		// make sure all threads always use the same currentBuffer (between two barriers)
		LOG_DEBUG(3, "sendReceive() past barrier2, buffer: " << thisBuffer);

		// now process buffers without any further thread-wide barriers
		// .. master will block until all out buffers are ready
		_prepBuffersByThread(out);

		assert( last.isUnused() || last.isReadyIn() || last.isBuildingIn() );

		// spin lock this thread until its last buffer is ready (last round's communication has finished)
		waitTime = MPI_Wtime();
		iterations = 1;
		while ( last.isBuildingIn() )
			WAIT_AND_WARN(iterations, "sendReceive() waiting for last buffer to finish building (" << TransmitBuffer::READY_IN << "): buffer: " << thisBuffer << " " << last.toString()); // spin lock
		this->threadWait( MPI_Wtime() - waitTime );

		long numReceived = 0;
		if (last.isReadyIn())
			numReceived += _processBuffersByThread(last);

		if (out.areAllFinal()) {
			// wait for all to process the last buffer, update checkpoints, then check it
			LOG_DEBUG_OPTIONAL(2, true, "sendReceive() Starting areAllFinal barrier on buffer: " << thisBuffer << " " << threadsSending);
			assert(out.isReadyOut());
			assert(last.isUnused());

#pragma omp barrier

			if (out.getBuildSize() == 0 && this->reachedCheckpoint(this->getNumThreads())) {
				if (omp_get_thread_num() == 0) {
					out.reset();
					in.reset();
				}
				LOG_DEBUG_OPTIONAL(1, omp_get_thread_num() == 0, "sendReceive() areAllFinal() Found checkpoint. master stopping on buffer: " << thisBuffer << " " << threadsSending);
				return numReceived;
			} else {
				LOG_DEBUG_OPTIONAL(2, true, "sendReceive() areAllFinal() did not reach checkpoint. buffer: " << thisBuffer << " " << threadsSending);
			}
		}

		if (omp_get_thread_num() == 0)
		{
			// spin lock master until all out buffers are ready (all threads have called _prepBuffersByThread )
			waitTime = MPI_Wtime();
			iterations = 1;
			while( ! out.areAllReadyOut() )
				WAIT_AND_WARN(iterations, "sendReceive() waiting for out buffer to finish building ( " << TransmitBuffer::READY_OUT << ") on all threads: buffer: " << thisBuffer << " " << out.toString()); // spin lock
			this->threadWait( MPI_Wtime() - waitTime );
			assert( in.areAllInState( TransmitBuffer::EMPTY_IN) );
			assert( out.areAllInState( TransmitBuffer::READY_OUT) );

			out.setTransmitSizes(out.areAllFinal(), this->_sendBytes, this->_recvBytes);

			LOG_DEBUG(4, "sendReceive(): Starting all2all on buffer: " << thisBuffer << " threadsSending: " << threadsSending);
			waitTime = MPI_Wtime();
			// mpi_alltoall
			TransmitBuffer::AllToAll(out, in, this->getWorld());
			this->transit(MPI_Wtime() - waitTime);
			this->newMessageDelivery();
			LOG_DEBUG(4, "sendReceive(): Finished all2all on buffer: " << thisBuffer << " threadsSending: " << threadsSending);
		}

		return numReceived;
	}
public:
	long sendReceive() {
		return sendReceive(false);
	}

	int getBytesInBuffer() const {
		int offset = 0;
		for(int rankDest = 0 ; rankDest < this->getWorldSize(); rankDest++) {
			for(int tagDest = 0 ; tagDest < (this->getNumThreads() * numTags); tagDest++) {
				assert(omp_get_thread_num() < this->getNumThreads());
				const BuildBuffer &bb = buildsTWT[omp_get_thread_num()][rankDest][tagDest];
				offset += bb.header.getOffset();
			}
		}
		return offset;
	}

	int processMessages(int sourceRank, MessageHeader *header) {
		assert(header->validate());

		Buffer buf = (Buffer) (header + 1);
		MessagePackage msgPkg(buf, header->getOffset(), sourceRank,
				header->getTag(), this);
		return this->processMessagePackage(msgPkg);
	}

	void finalize() {
		assert(omp_get_max_threads() == 1 || omp_in_parallel());
		LOG_DEBUG_OPTIONAL(2, true, "Entering finalize()");
		while (!this->reachedCheckpoint(this->getNumThreads())) {
			sendReceive(true);
		}
		this->syncPoint();
#pragma omp barrier
	}
	bool isReadyToSend(int offset, int trailingBytes) {
		assert(trailingBytes + this->getMessageSize() <= this->getBufferSize());
		return ( ( (offset >= this->getSoftMaxBufferSize()) & (threadsSending > this->getSoftNumThreads()) )
				| ((offset + trailingBytes + this->getMessageSize()) >= this->getBufferSize()) );
	}

	MessageClass *bufferMessage(int rankDest, int tagDest) {
		return bufferMessage(rankDest, tagDest, 0);
	}
	MessageClass *bufferMessage(int rankDest, int tagDest, int trailingBytes) {
		bool wasSent = false;
		long messages = 0;
		return bufferMessage(rankDest, tagDest, wasSent, messages,
				trailingBytes);
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg,
			int trailingBytes = 0) {
		char *buf = (char *) bufferMessage(rankDest, tagDest, trailingBytes);
		memcpy(buf, (char *) msg, this->getMessageSize() + trailingBytes);
	}

	MessageClass *bufferMessage(int rankDest, int tagDest, bool &wasSent,
			long &messages) {
		return bufferMessage(rankDest, tagDest, wasSent, messages, 0);
	}
	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest, bool &wasSent,
			long &messages, int trailingBytes) {
		int threadId = omp_get_thread_num();
		assert(threadId < this->getNumThreads() && rankDest < this->getWorldSize() && tagDest < this->numTags * this->getNumThreads());
		BuildBuffer &bb = buildsTWT[threadId][rankDest][tagDest];

		while (isReadyToSend(bb.header.getOffset(), trailingBytes))
			messages += sendReceive(false);

		assert(bb.buffer != NULL);
		MessageClass *buf =
				(MessageClass *) (bb.buffer + bb.header.getOffset());
		int msgSize = this->getMessageSize() + trailingBytes;
		bb.header.append(msgSize);
		this->newMessage();
		LOG_DEBUG(1, "bufferMessage(" << rankDest << ", " << tagDest << ", "
				<< wasSent << ", " << messages << ", " << trailingBytes
				<< "): " << (void*) buf << " size: " << msgSize);
		assert(((Buffer)buf) + msgSize < bb.buffer + this->getBufferSize());
		return buf;
	}
private:
	void _prepBuffersByThread(TransmitBuffer &out) {
		assert(omp_get_max_threads() == 1 || omp_in_parallel());
		int threadId = omp_get_thread_num();
		int numThreads = this->getNumThreads();
		int worldSize = this->getWorldSize();

		out.setThreadBuildingOut();

		double waitTime = MPI_Wtime();
#pragma omp critical
		{
			waitTime = MPI_Wtime() - waitTime;
			// allocate the out message headers
			for (int rankDest = 0; rankDest < worldSize; rankDest++) {
				for (int threadDest = 0; threadDest < (numThreads * numTags); threadDest++) {
					assert(threadId < this->getNumThreads() && rankDest < this->getWorldSize() && threadDest < this->numTags * this->getNumThreads());
					BuildBuffer &bb = buildsTWT[threadId][rankDest][threadDest];
					assert(bb.header.validate());
					out.setSize(rankDest, threadDest, threadId, bb);
				}
			}
		}
		this->threadWait(waitTime);

		// copy any headers and messages to out
		for (int rankDest = 0; rankDest < worldSize; rankDest++) {
			for (int threadDest = 0; threadDest < (numThreads * numTags); threadDest++) {
				assert(threadId < this->getNumThreads() && rankDest < this->getWorldSize() && threadDest < this->numTags * this->getNumThreads());
				BuildBuffer &bb = buildsTWT[threadId][rankDest][threadDest];
				assert(bb.header.validate());

				MessageHeader *outHeader = (MessageHeader *) out.getJumpBuffer(rankDest, threadDest, threadId);
				*outHeader = bb.header;
				if (bb.header.getOffset() > 0)
					memcpy(outHeader + 1, bb.buffer, bb.header.getOffset());

				LOG_DEBUG(5, "_prepBuffersByThread(): Prepared for rank " << rankDest << " threadDest " << threadDest << " " << outHeader->toString());
				// reset the counter
				bb.header.resetOffset();
			}
		}
		out.setThreadReadyOut();
	}

	long _processBuffersByThread(TransmitBuffer &in) {
		assert(omp_get_max_threads() == 1 || omp_in_parallel());
		int threadId = omp_get_thread_num();
		int numThreads = this->getNumThreads();
		int worldSize = this->getWorldSize();
		long numReceived = 0;

		in.setThreadDrainingIn();

		for (int sourceRank = 0; sourceRank < worldSize; sourceRank++) {
			int transmitSize = in.getSize(sourceRank);
			int size = in.getDataSize(sourceRank);
			assert(size >= 0);
			LOG_DEBUG(3, "_processBuffersByThread(): Received (" << sourceRank << "):"
					<< size << " bytes, " << transmitSize << " maxTransmitSize");
			if (size == 0) {
				// no MessageHeaders - set checkpoint
				this->checkpoint();
				continue;
			}
			int tagsProcessed = 0;
			Buffer buf = in.getBuffer(sourceRank);
			Buffer begin = buf;
			Buffer end = buf + size;
			while (begin != end) {
				assert(begin < end);
				MessageHeader *header = (MessageHeader*) begin;
				LOG_DEBUG(5, "_processBuffersByThread(): source " << sourceRank << " header " << (void*) begin << " "
						<< header->toString());

				assert(header->validate());

				if (threadId == (header->getTag() % numThreads)) {
					if (header->getOffset() > 0)
						numReceived += processMessages(sourceRank, header);
					tagsProcessed++;
					LOG_DEBUG(5, "_processBuffersByThread(): source " << sourceRank << " processed " << header->getTag() << " at " << header->getOffset() << " " << (void*) header);
				}
				begin = ((Buffer) (header)) + sizeof(MessageHeader) + header->getOffset();
			}
			if (tagsProcessed != numTags * numThreads)
				LOG_WARN(1, "_processBuffersByThread(): data count problem! source " << sourceRank << " processed " << tagsProcessed << " out of " << numTags);
			assert(tagsProcessed == numTags * numThreads);
		}
		in.setThreadFinishedIn();

		return numReceived;
	}

};

template<typename C, typename CProcessor >
class MPIRecvMessageBuffer: public MPIMessageBuffer<C, CProcessor > {
public:
	typedef MPIMessageBuffer<C, CProcessor > BufferBase;

	typedef typename BufferBase::Buffer Buffer;
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
	MPIRecvMessageBuffer(mpi::communicator &world, int messageSize, int tag =
			mpi::any_tag, MessageClassProcessor processor =
					MessageClassProcessor()) :
						BufferBase(world, messageSize, processor), _tag(tag), _isProcessing(
								false) {
		_recvBuffers = new Buffer[this->getWorld().size()];
		_requests = new MPIOptionalRequest[this->getWorld().size()];
		_requestAttempts.resize(this->getWorld().size(), 0);
		for (int destRank = 0; destRank < this->getWorld().size(); destRank++) {
			_recvBuffers[destRank] = NULL;
			_requests[destRank] = irecv(destRank);
		}
	}
	~MPIRecvMessageBuffer() {
		assert(_MessagePackageQueue.empty());
		cancelAllRequests();
		for (int i = 0; i < this->getWorld().size(); i++) {
			Buffer &buf = _recvBuffers[i];
			if (buf != NULL)
				delete [] buf;
			buf = NULL;
		}
		delete [] _recvBuffers;
		delete [] _requests;
	}
	void cancelAllRequests() {
		for (int source = 0; source < this->getWorld().size(); source++) {
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
			bool retry = _RETRY_MESSAGES && _requestAttempts[rankSource]
			                                                 > _RETRY_THRESHOLD;
			{
				MPIOptionalStatus = MPIOptionalRequest.get().test();
				if (retry && !MPIOptionalStatus) {
					LOG_WARN(1,
							"Canceling pending request that looks to be stuck tag: "
							<< _tag);
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

		LOG_DEBUG(4, _tag << ": receiving all messages for tag. cp: "
				<< this->getNumCheckpoints() << " msgCount: "
				<< this->getNumDeliveries());
		bool mayHavePending = true;
		while (mayHavePending) {
			mayHavePending = false;
			for (int rankSource = 0; rankSource < this->getWorld().size(); rankSource++) {
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
		long iterations = 1;
		long messages;
		do {
			messages = 0;
			messages += processQueue();
			messages += this->receiveAllIncomingMessages();
			messages += this->flushAll();
			if (this->reachedCheckpoint(checkpointFactor) && messages != 0) {
				LOG_DEBUG_OPTIONAL(
						3,
						true,
						"Recv " << _tag
						<< ": Achieved checkpoint but more messages are pending");
			}
			if (messages == 0)
				WAIT_AND_WARN(iterations, "_finalize(" << checkpointFactor << ") with messages: " << messages << " checkpoint: " << this->getNumCheckpoints());
		} while (!this->reachedCheckpoint(checkpointFactor));
		return messages;
	}
public:
	void finalize(int checkpointFactor = 1) {
		LOG_DEBUG_OPTIONAL(2, true, "Recv " << _tag
				<< ": Entering finalize checkpoint: "
				<< this->getNumCheckpoints() << " out of " << (checkpointFactor
						* this->getWorld().size()));

		long messages = 0;
		do {
			messages = _finalize(checkpointFactor);
			LOG_DEBUG_OPTIONAL(3, true, "waiting for message to be 0: "
					<< messages);
		} while (messages != 0);

		LOG_DEBUG_OPTIONAL(3, true, "Recv " << _tag
				<< ": Finished finalize checkpoint: "
				<< this->getNumCheckpoints());
		this->resetCheckpoints();
		this->syncPoint();
	}

private:
	MPIOptionalRequest irecv(int sourceRank) {
		MPIOptionalRequest oreq;
		LOG_DEBUG(5, "Starting irecv for " << sourceRank << "," << _tag);
		oreq = this->getWorld().irecv(sourceRank, _tag, getBuffer(sourceRank),
				this->getBufferSize());
		assert(!!oreq);
		_requestAttempts[sourceRank] = 0;
		return oreq;
	}
	bool queueMessage(mpi::status &status, int rankSource) {

		bool wasMessage = false;
		int source = status.source();
		int tag = status.tag();
		int size = status.count<char> ().get();
		if (status.cancelled()) {
			LOG_WARN(1, _tag << ": request was successfully canceled from "
					<< source << " size " << size);
		} else {

			assert( rankSource == source );
			assert( _tag == tag );

			Buffer &bufLoc = getBuffer(source);
			Buffer buf = bufLoc;
			bufLoc = NULL;
			MessagePackage MessagePackage(buf, size, source, tag, this);

			_MessagePackageQueue.push_back(MessagePackage);

			this->newMessageDelivery();
			LOG_DEBUG(4, _tag << ": received delivery "
					<< this->getNumDeliveries() << " from " << source << ","
					<< tag << " size " << size << " probe attempts: "
					<< _requestAttempts[rankSource]);

			wasMessage = true;
		}
		return wasMessage;

	}
	long processQueue() {
		long processCount = 0;
		// code to allow recursive message generation
		if (!_isProcessing) {
			_isProcessing = true;

			while (!_MessagePackageQueue.empty()) {
				MessagePackage MessagePackage = _MessagePackageQueue.front();

				_MessagePackageQueue.pop_front();

				if (MessagePackage.size == 0) {
					this->checkpoint();
					LOG_DEBUG(3, _tag << ": got checkpoint from "
							<< MessagePackage.source << "/"
							<< MessagePackage.tag << ": "
							<< this->getNumCheckpoints());
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

template<typename C, typename CProcessor >
class MPISendMessageBuffer: public MPIMessageBuffer<C, CProcessor > {
public:
	typedef MPIMessageBuffer<C, CProcessor > BufferBase;
	typedef MPIRecvMessageBuffer<C, CProcessor > RecvBuffer;
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

		SentBuffer(Buffer &_buffer, int _rank, int _tag, int _size, int _id) :
			request(), buffer(_buffer), destRank(_rank), destTag(_tag), size(
					_size), deliveryNum(_id), pollCount(0) {
		}
		SentBuffer() :
			request(), buffer(NULL), destRank(mpi::any_source), destTag(
					mpi::any_tag), size(0), deliveryNum(0), pollCount(0) {
		}
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
		bool operator()(const SentBuffer &sent) {
			return (sent.buffer == NULL);
		}
	};
	typedef std::list<SentBuffer> SentBuffers;
	typedef typename SentBuffers::iterator SentBuffersIterator;
	typedef C MessageClass;
	typedef CProcessor MessageClassProcessor;
	typedef std::vector<MPIMessageBufferBase*> RecvBufferCallbackVector;

protected:
	Buffer *_sendBuffers;
	int *_offsets;
	SentBuffers _sentBuffers;

public:
	MPISendMessageBuffer(mpi::communicator &world, int messageSize,
			MessageClassProcessor processor = MessageClassProcessor()) :
				BufferBase(world, messageSize, processor) {
		_sendBuffers = new Buffer[world.size()];
		_offsets = new int[this->getWorld().size()];
		for (int destRank = 0; destRank < this->getWorld().size(); destRank++) {
			_sendBuffers[destRank] = NULL;
			_offsets[destRank] = 0;
		}
	}
	~MPISendMessageBuffer() {
		for (int i = 0; i < this->getWorld().size(); i++) {
			Buffer &buf = _sendBuffers[i];
			if (buf != NULL)
				delete[] buf;
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

	MessageClass *bufferMessage(int rankDest, int tagDest, int trailingBytes =
			0) {
		bool wasSent;
		long messages = 0;
		return bufferMessage(rankDest, tagDest, wasSent, messages,
				trailingBytes);
	}
	// returns a pointer to the next message.  User can use this to create message
	MessageClass *bufferMessage(int rankDest, int tagDest, bool &wasSent,
			long &messages, int trailingBytes = 0) {
		int &offset = _offsets[rankDest];
		wasSent = flushIfFull(rankDest, tagDest, messages, trailingBytes);
		Buffer &buffStart = getBuffer(rankDest);
		assert(buffStart != NULL);
		MessageClass *buf = (MessageClass *) (buffStart + offset);
		offset += this->getMessageSize() + trailingBytes;
		this->newMessage();

		return buf;
	}
	// copies msg as the next message in the buffer
	void bufferMessage(int rankDest, int tagDest, MessageClass *msg,
			int trailingBytes = 0) {
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
	bool flushIfFull(int rankDest, int tagDest, long &messages,
			int trailingBytes) {
		bool wasSent = false;
		int &offset = _offsets[rankDest];
		while (offset != 0 && offset + this->getMessageSize() + trailingBytes
				>= this->getSoftMaxBufferSize()) {
			messages += flushMessageBuffer(rankDest, tagDest);
			wasSent = true;
		}
		return wasSent;
	}
	long flushMessageBuffer(int rankDest, int tagDest, bool sendZeroMessage =
			false) {
		long messages = 0;
		long newMessages = 0;

		int &offsetLocation = _offsets[rankDest];
		bool waitForSend = sendZeroMessage == true;

		long iterations = 1;
		messages += checkSentBuffers(waitForSend);
		while (sendZeroMessage && (offsetLocation > 0 || newMessages != 0)) {
			// flush all messages until there is nothing in the buffer
			newMessages = 0;
			WAIT_AND_WARN(iterations, "flushMessageBuffer(" << rankDest << ", " << tagDest << ", " << sendZeroMessage << ")");
			newMessages = flushMessageBuffer(rankDest, tagDest, false);
			newMessages += checkSentBuffers(waitForSend);
			messages += newMessages;
		}

		if (offsetLocation > 0 || sendZeroMessage) {
			Buffer &bufferLocation = getBuffer(rankDest);

			Buffer buffer = bufferLocation;
			bufferLocation = NULL;
			int offset = offsetLocation;
			offsetLocation = 0;

			this->newMessageDelivery();
			SentBuffer sent(buffer, rankDest, tagDest, offset,
					this->getNumDeliveries());
			LOG_DEBUG(3, "sending message to " << rankDest << ", " << tagDest
					<< " size " << offset);
			sent.request = this->getWorld().isend(rankDest, tagDest, buffer,
					offset);

			messages++;

			recordSentBuffer(sent);

		}

		// do not let the sent queue get too large
		iterations = 0;
		newMessages = 0;
		while (newMessages == 0 && _sentBuffers.size()
				>= (size_t) BufferBase::BUFFER_QUEUE_SOFT_LIMIT
				* this->getWorld().size()) {
			if (newMessages != 0)
				WAIT_AND_WARN(iterations, "flushMessageBuffer(" << rankDest << ", " << tagDest << ", " << sendZeroMessage << ") in checkSentBuffers loop");
			newMessages = checkSentBuffers(true);
			messages += newMessages;
		}

		return messages;
	}

	void recordSentBuffer(SentBuffer &sent) {
		_sentBuffers.push_back(sent);
	}
	void checkSent(SentBuffer &sent) {

		MPIOptionalStatus optStatus;
		{
			optStatus = sent.request.test();

			if (_RETRY_MESSAGES && (!optStatus) && ++sent.pollCount
					> _RETRY_THRESHOLD) {
				sent.request.cancel();
				if (!sent.request.test()) {
					sent.request = this->getWorld().isend(sent.destRank,
							sent.destTag, sent.buffer, sent.size);
					sent.pollCount = 0;
					LOG_WARN(1, "Canceled and retried pending message to "
							<< sent.destRank << ", " << sent.destTag
							<< " size " << sent.size << " deliveryCount "
							<< sent.deliveryNum);
				}
				optStatus = sent.request.test();
			}
		}

		if (!!optStatus) {

			mpi::status status = optStatus.get();
			// hmmm looks like error is sometimes populated with junk...
			//if (status.error() > 0)
			//	LOG_WARN(1, "sending message returned an error: " << status.error());
			LOG_DEBUG(4, "finished sending message to " << sent.destRank
					<< ", " << sent.destTag << " size " << sent.size
					<< " deliveryCount " << sent.deliveryNum << " status.error(): " << status.error());

			this->returnBuffer(sent.buffer);
			sent.reset();

		}

	}
	long checkSentBuffers(bool wait = false) {
		long messages = 0;
		long iterations = 0;
		while (wait || iterations++ == 0) {
			messages += this->receiveAllIncomingMessages(wait);
			for (SentBuffersIterator it = _sentBuffers.begin(); it
			!= _sentBuffers.end(); it++) {
				SentBuffer &sent = *it;
				checkSent(sent);
			}
			_sentBuffers.remove_if(SentBuffer());

			if (wait) {
				WAIT_AND_WARN(iterations, "checkSentBuffers()" );
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
		for (int rankDest = 0; rankDest < this->getWorld().size(); rankDest++)
			messages += flushMessageBuffer(rankDest, tagDest, sendZeroMessage);
		messages += this->flushAll();
		return messages;
	}
	void flushAllMessagesUntilEmpty(int tagDest) {
		long messages = 0;
		long iterations = 0;
		while (iterations == 0 || messages != 0) {
			if (iterations > 0)
				WAIT_AND_WARN(iterations, "flushAllMessageBuffersUntilEmpty(" << tagDest << ")");
			messages = flushAllMessageBuffers(tagDest);
			messages += checkSentBuffers(true);
		}

	}
	void finalize(int tagDest) {
		LOG_DEBUG_OPTIONAL(2, true, "Send " << tagDest
				<< ": entering finalize()");

		// first clear the buffer and in-flight messages
		flushAllMessagesUntilEmpty(tagDest);
		receiveAllIncomingMessages();

		LOG_DEBUG_OPTIONAL(3, true, "Send " << tagDest
				<< ": entering finalize stage2()");
		// send zero-message as checkpoint signal to stop
		flushAllMessageBuffers(tagDest, true);

		LOG_DEBUG_OPTIONAL(3, true, "Send " << tagDest
				<< ": entering finalize stage3()");
		// now continue to flush until there is nothing left in the buffer;
		flushAllMessagesUntilEmpty(tagDest);

		LOG_DEBUG_OPTIONAL(3, true, "Send " << tagDest
				<< ": finished finalize()");
		this->syncPoint();
	}
	long processPending() {
		return checkSentBuffers(false);
	}
};

#endif /* MPIBUFFER_H_ */
