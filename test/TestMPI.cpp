/*
 * TestMPI.cpp
 *
 *  Created on: Nov 12, 2010
 *      Author: regan
 */


#include <stdio.h>
#include <omp.h>
#include "config.h"
#include "Options.h"
#include "MPIBuffer.h"
#include "Log.h"

#include <boost/date_time/posix_time/posix_time.hpp>

#ifndef ENABLE_MPI
#error "mpi is required for this library"
#endif

class TextMessage;
typedef MPIMessageBufferBase< TextMessage > TextMessageBuffersBase;
typedef MPIRecvMessageBuffer< TextMessage > RecvTextMessageBufferBase;
typedef MPISendMessageBuffer< TextMessage > SendTextMessageBufferBase;

class RecvTextMessageBuffer : public RecvTextMessageBufferBase
{
public:
	RecvTextMessageBuffer(mpi::communicator &world, int messageSize, int tag) : RecvTextMessageBufferBase(world, messageSize, tag) {}
	virtual void processMessage(const char *text, int length) {
		LOG_DEBUG(4, "Recv: " << std::string(text, length));
	}
};
class SendTextMessageBuffer : public SendTextMessageBufferBase
{
public:
	SendTextMessageBuffer(mpi::communicator &world, int messageSize) : SendTextMessageBufferBase(world, messageSize) {}
};
class TextMessage {
public:
	int length;
	void set(const char *_data, int _length) {
		memcpy(getText(), _data, _length);
		LOG_DEBUG(5, "TextMessage::set(): " << std::string(getText(), _length));
		length = _length;
	}
	// THIS IS DANGEROUS unless allocated extra bytes!
	char *getText() {
		return (char*) (((char*)this)+sizeof(*this));
	}
	int process(RecvTextMessageBufferBase *_bufferCallback) {
		RecvTextMessageBuffer *bufferCallback = (RecvTextMessageBuffer *) _bufferCallback;
		bufferCallback->processMessage(getText(), length);
		LOG_DEBUG(4, "process(): length: " << length << " recvSize: " << bufferCallback->getRecvSize() << " msg#: " << bufferCallback->getNumMessages() << " deliver#: " << bufferCallback->getNumDeliveries() << " " << std::string(getText(), length))
		return length;
	}
};

int main(int argc, char **argv)
{
  int provided;
  provided = MPI::Init_thread(MPI_THREAD_MULTIPLE);
  mpi::environment env(argc, argv);
  mpi::communicator world;

  Logger::setWorld(&world);

  Options::parseOpts(argc, argv);

  printf("Hello World from Node %d\n", world.rank());

  int numThreads = omp_get_max_threads();
  RecvTextMessageBuffer *recv[numThreads];
  SendTextMessageBuffer *send[numThreads];


#pragma omp parallel num_threads(numThreads)
  {
	  int threadId = omp_get_thread_num();
	  recv[threadId] = new RecvTextMessageBuffer(world, sizeof(TextMessage), threadId);
	  send[threadId] = new SendTextMessageBuffer(world, sizeof(TextMessage));
	  send[threadId]->addReceiveAllCallback(recv[threadId]);
  }

  int spamMax = TextMessageBuffersBase::MESSAGE_BUFFER_SIZE;
  char spam[spamMax];

  #pragma omp parallel for
  for(int i = 0 ; i < spamMax; i++) {
	  spam[i] = 'a' + (i%26);
  }

  int mb = 1024;
  LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Send/Recv " << (mb*numThreads) << "MB per node");
  world.barrier();
  boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();

  #pragma omp parallel num_threads(numThreads)
  {
	  int threadId = omp_get_thread_num();
	  int msgSize = 16*1024 - sizeof(TextMessage);
	  int msgPerMb = 1024*1024 / msgSize;
	  int destRank = threadId;
	  for(int i = 0; i < msgPerMb * mb; i++) {
		  while (++destRank % world.size() == world.rank());
		  send[threadId]->bufferMessage(destRank % world.size(), threadId, msgSize)->set(spam, msgSize);
	  }
	  send[threadId]->flushAllMessageBuffers(threadId);
	  send[threadId]->checkSentBuffers(true);
	  LOG_VERBOSE_OPTIONAL(1, true, "sent " << send[threadId]->getNumDeliveries());
  }

  std::stringstream ss;
  #pragma omp parallel num_threads(numThreads)
  {
	  int threadId = omp_get_thread_num();
	  send[threadId]->finalize(threadId);
	  recv[threadId]->finalize();

	  #pragma omp critical
	  {
		  ss << threadId << " " << send[threadId]->getNumMessages()
		     << " " << recv[threadId]->getNumMessages() << std::endl;
	  }
	  delete recv[threadId];
	  delete send[threadId];
  }
  boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration elapsed = end - start;

  std::string s = ss.str();
  LOG_VERBOSE_OPTIONAL(1, true, s);
  double rate = (double) mb*numThreads / (double) elapsed.total_milliseconds() * 1000.0;
  LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, rate << " MB/s " << rate * 8 << " Mbit/s " << elapsed.total_milliseconds() << "ms");

  MPI::Finalize();
  return 0;
}

