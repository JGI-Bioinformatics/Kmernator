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
};

class TextMessageProcessor;
typedef MPIMessageBuffer< TextMessage, TextMessageProcessor > TextMessageBufferBase;
class TextMessageProcessor {
public:
	long process(TextMessage *msg, MessagePackage &msgPkg) {
		LOG_DEBUG(4, "TextMessageProcessor::process(): length: " << msg->length << " recvSize: " << msgPkg.size << " msg#: " << msgPkg.bufferCallback->getNumMessages() << " deliver#: " << msgPkg.bufferCallback->getNumDeliveries() << " " << std::string(msg->getText(), msg->length))
		return msg->length;
	}
};
typedef MPIRecvMessageBuffer< TextMessage, TextMessageProcessor > RecvTextMessageBufferBase;
typedef MPISendMessageBuffer< TextMessage, TextMessageProcessor > SendTextMessageBufferBase;
typedef MPIAllToAllMessageBuffer< TextMessage, TextMessageProcessor > A2ABufferBase;

class RecvTextMessageBuffer : public RecvTextMessageBufferBase
{
public:
	RecvTextMessageBuffer(mpi::communicator &world, int messageSize, int tag) : RecvTextMessageBufferBase(world, messageSize, tag) {}
	virtual void processTextMessage(const char *text, int length) {
		LOG_DEBUG(4, "Recv: " << std::string(text, length));
	}
};
class SendTextMessageBuffer : public SendTextMessageBufferBase
{
public:
	SendTextMessageBuffer(mpi::communicator &world, int messageSize) : SendTextMessageBufferBase(world, messageSize) {}
};

class _MPITestOptions : public _MPIOptions {
public:
	void _resetDefaults() {}
	void _setOptions(po::options_description &desc,
				po::positional_options_description &p) {
		GeneralOptions::getOptions()._setOptions(desc,p);
		_MPIOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::getOptions()._parseOptions(vm);
		ret &= _MPIOptions::_parseOptions(vm);
		return ret;
	}
};

typedef OptionsBaseTemplate< _MPITestOptions > MPITestOptions;

int main(int argc, char **argv)
{
  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  mpi::environment env(argc, argv);
  mpi::communicator world;

  MPITestOptions::parseOpts(argc,argv);
  Logger::setWorld(&world);

  printf("Hello World from Node %d, mpi support %d\n", world.rank(), provided);

  int numThreads = omp_get_max_threads();
  int stdNumThreads = numThreads;
  if (!MPIMessageBufferBase::isThreadSafe())
	  stdNumThreads = 1;


  A2ABufferBase a2a(world, sizeof(TextMessage));

  RecvTextMessageBuffer *recv[numThreads];
  SendTextMessageBuffer *send[numThreads];

#pragma omp parallel num_threads(stdNumThreads)
  {
	  //int threadId = omp_get_thread_num();
	  //recv[threadId] = new RecvTextMessageBuffer(world, sizeof(TextMessage), threadId);
	  //send[threadId] = new SendTextMessageBuffer(world, sizeof(TextMessage));
	  //send[threadId]->addReceiveAllCallback(recv[threadId]);
  }
  int spamMax = MPIOptions::getOptions().getTotalBufferSize() / world.size() / 20;

  char spam[spamMax];

  #pragma omp parallel for
  for(int i = 0 ; i < spamMax; i++) {
	  spam[i] = 'a' + (i%26);
  }

  int mb = 1;
  LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Send/Recv " << (mb) << "MB per rank");
  world.barrier();
  boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();
  int msgSize = 16*1024 - sizeof(TextMessage);
  int msgPerMb = 1024*1024 / msgSize;

  #pragma omp parallel num_threads(stdNumThreads)
  {
	  int threadId = omp_get_thread_num();
	  for(int i = 0; i < msgPerMb * mb / stdNumThreads / world.size(); i++) {
		  for(int w=0; w < world.size() ; w++) {
			  //send[threadId]->bufferMessage(w, threadId, msgSize)->set(spam, msgSize);
		  }
	  }
	 // send[threadId]->finalize(threadId);
	  //recv[threadId]->finalize();

	 // delete recv[threadId];
	 // delete send[threadId];
  }

  boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration elapsed = end - start;

  double rate = (double) mb  / (double) elapsed.total_milliseconds() * 1000.0;
  LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "MPIMessagBuffer: " << rate << " MB/s " << rate * 8 << " Mbit/s " << elapsed.total_milliseconds() << "ms");


  world.barrier();
  start = boost::posix_time::microsec_clock::local_time();

  #pragma omp parallel num_threads(numThreads)
  {
	  int threadId = omp_get_thread_num();
	  for(int i = 0; i < msgPerMb * mb / numThreads / world.size(); i++) {
		  for(int w=0; w < world.size() ; w++) {
			  a2a.bufferMessage(w, threadId, msgSize)->set(spam, msgSize);
		  }
	  }
	  a2a.finalize();
  }

  end = boost::posix_time::microsec_clock::local_time();
  elapsed = end - start;
  rate = (double) mb / (double) elapsed.total_milliseconds() * 1000.0;
  LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "MPI all2all messagebuffer: " << rate << " MB/s " << rate * 8 << " Mbit/s " << elapsed.total_milliseconds() << "ms");

  char in[spamMax * world.size()], out[spamMax * world.size()];
  int inSize[world.size()], outSize[world.size()], inDisp[world.size()], outDisp[world.size()];

  world.barrier();
  start = boost::posix_time::microsec_clock::local_time();

  for(int i = 0; i < mb*1024*1024 / spamMax / world.size(); i++) {
	  for(int w = 0 ; w < world.size(); w++) {
		  inSize[w] = outSize[w] = spamMax;
		  inDisp[w] = outDisp[w] = spamMax * w;
		  memcpy(out + outDisp[w], spam, spamMax);
	  }
	  LOG_DEBUG_OPTIONAL(2, world.rank() == 0, "All2All: " << i);
	  MPI_Alltoallv(out, outSize, outDisp, MPI_BYTE, in, inSize, inDisp, MPI_BYTE, world);
	  for(int w = 0 ; w < world.size(); w++) {
		  if (memcmp(out + outDisp[w], in + inDisp[w], inSize[w]) != 0)
			  LOG_ERROR(1, "invalid transfer: " << w);
	  }
  }

  end = boost::posix_time::microsec_clock::local_time();
  elapsed = end - start;
  rate = (double) mb / (double) elapsed.total_milliseconds() * 1000.0;
  LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "MPI all2all: " << rate << " MB/s " << rate * 8 << " Mbit/s " << elapsed.total_milliseconds() << "ms");


  MPI_Finalize();
  return 0;
}

