/*
 * TestMPI.cpp
 *
 *  Created on: Nov 12, 2010
 *      Author: regan
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


#include <stdio.h>
#include <algorithm>

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
		LOG_DEBUG(1, "TextMessage::set(" << (void *)this << " to " << (void*) (getText() + _length)  << ", " << _length << "): " << std::string(getText(), std::min(_length, 25)) );
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
		LOG_DEBUG(4, "TextMessageProcessor::process(): length: " << msg->length << " recvSize: " << msgPkg.size << " msg#: " << msgPkg.bufferCallback->getNumMessages() << " deliver#: " << msgPkg.bufferCallback->getNumDeliveries() << " " << std::string(msg->getText(), std::min(msg->length, 25)))
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

class _MPITestOptions : public OptionsBaseInterface {
public:
	void _resetDefaults() {
		GeneralOptions::getOptions().getDebug() = 1;
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {
		GeneralOptions::getOptions()._setOptions(desc,p);
		MPIOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::getOptions()._parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);
		return ret;
	}
};

typedef OptionsBaseTemplate< _MPITestOptions > MPITestOptions;

int main(int argc, char **argv)
{
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
	mpi::environment env(argc, argv);
	mpi::communicator world;

	MPITestOptions::parseOpts(argc,argv);
	Logger::setWorld(&world);

	printf("Hello World from Node %d, mpi support %d\n", world.rank(), provided);

	if (provided < MPI_THREAD_FUNNELED)
		omp_set_num_threads(1);

	int numThreads = omp_get_max_threads();

	char *spam = NULL;

{
	A2ABufferBase a2a(world, sizeof(TextMessage));

	long unsigned int spamMax = a2a.getBufferSize() - sizeof(TextMessage);
	int msgSize = std::min(16*1024 - sizeof(TextMessage), spamMax);
	int msgPerMb = 1024*1024 / msgSize;
	assert(msgPerMb > 0);

	LOG_VERBOSE_OPTIONAL(1, world.rank()==0, "Sending messages of size: " << msgSize);
	spam = new char[spamMax];

#pragma omp parallel for
	for(int i = 0 ; i < (int) spamMax; i++) {
		spam[i] = 'a' + (i%26);
	}

	int mb = 1;
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Send/Recv " << (mb) << "MB per rank");
	world.barrier();
	boost::posix_time::ptime start = boost::posix_time::microsec_clock::local_time();

	boost::posix_time::ptime end = boost::posix_time::microsec_clock::local_time();
	boost::posix_time::time_duration elapsed = end - start;

	double rate = (double) mb  / (double) elapsed.total_milliseconds() * 1000.0;

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

	char *in = new char[spamMax * world.size()];
	char *out = new char[spamMax * world.size()];
	int inSize[world.size()], outSize[world.size()], inDisp[world.size()], outDisp[world.size()];

	world.barrier();
	start = boost::posix_time::microsec_clock::local_time();

	for(int i = 0; i < (int) (mb*1024*1024 / spamMax / world.size()); i++) {
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

	delete [] in;
	delete [] out;
}
	delete [] spam;

	MPI_Finalize();

	return 0;
}

