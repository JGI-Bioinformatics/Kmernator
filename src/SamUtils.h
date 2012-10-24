/* Some code was incorporated from samtools 0.1.18 */
/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
 */

/*
 * SamUtils.h
 *
 *  Created on: Sep 21, 2012
 *      Author: regan
 */

#ifndef SAMUTILS_H_
#define SAMUTILS_H_

#include "mpi.h"
#include <stdint.h>
#include <string>
#include <algorithm>

#include "sam.h"
#include "bam.h"
#include "Log.h"
#include "BgzfStream.h"
#include "MemoryBufferStream.h"

namespace bam_endian {
#include "bam_endian.h"
static void swap_endian_data(const bam1_core_t *c, int data_len, uint8_t *data)
{
	uint8_t *s;
	int32_t i;
	uint32_t *cigar = (uint32_t*)(data + c->l_qname);
	s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
	for (i = 0; i < c->n_cigar; ++i) bam_swap_endian_4p(&cigar[i]);
	while (s < data + data_len) {
		uint8_t type;
		s += 2; // skip key
		type = toupper(*s); ++s; // skip type
		if (type == 'C' || type == 'A') ++s;
		else if (type == 'S') { bam_swap_endian_2p(s); s += 2; }
		else if (type == 'I' || type == 'F') { bam_swap_endian_4p(s); s += 4; }
		else if (type == 'D') { bam_swap_endian_8p(s); s += 8; }
		else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
		else if (type == 'B') {
			int32_t n, Bsize = bam_aux_type2size(*s);
			memcpy(&n, s + 1, 4);
			if (1 == Bsize) {
			} else if (2 == Bsize) {
				for (i = 0; i < n; i += 2)
					bam_swap_endian_2p(s + 5 + i);
			} else if (4 == Bsize) {
				for (i = 0; i < n; i += 4)
					bam_swap_endian_4p(s + 5 + i);
			}
			bam_swap_endian_4p(s+1);
		}
	}
}

}

class BamStreamUtils {
public:
	typedef std::vector< bam1_t * > BamVector;
	typedef std::vector< bam1_core_t > BamCoreVector;
	static std::ostream &writeBamHeaderPart1(std::ostream &os, bam_header_t *header) {
		return writeBamHeaderPart1(os, header->l_text, header->text, header->n_targets);
	}
	static std::ostream &writeBamHeaderPart1(std::ostream &os, int32_t l_text, char *text, int32_t n_targets) {
		int x;
		// write "BAM1"
		os.write("BAM\001", 4);
		// write plain text and the number of reference sequences
		if (bam_is_be) {
			x = bam_endian::bam_swap_endian_4(l_text);
			os.write((char*) &x, 4);
			if (l_text)
				os.write(text, l_text);
			x = bam_endian::bam_swap_endian_4(n_targets);
			os.write((char*) &x, 4);
		} else {
			os.write((char*) &l_text, 4);
			if (l_text)
				os.write((char*) text, l_text);
			os.write((char*) &n_targets, 4);
		}
		return os;

	}
	static std::ostream &writeBamHeaderPart2(std::ostream &os, bam_header_t *header) {
		int32_t i, name_len, x;
		// write sequence names and lengths
		for (i = 0; i != header->n_targets; ++i) {
			char *p = header->target_name[i];
			name_len = strlen(p) + 1;
			if (bam_is_be) {
				x = bam_endian::bam_swap_endian_4(name_len);
				os.write((char*) &x, 4);
			} else
				os.write((char*) &name_len, 4);
			os.write((char*) p, name_len);
			if (bam_is_be) {
				x = bam_endian::bam_swap_endian_4(header->target_len[i]);
				os.write((char*) &x, 4);
			} else
				os.write((char*) &header->target_len[i], 4);
		}
		os.flush();
		return os;
	}
	static std::ostream &writeBam(std::ostream &os, bam1_t &bam) {
		const bam1_core_t *c = &bam.core;
		int &data_len = bam.data_len;
		uint8_t *data = bam.data;
		uint32_t x[8], block_len = data_len + BAM_CORE_SIZE, y;
		int i;
		assert(BAM_CORE_SIZE == 32);
		x[0] = c->tid;
		x[1] = c->pos;
		x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | c->l_qname;
		x[3] = (uint32_t)c->flag<<16 | c->n_cigar;
		x[4] = c->l_qseq;
		x[5] = c->mtid;
		x[6] = c->mpos;
		x[7] = c->isize;

		if (bam_is_be) {
			for (i = 0; i < 8; ++i) bam_endian::bam_swap_endian_4p(x + i);
			y = block_len;
			os.write((char*) bam_endian::bam_swap_endian_4p(&y), 4);
			bam_endian::swap_endian_data(c, data_len, data);
		} else
			os.write((char*) &block_len, 4);

		os.write((char*) x, BAM_CORE_SIZE);
		os.write((char*) data, data_len);
		if (bam_is_be) bam_endian::swap_endian_data(c, data_len, data);
		return os;
	}

	static samfile_t *openSamOrBam(std::string fileName) {
		samfile_t *in = samopen(fileName.c_str(), "rb", 0);
		if (in == NULL || in->header == NULL) {
			//printf("Checking if %s is a SAM formatted file\n", fileName);
			in = samopen(fileName.c_str(), "r", 0);
			if (in == NULL || in->header == NULL) {
				//printf("Error! Failed to open BAM/SAM file %s\n", fileName);
				exit(1);
			}
		}
		return in;
	}
	static void closeSamOrBam(samfile_t *fh) {
		samclose(fh);
	}

	static int readNextBam(samfile_t *ins, bam1_t *read) {
		assert(ins != NULL);
		assert(read != NULL);

		int bytesRead = samread(ins, read);
		return bytesRead;
	}

	static void readBamFile(samfile_t *ins, BamVector &reads) {
		while (true) {
			bam1_t *bam = bam_init1();
			if (readNextBam(ins, bam) >= 0) {
				reads.push_back(bam);
			} else {
				bam_destroy1(bam);
				break;
			}
		}
	}
	static std::string getPairTag(bam1_t *bam) {
		std::string s;
		if ((bam->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
			s = std::string((bam->core.flag & BAM_FREAD1) ? "/1" : "/2" );
		}
		return s;
	}
	static std::string getSequence(bam1_t *bam) {

		if (bam->core.l_qseq) {
			std::stringstream ss;
			uint8_t *s = bam1_seq(bam);
			for (int32_t i = 0; i < bam->core.l_qseq; ++i)
				ss << bam_nt16_rev_table[bam1_seqi(s, i)];
			return ss.str();
		} else {
			return std::string('N', 1);
		}

	}
	static std::string getQualFasta(bam1_t *bam, int offset = 33) {
		uint8_t *t = bam1_qual(bam);
		int32_t len = bam->core.l_qseq;
		if (len == 0)
			return std::string((char) offset, 1);
		if (t[0] == 0xff) {
			return std::string((char) 60+offset, len);
		} else {
			std::stringstream ss;
			for (int32_t i = 0; i < len; ++i)
				ss << (char) t[i] + offset;
			return ss.str();
		}
	}
	static std::ostream &writeFastq(std::ostream &os, bam1_t *bam, int offset = 33) {
		os << '@' << bam1_qname(bam) << getPairTag(bam) << "\n"
				<< getSequence(bam) << "\n+\n"
				<< getQualFasta(bam, offset) << "\n";
		return os;
	}
};

// define write operator for a bam record
std::ostream &operator<<(std::ostream& os, bam1_t &bam) {
	return BamStreamUtils::writeBam(os, bam);
}

// define write operator for a bam header
std::ostream &operator<<(std::ostream& os, bam_header_t &_header) {
	bam_header_t *header = &_header;
	BamStreamUtils::writeBamHeaderPart1(os, header);
	BamStreamUtils::writeBamHeaderPart2(os, header);
	return os;
}
typedef BamStreamUtils::BamVector BamVector;
typedef BamStreamUtils::BamCoreVector BamCoreVector;

class SamUtils {
public:
	static const int READS_PER_BATCH = 32768;
	typedef std::vector< int > IntVector;
	typedef bam1_t *bam1_p;
	typedef bam1_p *bam1_pp;
	typedef std::pair< bam1_p, bam1_pp > bam1_p_pair;
	typedef std::pair< bam1_p, bam1_p > BamPair;
	typedef bam1_core_t *bam1_core_p;

	static int roundup(int val, int factor) {
		int mod = val % factor;
		if (mod > 0)
			val += (factor - mod);
		return val;
	}
	static int *packVector(IntVector &in) {
		int *packed = (int*) calloc(in.size(), sizeof(int));
		for(int i = 0; i < (int) in.size(); i++)
			packed[i] = in[i];
		return packed;
	}
	static IntVector unpackVector(int *in, int size) {
		IntVector unpacked;
		unpacked.reserve(size);
		for(int i = 0; i < size ; i++)
			unpacked.push_back(in[i]);
		return unpacked;
	}

	class SortByPosition {

	public:
		inline bool operator()(const bam1_p_pair a, const bam1_p_pair b) const {
			return cmpBam(a.first, b.first);
		}
		inline bool operator()(const bam1_pp a, const bam1_pp b) const {
			return cmpBam(*a, *b) < 0;
		}
		inline bool operator()(const bam1_p a, const bam1_p b) const {
			return cmpBam(a,b) < 0;
		}
		inline bool operator()(const bam1_core_t &a, const bam1_core_t &b) const {
			return cmpBamCore(a,b) < 0;
		}
		inline bool operator()(const bam1_core_p a, const bam1_core_p b) const {
			return cmpBamCore(*a, *b) < 0;
		}
		static inline int cmpBam(const bam1_p a, const bam1_p b)
		{
			int _cmp = cmpBamCore(a->core, b->core);
			if (_cmp == 0) {
				if (a->data_len == 0) { // let empty sort first
					return -1;
				} else if (b->data_len == 0) {
					return 1;
				} else {
					_cmp = strcmp(bam1_qname(a), bam1_qname(b));
					if (_cmp == 0)
						return isMapped(a) ? -1 : (isMapped(b) ? 1 : (isRead1(a) ? -1 : 1));
					else
						return _cmp;
				}
			} else
				return _cmp;
		}
		static inline int cmpBamCore(const bam1_core_t &a, const bam1_core_t &b) {
			uint64_t ax = ((uint64_t)a.tid<<32|(a.pos+1));
			uint64_t bx = ((uint64_t)b.tid<<32|(b.pos+1));
			if (ax < bx)
				return -1;
			else if (ax == bx)
				return 0;
			else
				return 1;
		}

	};

	static long concatenateOutput(const MPI_Comm &comm, MPI_File &ourFile, long long int myLength, std::istream &data) {
		LOG_VERBOSE_OPTIONAL(1, true, "concatenateOutput(): writing: " << myLength);
		int rank,size;
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &size);

		// exchange sizes
		long long int *ourSizes = (long long int*) calloc(size, sizeof(myLength));
		ourSizes[rank] = myLength;

		if (MPI_SUCCESS != MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ourSizes, 1, MPI_LONG_LONG_INT, comm))
			LOG_THROW("MPI_Allgather() failed");

		long long int totalSize = 0, myOffset = 0, myEnd = 0;;

		if (rank == 0) {
			if (MPI_SUCCESS != MPI_File_get_size(ourFile, &totalSize))
				LOG_THROW("MPI_File_get_size() failed");
		}
		MPI_Bcast(&totalSize, 1, MPI_LONG_LONG_INT, 0, comm);
		LOG_DEBUG_OPTIONAL(2, rank == 0, "ourFile is already: " << totalSize);

		for(int i = 0; i < size; i++) {
			if (i == rank)
				myOffset = totalSize;
			totalSize += ourSizes[i];
			if (i == rank)
				myEnd = totalSize;
		}
		free(ourSizes);
		// MPI_FILE_WRITE
		LOG_DEBUG(2, "myOffset: " << myOffset << " total size will be: " << totalSize << " MySize: " << myLength << " myEnd: " << myEnd);

		if (MPI_SUCCESS != MPI_File_set_size(ourFile, totalSize))
			LOG_THROW("MPI_File_set_size failed");

		MPI_Offset totalWritten = 0;
		int bufSize = 16 * 1024 * 1024; // 16 MB chunks
		char *buf = (char*) calloc(bufSize, 1);
		while(!data.eof() & !data.fail()) {
			data.read(buf, bufSize);
			int count = data.gcount();
			MPI_Status status;
			if (MPI_SUCCESS != MPI_File_write_at(ourFile, myOffset, buf, count, MPI_BYTE, &status))
				LOG_THROW("MPI_File_write_at() failed");

			int writeCount;
			MPI_Get_count(&status, MPI_BYTE, &writeCount);
			if (writeCount != count)
				LOG_THROW("writeCount " << writeCount << " != count " << count);
			myOffset += count;
			totalWritten += count;
			LOG_DEBUG(3, "concatenateOutput(): wrote " << count);
		}
		LOG_DEBUG_OPTIONAL(2, true, "concatenateOutput: completed writing: " << totalWritten << " myOffset: " << myOffset);
		assert(myOffset == myEnd);
		free(buf);
		return totalSize;
	}

	static long writePartialSortedBamVector(const MPI_Comm &comm, MPI_File &ourFile, BamVector &reads, IntVector &sortedCounts, bam_header_t *header = NULL) {
		LOG_VERBOSE_OPTIONAL(1, true, "writePartialSortedBamVector(): with numReads: " << reads.size());
		bool destroyBam = true;
		int rank,size;
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &size);
		assert((int) sortedCounts.size() == size);

		// prepare the heap
		std::vector< bam1_p_pair > heap;
		heap.reserve(size + 1);
		long totalCount = 0;
		for(int i = 0; i < size; i++) {
			if (sortedCounts[i] > 0 && reads[totalCount] != NULL) {
				heap.push_back( bam1_p_pair(reads[totalCount], &reads[totalCount]));
				reads[totalCount] = NULL;
			}
			totalCount += sortedCounts[i];
		}
		assert(totalCount == (int) reads.size());
		reads.push_back(NULL); // terminate condition

		SortByPosition sbp;
		std::make_heap(heap.begin(), heap.end(), sbp);

		// output while processing the heap
		MemoryBuffer myBams;
		{
			MemoryBuffer::ostream os(myBams);
			bgzf_ostream bgzfo(os, rank == size-1);
			if (header != NULL && rank == 0) {
				LOG_DEBUG_OPTIONAL(1, true, "Writing header");
				bgzfo << *header;
			}
			while(!heap.empty()) {
				bam1_p_pair bppair = heap.front();
				bam1_pp bampp = bppair.second;
				bam1_t *bam = bppair.first;
				LOG_DEBUG_OPTIONAL(3, true, bam1_qname(bam));;
				bgzfo << *bam;
				if (destroyBam) {
					bam_destroy1(bam);
				}
				std::pop_heap(heap.begin(), heap.end(), sbp); heap.pop_back();
				bampp++;
				if (*bampp != NULL) {
					bppair = bam1_p_pair(*bampp, bampp);
					*bampp = NULL;
					heap.push_back(bppair); std::push_heap(heap.begin(), heap.end(), sbp);
				}
			}
		}
		long count = 0;
		{
			long long int myLength = myBams.tellp();

			MemoryBuffer::istream is(myBams);
			count = concatenateOutput(comm, ourFile, myLength, is);
		}
		if (! destroyBam )
			reads.pop_back(); // remove terminate condition

		return count;
	}
	static long writeBamVector(const MPI_Comm &comm, MPI_File &ourFile, BamVector &reads, bam_header_t *header = NULL, bool destroyBam = true) {
		LOG_VERBOSE_OPTIONAL(1, true, "writeBamVector(): with numReads: " << reads.size());
		int rank,size;
		MPI_Comm_rank(comm, &rank);
		MPI_Comm_size(comm, &size);
		MemoryBuffer myBams;
		{
			MemoryBuffer::ostream os(myBams);
			bgzf_ostream bgzfo(os, rank == size-1);
			if (header != NULL && rank == 0) {
				LOG_DEBUG_OPTIONAL(1, true, "Writing header");
				bgzfo << *header;
			}
			for(BamVector::iterator it = reads.begin(); it != reads.end(); it++) {
				bam1_t *bam = *it;
				bgzfo << *bam;
				LOG_DEBUG_OPTIONAL(4, true, bam1_qname(bam));
				if (destroyBam) {
					bam_destroy1(bam);
					*it = NULL;
				}
			}
		}
		long count = 0;
		{
			long long int myLength = myBams.tellp();

			MemoryBuffer::istream is(myBams);
			count = concatenateOutput(comm, ourFile, myLength, is);
		}
		return count;
	}
	static long writePartialSortedBamVector(const MPI_Comm &comm, std::string ourFileName, BamVector &reads, IntVector &sortedCounts, bam_header_t *header = NULL) {
		MPI_File outfh;
		if (MPI_SUCCESS != MPI_File_open(comm, (char*) ourFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfh))
			LOG_THROW("Could not open/write: " << ourFileName);

		long val = writePartialSortedBamVector(comm, outfh,reads, sortedCounts, header);

		if (MPI_SUCCESS != MPI_File_close(&outfh))
			LOG_THROW("Could not close: " << ourFileName);
		return val;
	}
	static long writeBamVector(const MPI_Comm &comm, std::string ourFileName, BamVector &reads, bam_header_t *header = NULL, bool destroyBam = true) {
		MPI_File outfh;
		if (MPI_SUCCESS != MPI_File_open(comm, (char*) ourFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfh))
			LOG_THROW("Could not open/write: " << ourFileName);

		long val = writeBamVector(comm, outfh,reads, header,destroyBam);

		if (MPI_SUCCESS != MPI_File_close(&outfh))
			LOG_THROW("Could not close: " << ourFileName);
		return val;
	}

	class TransferStats {
	public:
		IntVector recvRemaining, globalRemaining, sent, remainingToSend, received;
		int myTotalRemainingToSend;
		TransferStats(int size) {
			clear(size);
		}
		void clear(int size) {
			received.resize(size, 0);
			recvRemaining.resize(size, 0);
			globalRemaining.resize(size, 0);
			sent.resize(size, 0);
			remainingToSend.resize(size, 0);
			myTotalRemainingToSend = 0;
		}
		void setMySending(int rank, int _sent, int _remainingToSend) {
			assert(_remainingToSend >= 0);
			assert(_sent >= 0);
			myTotalRemainingToSend += _remainingToSend;
			remainingToSend[rank] = _remainingToSend;
			sent[rank] = _sent;
		}
		void *setSendingBuffer(void *_ptr, int rank) {
			int *ptr = (int*) _ptr;
			*(ptr++) = myTotalRemainingToSend;
			*(ptr++) = remainingToSend[rank];
			*(ptr++) = sent[rank];
			return ptr;
		}
		void setRemainingFromRank(int rank, int _received, int _recvRemaining, int _totalRemaining) {
			received[rank] = _received;
			recvRemaining[rank] = _recvRemaining;
			globalRemaining[rank] = _totalRemaining;
		}
		void *setRemainingFromRank(void *_ptr, int rank) {
			int *ptr = (int*) _ptr;
			globalRemaining[rank] = *(ptr++);
			recvRemaining[rank] = *(ptr++);
			received[rank] = *(ptr++);
			return ptr;
		}
		int getMyTotalRemainingToSend() const {
			return myTotalRemainingToSend;
		}
		int getMyRemainingToSendToRank(int rank) const {
			return remainingToSend[rank];
		}
		int getMySendingToRank(int rank) const {
			return sent[rank];
		}
		int getRemainingFromRank(int rank) const {
			return recvRemaining[rank];
		}
		int getReceivedFromRank(int rank) const {
			return received[rank];
		}
		void *setSendingBuffer(int *ptr, int rank) {
			*(ptr++) = getMyTotalRemainingToSend();
			*(ptr++) = getMyRemainingToSendToRank(rank);
			*(ptr++) = getMySendingToRank(rank);
			return ptr;
		}
		long getGlobalRemaining() const {
			long count = 0;
			for(IntVector::const_iterator it = globalRemaining.begin(); it != globalRemaining.end(); it++)
				count += *it;
			return count;
		}
		bool isDone() const {
			return getGlobalRemaining() == 0;
		}
	};

	class MPIReadExchanger {
	public:
		MPIReadExchanger(const MPI_Comm &_comm)
		: buf1(NULL), buf2(NULL)  {
			myComm = mpi::communicator(_comm, mpi::comm_duplicate);
		}
		~MPIReadExchanger() {
			if (buf1 != NULL) free(buf1);
			if (buf2 != NULL) free(buf2);
		}

		// sendCounts must be known a priori, recvCounts can be specified or be empty if it is unknown by this rank.
		// This is expected to be run until TransferStats.isDone() and will make forward progress populating recvReads and incrementing sendOffsets & recvOffsets
		// _recvCounts will have the final number of reads from each rank after TransferStats.isDone()
		TransferStats exchangeReads(const BamVector &sendReads, IntVector &sendOffsets, IntVector &_sendCounts, BamVector &recvReads, IntVector &recvOffsets, IntVector &_recvCounts, bool copyDataIfUnmapped = true) {
			int rank,size;
			MPI_Comm_rank(myComm, &rank);
			MPI_Comm_size(myComm, &size);
			TransferStats stats(size);

			int maxReadsPerRank = (READS_PER_BATCH + size - 1)/ size;
			IntVector sendBytes(size, 0), sendByteDispl(size, 0), recvBytes(size,0), recvByteDispl(size, 0);
			IntVector sendCounts = _sendCounts, recvCounts = _recvCounts;
			IntVector sendCountDispl(size, 0), recvCountDispl(size, 0);

			int totalToReceive = 0;
			if (_recvCounts.empty()) {

				// Determine the total recvCounts automatically
				LOG_DEBUG_OPTIONAL(3, true, "MPIReadExchanger::exchangeReads(): receiveCounts unknown, exchanging sendCounts");

				_recvCounts.resize(size, 0);
				int *_s = packVector(_sendCounts);
				int *_r = (int*) calloc(size, sizeof(int));
				if (MPI_SUCCESS != MPI_Alltoall(_s, 1, MPI_INT, _r, 1, MPI_INT, myComm))
					LOG_THROW("MPI_Alltoall() failed: ");
				free(_s);
				_recvCounts = unpackVector(_r, size); free(_r);

				recvCounts = _recvCounts;
			}
			for(int i = 0; i < (int) _recvCounts.size(); i++)
				totalToReceive += _recvCounts[i];
			if ((int) recvReads.size() < totalToReceive)
				initBamVector(recvReads, totalToReceive);

			if (Log::isDebug(3)) {
				std::ostringstream ss;
				ss << "_sendCounts ";
				for(int i = 0; i < size; i++)
					ss << _sendCounts[i] << " - " << sendOffsets[i] << ", ";
				ss << ". _recvCounts ";
				for(int i = 0; i < size; i++)
					ss << _recvCounts[i] << " - " << recvOffsets[i] << ", ";

				std::string s = ss.str();
				LOG_DEBUG(3, "exchangeReads(): " << rank << "of" << size << " total blocks " << s);
			}
			assert((int) _sendCounts.size() == size);
			assert((int) _recvCounts.size() == size);
			assert((int) sendOffsets.size() == size);
			assert((int) recvOffsets.size() == size);

			bam1_t *tmpBam;
			int sendReadsIdx, recvReadsIdx;
			int totalSendBytes = 0, totalRecvBytes = 0;
			int totalSendCounts = 0, totalRecvCounts = 0;
			int offsetBytes = sizeof(int) * 3;

			long totalRemainingToSend = 0, totalRemainingToRecv = 0;
			for(int i = 0; i < size; i++) {
				sendCounts[i] = _sendCounts[i] - sendOffsets[i];
				totalRemainingToSend += sendCounts[i];

				if (sendCounts[i] > maxReadsPerRank) {
					if (i != rank)
						sendCounts[i] = maxReadsPerRank;
				}
				sendCountDispl[i] = totalSendCounts + sendOffsets[i];
				totalSendCounts += _sendCounts[i];
				stats.setMySending(i, sendCounts[i], _sendCounts[i] - sendOffsets[i] - sendCounts[i]);

				recvCounts[i] = _recvCounts[i] - recvOffsets[i];
				totalRemainingToRecv += recvCounts[i];
				if (recvCounts[i] > maxReadsPerRank) {
					if (i != rank)
						recvCounts[i] = maxReadsPerRank;
				}
				recvCountDispl[i] = totalRecvCounts + recvOffsets[i];
				totalRecvCounts += _recvCounts[i];

				sendByteDispl[i] = totalSendBytes;
				recvByteDispl[i] = totalRecvBytes;

				if (i == rank) {
					sendBytes[i] = 0;
					recvBytes[i] = 0;
				} else {
					sendBytes[i] += offsetBytes;
					recvBytes[i] += offsetBytes;
					sendBytes[i] += sendCounts[i] * sizeof(bam1_t);

					recvBytes[i] += recvCounts[i] * sizeof(bam1_t);
				}
				totalSendBytes += sendBytes[i];
				totalRecvBytes += recvBytes[i];
			}
			LOG_DEBUG_OPTIONAL(2, true, "MPIReadExchanger::exchangeReads(): remaining to send: " << totalRemainingToSend << " recv: " << totalRemainingToRecv);

			if (Log::isDebug(3)) {
				std::ostringstream ss;
				ss << "sendCounts ";
				for(int i = 0; i < size; i++)
					ss << sendCounts[i] << ", ";
				ss << ". recvCounts ";
				for(int i = 0; i < size; i++)
					ss << recvCounts[i] << ", ";

				std::string s = ss.str();
				LOG_DEBUG(3, "exchangeReads(): exchange blocks " << s);
			}

			IntVector sendCounts2(size, 0);
			IntVector recvCounts2(size, 0);

			// buffer concats 3 ints (totalRemaining, remainingToRank, sentToRank) then the bam1_t structures for each rank
			buf1 = (char*) realloc(buf1, totalSendBytes);
			buf2 = (char*) realloc(buf2, totalRecvBytes);

			tmpBam = (bam1_t*) buf1;

			for(int i = 0; i < size; i++) {
				sendReadsIdx = sendCountDispl[i];
				if (i == rank) {
					recvReadsIdx = recvCountDispl[i];
					LOG_DEBUG_OPTIONAL(3, true, "Copying " << sendCounts[i] << " from sendIdx " << sendReadsIdx << " to recvIdx: " << recvReadsIdx);
					for(int j = 0 ; j < sendCounts[i]; j++) {
						assert(sendReadsIdx < (int) sendReads.size());
						assert(recvReadsIdx < (int) recvReads.size());
						assert(sendReads[sendReadsIdx]->data != NULL);
						assert(recvReads[recvReadsIdx]->data == NULL);
						bam_copy1(recvReads[recvReadsIdx++], sendReads[sendReadsIdx++]);
					}
					continue;
				}

				LOG_DEBUG_OPTIONAL(3, true, "Preping bytes: " << i << " scount: " << sendCountDispl[i] << " soff: " << sendOffsets[i] << " idx: " << sendReadsIdx << " sendReads.size() " << sendReads.size() << " sendCounts: " << sendCounts[i]);

				tmpBam = (bam1_t*) stats.setSendingBuffer(tmpBam, i);

				for(int j = 0; j < sendCounts[i]; j++) {
					assert(sendReadsIdx < (int) sendReads.size());
					memcpy(tmpBam, sendReads[sendReadsIdx], sizeof(bam1_t));
					int dataLen = (copyDataIfUnmapped | isMapped(tmpBam)) ? tmpBam->data_len : 0;
					sendCounts2[i] += dataLen;
					LOG_DEBUG_OPTIONAL(3, true, "will send to rank: " << i << " j: " << j << " sendidx: " << sendReadsIdx << " dataLen: " << dataLen);
					tmpBam++; sendReadsIdx++;
				}
			}

			if (MPI_SUCCESS != MPI_Alltoallv(buf1, &sendBytes[0], &sendByteDispl[0], MPI_BYTE, buf2, &recvBytes[0], &recvByteDispl[0], MPI_BYTE, myComm))
				LOG_THROW("MPI_Alltoallv() failed: ");

			for(int i = 0; i < size; i++) {
				if (i == rank) {
					stats.setRemainingFromRank(i, sendCounts[i], 0, stats.getMyTotalRemainingToSend());
				} else {
					stats.setRemainingFromRank(buf2 + recvByteDispl[i], i);
				}
				assert( recvCounts[i] == stats.getReceivedFromRank(i) );
			}

			totalSendBytes = totalRecvBytes = 0;
			tmpBam = (bam1_t*) buf2;
			for(int i = 0; i < size; i++) {
				sendReadsIdx = sendCountDispl[i];
				recvReadsIdx = recvCountDispl[i];
				if (i == rank) {
					assert(recvBytes[i] == 0);
					assert(sendBytes[i] == 0);
					// copy already happened
				} else {
					tmpBam = (bam1_t*) (((char*)tmpBam) + offsetBytes);
					int dataLen = 0;
					int recvCount = stats.getReceivedFromRank(i);
					for(int j = 0; j < recvCounts[i]; j++) {
						assert(recvReadsIdx < (int) recvReads.size());
						if (j < recvCount) {
							dataLen = (copyDataIfUnmapped | isMapped(tmpBam)) ? tmpBam->data_len : 0;
							recvCounts2[i] += dataLen;

							clearData(tmpBam);
							bam_copy1(recvReads[recvReadsIdx], tmpBam);
							recvReads[recvReadsIdx]->data_len = dataLen; // copy it back for now...
						}
						LOG_DEBUG_OPTIONAL(3, true, "will recv from rank: " << i << " j: " << j << " recvIdx: " << recvReadsIdx << " dataLen: " << dataLen << " isMapped: " << isMapped(tmpBam) << " recvCount: " << recvCount);
						tmpBam++;
						recvReadsIdx++;
					}
				}
				sendByteDispl[i] = totalSendBytes;
				totalSendBytes += sendCounts2[i];
				recvByteDispl[i] = totalRecvBytes;
				totalRecvBytes += recvCounts2[i];
				LOG_DEBUG(3, "rank: bytes " << i << " sendD " << sendByteDispl[i] << " totalSend " << totalSendBytes << " recvD " << recvByteDispl[i] << " totalRecv " << totalRecvBytes);
			}

			buf1 = (char*) realloc(buf1, totalSendBytes);
			char *ptr = buf1;
			int countBytes = 0;
			for(int i = 0; i < size; i++) {
				if (i == rank) {
					continue;
				}
				sendReadsIdx = sendCountDispl[i];
				for(int j = 0; j < sendCounts[i]; j++) {
					assert(sendReadsIdx < (int) sendReads.size());
					if (copyDataIfUnmapped | isMapped(sendReads[sendReadsIdx])) {
						LOG_DEBUG_OPTIONAL(3, true, "packing data for rank: " << i << " sendidx: " << sendReadsIdx << " j " << j << " storing: " << sendReads[sendReadsIdx]->data_len << " countBytes: " << countBytes << " totalSend " << totalSendBytes << " " << sendReads[sendReadsIdx]->data << " " << sendReads[sendReadsIdx]->core.pos);
						memcpy(ptr, sendReads[sendReadsIdx]->data, sendReads[sendReadsIdx]->data_len);
						ptr += sendReads[sendReadsIdx]->data_len;
						countBytes += sendReads[sendReadsIdx]->data_len;
					}
					sendReadsIdx++;
				}
			}
			buf2 = (char*) realloc(buf2, totalRecvBytes);
			if (MPI_SUCCESS != MPI_Alltoallv(buf1, &sendCounts2[0], &sendByteDispl[0], MPI_BYTE, buf2, &recvCounts2[0], &recvByteDispl[0], MPI_BYTE, myComm))
				LOG_THROW("MPI_Alltoallv() failed: ");

			ptr = buf2;
			int ptrOffset = 0;
			for(int i = 0; i < size; i++) {
				if (i == rank) {
					assert(recvCounts2[i] == 0);
					continue;
				}
				recvReadsIdx = recvCountDispl[i];
				for(int j = 0; j < recvCounts[i]; j++) {
					assert(recvReadsIdx < (int) recvReads.size());
					int dataLen = recvReads[recvReadsIdx]->data_len;
					LOG_DEBUG_OPTIONAL(3, ptrOffset < totalRecvBytes, "expanding data for rank " << i << " j " << j << " recvIdx: " << recvReadsIdx << " datalen: " << dataLen << " ptrOffset: " << ptrOffset << " totalRecvBytes: " << totalRecvBytes << " " << ptr << " " << recvReads[recvReadsIdx]->core.pos);
					if (dataLen > 0) {
						assert(ptrOffset + dataLen <= totalRecvBytes);
						assert(copyDataIfUnmapped | isMapped(recvReads[recvReadsIdx]));
						int m_data = roundup(dataLen,8);
						if (m_data > recvReads[recvReadsIdx]->m_data) {
							recvReads[recvReadsIdx]->m_data = m_data;
							recvReads[recvReadsIdx]->data = (uint8_t*) realloc(recvReads[recvReadsIdx]->data, m_data);
						}
						memcpy(recvReads[recvReadsIdx]->data, ptr, dataLen);
					}
					ptr += dataLen;
					ptrOffset += dataLen;
					recvReadsIdx++;
				}
			}

			for(int i = 0; i < size; i++) {
				sendOffsets[i] += sendCounts[i];
				recvOffsets[i] += recvCounts[i];
			}

			LOG_DEBUG_OPTIONAL(2, rank == 0, "exchangeReads() Finished.  Remaining: " << stats.getGlobalRemaining());
			return stats;
		}
	private:
		mpi::communicator myComm;
		char *buf1, *buf2;
	};

	class MPIMergeSam {
		// This class expects to be run with a small MPI communicator with the same number
		// of ranks as there are version of the inputFile.

		// Merges 2+ streams of the same read aligned to separate partitions of the index
	public:

		MPIMergeSam(const MPI_Comm &_comm, std::string _inputFile, BamVector &_reads)
		: inputFile(_inputFile), myGlobalHeaderOffset(0), myReads(_reads), readExchanger(_comm){
			myComm = mpi::communicator(_comm, mpi::comm_duplicate),
			fh = BamStreamUtils::openSamOrBam(inputFile);
			setHeaderCounts();
			pickBestAlignments();
		}
		~MPIMergeSam() {
			LOG_DEBUG_OPTIONAL(2, true, "~MPIMergeSam()");
			samclose(fh);
			free(headerCounts);
			destroyBamVector(tmpReads);
		}
		bam_header_t *getHeader() {
			return fh->header;
		}
		void syncAndPick(BamVector &reads) {
			BamVector pickedReads = syncAndPick(reads, reads.size());
			destroyBamVector(reads);
			reads.swap(pickedReads);
		}
		// all processes write the combined BGZF compressed header to outputFile
		// returns the fileoffset for the end of the header
		long outputMergedHeader(std::string ourFileName) {
			LOG_VERBOSE_OPTIONAL(1, true, "outputMergedHeader()");
			MPI_File outfh;
			if (MPI_SUCCESS != MPI_File_open(myComm, (char*) ourFileName.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outfh))
				LOG_THROW("Could not open/write: " << ourFileName);

			long val = outputMergedHeader(outfh);

			if (MPI_SUCCESS != MPI_File_close(&outfh))
				LOG_THROW("Could not close: " << ourFileName);
			LOG_DEBUG_OPTIONAL(2, true, "Finished outputMergedHeader(): " << val);
			return val;

		}
		long outputMergedHeader(MPI_File &ourFile) {
			int rank,size;
			MPI_Comm_rank(myComm, &rank);
			MPI_Comm_size(myComm, &size);

			bam_header_t *header = fh->header;
			// fix @SQ lines within plaintext header->text
			{
				int root = 0;
				//std::string s(header->text, header->l_text);
				char *tmp, *endSQ = NULL, *firstSQ = strstr(header->text, "@SQ");
				if (firstSQ == NULL) LOG_THROW("Could not find @SQ in header: " << header->text);
				tmp = firstSQ + 1; 
				while ((tmp = strstr(tmp, "@SQ")) != NULL) {
					endSQ = tmp;
					tmp++;
				}
				if (endSQ == NULL) LOG_THROW("Could not find last @SQ record");
				tmp = strstr(endSQ + 1, "@");
				if (tmp == NULL) LOG_THROW("Could not find any header after @SQ section");
				endSQ = tmp;

				int myCount = endSQ - firstSQ;
				int *counts = NULL;
				int *displs = NULL;
				char *buf = NULL;
				if (rank == root) {
					counts = (int*) calloc(size, sizeof(myCount));
					displs = (int*) calloc(size, sizeof(myCount));
				}
				LOG_DEBUG(2, "Sending " << myCount << " bytes of @SQ text to root: " << (firstSQ - header->text) << " to " << (endSQ - header->text));

				if (MPI_SUCCESS != MPI_Gather(&myCount, 1, MPI_INT, counts, 1, MPI_INT, root, myComm))
					LOG_THROW("MPI_Gather() failed");

				int total = 0;
				if (rank == root) {
					for(int i = 0; i < size; i++) {
						displs[i] = total;
						total += counts[i];
					}
					buf = (char*) calloc(total, 1);
				}

				LOG_DEBUG(2, "Sending @SQ lines:\n" << std::string(firstSQ, myCount));
				if (MPI_SUCCESS != MPI_Gatherv(firstSQ, myCount, MPI_BYTE, buf, counts, displs, MPI_BYTE, root, myComm))
					LOG_THROW("MPI_Gatherv failed");

				if (rank == root) {
					free(counts);
					free(displs);
					std::string s = std::string(header->text, firstSQ - header->text);
					s += std::string(buf, total);
					free(buf);
					s += std::string(endSQ, header->l_text - (endSQ - header->text));
					header->text = (char*) realloc(header->text, s.size());
					header->l_text = s.size();
					s.copy(header->text, s.size());
					LOG_DEBUG_OPTIONAL(2, true, "New header is " << s.length() << "\n" << s);
				}
			}

			// build (bgzf compressed) header from own reads
			std::stringstream myHeaderPart;
			{
				bgzf_ostream bgzfo(myHeaderPart, false);
				if (rank == 0) {
					// set total contigs
					int32_t totalContigs = 0;
					for(int i = 0 ; i < size; i++)
						totalContigs += headerCounts[i];
					LOG_DEBUG_OPTIONAL(2, true, "Writing header.  textSize: " << header->l_text << " contigs: " << totalContigs);
					BamStreamUtils::writeBamHeaderPart1(bgzfo, header->l_text, header->text, totalContigs);
				}
				BamStreamUtils::writeBamHeaderPart2(bgzfo, header);
			}
			bam_header_destroy(header);
			fh->header = NULL;

			int mySize = myHeaderPart.tellg();
			LOG_DEBUG_OPTIONAL(2, true, "my bgzf compressed header size: " << mySize);

			return concatenateOutput(myComm, ourFile, myHeaderPart.tellp(), myHeaderPart);
		}

		int getMyGlobalHeaderOffset() const {
			return myGlobalHeaderOffset;
		}
	protected:
		void setHeaderCounts() {
			int rank,size;
			MPI_Comm_rank(myComm, &rank);
			MPI_Comm_size(myComm, &size);
			headerCounts = (int*) calloc(size, sizeof(int));
			headerCounts[rank] = fh->header->n_targets;

			if (MPI_SUCCESS != MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, headerCounts, 1, MPI_INT, myComm))
				LOG_THROW("MPI_Allgather() failed");

			for(int i = 0; i < rank; i++)
				myGlobalHeaderOffset += headerCounts[i];
			int totalHeaderCounts = myGlobalHeaderOffset;
			for(int i = rank; i < size; i++)
				totalHeaderCounts += headerCounts[i];
			LOG_DEBUG_OPTIONAL(1, true, "setHeaderCounts(): globalHeaderOffset: " << myGlobalHeaderOffset << " my n_targets: " << fh->header->n_targets << " totalContigs: " << totalHeaderCounts);
		}
		void pickBestAlignments() {
			LOG_DEBUG_OPTIONAL(2, true, "pickBestAlignment()");
			BamVector batchReads;
			initBamVector(batchReads, READS_PER_BATCH);

			long totalBytes = 0;
			int count = 0;
			while (count < READS_PER_BATCH) {
				bam1_t *tmpBam = batchReads[count];
				int bytesRead = BamStreamUtils::readNextBam(fh, tmpBam);
				if (bytesRead < 0)
					break;
				fixtids(tmpBam);
				totalBytes += bytesRead;
				count++;
				if (count == READS_PER_BATCH) {
					BamVector v = syncAndPick(batchReads, count);
					myReads.insert(myReads.end(), v.begin(), v.end());
					count = 0;
				}
			}
			if (count > 0) {
				BamVector v = syncAndPick(batchReads, count);
				myReads.insert(myReads.end(), v.begin(), v.end());
			}

			destroyBamVector(batchReads);
		}
		void fixtids(bam1_t *tmpBam) {
			bam1_core_t &c = tmpBam->core;
			if (c.tid >= 0)
				c.tid += myGlobalHeaderOffset;
			if (c.mtid >= 0)
				c.mtid += myGlobalHeaderOffset;
		}
		void fixMates(bam1_t *pre, bam1_t *cur) {
			// derived from samtools 0.1.18 bam_mate.c
			if (strcmp(bam1_qname(cur), bam1_qname(pre)) == 0) { // identical pair name
				cur->core.mtid = pre->core.tid; cur->core.mpos = pre->core.pos;
				pre->core.mtid = cur->core.tid; pre->core.mpos = cur->core.pos;

				// fix FMUNMAP & FMREVERSE for both reads
				// Added to ensure each read could be derived from separate index partitions
				// where one mate has no information on how the other mate may or may not map
				if ((pre->core.flag & BAM_FUNMAP) == BAM_FUNMAP)
					cur->core.flag |= BAM_FMUNMAP;
				else
					cur->core.flag &= ~BAM_FMUNMAP;

				if ((pre->core.flag & BAM_FREVERSE) == BAM_FREVERSE)
					cur->core.flag |= BAM_FMREVERSE;
				else
					cur->core.flag &= ~BAM_FMREVERSE;

				if ((cur->core.flag & BAM_FUNMAP) == BAM_FUNMAP)
					pre->core.flag |= BAM_FMUNMAP;
				else
					pre->core.flag &= ~BAM_FMUNMAP;

				if ((cur->core.flag & BAM_FREVERSE) == BAM_FREVERSE)
					pre->core.flag |= BAM_FMREVERSE;
				else
					pre->core.flag &= ~BAM_FMREVERSE;

				if (pre->core.tid == cur->core.tid && !(cur->core.flag&(BAM_FUNMAP|BAM_FMUNMAP))
						&& !(pre->core.flag&(BAM_FUNMAP|BAM_FMUNMAP)))
				{
					uint32_t cur5, pre5;
					cur5 = (cur->core.flag&BAM_FREVERSE)? bam_calend(&cur->core, bam1_cigar(cur)) : cur->core.pos;
					pre5 = (pre->core.flag&BAM_FREVERSE)? bam_calend(&pre->core, bam1_cigar(pre)) : pre->core.pos;
					cur->core.isize = pre5 - cur5; pre->core.isize = cur5 - pre5;
				} else cur->core.isize = pre->core.isize = 0;
				if (pre->core.flag&BAM_FREVERSE) cur->core.flag |= BAM_FMREVERSE;
				else cur->core.flag &= ~BAM_FMREVERSE;
				if (cur->core.flag&BAM_FREVERSE) pre->core.flag |= BAM_FMREVERSE;
				else pre->core.flag &= ~BAM_FMREVERSE;
				if (cur->core.flag & BAM_FUNMAP) { pre->core.flag |= BAM_FMUNMAP; pre->core.flag &= ~BAM_FPROPER_PAIR; }
				if (pre->core.flag & BAM_FUNMAP) { cur->core.flag |= BAM_FMUNMAP; cur->core.flag &= ~BAM_FPROPER_PAIR; }

			} else { // unpaired or singleton
				pre->core.mtid = -1; pre->core.mpos = -1; pre->core.isize = 0;
				if (pre->core.flag & BAM_FPAIRED) {
					pre->core.flag |= BAM_FMUNMAP;
					pre->core.flag &= ~BAM_FMREVERSE & ~BAM_FPROPER_PAIR;
				}
			}
		}
		void fixMates(BamVector &reads) {
			bool hasPrev = false;
			for(int i = 0; i < (int) reads.size(); i++) {
				if (isPaired(reads[i])) {
					if (hasPrev) {
						fixMates(reads[i-1], reads[i]);
						hasPrev = false;
					} else {
						hasPrev = true;
					}
				} else {
					hasPrev = false;
				}
			}
		}
		BamVector syncAndPick(const BamVector &batchReads, int batchSize) {
			LOG_DEBUG(3, "syncAndPick(): " << batchReads.size() << " batchSize: " << batchSize);
			int rank,size;
			MPI_Comm_rank(myComm, &rank);
			MPI_Comm_size(myComm, &size);
			assert(batchSize <= (int) batchReads.size());

			BamVector pickedReads;
			IntVector sendCounts(size,0), sendOffsets(size,0), recvCounts, recvOffsets(size,0);
			int blockSize = (batchSize + size - 1) / size;
			if ((blockSize & 0x1) == 1)
				blockSize++; // make sure pairs are kept together
			int count = 0;
			int myCount = 0;
			for(int i = 0; i < size; i++) {
				int c = std::min(blockSize, (int) (batchSize - count));
				sendCounts[i] = c;
				if (c > 0) {
					assert(batchReads[count] != NULL);
					assert(batchReads[count]->data != NULL);
					assert(batchReads[count+c-1] != NULL);
					assert(batchReads[count+c-1]->data != NULL);
				}
				if (i == rank) {
					myCount = c;
					recvCounts.resize(size, myCount);
					tmpReads.reserve(myCount * size);
				}
				count += c;
			}
			assert(count == batchSize);
			while (! readExchanger.exchangeReads(batchReads, sendOffsets, sendCounts, tmpReads, recvOffsets, recvCounts, true).isDone() );
			assert((int) tmpReads.size() >= myCount * size);
			pickedReads.reserve(myCount);

			BamVector reads1, reads2;
			reads1.resize(size, NULL);
			reads2.resize(size, NULL);
			for(int i = 0; i < myCount; i++) {
				bool _isPaired = false;
				for(int j = 0; j < size; j++) {
					bam1_t *x =  tmpReads[ j * myCount + i];
					if (isRead1(x))
						reads1[j] = x;
					else
						reads2[j] = x;
					if (!_isPaired && isPaired(x))
						_isPaired = true;
				}
				if (_isPaired) {
					i++;
					for(int j = 0; j < size; j++) {
						bam1_t *x =  tmpReads[ j * myCount + i];
						if (isRead1(x))
							reads1[j] = x;
						else
							reads2[j] = x;
					}
					BamPair best = pickBestReadPair(reads1, reads2);
					assert(strcmp(bam1_qname(best.first), bam1_qname(best.second)) == 0);

					pickedReads.push_back( bam_dup1(best.first) );
					pickedReads.push_back( bam_dup1(best.second) );
				} else {
					pickedReads.push_back( bam_dup1( pickBestRead(reads1) ) );
				}
			}

			fixMates(pickedReads);
			return pickedReads;

		}
	private:
		mpi::communicator myComm;
		std::string inputFile;
		int myGlobalHeaderOffset, *headerCounts;
		samfile_t *fh;
		BamVector tmpReads;
		BamVector &myReads;
		MPIReadExchanger readExchanger;

	};

	class MPISortBam {
		// A single, coordinated sorted outputFile will be constructed.
	public:

		MPISortBam(const MPI_Comm &_comm, BamVector &_reads)
		: myReads(_reads), readExchanger(_comm) {
			myComm = mpi::communicator(_comm, mpi::comm_duplicate);
			sortGlobal();
		}
		MPISortBam(const MPI_Comm &_comm, BamVector &_reads, std::string outputFileName, bam_header_t *header)
		: myReads(_reads), readExchanger(_comm) {
			myComm = mpi::communicator(_comm, mpi::comm_duplicate),
			sortGlobal(outputFileName, header);
		}
		~MPISortBam() {
		}

		static void copyBamCore(bam1_core_t *dest, const bam1_core_t *src) {
			memcpy(dest, src, sizeof(bam1_core_t));
		}

		static BamCoreVector calculateGlobalPartitions(BamVector &sortedReads, const MPI_Comm &mpicomm, int granularity = -1) {
			LOG_DEBUG_OPTIONAL(2, true, "calculateGlobalPartitions(): " << sortedReads.size());
			int myRank, ourSize;
			MPI_Comm_rank(mpicomm, &myRank);
			MPI_Comm_size(mpicomm, &ourSize);
			if (ourSize == 1)
				return BamCoreVector();

			if (granularity == -1)
				granularity = (ourSize-1)*11;

			// find first unmapped read
			long myUnmapped, myTotal = sortedReads.size();
			if (sortedReads.empty()) {
				myUnmapped = 0;
			} else {
				BamVector::const_iterator it = std::lower_bound(sortedReads.begin(), sortedReads.end(), nullRead(), SortByPosition());
				myUnmapped = sortedReads.end() - it;
			}

			long *ourDataSizes = (long*) calloc(ourSize*2, sizeof(myUnmapped));

			ourDataSizes[myRank*2] = myTotal;
			ourDataSizes[myRank*2+1] = myUnmapped;
			long myMapped = myTotal - myUnmapped;
			LOG_DEBUG_OPTIONAL(2, true, "myTotal: " << ourDataSizes[myRank*2] << " myunmapped: " << ourDataSizes[myRank*2+1] << " myMapped: " << myMapped);

			if (MPI_SUCCESS != MPI_Allgather(MPI_IN_PLACE, 0, MPI_LONG, ourDataSizes, 2, MPI_LONG, mpicomm))
				LOG_THROW("MPI_Allgather() failed");

			long ourTotalSize = 0, ourTotalUnmapped = 0;
			for(int i = 0; i < ourSize; i++) {
				ourTotalSize += ourDataSizes[i*2];
				ourTotalUnmapped += ourDataSizes[i*2+1];
			}

			int numSamples = 0, mySampleOffset = 0;
			int *ourGranularity = (int*) calloc(ourSize, sizeof(granularity));
			for (int i = 0; i < ourSize; i++) {
				ourGranularity[i] = (ourSize * ourDataSizes[i*2] * granularity + ourTotalSize - 1) / ourTotalSize ;
				ourGranularity[i] = roundup(ourGranularity[i],(ourSize-1)); // add extra samples to make a multiple of (ourSize-1)
				numSamples += ourGranularity[i];
				if (i < myRank)
					mySampleOffset = numSamples;
				LOG_DEBUG_OPTIONAL(3, true, "granularity for rank: " << i << ourGranularity[i] << " ourSizes[i] " << ourDataSizes[i]);
			}
			free(ourDataSizes);
			LOG_DEBUG_OPTIONAL(2, true, "myGranularity: "<< ourGranularity[myRank] << " / " << numSamples << " myTotal: " << myTotal << " ourTotal: " << ourTotalSize << " ourUnmapped: " << ourTotalUnmapped);

			bam1_core_t *bamCoreBuf = (bam1_core_t*) calloc(numSamples, sizeof(bam1_core_t));
			long partitionIdx[ ourGranularity[myRank] ];
			if (ourGranularity[myRank] > 0) {
				long blockSize = myTotal / (ourGranularity[myRank]+1);

				// find rank partitions of rank-tile break points
				// store a buffer containing all the percentile ranking cores

				for (int i = 0 ; i < ourGranularity[myRank]; i++) {
					long idx = std::min((long) (myTotal-1), (long)blockSize*(i+1));
					copyBamCore((bamCoreBuf + mySampleOffset + i), (idx >= 0) ? &(sortedReads[idx]->core) : nullCore());
					partitionIdx[i] = idx;
				}
			}
			if (Log::isDebug(2)) {
				std::stringstream ss;
				ss << "My partitions: ";
				for(int i = mySampleOffset; i < mySampleOffset + ourGranularity[myRank]; i++)
					ss << bamCoreBuf[i].tid << "," << bamCoreBuf[i].pos << " (" << partitionIdx[i-mySampleOffset] << ") ";
				ss << " unmapped: " << myUnmapped;
				std::string s = ss.str();
				LOG_DEBUG_OPTIONAL(2, true, s);
			}
			// fix ourGranularity to be bytes sent

			int totaldispl = 0;
			int *displ = (int*) calloc(ourSize, sizeof(totaldispl));
			for(int i = 0; i < ourSize; i++) {
				ourGranularity[i] *= sizeof(bam1_core_t);
				displ[i] = totaldispl;
				totaldispl += ourGranularity[i];
			}

			if (MPI_SUCCESS != MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_BYTE, bamCoreBuf, ourGranularity, displ, MPI_BYTE, mpicomm))
				LOG_THROW("MPI_Allgather() failed");
			free(ourGranularity);
			free(displ);

			BamCoreVector partitions;
			partitions.resize(ourSize-1);

			bam1_core_t **tmp = (bam1_core_t**) calloc(numSamples, sizeof(bam1_core_t*));
			for(int i = 0; i < numSamples; i++) {
				*(tmp+i) = bamCoreBuf + i;
			}
			std::sort(tmp, tmp+numSamples, SortByPosition());

			int blockSize = numSamples / ourSize;
			if (blockSize == 0)
				blockSize = 1;
			for(int i = 0; i < ourSize - 1; i++) {
				int idx = (i+1) * blockSize - 1;
				bam1_core_t *b;
				if (idx >= numSamples)
					b = nullCore();
				else
					b = *(tmp+idx);
				LOG_DEBUG_OPTIONAL(3, true, "Chose partition for i: " << i << " idx: " << idx << " blocksize: " << blockSize << " numSamples: " << numSamples);
				copyBamCore(&partitions[i], b);
			}

			free(tmp);
			free(bamCoreBuf);
			if (Log::isDebug(2)) {
				std::stringstream ss;
				ss << "partitions: ";
				for(int i = 0; i < ourSize - 1; i++)
					ss << partitions[i].tid << "," << partitions[i].pos << " ";
				std::string s = ss.str();
				LOG_DEBUG_OPTIONAL(2, myRank==0, s);
			}

			return partitions;
		}

		void sortGlobal() {
			_sortGlobal();

			// sortLocal
			sortLocal(myReads);
		}
		void sortGlobal(std::string outputFileName, bam_header_t *header) {
			IntVector sortedCounts = _sortGlobal();
			writePartialSortedBamVector(myComm, outputFileName, myReads, sortedCounts, header);
		}

	protected:
		IntVector _sortGlobal() {
			LOG_DEBUG_OPTIONAL(1, true, "sortGlobal(): " << myReads.size());
			int myRank, ourSize;
			MPI_Comm_rank(myComm, &myRank);
			MPI_Comm_size(myComm, &ourSize);

			sortLocal(myReads);

			if (ourSize == 1) {
				IntVector sortedCounts;
				sortedCounts.push_back(myReads.size());
				return sortedCounts;
			}

			BamCoreVector partitions = calculateGlobalPartitions(myReads, myComm);

			// calculate full read counts
			IntVector sendCounts(ourSize, 0);
			int totalSendCounts = 0;

			bam1_t *b = bam_init1();
			BamVector::iterator start = myReads.begin(), it;
			for(int i = 0 ; i < ourSize; i++) {
				if (i == ourSize - 1) {
					sendCounts[i] = myReads.end() - start;
				} else {
					if (start != myReads.end()) {
						if (partitions[i].tid >= 0) {
							memcpy(&(b->core), &partitions[i], sizeof(partitions[i]));
							it = std::upper_bound(start, myReads.end(), b, SortByPosition());
							sendCounts[i] = (it - start);
						} else {
							it = std::lower_bound(start, myReads.end(), nullRead(), SortByPosition());
							if (it > start) {
								sendCounts[i] = (it - start);
								start = it;
							}
							int remainingRanks = ourSize - i;
							int remainingReads = (myReads.end() - start + remainingRanks - 1) / remainingRanks;
							sendCounts[i] += remainingReads;
							it = start + remainingReads;
						}
					} else {
						sendCounts[i] = 0;
					}
				}
				totalSendCounts += sendCounts[i];
				start = it;
			}
			bam_destroy1(b);

			assert(totalSendCounts == (int) myReads.size());

			if (Log::isDebug(2)) {
				std::stringstream ss;
				ss << "My sorted exchange: " << std::endl;
				int x = 0;
				for(int i = 0; i < ourSize; i++) {
					ss << "torank: " << i << " " << sendCounts[i];
					if (sendCounts[i] > 0 ) {
						ss << " first: " << myReads[x]->data << " " << myReads[x]->core.pos;
						x += sendCounts[i];
						ss << " last: " << myReads[x-1]->data << " " << myReads[x-1]->core.pos;
					}
					ss << std::endl;
				}
				std::string s = ss.str();
				LOG_DEBUG_OPTIONAL(2, true, s);
			}
			// exchange reads (in batches)
			IntVector sendOffsets(ourSize, 0), recvOffsets(ourSize, 0), recvCounts;
			BamVector globalReads;
			while (! readExchanger.exchangeReads(myReads, sendOffsets, sendCounts, globalReads, recvOffsets, recvCounts, true).isDone() );
			std::swap(myReads, globalReads);
			destroyBamVector(globalReads);
			LOG_DEBUG_OPTIONAL(1, true, "sortGlobal() finished: " << myReads.size());
			return recvCounts;
		}

	private:
		mpi::communicator myComm;
		BamVector &myReads;
		MPIReadExchanger readExchanger;

	};

	// in-place sorts a BamVector
	static void sortLocal(BamVector &reads) {
		LOG_DEBUG_OPTIONAL(1, true, "sortLocal(): " << reads.size());
		std::sort(reads.begin(), reads.end(), SortByPosition());
		LOG_DEBUG_OPTIONAL(2, true, "sortLocal(): finished");
	}

	// takes all unmapped reads in reads and puts them in unmappedReads
	// migrated reads will have NULL in reads (use collapseVector to fix)
	static void splitUnmapped(BamVector &reads, BamVector &unmappedReads, bool onlyUnmappedPairs = false) {
		if (onlyUnmappedPairs) {
			int testPairedOrdering = 0;
			const int MAX_TEST = 1000;
			for(long i = 0; i < (long) reads.size()-1; i+=2) {
				if (testPairedOrdering++ < MAX_TEST) {
					if (strcmp(bam1_qname(reads[i]), bam1_qname(reads[i+1])) != 0) {
						LOG_WARN(1, "splitUnmapped(): Detected bam not sorted by name.  Aborting by onlyUnmappedPairs");
						break;
					}
				}
				if ((!isMapped(reads[i])) && (!isMapped(reads[i+1]))) {
					unmappedReads.push_back(reads[i]);
					unmappedReads.push_back(reads[i+1]);
					reads[i] = NULL;
					reads[i+1] = NULL;
				}
			}
		} else {
			for(long i = 0; i < (long) reads.size(); i++) {
				if (reads[i] != NULL && !isMapped(reads[i])) {
					unmappedReads.push_back(reads[i]);
					reads[i] = NULL;
				}
			}
		}
	}

	// removes NULL entries from reads
	static void collapseVector(BamVector &reads) {
		long countBams = 0;
		for(long i = 0 ; i < (long) reads.size(); i++) {
			if (reads[i] != NULL) {
				if (i != countBams)
					reads[countBams] = reads[i];
				countBams++;
			}
		}
		reads.resize(countBams);
	}
	static void dumpBamVector(int debugLevel, BamVector &reads, std::string msg) {
		if (Log::isDebug(debugLevel)) {
			std::stringstream ss;
			ss << msg << std::endl;
			int count = 0;
			for(BamVector::iterator it = reads.begin(); it != reads.end(); it++) {
				if (*it == NULL)
					ss << "null" << std::endl;
				else
					ss << (*it)->core.tid << "," << (*it)->core.pos << " " << (*it)->data << " " << (void*)(*it)->data << std::endl;
			    if (count++ > 100)
			    	break;
			}
			std::string s = ss.str();
			LOG_DEBUG_OPTIONAL(debugLevel, true, s);
		}
	}


	static bam1_core_t *nullCore() {
		static bam1_core_t *x = NULL;
		if (x == NULL) {
			x = (bam1_core_t*) calloc(sizeof(bam1_core_t), 1);
			x->tid=-1;
			x->pos=-1;
			x->bin=0;
			x->qual=0;
			x->l_qname=0;
			x->n_cigar=0;
			x->flag=4;
			x->l_qseq=0;
			x->mtid=-1;
			x->mpos=-1;
			x->isize=0;
		}
		return x;
	}
	static bam1_t *nullRead() {
		static bam1_t *x = NULL;
		if (x == NULL) {
			x = (bam1_t*) calloc(1, sizeof(bam1_t));
			memcpy(&x->core, nullCore(), sizeof(bam1_core_t));
		}
		return x;
	}

	static bam1_t *pickBestRead( BamVector &reads, bool forceSingle = false ) {
		assert(reads.size() > 0);
		bam1_t *best = reads[0];
		for(int j = 1; j < (int) reads.size() ; j++) {
			best = isBetter(best, reads[j], forceSingle);
		}
		return best;
	}

	static BamPair pickBestReadPair( BamVector &reads1, BamVector &reads2 ) {
		assert(reads1.size() > 0);
		assert(reads2.size() == reads1.size());
		BamPair best(reads1[0], reads2[0]);
		// assess reads as Pairs first
		for(int j = 1; j < (int) reads1.size(); j++) {
			best = isBetter(best, BamPair(reads1[j], reads2[j]));
		}
		if ((isProperPair(best.first) && isProperPair(best.second)) || best.first->core.tid == best.second->core.tid)
			return best;
		// no pairing is proveably best.  Get best individual reads (i.e. chimers)
		best.first = pickBestRead(reads1, true);
		best.second= pickBestRead(reads2, true);
		return best;
	}

	static BamPair isBetter(BamPair a, BamPair b) {
		static char *unmapped = (char*) "unmapped";
		int ranka = getRankedState(a.first) + getRankedState(a.second);
		int rankb = getRankedState(b.first) + getRankedState(b.second);
		if (ranka == rankb) {
			// Dive deeper into AS (Alignment Score generated by aligner)
			// less reliable as partitioned aligner did not see all the data..
			ranka += bam_aux2i(bam_aux_get(a.first, "AS"));
			ranka += bam_aux2i(bam_aux_get(a.second, "AS"));
			rankb += bam_aux2i(bam_aux_get(b.first, "AS"));
			rankb += bam_aux2i(bam_aux_get(b.second, "AS"));
		}
		BamPair *better = ranka < rankb ? &a : &b;
		if (Log::isDebug(2)) {
			char *a1str = a.first->data != NULL ? bam_format1(NULL, a.first) : unmapped;
			char *a2str = a.second->data != NULL ? bam_format1(NULL, a.second) : unmapped;
			char *b1str = b.first->data != NULL ? bam_format1(NULL, b.first) : unmapped;
			char *b2str = b.second->data != NULL ? bam_format1(NULL, b.second) : unmapped;
			LOG_DEBUG_OPTIONAL(2, true, "Selected " << (better == &a ? "a " : "b ") << ranka << " " << rankb  << "\na1 " << a1str << "\na2 " << a2str << "\nb1 " << b1str << "\nb2 " << b2str);
			if (a1str != unmapped)
				free(a1str);
			if (b1str != unmapped)
				free(b1str);
			if (a2str != unmapped)
				free(a2str);
			if (b2str != unmapped)
				free(b2str);
		}
		return *better;
	}

	static bam1_t *isBetter(bam1_t *a, bam1_t *b, bool forceSingle = false) {
		static char *unmapped = (char*) "unmapped";
		int ranka = getRankedState(a, forceSingle);
		int rankb = getRankedState(b, forceSingle);
		if (ranka == rankb) {
			// Dive deeper into AS (Alignment Score generated by aligner)
			ranka += bam_aux2i(bam_aux_get(a, "AS"));
			rankb += bam_aux2i(bam_aux_get(b, "AS"));
		}

		bam1_t *better = ranka < rankb ? a : b;
		if (Log::isDebug(2)) {
			char *astr = a->data != NULL ? bam_format1(NULL, a) : unmapped;
			char *bstr = b->data != NULL ? bam_format1(NULL, b) : unmapped;
			LOG_DEBUG_OPTIONAL(2, true, "Selected " << (better == a ? "a " : "b ") << ranka << " " << rankb << "\n" << astr << "\n" << bstr);
			if (astr != unmapped)
				free(astr);
			if (bstr != unmapped)
				free(bstr);
		}
		return better;
	}

	static inline bool isPaired(const bam1_t *b) {
		return (b->core.flag & BAM_FPAIRED) == BAM_FPAIRED;
	}
	static inline bool isProperPair(const bam1_t *b) {
		return (b->core.flag & BAM_FPROPER_PAIR) == BAM_FPROPER_PAIR;
	}
	static inline bool isUnMapped(const bam1_t *b) {
		return (b->core.flag & BAM_FUNMAP) == BAM_FUNMAP;
	}
	static inline bool isMapped(const bam1_t *b) {
		return (b->core.flag & BAM_FUNMAP) != BAM_FUNMAP;
	}

	static inline bool isMateUnMapped(const bam1_t *b) {
		return (b->core.flag & BAM_FMUNMAP) == BAM_FMUNMAP;
	}
	static inline bool isMateMapped(const bam1_t *b) {
		return (b->core.flag & BAM_FMUNMAP) != BAM_FMUNMAP;
	}

	static inline bool isReverse(const bam1_t *b) {
		return (b->core.flag & BAM_FREVERSE) == BAM_FREVERSE;
	}
	static inline bool isMateReverse(const bam1_t *b) {
		return (b->core.flag & BAM_FMREVERSE) == BAM_FMREVERSE;
	}
	static inline bool isRead1(const bam1_t *b) {
		return (b->core.flag & BAM_FREAD1) == BAM_FREAD1;
	}
	static inline bool isRead2(const bam1_t *b) {
		return (b->core.flag & BAM_FREAD2) == BAM_FREAD2;
	}
	static inline bool isSecondary(const bam1_t *b) {
		return (b->core.flag & BAM_FSECONDARY) == BAM_FSECONDARY;
	}
	static inline bool isQCFail(const bam1_t *b) {
		return (b->core.flag & BAM_FQCFAIL) == BAM_FQCFAIL;
	}
	static inline bool isDup(const bam1_t *b) {
		return (b->core.flag & BAM_FDUP) == BAM_FDUP;
	}

	// getRankedState used for comparison whether one alignment is "better" than another.  Lower is better
	static inline int getRankedState(const bam1_t *bam, bool forceSingle = false) {
		if (isPaired(bam) && ! forceSingle)
			return getRankedPairedState(bam);
		else
			return getRankedSingleState(bam);
	}
	static inline int getRankedSingleState(const bam1_t *bam) {
		// %aligned + 100 - 100 * alignlen / seq_len (0-6th bits)
		// qual     + (31-qual & 31)<<7 (7-12th bits)
		// NM       + (31-NM & 31)<<14 (14-19th bits)
		// unmapped + 1<<20 (20th bit)
		// nodata   + 1<<30 (30th bit)
		int rank = 0;
		if (bam->data == NULL)
			rank += 1<<30;
		else if (isUnMapped(bam)) {
			rank += 1<<20;
		} else {
			int nm = bam_aux2i(bam_aux_get(bam, "NM"));
			rank += (nm & 31) << 14;
			int qual = bam->core.qual;
			if (qual <= 31)
				rank += (31 - (qual & 31)) << 7;
			if (bam->core.l_qseq > 0)
				rank += 100 - (100 * (bam_calend(&(bam->core), bam1_cigar(bam)) - bam->core.pos)) / bam->core.l_qseq;
		}
		return rank;
	}
	static inline int getRankedPairedState(const bam1_t *bam) {
		// isMapped - getRankedSingleState (0-20th bits)
		// populate every other bit to allow room to add both mates
		// isMateUnmapped 21+ bits
		// isMateMapped 1<<21 (21th bit)
		// isProperPair 1<<23
		// both mapped to same contig 1<<25
		int rank = getRankedSingleState(bam);
		if (isMateUnMapped(bam)) {
			assert(!isProperPair(bam));
			rank += 255<<21;
		} else {
			if (!isProperPair(bam)) {
				rank += 1<<23;
			}
			if (bam->core.tid != bam->core.mtid) {
				rank += 1<<25;
			}
		}
		return rank;
	}

	// only use this when memory is not allocated!!!
	static void clearData(bam1_t *bam) {
		bam->data_len = 0; bam->m_data = 0; bam->data = NULL;
	}

	static char *copyBAMCore(char *dest, const bam1_t *read) {
		memcpy(dest, &read->core, sizeof(bam1_core_t));
		return dest + sizeof(bam1_core_t);
	}
	static int getBamDataSize(const bam1_t *read) {
		return sizeof(int)*2 + read->data_len;
	}
	static char *copyBAMData(char *dest, const bam1_t *read) {
		*(((int*)dest)) = read->l_aux;
		*(((int*)dest)+1) = read->data_len;
		dest += sizeof(int)*2;
		memcpy(dest, read->data, read->data_len);
		return dest + read->data_len;
	}
	static void restoreBAM(bam1_t *read, const bam1_core_t *core, const char *dataCopy) {
		memcpy(&read->core, core, sizeof(bam1_core_t));
		int *ptr = (int*) dataCopy;
		read->l_aux = *(ptr++);
		read->data_len = *(ptr++);
		if (read->m_data < read->data_len) {
			int m_data = roundup(read->data_len, 8);
			read->m_data =m_data;
			read->data = (uint8_t*) realloc(read->data, m_data);
		}
		memcpy(read->data, (char*) ptr, read->data_len);
	}

	static void initBamVector(BamVector &reads, int newSize) {
		reads.reserve(newSize+1); // add extra for option to add extra NULL at end for heapsort
		for(int i = reads.size() ; i < newSize ; i++)
			reads.push_back(bam_init1());
	}
	static void destroyBamVector(BamVector &reads) {
		for(BamVector::iterator it = reads.begin(); it != reads.end(); it++)
			if (*it != NULL)
				bam_destroy1(*it);
		reads.clear();
	}

};

#endif /* SAMUTILS_H_ */
