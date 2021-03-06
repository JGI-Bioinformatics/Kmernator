/*
 * MatcherInterface.h
 *
 *  Created on: Nov 1, 2011
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

#ifndef MATCHERINTERFACE_H_
#define MATCHERINTERFACE_H_

#include <vector>
#include "boost/unordered_set.hpp"
#include "boost/random/uniform_01.hpp"

#include "ReadSet.h"
#include "Utils.h"
#include "Options.h"
#include "MPIBuffer.h"
#include "KmerAlign.h"

class _MatcherInterfaceOptions  : public OptionsBaseInterface {
public:
	_MatcherInterfaceOptions() : maxReadMatches(450), maxReadDepthMatches(0), includeMate(true), minOverlap(51), minIdentity(0.986), returnOverlapOnly(true) {}
	virtual ~_MatcherInterfaceOptions() {}

	int &getMaxReadMatches() {
		return maxReadMatches;
	}
	int &getMaxReadDepthMatches() {
		return maxReadDepthMatches;
	}
	bool &getIncludeMate() {
		return includeMate;
	}
	void setIncludeMate(bool v) {
		if (v)
			includeMate = 1;
		else
			includeMate = 0;
	}
	int &getMinOverlap() {
		return minOverlap;
	}
	float &getMinIdentity() {
		return minIdentity;
	}
	bool &getReturnOverlapOnly() {
		return returnOverlapOnly;
	}

	// use to set/overrided any defaults on options that are stored persistently
	void _resetDefaults() {}
	// use to set the description of all options
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("Matching Options");
		opts.add_options()
						("max-read-matches", po::value<int>()->default_value(maxReadMatches), "maximum number of (randomly sampled) reads to return for matching. '0' disables.")

						("max-read-depth-matches", po::value<int>()->default_value(maxReadDepthMatches), "maximum number of (randomly sampled) reads per query length to return for matching. '0' disables.")

						("include-mate", po::value<bool>()->default_value(includeMate), "true - include mates, false - do not")

						("min-match-overlap", po::value<int>()->default_value(minOverlap), "The minimum amount of overlap for a matching read")

						("min-identity-fraction", po::value<float>()->default_value(minIdentity), "The minimum fraction identity for a matching read (to the end)")

						("return-overlap-only", po::value<bool>()->default_value(returnOverlapOnly), "If set true, only overlapping reads (or unaligned mates) will be returned")

						;
		desc.add(opts);
	}
	// use to post-process options, returning true if everything is okay
	bool _parseOptions(po::variables_map &vm) {
		setOpt("max-read-matches", maxReadMatches);
		setOpt("max-read-depth-matches", maxReadDepthMatches);\
		setOpt("include-mate", includeMate);
		setOpt("min-match-overlap", minOverlap);
		setOpt("min-identity-fraction", minIdentity);
		setOpt("return-overlap-only", returnOverlapOnly);

		return true;
	}
protected:
	int maxReadMatches;
	int maxReadDepthMatches;
	bool includeMate;
	int minOverlap;
	float minIdentity;
	bool returnOverlapOnly;

};
typedef OptionsBaseTemplate< _MatcherInterfaceOptions > MatcherInterfaceOptions;


class MatcherInterface : public Timer {
public:
	typedef boost::unordered_set< ReadSet::ReadSetSizeType > ReadIdxSet;
	typedef ReadIdxSet MatchHitSet;
	typedef std::vector< MatchHitSet > MatchResults;
	typedef std::vector< ReadSet::ReadSetSizeType > MatchHitVector;
	typedef ReadSet::ReadSetVector MatchReadResults;

	MatcherInterface(mpi::communicator &world, const ReadSet &target)
	: _world(world, mpi::comm_duplicate), _target(target), globalQueryFile(), rmGlobalQueryFile(false) {
		assert(_target.isGlobal() && _target.getGlobalSize() > 0);
		assert(!omp_in_parallel());
		FloatRand::getInstance(omp_get_max_threads());	// initialize random number thread-safe generators
	}
	~MatcherInterface() {
		cleanGlobalQueryFile();
	}
	void cleanGlobalQueryFile() {
		if (rmGlobalQueryFile) {
			unlink(globalQueryFile.c_str());
			globalQueryFile.clear();
		}
		rmGlobalQueryFile = false;
		_world.barrier();
	}

	// returns a list of sets of local reads (i.e. target reads idxs in their globalReadIdx space)
	// to a query, which should be a global ReadSet
	MatchResults matchLocal(const ReadSet &query) {
		assert(query.isGlobal());
		cleanGlobalQueryFile();
		std::string queryFile = UniqueName::generateUniqueGlobalName( GeneralOptions::getOptions().getTmpDir() + "/MatcherInterface-" );
		globalQueryFile = DistributedOfstreamMap::writeGlobalReadSet(_world, query, queryFile, "", FormatOutput::Fasta(), false);
		rmGlobalQueryFile = true;

		LOG_DEBUG(3, "Running matchLocal on " << globalQueryFile);
		MatchResults results = this->matchLocalImpl(globalQueryFile);
		assert(results.size() == query.getGlobalSize());
		return results;
	}

	// returns a list of sets of local reads (i.e. target reads idxs in their globalReadIdx space)
	// to a query, in its localReadIdxSpace space (should be full global copy, not a global ReadSet)
	MatchResults matchLocal(std::string queryFile) {
		globalQueryFile = queryFile;
		rmGlobalQueryFile = false;
		MatchResults results = this->matchLocalImpl(queryFile);
		return results;
	}

	virtual MatchResults matchLocalImpl(std::string(queryFile)) {
		LOG_THROW("Base class must implement the matchLocalImpl(string) method!");
	}

	// returns a ReadSetVector of (possibly copied Read) results over the target global reads (i.e. full set)
	// matching the local reads in the query global ReadSet (i.e. partition of the global)
	// queryFile is expected to contain the entire query set and will be used in matchLocal
	MatchReadResults match(const ReadSet &query, std::string queryFile) {
		assert(query.isGlobal() && query.getGlobalSize() > 0);
		MatchResults matchResults = this->matchLocal(queryFile);
		matchResults.resize(query.getGlobalSize(), MatchHitSet());
		return convertLocalMatchesToGlobalReads(query, matchResults);
	}

	// returns a ReadSetVector of (possibly copied Read) results over the target global reads (i.e. full set)
	// matching the local reads in the query global ReadSet (i.e. partition of the global)
	MatchReadResults match(const ReadSet &query) {
		assert(query.isGlobal() && query.getGlobalSize() > 0);
		MatchResults matchResults = this->matchLocal(query);
		assert(matchResults.size() == query.getGlobalSize());
		return convertLocalMatchesToGlobalReads(query, matchResults);
	}

	static void debuglog(int level, std::string label, const MatchResults &matchResults) {
		if (Log::isDebug(level)) {
			std::stringstream ss;
			ss << label << ":" << std::endl;
			for(int i = 0 ; i < (int) matchResults.size(); i++) {
				ss << i << ":";
				for(MatchHitSet::iterator it = matchResults[i].begin(); it != matchResults[i].end(); it++)
					ss << " " << *it;
				ss << std::endl;
			}
			std::string s = ss.str();
			LOG_DEBUG(level, s);
		}
	}
	static void debuglog(int level, std::string label, const MatchReadResults &matchReadResults) {
		if (Log::isDebug(level)) {
			std::stringstream ss;
			ss << label << ":" << std::endl;
			for(int i = 0 ; i < (int) matchReadResults.size(); i++) {
				ss << i << ":";
				for(ReadSet::ReadSetSizeType j = 0; j < matchReadResults[i].getSize(); j++)
					ss << " " << matchReadResults[i].getRead(j).getName();
				ss << std::endl;
			}
			std::string s = ss.str();
			LOG_DEBUG(level, s);
		}
	}

	// takes a query global ReadSet
	// a MatchResults of globalReadIdx
	// and returns the global ReadSetVector (possibly copied reads) for the local reads within the query
	MatchReadResults convertLocalMatchesToGlobalReads(const ReadSet &query, MatchResults &matchResults) {
		LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "convertLocalMatches():");
		assert(query.isGlobal() && query.getGlobalSize() > 0);
		assert(matchResults.size() == query.getGlobalSize());
		debuglog(3, "MatcherInterface::converLocalMatchesToGlobalReads(): LocalMatches", matchResults);
		MatchReadResults localReads = getLocalReads(matchResults);
		assert(localReads.size() == query.getGlobalSize());
		debuglog(3, "MatcherInterface::converLocalMatchesToGlobalReads(): LocalReads", localReads);
		MatchReadResults globalReads = exchangeGlobalReads(query, localReads);
		assert(globalReads.size() == query.getSize());
		debuglog(3, "MatcherInterface::converLocalMatchesToGlobalReads(): GlobalReads", globalReads);
		recordTime("returnMatch", MPI_Wtime());
		return globalReads;
	}

	// randomly downsizes overly full MatchHitSets in the MatchResults
	void sampleMatches(MatchResults &matchResults, ReadSet::ReadSetSizeType maxMatches = MatcherInterfaceOptions::getOptions().getMaxReadMatches()) {

		if (maxMatches <= 0)
			return;
#pragma omp parallel for
		for(long i = 0 ; i < (long) matchResults.size(); i++) {
			MatchHitSet &mhs = matchResults[i];
			if (mhs.size() > maxMatches) {
				MatchHitVector mhv(mhs.begin(), mhs.end());
				ReadSet::ReadSetSizeType oldSize = mhv.size();
				mhs.clear();
				Random<ReadSet::ReadSetSizeType>::Set sampledSet = Random<ReadSet::ReadSetSizeType>::sample(oldSize, maxMatches);
				for(Random<ReadSet::ReadSetSizeType>::SetIterator it = sampledSet.begin(); it != sampledSet.end(); it++) {
					mhs.insert(mhv[*it]);
				}
				LOG_DEBUG_OPTIONAL(2, true, "Reduced " << i << " from " << oldSize << " to " << mhs.size());
			}
		}
		return;
	}

	MatchHitSet sampleMatches(MatchHitSet &all, float fraction) {
		FloatRand &myRand = FloatRand::getInstance();
		MatchHitSet keep;
		for(MatchHitSet::iterator it = all.begin(); it != all.end(); it++)
			if (myRand.getRand() <= fraction)
				keep.insert(*it);
		return keep;
	}

	// returns a ReadSetVector of reads that are local (i.e. targets in globalReadIdx space local to the node)
	// exchanges globalReadSetIds with the responsible nodes, if needed
	MatchReadResults getLocalReads(MatchResults &globalMatchResults) {

		MatchResults matchResults = exchangeGlobalReadIdxs(globalMatchResults);
		debuglog(3, "MatcherInterface::getLocalReads(): after exchange LocalMatches", matchResults);
		assert(matchResults.size() == globalMatchResults.size());

		bool isPaired = MatcherInterfaceOptions::getOptions().getIncludeMate();
		bool screenForOverhang = MatcherInterfaceOptions::getOptions().getReturnOverlapOnly();
		bool screenNow = screenForOverhang && !globalQueryFile.empty();
		ReadSetStream *rss = NULL;
		if (screenNow) {
			rss = new ReadSetStream(globalQueryFile);
		} else {
			sampleMatches(matchResults, MatcherInterfaceOptions::getOptions().getMaxReadMatches() * 20);
		}

		long totalMatches = 0, maxMatches = 0;

		int myRank = _world.rank();
		MatchReadResults matchReadResults(matchResults.size(), ReadSet());
		for(int i = 0; i < (int) matchResults.size(); i++) {
			MatchHitVector mhv(matchResults[i].begin(), matchResults[i].end());
			std::sort(mhv.begin(), mhv.end());
			matchResults[i].clear();
			for(MatchHitVector::iterator it = mhv.begin(); it != mhv.end(); it++) {
				assert(getTarget().isLocalRead( *it ));
				ReadSet::ReadSetSizeType localReadIdx = getTarget().getLocalReadIdx(myRank, *it);
				const Read read = getTarget().getRead( localReadIdx );
				if (! read.isDiscarded() )
					matchReadResults[i].append( read );
			}
			if (screenNow) {
				if (isPaired)
					matchReadResults[i].identifyPairs();
				if (rss->hasNext()) {
					Read contig = rss->getRead();
					matchReadResults[i] = screenAlignmentsForOverhang(contig, matchReadResults[i], isPaired);
				} else {
					LOG_WARN(1, "mismatch between globalQuery file and number of reads");
					assert(false); // should not get here
				}
			}
			long size = matchReadResults[i].getSize();
			totalMatches += size;
			if (maxMatches < size)
				maxMatches = size;
		}

		if (screenNow)
			delete rss;

		recordTime("returnLocalReads", MPI_Wtime());
		LOG_DEBUG_GATHER(1, "getLocalReads(): finished. totalMatches: " << totalMatches << ", maxMatches: " << maxMatches << ". " << MemoryUtils::getMemoryUsage());

		return matchReadResults;
	}
	// transfers the proper matching localReadSetVector to
	// a globalReadSetVector to the node controlling the global reads over the query
	// reads will be copied from the localReadSet to the destination node
	// the returning matches are for the local reads in the query, global reads from the Target
	// localReadSetVector will be empty after this routine is called
	MatchReadResults exchangeGlobalReads(const ReadSet &query, MatchReadResults &localReadSetVector) {
		assert(query.isGlobal() && query.getGlobalSize() > 0);
		assert(localReadSetVector.size() == query.getGlobalSize());
		MatchReadResults globalReadSetVector, tmpReadSetVector;
		globalReadSetVector.resize(query.getSize(), ReadSet());
		tmpReadSetVector.resize(query.getGlobalSize(), ReadSet());

		bool isPaired = MatcherInterfaceOptions::getOptions().getIncludeMate();
		bool screenForOverlap = MatcherInterfaceOptions::getOptions().getReturnOverlapOnly() && globalQueryFile.empty();
		LOG_DEBUG_GATHER(1, "exchangeGlobalReads(): will " << (screenForOverlap? " " : "NOT ") << " screenForOverhang. " << MemoryUtils::getMemoryUsage());

		int myRank = _world.rank();

		int sendBytes[_world.size()], recvBytes[_world.size()],
		sendDisp[_world.size()], recvDisp[_world.size()];

		ReadIdxSet sendingGlobalContigIdx;
		ReadSet::ReadSetSizeType globalSize = query.getGlobalSize();
		ReadSet::ReadSetSizeType sendingReads;

		ReadSet empty;
		int emptyReadSetSize = empty.getStoreSize();
		int emptyBase = emptyReadSetSize * globalSize;
		// maximum size for a single ReadSet per rank (1/2 the full buffer)
		int maxRankTransmitSize1 = ( MPIOptions::getOptions().getTotalBufferSize() / _world.size() / 2);
		// minimum size
		if (maxRankTransmitSize1 < 4096)
			maxRankTransmitSize1 = 4096;
		// maximum size for a round of communication (2 readSets + all emptySet signals)
		int maxRankTransmitSize = maxRankTransmitSize1*2 + emptyBase + 64;
		maxRankTransmitSize += 1024 - (maxRankTransmitSize % 1024); // make a multiple of 1kb
		long maxBuffer = (maxRankTransmitSize) * _world.size() + 1024;
		bool isDone = false;
		long iteration = 0;
		while (! isDone ) {
			iteration++;
			long localReadSets = 0;
			LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "exchangeGlobalReads: " << iteration << ". Maximum transmit buffer is: " << maxRankTransmitSize1 << ", " << maxRankTransmitSize << " / " << maxBuffer);
			isDone = true;
			int totalSend = 0, totalRecv = 0; sendingReads = 0;
			for (int rank = 0; rank < _world.size(); rank++) {
				sendBytes[rank] = 1; // artificially pad by one byte to ensure negative 'isDone' flag can be sent
				recvBytes[rank] = 0;
			}

			for (ReadSet::ReadSetSizeType globalContigIdx = 0; globalContigIdx < globalSize; globalContigIdx++) {
				int rank;
				ReadSet::ReadSetSizeType rankReadIdx;
				query.getRankReadForGlobalReadIdx(globalContigIdx, rank, rankReadIdx);
				if (rank == myRank) {
					// do not encode and send reads to self
					assert(query.isLocalRead(globalContigIdx));
					ReadSet::ReadSetSizeType localIdx = query.getLocalReadIdx(globalContigIdx);
					globalReadSetVector[localIdx].append(localReadSetVector[globalContigIdx]);
					localReadSetVector[globalContigIdx].clear();
					tmpReadSetVector[globalContigIdx].clear();
					localReadSets++;
					continue;
				}
				int sendByteCount = localReadSetVector[globalContigIdx].getStoreSize();
				while (sendByteCount >= maxRankTransmitSize1) {
					// reduce the transmission read set into a smaller chunk
					ReadSet::ReadSetSizeType size = localReadSetVector[globalContigIdx].getSize();
					ReadSet::ReadSetSizeType maxSend = size / 2;
					if (maxSend <= 0) {
						LOG_WARN(1, "Could not send any reads for " << globalContigIdx << " as buffer size is too small");
						localReadSetVector[globalContigIdx].clear();
						tmpReadSetVector[globalContigIdx].clear();
					} else {
						LOG_DEBUG_OPTIONAL(1, true, "globalContig: " << globalContigIdx << " with " << localReadSetVector[globalContigIdx].getSize() << " reads is too big to transmit, reducing to " << maxSend);
						// temporarily store in globalReadSetVector
						tmpReadSetVector[globalContigIdx].append(localReadSetVector[globalContigIdx].truncate(maxSend));
					}
					sendByteCount = localReadSetVector[globalContigIdx].getStoreSize();
				}
				if (sendByteCount <= emptyReadSetSize || (int) (sendBytes[rank] + sendByteCount + (emptyReadSetSize*(globalSize - globalContigIdx - 1))) < maxRankTransmitSize) {
					// enough room for this plus all remaining sets if they are empty
					sendingGlobalContigIdx.insert(globalContigIdx);
					sendingReads += localReadSetVector[globalContigIdx].getSize();
				} else {
					// no room for this set this round
					// keep reads in localReadSetVector[globalContigIdx] for next round
					isDone = false;
					long requiredSendByteCount = sendByteCount;
					sendByteCount = empty.getStoreSize();
					ReadIdxSet::iterator readIdIt = sendingGlobalContigIdx.find(globalContigIdx);
					if (readIdIt != sendingGlobalContigIdx.end())
						sendingGlobalContigIdx.erase(readIdIt);
					LOG_DEBUG_OPTIONAL(3, true, "Skipping globalContigIdx this round(" << iteration << "): " << globalContigIdx << " need " << requiredSendByteCount << " prepared: " << sendBytes[rank] << " remaining: " << maxRankTransmitSize - sendBytes[rank]);
				}
				if (tmpReadSetVector[globalContigIdx].getSize() > 0) {
					// some reads are waiting for a future communication
					isDone = false;
				}
				sendBytes[rank] += sendByteCount;
				assert(sendBytes[rank] <= maxBuffer / _world.size() && sendBytes[rank] >= 0); // i.e. no overflow
				LOG_DEBUG(4, "GlobalContig: " << globalContigIdx << " sending " << localReadSetVector[ globalContigIdx ].getSize() << " reads / "<< sendByteCount << " bytes to " << rank << " total: " << sendBytes[rank]);
			}
			for (int rank = 0; rank < _world.size(); rank++) {
				sendDisp[rank] = totalSend;
				assert(sendDisp[rank] >= 0); // i.e. no overflow
				totalSend += sendBytes[rank] - 1; // do not count the artificial padding
				if (! isDone ) {
					// signal others this rank is not done
					sendBytes[rank] = 0 - sendBytes[rank];
				}
				assert(totalSend >= 0); // i.e. no overflow
				LOG_DEBUG_OPTIONAL(2, true, "Sending to rank " << rank << " sendBytes " << sendBytes[rank] << " sendDisp[] " << sendDisp[rank] << ". 0 == recvBytes[] "<< recvBytes[rank]);
			}
			LOG_DEBUG_OPTIONAL(1, true, "iteration " << iteration << " Sending bytes: " << totalSend << " readSets: " << sendingGlobalContigIdx.size() + localReadSets << " sending Reads: " << sendingReads << (isDone ? " isDone" : ""));

			MPI_Alltoall(sendBytes, 1, MPI_INT, recvBytes, 1, MPI_INT, _world);

			for (int rank = 0; rank < _world.size(); rank++) {

				if (sendBytes[rank] < 0) {
					// restore proper send bytes as they will be used for transmit sizes
					sendBytes[rank] = 0 - sendBytes[rank];
				}
				sendBytes[rank] -= 1;  // remove the artificially padded byte

				recvDisp[rank] = totalRecv;
				if (recvBytes[rank] < 0) {
					// intercept notDone signal
					isDone = false;
					recvBytes[rank] = 0 - recvBytes[rank];
				}
				recvBytes[rank] -= 1; // remove the artificially padded byte

				totalRecv += recvBytes[rank];
				assert(totalRecv >= 0); // i.e. no overflow
				LOG_DEBUG_OPTIONAL(3, true, "to/from rank " << rank << ": sendDisp[] " << sendDisp[rank] << " sendBytes[] " << sendBytes[rank] << " recvDisp[] " << recvDisp[rank] << " recvBytes[] " << recvBytes[rank] );
			}
			LOG_DEBUG_OPTIONAL(2, true, "Sending " << totalSend << " Receiving " << totalRecv);
			if (totalSend < 0 || totalRecv < 0)
				LOG_THROW("int overflow: too many reads to send / receive in one message packet");

			char *sendBuf = (char*) (malloc(totalSend == 0 ? 8 : totalSend));
			if (sendBuf == NULL)
				LOG_THROW("MatcherInterface::exchangeGlobalReads(): Could not allocate totalSend bytes! " << totalSend);

			char *tmp = sendBuf;
			for (ReadSet::ReadSetSizeType globalContigIdx = 0; globalContigIdx < globalSize; globalContigIdx++) {
				int storeSize = 0;
				if (query.isLocalRead(globalContigIdx)) {
					// do not encode and send reads to self
					assert(localReadSetVector[globalContigIdx].getSize() == 0);
					assert(tmpReadSetVector[globalContigIdx].getSize() == 0);
				} else {
					if (sendingGlobalContigIdx.find(globalContigIdx) != sendingGlobalContigIdx.end()) {
						storeSize = localReadSetVector[globalContigIdx].store(tmp);
						localReadSetVector[globalContigIdx].clear();
						// swap back any truncated reads for next run
						tmpReadSetVector[globalContigIdx].swap(localReadSetVector[globalContigIdx]);
					} else {
						storeSize = empty.store(tmp);
						LOG_DEBUG_OPTIONAL(3, true, "Not sending " << globalContigIdx << " with " << localReadSetVector[globalContigIdx].getSize() << " and " << tmpReadSetVector[globalContigIdx].getSize());
					}
				}
				tmp += storeSize;
				LOG_DEBUG(4, "GlobalContig: " << globalContigIdx << " prepared send " << storeSize << " bytes");
			}
			assert(tmp == sendBuf + totalSend);

			char *recvBuf = (char*) (malloc(totalRecv == 0 ? 8 : totalRecv));
			if (recvBuf == NULL)
				LOG_THROW("Could not allocate totalRecv bytes! " << totalRecv);

			MPI_Alltoallv(sendBuf, sendBytes, sendDisp, MPI_BYTE, recvBuf,
					recvBytes, recvDisp, MPI_BYTE, _world);

			free(sendBuf);
			// set localReadSetVector to hold read sets for local set of query
			for (int rank = 0; rank < _world.size(); rank++) {
				LOG_DEBUG(3, "restoring for rank " << rank << " recvDisp[rank] " << recvDisp[rank] << " recvBytes[rank] " << recvBytes[rank] << " totalRecv " << totalRecv);
				tmp = recvBuf + recvDisp[rank];
				if (rank == myRank) {
					// no reads will be sent to self
					assert(recvBytes[rank] == 0);
					continue;
				}
				ReadSet::ReadSetSizeType globalContigIdx = query.getGlobalOffset();
				while (tmp != recvBuf + recvDisp[rank] + recvBytes[rank]) {
					assert(query.isLocalRead(globalContigIdx));
					assert(globalContigIdx < query.getGlobalSize());
					assert(tmp < (recvBuf + recvDisp[rank] + recvBytes[rank]));
					assert(tmp < (recvBuf + totalRecv));
					assert(globalReadSetVector.size() == query.getSize());
					ReadSet rs;
					tmp = (char*) rs.restore(tmp);
					ReadSet::ReadSetSizeType localContigIdx = query.getLocalReadIdx(globalContigIdx);
					globalReadSetVector[localContigIdx].append(rs);
					LOG_DEBUG_OPTIONAL(4, true, "Restored from rank " << rank << " ReadSet " << rs.getSize() << " totaling " << globalReadSetVector[localContigIdx].getSize() << " " << globalReadSetVector.size() << " restore " << globalContigIdx << " is local: " << query.isLocalRead(rank,globalContigIdx));
					LOG_DEBUG_OPTIONAL(3, rs.getSize() > 0, query.getRead( localContigIdx ).getName() << " (" << globalContigIdx << " / " << localContigIdx << ") has " << globalReadSetVector[localContigIdx].getSize() << " reads in the pool");
					globalContigIdx++;
				}
				assert(tmp == recvBuf + recvDisp[rank] + recvBytes[rank]);
				LOG_DEBUG(3, "finished restoring from rank " << rank);
			}
			assert(tmp == (recvBuf + totalRecv));
			free(recvBuf);
			recordTime("EGR" + boost::lexical_cast<std::string>(iteration), MPI_Wtime());
		}
		localReadSetVector.clear();

		ReadSet::ReadSetSizeType maxReadMatches = MatcherInterfaceOptions::getOptions().getMaxReadMatches();
		ReadSet::ReadSetSizeType maxReadDepth = MatcherInterfaceOptions::getOptions().getMaxReadDepthMatches();

		#pragma omp parallel for
		for (long localContigIdx = 0; localContigIdx < (long) globalReadSetVector.size(); localContigIdx++) {
			LOG_DEBUG_OPTIONAL(2, true, "exchangeGlobalReads(): contig " << localContigIdx << ", " << globalReadSetVector[localContigIdx].getSize() << " reads");
			if (screenForOverlap) {
				if (isPaired)
					globalReadSetVector[localContigIdx].identifyPairs();
				globalReadSetVector[localContigIdx] = screenAlignmentsForOverhang(query.getRead(localContigIdx), globalReadSetVector[localContigIdx], isPaired);
			}
			ReadSet::ReadSetSizeType maxReads = std::max(maxReadMatches, maxReadDepth * query.getRead(localContigIdx).getLength() / (getTarget().getSize() > 0 ? getTarget().getAvgSequenceLength() : 76) );
			if (maxReads > 0) {
				ReadSet::ReadSetSizeType numReads = globalReadSetVector[localContigIdx].getSize();
				if (numReads > maxReads) {
					LOG_DEBUG_OPTIONAL(2, true, "for " << localContigIdx << " sampled from " << numReads << " to " << maxReads);
					globalReadSetVector[localContigIdx] = globalReadSetVector[localContigIdx].randomlySample(maxReads);
				}
			}
		}

		LOG_DEBUG_GATHER(1, "exchangeGlobalReads(): Done " << iteration << ". " << MemoryUtils::getMemoryUsage());

		return globalReadSetVector;
	}

	bool isPassingRead(KmerAlign &kalign, const Read &read, SequenceLengthType minOverlap) {
		Alignment bestAlignment;
		return isPassingRead(kalign, read, bestAlignment, minOverlap);
	}
	bool isPassingRead(KmerAlign &kalign, const Read &read, Alignment &bestAlignment, SequenceLengthType minOverlap) {
		bestAlignment = kalign.getAlignment(read);
		float minIdentity = MatcherInterfaceOptions::getOptions().getMinIdentity();

		if (bestAlignment.getOverlap() < minOverlap || bestAlignment.getIdentity() < minIdentity)
			return false;
		if (bestAlignment.isAtEnd(kalign.getTarget(), read))
			return true;
		return false;
	}
	// screen matches for those that pass over the edge of the contig
	// or, if paired, are full match and have pair with no alignment
	ReadSet screenAlignmentsForOverhang(const Read &contig, ReadSet &matches, bool isPaired) {
		LOG_DEBUG(3, "screenAlignmentsForOverhang() on " << contig.toString() << " with " << matches.getSize() << " prospective matches");
		ReadSet screenedMatches;
		SequenceLengthType minOverlap = MatcherInterfaceOptions::getOptions().getMinOverlap();
		KmerAlign kalign(contig);
		LOG_DEBUG_OPTIONAL(2, true, "screenAlignmentsForOverhang(): contig " << contig.getName() << ", " << matches.getSize() << " reads, " << matches.getPairSize() << " pairs");
		if (isPaired && matches.hasPairs()) {
			for(ReadSet::ReadSetSizeType matchidx = 0; matchidx < matches.getPairSize(); matchidx++) {
				ReadSet::Pair &pair = matches.getPair(matchidx);
				Alignment aln1, aln2;
				bool r1 = false, r2 = false, p1 = false, p2 = false, end1 = false, end2 = false;
				r1 = matches.isValidRead(pair.read1);
				r2 = matches.isValidRead(pair.read2);
				Read read1, read2;
				if (r1) {
					read1 = matches.getRead(pair.read1);
					SequenceLengthType r1len = read1.getLength();
					p1 = isPassingRead(kalign, read1, aln1, minOverlap);
					end1 = aln1.targetAln.isAtEnd(contig);
					if (p1 & !end1) // allow perfect full length matches up to half a readlength away in the contig
						end1 |= (aln1.targetAln.isAtEnd(contig, r1len/2) & (aln1.getIdentity() == 1.0) & (aln1.getOverlap() == r1len));
				}
				if (r2) {	
					read2 = matches.getRead(pair.read2);
					SequenceLengthType r2len = read2.getLength();
					p2 = isPassingRead(kalign, read2, aln2, minOverlap);
					end2 = aln2.targetAln.isAtEnd(contig);
					if (p2 & !end2) // allow perfect full length matches up to half a readlength away in the contig
						end2 |= (aln2.targetAln.isAtEnd(contig, r2len/2) & (aln2.getIdentity() == 1.0) & (aln2.getOverlap() == r2len));
				}

				bool add1 = false, add2 = false;
				if (p1) {
					// only include read1 if it overlaps the end
					if (end1)
						add1 = true;
					if (r2 && ((p2 & end2) || aln2.targetAln.getOverlap() < minOverlap))
						add2 = true;
				} else if (p2) {
					if (r1 && ((p1 & end1) || aln1.targetAln.getOverlap() < minOverlap))
						add1 = true;
					// only include read2 if it overlaps the end
					if (end2)
						add2 = true;
				}
				LOG_DEBUG(4, "screenAlignmentsForOverhang() (" << add1 << "," << add2 << ") " << p1 << "/" << end1 << " " << aln1.toString() << " pair:" << p2 << "/" << end2 << " " << aln2.toString() << read1.toString() << read2.toString());

				if (add1)
					screenedMatches.append(read1);
				if (add2)
					screenedMatches.append(read2);
			}
		} else {
			for(ReadSet::ReadSetSizeType matchidx = 0; matchidx < matches.getSize(); matchidx++) {
				const Read &match = matches.getRead(matchidx);
				Alignment aln;
				if (isPassingRead(kalign, match, aln, minOverlap))
					if (aln.targetAln.isAtEnd(contig))
						screenedMatches.append(match);
			}
		}
		LOG_DEBUG_OPTIONAL(2, true, "screenedMatches: " << screenedMatches.getSize());
		return screenedMatches;
	}

	// const version of exchangeGlobalReads
	//MatchReadResults exchangeGlobalReads(const ReadSet &query, const MatchReadResults &localReadSetVector) {
	//	MatchReadResults copy = localReadSetVector;
	//	return exchangeGlobalReads(query, copy);
	//}
	bool areAllReadsLocal(const MatchResults &matchResults) {
		int myRank = _world.rank();
		for(MatchResults::const_iterator it = matchResults.begin(); it != matchResults.end(); it++) {
			for(MatchHitSet::iterator it2 = it->begin(); it2 != it->end(); it2++) {
				if (!getTarget().isLocalRead(myRank, *it2)) {
					return false;
				}
			}
		}
		return true;
	}

	// exchanges globalReadIdx matches with responsible nodes
	// so all nodes end up with reads in their own globalIndex ranges
	MatchResults exchangeGlobalReadIdxs(MatchResults &globalMatchResults) {
		int numMatchHitSets = globalMatchResults.size();
		int numRanks = _world.size();
		// sort globalReadIdxs by rank, matchHitSet
		std::vector< std::vector< MatchHitVector > > rankHitGlobalIndexes(numRanks);
		for(long i = 0; i < numRanks; i++)
			rankHitGlobalIndexes[i].resize(numMatchHitSets);
		for(long i = 0; i < numMatchHitSets; i++) {
			for(MatchHitSet::const_iterator it = globalMatchResults[i].begin(); it != globalMatchResults[i].end(); it++) {
				int destRank;
				ReadSet::ReadSetSizeType destRankLocalIdx;
				getTarget().getRankReadForGlobalReadIdx(*it, destRank, destRankLocalIdx);
				rankHitGlobalIndexes[destRank][i].push_back(*it);
			}
		}
		// send each rank the count of globalReadIdx that will send alltoall
		int *sendCounts = new int[ numRanks * numMatchHitSets ];
		int *sendRankCount = new int[numRanks];
		int *sendRankDispl = new int[numRanks];
		long totalSendCount = 0;
		long maxSendCount = 0;

		int *recvCounts = new int[ numRanks * numMatchHitSets ];
		int *recvRankCount = new int[numRanks];
		int *recvRankDispl = new int[numRanks];
		long totalRecvCount = 0;

		for(int i = 0; i < numRanks; i++) {
			sendRankCount[i] = 0;
			sendRankDispl[i] = 0;
			for(int j = 0 ; j < numMatchHitSets ; j++) {
				long pos = i * numMatchHitSets + j;
				sendCounts[pos] = rankHitGlobalIndexes[i][j].size();
				sendRankCount[i] += sendCounts[pos];
				totalSendCount += sendCounts[pos];
				if (maxSendCount < sendCounts[pos])
					maxSendCount = sendCounts[pos];
			}
			if(i > 0) {
				sendRankDispl[i] = sendRankDispl[i-1] + sendRankCount[i-1];
			}
		}

		LOG_DEBUG_GATHER(1, "exchangeGlobalReadIdxs(): sending: " << totalSendCount << " readIndexes, max: " << maxSendCount<< ". " << MemoryUtils::getMemoryUsage());

		MPI_Alltoall(sendCounts, numMatchHitSets, MPI_INT, recvCounts, numMatchHitSets, MPI_INT, _world);
		// Prepare accounting arrays for alltoallv
		for(int i = 0; i < numRanks; i++) {
			recvRankCount[i] = 0;
			recvRankDispl[i] = 0;
			for(int j = 0; j < numMatchHitSets ; j++) {
				long pos = i * numMatchHitSets + j;
				recvRankCount[i] += recvCounts[pos];
				totalRecvCount += recvCounts[pos];
			}
			if(i > 0) {
				recvRankDispl[i] = recvRankDispl[i-1] + recvRankCount[i-1];
			}
		}

		// build sending buffer
		ReadSet::ReadSetSizeType *sendBuf = new ReadSet::ReadSetSizeType[totalSendCount];
		ReadSet::ReadSetSizeType *tmp;
		for(int i = 0; i < numRanks; i++) {
			tmp = sendBuf + sendRankDispl[i];
			for(int j = 0; j < numMatchHitSets ; j++) {
				long pos = i * numMatchHitSets + j;
				//std::set<ReadSet::ReadSetSizeType> test(rankHitGlobalIndexes[i][j].begin(), rankHitGlobalIndexes[i][j].end());
				//LOG_DEBUG_OPTIONAL(1, test.size() != sendCounts[pos], "buildingSendBuffer mismatch problem: " << i << ", " << j << " size " << test.size() << " vs. " << sendCounts[pos]);
				for(int k = 0; k < sendCounts[pos]; k++) {
					//LOG_DEBUG_OPTIONAL(1, test.size() != sendCounts[pos], "buildingSendBuffer: " << rankHitGlobalIndexes[i][j][k]);
					*tmp = rankHitGlobalIndexes[i][j][k];
					tmp++;
				}
			}
			assert(tmp == sendBuf + sendRankDispl[i] + sendRankCount[i]);
		}
		// free some memory before using more
		delete [] sendCounts;
		rankHitGlobalIndexes.clear();
		// allocate receiving buffer
		ReadSet::ReadSetSizeType *recvBuf = new ReadSet::ReadSetSizeType[totalRecvCount];

		// exchange globalReadIdx to rank that owns it
		MPI_Alltoallv(sendBuf, sendRankCount, sendRankDispl, MPIReadSetSizeType, recvBuf, recvRankCount, recvRankDispl, MPIReadSetSizeType, _world);
		delete [] sendBuf;
		delete [] sendRankCount;
		delete [] sendRankDispl;

		bool includeMates = getTarget().hasPairs() && MatcherInterfaceOptions::getOptions().getIncludeMate();

		// consolidate into localMatchResults
		MatchResults localMatchResults(globalMatchResults.size());
		for(int i = 0; i < numRanks; i++) {
			tmp = recvBuf + recvRankDispl[i];
			for(int j = 0; j < numMatchHitSets ; j++) {
				long pos = i * numMatchHitSets + j;
				for(int k = 0; k < recvCounts[pos]; k++) {
					ReadSet::ReadSetSizeType globalReadIdx = *tmp;
					ReadSet::ReadSetSizeType localReadIdx = getTarget().getLocalReadIdx( globalReadIdx );
					LOG_DEBUG(4, "exchangeGlobalReadIdxs(): Adding result: " << j << " " << globalReadIdx << " " << getTarget().getRead(localReadIdx).getName());
					localMatchResults[j].insert( globalReadIdx );
					if (includeMates) {
						ReadSet::ReadSetSizeType localReadPairIdx = getTarget().getLocalPairIdx( localReadIdx );
						if ( getTarget().isValidRead(localReadPairIdx) ) {
							ReadSet::ReadSetSizeType globalReadPairIdx = getTarget().getGlobalReadIdx( localReadPairIdx );
							LOG_DEBUG(4, "exchangeGlobalReadIdxs(): Adding mate for result: " << j << " " << globalReadPairIdx << " " << getTarget().getRead(localReadIdx).getName() << ": " << getTarget().getRead(localReadPairIdx).getName() );
							localMatchResults[j].insert(globalReadPairIdx);
						}
					}
					tmp++;
				}
			}
			assert(tmp == recvBuf + recvRankDispl[i] + recvRankCount[i]);
		}
		long totalLocalReads = 0;
		long maxLocalReads = 0;
		for(int j = 0; j < numMatchHitSets ; j++) {
			long size = localMatchResults[j].size();
			totalLocalReads += size;
			if (maxLocalReads < size)
				maxLocalReads = size;
		}

		delete [] recvRankDispl;
		delete [] recvRankCount;
		delete [] recvCounts;
		delete [] recvBuf;
		recordTime("EGRI", MPI_Wtime());
		LOG_DEBUG_GATHER(1, "exchangeGlobalReadIdxs(): totalLocalReads: " << totalLocalReads << " max: " << maxLocalReads<< ". " << MemoryUtils::getMemoryUsage());

		assert(areAllReadsLocal(localMatchResults));
		return localMatchResults;
	}

	mpi::communicator &getWorld() { return _world; }
	const ReadSet &getTarget() const { return _target; }

protected:
	mpi::communicator _world;
	const ReadSet &_target;
	std::string globalQueryFile;
	bool rmGlobalQueryFile;
};

#endif /* MATCHERINTERFACE_H_ */
