/*
 * MatcherInterface.h
 *
 *  Created on: Nov 1, 2011
 *      Author: regan
 */

#ifndef MATCHERINTERFACE_H_
#define MATCHERINTERFACE_H_

#include <vector>
#include "boost/unordered_set.hpp"

#include "ReadSet.h"
#include "Utils.h"
#include "Options.h"
#include "MPIBuffer.h"

class _MatcherInterfaceOptions  : public OptionsBaseInterface {
public:
	_MatcherInterfaceOptions() : maxReadMatches(500) {}
	virtual ~_MatcherInterfaceOptions() {}

	int &getMaxReadMatches() {
		return maxReadMatches;
	}
	// use to set/overrided any defaults on options that are stored persistently
	void _resetDefaults() {}
	// use to set the description of all options
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("Matching Options");
		opts.add_options()
		("max-read-matches", po::value<int>()->default_value(maxReadMatches),
						"maximum number of (randomly sampled) reads to return for matching. '0' disables.")
						;
		desc.add(opts);
	}
	// use to post-process options, returning true if everything is okay
	bool _parseOptions(po::variables_map &vm) {
		setOpt<int>("max-read-matches", maxReadMatches);
		return true;
	}
protected:
	int maxReadMatches;
};
typedef OptionsBaseTemplate< _MatcherInterfaceOptions > MatcherInterfaceOptions;


class MatcherInterface : public Timer {
public:
	typedef std::set< ReadSet::ReadSetSizeType > MatchHitSet;
	typedef std::vector< ReadSet::ReadSetSizeType > MatchHitVector;
	typedef std::vector< MatchHitSet > MatchResults;
	typedef ReadSet::ReadSetVector MatchReadResults;
	typedef boost::unordered_set<ReadSet::ReadSetSizeType> ReadIdxSet;

	MatcherInterface(mpi::communicator &world, const ReadSet &target, bool returnPairedMatches = true)
	: _world(world), _target(target), _returnPairedMatches(returnPairedMatches) {
		assert(_target.isGlobal() && _target.getGlobalSize() > 0);
	}

	// returns a list of sets of local reads (i.e. target reads idxs in their globalReadIdx space)
	// to a query, which should be a global ReadSet
	virtual MatchResults matchLocal(const ReadSet &query) {
		assert(query.isGlobal());
		std::string queryFile = UniqueName::generateUniqueGlobalName( GeneralOptions::getOptions().getTmpDir() + "/MatcherInterface-" );
		queryFile = DistributedOfstreamMap::writeGlobalReadSet(_world, query, queryFile, "", FormatOutput::Fasta());
		LOG_DEBUG(3, "Running matchLocal on " << queryFile);
		MatchResults results = this->matchLocal(queryFile);
		assert(results.size() == query.getGlobalSize());
		LOG_DEBUG(3, "Waiting for rest of world to finish");
		_world.barrier();
		if (_world.rank() == 0)
			unlink(queryFile.c_str());
		return results;
	}

	// returns a list of sets of local reads (i.e. target reads idxs in their globalReadIdx space)
	// to a query, in its localReadIdxSpace space (should be full global copy, not a global ReadSet)
	virtual MatchResults matchLocal(std::string queryFile) {
		LOG_THROW("You must implement this function for your specific MatcherInterface");
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
		assert(query.isGlobal() && query.getGlobalSize() > 0);
		assert(matchResults.size() == query.getGlobalSize());
		debuglog(3, "LocalMatches", matchResults);
		MatchReadResults localReads = getLocalReads(matchResults);
		assert(localReads.size() == query.getGlobalSize());
		debuglog(3, "LocalReads", localReads);
		MatchReadResults globalReads = exchangeGlobalReads(query, localReads);
		assert(globalReads.size() == query.getSize());
		debuglog(3, "GlobalReads", globalReads);
		recordTime("returnMatch", MPI_Wtime());
		return globalReads;
	}

	// randomly downsizes overly full MatchHitSets in the MatchResults
	void sampleMatches(MatchResults &matchResults) {
		ReadSet::ReadSetSizeType maxMatches = MatcherInterfaceOptions::getOptions().getMaxReadMatches();

		if (maxMatches <= 0)
			return;

		for(ReadSet::ReadSetSizeType i = 0 ; i < matchResults.size(); i++) {
			MatchHitSet &mhs = matchResults[i];
			if (mhs.size() > maxMatches) {
				MatchHitVector mhv(mhs.begin(), mhs.end());
				ReadSet::ReadSetSizeType oldSize = mhv.size();
				mhs.clear();
				Random<ReadSet::ReadSetSizeType>::Set sampledSet = Random<ReadSet::ReadSetSizeType>::sample(oldSize, maxMatches);
				for(Random<ReadSet::ReadSetSizeType>::SetIterator it = sampledSet.begin(); it != sampledSet.end(); it++) {
					mhs.insert(mhv[*it]);
				}
				LOG_DEBUG(3, "Reduced " << i << " from " << oldSize << " to " << mhs.size());
			}
		}
		return;
	}
	// returns a ReadSetVector of reads that are local (i.e. targets in globalReadIdx space local to the node)
	// exchanges globalReadSetIds with the responsible nodes, if needed
	MatchReadResults getLocalReads(MatchResults &globalMatchResults) {
		MatchResults matchResults = exchangeGlobalReadIdxs(globalMatchResults);
		assert(matchResults.size() == globalMatchResults.size());

		int myRank = _world.rank();
		MatchReadResults matchReadResults(matchResults.size(), ReadSet());
		for(int i = 0; i < (int) matchResults.size(); i++) {
			for(MatchHitSet::iterator it = matchResults[i].begin(); it != matchResults[i].end(); it++) {
				assert(getTarget().isLocalRead( *it ));
				ReadSet::ReadSetSizeType localReadIdx = getTarget().getLocalReadIdx(myRank, *it);
				matchReadResults[i].append( getTarget().getRead( localReadIdx ));
			}
		}
		recordTime("returnLocalReads", MPI_Wtime());
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
		int myRank = _world.rank();

		int sendBytes[_world.size()], recvBytes[_world.size()],
				sendDisp[_world.size()], recvDisp[_world.size()];

		ReadIdxSet sendingGlobalContigIdx;

		ReadSet empty;
		int emptyBase = empty.getStoreSize() * query.getGlobalSize();
		int maxRankTransmitSize1 = ( MPIOptions::getOptions().getTotalBufferSize() / _world.size() / 2);
		if (maxRankTransmitSize1 < 10240)
			maxRankTransmitSize1 = 10240;
		int maxRankTransmitSize = maxRankTransmitSize1 + emptyBase / _world.size() + 1024;
		long maxBuffer = maxRankTransmitSize * _world.size();
		bool isDone = false;
		long iteration = 0;
		while (! isDone ) {
			iteration++;
			long localReadSets = 0;
			LOG_DEBUG_OPTIONAL(1, _world.rank() == 0, "exchangeGlobalReads: " << iteration << ". Maximum transmit buffer is: " << maxRankTransmitSize1 << ", " << maxRankTransmitSize << " / " << maxBuffer);
			isDone = true;
			int totalSend = 0, totalRecv = 0;
			for (int rank = 0; rank < _world.size(); rank++) {
				sendBytes[rank] = 1; // artificially pad by one byte to ensure negative 'isDone' flag can be sent
				recvBytes[rank] = 0;
			}

			for (ReadSet::ReadSetSizeType globalContigIdx = 0; globalContigIdx < query.getGlobalSize(); globalContigIdx++) {
				int rank;
				ReadSet::ReadSetSizeType rankReadIdx;
				query.getRankReadForGlobalReadIdx(globalContigIdx, rank, rankReadIdx);
				if (rank == myRank) {
					// do not encode and send reads to self
					assert(query.isLocalRead(globalContigIdx));
					ReadSet::ReadSetSizeType localIdx = query.getLocalReadIdx(globalContigIdx);
					globalReadSetVector[localIdx].append(localReadSetVector[globalContigIdx]);
					localReadSets++;
					continue;
				}
				int sendByteCount = localReadSetVector[globalContigIdx].getStoreSize();
				while (sendByteCount >= maxRankTransmitSize1) {
					ReadSet::ReadSetSizeType size = localReadSetVector[globalContigIdx].getSize();
					ReadSet::ReadSetSizeType maxSend = size * maxRankTransmitSize1 * 3 / sendByteCount / 2;
					if (maxSend <= 0) {
						LOG_WARN(1, "Could not send any reads for " << globalContigIdx << " as buffer size is too small");
						localReadSetVector[globalContigIdx].clear();
						tmpReadSetVector[globalContigIdx].clear();
					} else {
						LOG_DEBUG(4, "globalContig: " << globalContigIdx << " with " << localReadSetVector[globalContigIdx].getSize() << " reads is too big to transmit, reducing to " << maxSend);
						// temporarily store in globalReadSetVector
						tmpReadSetVector[globalContigIdx].append(localReadSetVector[globalContigIdx].truncate(maxSend));
					}
					sendByteCount = localReadSetVector[globalContigIdx].getStoreSize();
				}
				if (sendBytes[rank] + sendByteCount < maxRankTransmitSize1) {
					sendingGlobalContigIdx.insert(globalContigIdx);
				} else {
					isDone = false;
					long requiredSendByteCount = sendByteCount;
					sendByteCount = empty.getStoreSize();
					ReadIdxSet::iterator readIdIt = sendingGlobalContigIdx.find(globalContigIdx);
					if (readIdIt != sendingGlobalContigIdx.end())
						sendingGlobalContigIdx.erase(readIdIt);
					LOG_DEBUG(3, "Skipping globalContigIdx this round(" << iteration << "): " << globalContigIdx << " need " << requiredSendByteCount << " prepared: " << sendBytes[rank] << " remaining: " << maxRankTransmitSize - sendBytes[rank]);
				}
				if (tmpReadSetVector[globalContigIdx].getSize() > 0) {
					// some reads are waiting for a future communication
					isDone = false;
				}
				sendBytes[rank] += sendByteCount;
				assert(sendBytes[rank] <= maxRankTransmitSize + emptyBase && sendBytes[rank] >= 0); // i.e. no overflow
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
			LOG_DEBUG_OPTIONAL(1, true, "iteration " << iteration << " Sending bytes: " << totalSend << " readSets: " << sendingGlobalContigIdx.size() + localReadSets << (isDone ? " isDone" : ""));

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
				LOG_THROW("Could not allocate totalSend bytes! " << totalSend);

			char *tmp = sendBuf;
			for (ReadSet::ReadSetSizeType globalContigIdx = 0; globalContigIdx
			< query.getGlobalSize(); globalContigIdx++) {
				int storeSize = 0;
				if (query.isLocalRead(globalContigIdx)) {
					// do not encode and send reads to self
				} else {
					if (sendingGlobalContigIdx.find(globalContigIdx) != sendingGlobalContigIdx.end()) {
						storeSize = localReadSetVector[globalContigIdx].store(tmp);
						localReadSetVector[globalContigIdx].clear();
						// swap back any truncated reads for next run
						tmpReadSetVector[globalContigIdx].swap(localReadSetVector[globalContigIdx]);
					} else {
						storeSize = empty.store(tmp);
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

		ReadSet::ReadSetSizeType maxReads = MatcherInterfaceOptions::getOptions().getMaxReadMatches();
		if (maxReads > 0) {
			for (ReadSet::ReadSetSizeType localContigIdx = 0; localContigIdx
						< globalReadSetVector.size(); localContigIdx++) {
				ReadSet::ReadSetSizeType numReads = globalReadSetVector[localContigIdx].getSize();
				if (numReads > maxReads) {
					LOG_DEBUG_OPTIONAL(2, true, "for " << localContigIdx << " sampled from " << numReads << " to " << maxReads);
					globalReadSetVector[localContigIdx] = globalReadSetVector[localContigIdx].randomlySample(maxReads);
				}
			}
		}
		LOG_DEBUG_OPTIONAL(1, true, "exchangeGlobalReads(): Done " << iteration);

		return globalReadSetVector;
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
		std::vector< std::vector< std::vector< ReadSet::ReadSetSizeType > > > rankHitGlobalIndexes(numRanks);
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
			}
			if(i > 0) {
				sendRankDispl[i] = sendRankDispl[i-1] + sendRankCount[i-1];
			}
		}

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
				for(int k = 0; k < sendCounts[pos]; k++) {
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

		// consolidate into localMatchResults
		MatchResults localMatchResults(globalMatchResults.size());
		for(int i = 0; i < numRanks; i++) {
			tmp = recvBuf + recvRankDispl[i];
			for(int j = 0; j < numMatchHitSets ; j++) {
				long pos = i * numMatchHitSets + j;
				for(int k = 0; k < recvCounts[pos]; k++) {
					localMatchResults[j].insert( *tmp );
					tmp++;
				}
			}
			assert(tmp == recvBuf + recvRankDispl[i] + recvRankCount[i]);
		}
		sampleMatches(localMatchResults);

		delete [] recvRankDispl;
		delete [] recvRankCount;
		delete [] recvCounts;
		delete [] recvBuf;
		recordTime("EGRI", MPI_Wtime());

		assert(areAllReadsLocal(localMatchResults));
		return localMatchResults;
	}

	mpi::communicator &getWorld() { return _world; }
	const ReadSet &getTarget() const { return _target; }
	bool isReturnPairedMatches() const { return _returnPairedMatches; }

protected:
	mpi::communicator &_world;
	const ReadSet &_target;
	bool _returnPairedMatches;
};

#endif /* MATCHERINTERFACE_H_ */
