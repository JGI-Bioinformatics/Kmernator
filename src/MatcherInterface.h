/*
 * MatcherInterface.h
 *
 *  Created on: Nov 1, 2011
 *      Author: regan
 */

#ifndef MATCHERINTERFACE_H_
#define MATCHERINTERFACE_H_

#include <vector>
#include <set>

#include "ReadSet.h"
#include "Utils.h"
#include "Options.h"

class MatcherInterface : public Timer {
public:
	typedef std::set< ReadSet::ReadSetSizeType > MatchHitSet;
	typedef std::vector< MatchHitSet > MatchResults;
	typedef ReadSet::ReadSetVector MatchReadResults;

	MatcherInterface(mpi::communicator &world, const ReadSet &target, bool returnPairedMatches = true)
	: _world(world), _target(target), _returnPairedMatches(returnPairedMatches) {
	}

	// returns a list of sets of local reads (targets in the globalReadIdx space) that match the query
	virtual MatchResults matchLocal(const ReadSet &query) {
		std::string queryFile = UniqueName::generateUniqueName( GeneralOptions::getOptions().getTmpDir() + "/MatcherInterface-" );
		{
			OfstreamMap om(queryFile);
			query.writeAll(om.getOfstream(""));
		}
		MatchResults results = this->matchLocal(queryFile);
		assert(results.size() == query.getSize());
		unlink(queryFile.c_str());
		return results;
	}
	virtual MatchResults matchLocal(std::string queryFile) {
		throw;
	}

	// returns a ReadSetVector of global reads (copied from each local instance) matching the query
	virtual MatchReadResults match(const ReadSet &query, std::string queryFile) {
		recordTime("startMatch", MPI_Wtime());
		MatchResults matchResults;
		if (query.getGlobalSize() != query.getSize()){
			matchResults = this->matchLocal(queryFile);
		} else {
			matchResults = this->matchLocal(query);
		}
		MatchReadResults mrr = exchangeGlobalReads(query, getLocalReads(matchResults, query));
		recordTime("returnMatch", MPI_Wtime());
		return mrr;
	}
	virtual MatchReadResults match(const ReadSet &query) {
		if (query.getSize() != query.getGlobalSize())
			LOG_THROW("Can not run MatcherInterface::match(ReadSet&) on global ReadSet (yet)");
		recordTime("startMatch", MPI_Wtime());
		MatchResults matchResults = this->matchLocal(query);
		MatchReadResults mrr = getLocalReads(matchResults, query);
		recordTime("returnMatch", MPI_Wtime());
		return mrr;
	}

	// returns a ReadSetVector of reads that are local (targets in globalReadIdx space)
	MatchReadResults getLocalReads(MatchResults &matchResults, const ReadSet &query) {
		assert(query.getGlobalSize() == matchResults.size());
		recordTime("getLocalReads", MPI_Wtime());
		if (!areAllReadsLocal(matchResults)) {
			matchResults = exchangeGlobalReadIdxs(matchResults);
		}

		int myRank = _world.rank();
		MatchReadResults matchReadResults(matchResults.size(), ReadSet());
		for(int i = 0; i < (int) matchResults.size(); i++) {
			for(MatchHitSet::iterator it = matchResults[i].begin(); it != matchResults[i].end(); it++) {
				ReadSet::ReadSetSizeType localReadIdx = getTarget().getLocalReadIdx(myRank, *it);
				matchReadResults[i].append( getTarget().getRead( localReadIdx ));
			}
		}
		recordTime("returnLocalReads", MPI_Wtime());
		return matchReadResults;
	}
	// transfers the proper matching localRedSetVector to
	// a globalReadSetVector to the node controlling the global reads over the query
	// reads will be copied from the localReadSet to the destination node
	MatchReadResults exchangeGlobalReads(const ReadSet &query, const MatchReadResults &localReadSetVector) {
		assert(query.getSize() == localReadSetVector.size());
		MatchReadResults globalReadSetVector;

		int sendBytes[_world.size()], recvBytes[_world.size()],
				sendDisp[_world.size()], recvDisp[_world.size()];
		int totalSend = 0, totalRecv = 0;

		for (int rank = 0; rank < _world.size(); rank++) {
			sendBytes[rank] = 0;
			recvBytes[rank] = 0;
		}
		for (ReadSet::ReadSetSizeType globalContigIdx = 0; globalContigIdx
				< query.getGlobalSize(); globalContigIdx++) {
			int rank;
			ReadSet::ReadSetSizeType rankReadIdx;
			query.getRankReadForGlobalReadIdx(globalContigIdx, rank, rankReadIdx);
			int sendByteCount = localReadSetVector[globalContigIdx].getStoreSize();
			sendBytes[rank] += sendByteCount;
			LOG_DEBUG_OPTIONAL(3, true, "GlobalContig: " << globalContigIdx << " sending " << localReadSetVector[ globalContigIdx ].getSize() << " reads / "<< sendByteCount << " bytes to " << rank << " total: " << sendBytes[rank]);
		}
		for (int rank = 0; rank < _world.size(); rank++) {
			sendDisp[rank] = totalSend;
			totalSend += sendBytes[rank];
			LOG_DEBUG_OPTIONAL(3, true, "Sending to rank " << rank << " sendBytes " << sendBytes[rank] << " sendDisp[] " << sendDisp[rank] << ". 0 == recvBytes[] "<< recvBytes[rank]);
		}

		recordTime("EGR-prepareSizes",MPI_Wtime());;

		MPI_Alltoall(sendBytes, 1, MPI_INT, recvBytes, 1, MPI_INT, _world);

		recordTime("EGR-exchangeSizesTime",MPI_Wtime());

		for (int rank = 0; rank < _world.size(); rank++) {
			recvDisp[rank] = totalRecv;
			totalRecv += recvBytes[rank];
			LOG_DEBUG_OPTIONAL(3, true, "to/from rank " << rank << ": sendDisp[] " << sendDisp[rank] << " sendBytes[] " << sendBytes[rank] << " recvDisp[] " << recvDisp[rank] << " recvBytes[] " << recvBytes[rank] );
		}
		LOG_DEBUG_OPTIONAL(3, true, "Sending " << totalSend << " Receiving " << totalRecv);

		char *sendBuf = (char*) (malloc(totalSend == 0 ? 8 : totalSend));
		if (sendBuf == NULL)
			LOG_THROW("Could not allocate totalSend bytes! " << totalSend);
		char *tmp = sendBuf;
		for (ReadSet::ReadSetSizeType contigGlobalIndex = 0; contigGlobalIndex
				< query.getGlobalSize(); contigGlobalIndex++) {
			tmp += localReadSetVector[contigGlobalIndex].store(tmp);
		}


		char *recvBuf = (char*) (malloc(totalRecv == 0 ? 8 : totalRecv));
		if (recvBuf == NULL)
			LOG_THROW("Could not allocate totalRecv bytes! " << totalRecv);

		recordTime("EGR-prepareSendTime", MPI_Wtime());
		MPI_Alltoallv(sendBuf, sendBytes, sendDisp, MPI_BYTE, recvBuf,
				recvBytes, recvDisp, MPI_BYTE, _world);
		recordTime("EGR-exchangeReadsTime", MPI_Wtime());

		free(sendBuf);
		// set localReadSetVector to hold read sets for local set of query
		globalReadSetVector.resize(query.getSize(), ReadSet());
		for (int rank = 0; rank < _world.size(); rank++) {
			ReadSet::ReadSetSizeType contigIdx = 0;
			tmp = recvBuf + recvDisp[rank];
			while (tmp != recvBuf + recvDisp[rank] + recvBytes[rank]) {
				ReadSet rs;
				tmp = (char*) rs.restore(tmp);
				globalReadSetVector[contigIdx].append(rs);
				assert (query.isLocalRead(rank, contigIdx) || globalReadSetVector[contigIdx].getSize() == 0);
				LOG_DEBUG_OPTIONAL(3, true, "Restored from rank " << rank << " ReadSet " << rs.getSize() << " totaling " << globalReadSetVector[contigIdx].getSize());
				contigIdx++;
			}
		}
		free(recvBuf);

		return globalReadSetVector;
	}
	bool areAllReadsLocal(MatchResults &matchResults) {
		int myRank = _world.rank();
		for(MatchResults::iterator it = matchResults.begin(); it != matchResults.end(); it++) {
			for(MatchHitSet::iterator it2 = it->begin(); it2 != it->end(); it2++) {
				if (!getTarget().isLocalRead(myRank, *it2)) {
					return false;
				}
			}
		}
		return true;
	}
	MatchResults exchangeGlobalReadIdxs(const MatchResults &globalMatchResults) {
		int numMatchHitSets = globalMatchResults.size();
		int numRanks = _world.size();
		// sort globalReadIdxs by rank, matchHitSet
		std::vector< std::vector< std::vector< ReadSet::ReadSetSizeType > > > rankHitGlobalIndexes(numRanks);
		for(long i = 0; i < numRanks; i++)
			rankHitGlobalIndexes[i].resize(numMatchHitSets);
		for(long i = 0; i < numMatchHitSets; i++) {
			for(MatchHitSet::iterator it = globalMatchResults[i].begin(); it != globalMatchResults[i].end(); it++) {
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
		recordTime("EGRI-prepareSizes", MPI_Wtime());
		MPI_Alltoall(sendCounts, numMatchHitSets, MPI_INT, recvCounts, numMatchHitSets, MPI_INT, _world);
		recordTime("EGRI-exchangeSizes", MPI_Wtime());
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
		for(int i = 0; i < numRanks; i++) {
			for(int j = 0; j < numMatchHitSets ; j++) {
				long pos = i * numMatchHitSets + j;
				for(int k = 0; k < sendCounts[pos]; k++) {
					*(sendBuf + sendRankDispl[i] + k) = rankHitGlobalIndexes[i][j][k];
				}
			}
		}
		// free some memory before using more
		delete [] sendCounts;
		rankHitGlobalIndexes.clear();
		// allocate receiving buffer
		ReadSet::ReadSetSizeType *recvBuf = new ReadSet::ReadSetSizeType[totalRecvCount];
		recordTime("EGRI-PrepareReads", MPI_Wtime());
		// exchange globalReadIdx to rank that owns it
		MPI_Alltoallv(sendBuf, sendRankCount, sendRankDispl, MPIReadSetSizeType, recvBuf, recvRankCount, recvRankDispl, MPIReadSetSizeType, _world);
		recordTime("EGRI-ExchangeReads", MPI_Wtime());
		delete [] sendBuf;
		delete [] sendRankCount;
		delete [] sendRankDispl;

		// consolidate into localMatchResults
		MatchResults localMatchResults(globalMatchResults.size());
		for(int i = 0; i < numRanks; i++) {
			for(int j = 0; j < numMatchHitSets ; j++) {
				long pos = i * numMatchHitSets + j;
				for(int k = 0; k < recvCounts[pos]; k++) {
					localMatchResults[j].insert( *(recvBuf + recvRankDispl[i] + k) );
				}
			}
		}

		delete [] recvRankDispl;
		delete [] recvRankCount;
		delete [] recvCounts;
		delete [] recvBuf;

		assert(areAllReadsLocal(localMatchResults));
		return localMatchResults;
	}

	mpi::communicator &getWorld() { return _world; }
	const ReadSet &getTarget() { return _target; }
	bool isReturnPairedMatches() { return _returnPairedMatches; }

protected:
	mpi::communicator &_world;
	const ReadSet &_target;
	bool _returnPairedMatches;
};

#endif /* MATCHERINTERFACE_H_ */
