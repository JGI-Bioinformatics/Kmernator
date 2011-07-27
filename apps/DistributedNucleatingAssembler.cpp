//
// DistributedNucleatingAssembler.cpp
// Author: Rob Egan
//
// Copyright 2011 The Regents of the University of California.
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
// binaries display the following acknowledgment:  "This product
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

#include "config.h"
#include "ReadSet.h"
#include "Options.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "DuplicateFragmentFilter.h"
#include "ContigExtender.h"
#include "DistributedFunctions.h"
#include "Log.h"

using namespace std;
typedef TrackingDataMinimal4f DataType;
typedef KmerSpectrum<DataType, DataType> KS;


class DistributedNucleatingAssemblerOptions : public ContigExtenderOptions {
public:
	static std::string getVmatchOptions() {
		return getVarMap()["vmatch-options"].as<std::string>();
	}
	static std::string getVmatchIndexPath() {
		return getVarMap()["vmatch-index-path"].as<std::string>();
	}
	static int getMaxIterations() {
		return getVarMap()["max-iterations"].as<int>();
	}
	static bool parseOpts(int argc, char *argv[]) {

		getDesc().add_options()

		("max-iterations", po::value<int>()->default_value(1000),
				"the maximum number of rounds to extend the set of contigs")

		("vmatch-options", po::value<std::string>()->default_value("-d -p -seedlength 10 -l 50 -e 3"),
				"options with which to call vmatch")

		("vmatch-index-path", po::value<std::string>()->default_value("."),
				"top level directory under which to create the vmatch index directories for each rank")

				;

		bool ret = ContigExtenderOptions::parseOpts(argc, argv);

		if (Options::getOutputFile().empty()) {
			LOG_WARN(1, "You must specify an --output");
			ret = -1;
		}

		return ret;
	}
};


#include <stdlib.h>
#include "Utils.h"
class Vmatch {
public:
	class FieldsType {
	public:

		FieldsType(std::string &_line) : line(_line) {
			static const char *format = "%d %d %d %c %d %d %d %d %e %d %f";
			sscanf(line.c_str(), format, &subjectLength, &subjectNumber, &subjectPosition, &type, &queryLength, &queryNumber, &queryPosition, &distance, &eValue, &scoreValue, &percentIdentity);
		}

		int subjectLength;
		int subjectNumber;
		int subjectPosition;
		char type; // D or P;
		int queryLength;
		int queryNumber;
		int queryPosition;
		int distance;
		float eValue;
		int scoreValue;
		float percentIdentity;
		std::string line;
	};
	typedef std::vector<FieldsType> MatchResults;

private:
	std::string _indexName;
	MatchResults _results;
public:
	Vmatch(std::string indexName, ReadSet &inputs) : _indexName(indexName) {
		if (FileUtils::getFileSize(indexName + ".suf") > 0) {
			LOG_VERBOSE_OPTIONAL(1, true, "Vmatch(" << indexName << ", reads(" << inputs.getSize() << ")): Index already exists: " << indexName);
		} else {
			buildVmatchIndex(indexName, inputs);
		}
	}
	static void buildVmatchIndex(std::string indexName, ReadSet &inputs) {
		std::string inputFile = indexName + ".tmp.input";
		OfstreamMap ofm(inputFile, "");
		inputs.writeAll(ofm.getOfstream(""), FormatOutput::FASTA);
		ofm.clear();
		buildVmatchIndex(indexName, inputFile);
		LOG_DEBUG_OPTIONAL(1, true, "Removing temporary inputFile" << inputFile);
		unlink(inputFile.c_str());
	}
	static void buildVmatchIndex(std::string indexName, std::string inputFasta) {
		std::string cmd("mkvtree -dna -allout -pl -indexname " + indexName + " -db " + inputFasta);
		LOG_VERBOSE_OPTIONAL(1, true, "Building vmatch index " << indexName << " : " << cmd);
		int ret = system(cmd.c_str());
		if (ret != 0)
			LOG_THROW("mkvtree failed to build(" << ret << "): " << cmd);
	}
	MatchResults &match(std::string queryFile, std::string options = "") {
		_results.clear();
		std::string cmd("vmatch " + options + " -q " + queryFile + " " + _indexName);
		LOG_DEBUG_OPTIONAL(1, true, "Executing vmatch: " << cmd);
		IPipestream vmatchOutput(cmd);
		std::string line;
		while (!vmatchOutput.eof()) {
			getline(vmatchOutput, line);
			if (line.length() == 0)
				break;
			if (line[0] == '#')
				continue;
			_results.push_back(FieldsType(line));
		}
		LOG_DEBUG_OPTIONAL(1, true, "Vmatch::match(,): Found " << _results.size() << " results");
		return _results;
	}
};


int main(int argc, char *argv[]) {
	// do not apply artifact filtering by default
	Options::getSkipArtifactFilter() = 1;
	// override the default output format!
	Options::getFormatOutput() = 3;
	Options::getMmapInput() = 0;
	Options::getVerbosity() = 1;
	Options::getMaxThreads() = 1;

	int threadProvided;
	int threadRequest = omp_get_max_threads() == 1 ? MPI_THREAD_SINGLE : MPI_THREAD_FUNNELED;
	MPI_Init_thread(&argc, &argv, threadRequest, &threadProvided);
	mpi::environment env(argc, argv);
	mpi::communicator world;
	MPI_Comm_set_errhandler( world, MPI::ERRORS_THROW_EXCEPTIONS );

	try {
		Logger::setWorld(&world);

		if (!DistributedNucleatingAssemblerOptions::parseOpts(argc, argv))
			throw std::invalid_argument("Please fix the command line arguments");

		if (Options::getGatheredLogs())
			Logger::setWorld(&world, Options::getDebug() >= 2);

	} catch (...) {
		std::cerr << DistributedNucleatingAssemblerOptions::getDesc() << std::endl;
		std::cerr << std::endl << "Please fix the options and/or MPI environment" << std::endl;
		exit(1);
	}
	world.barrier();

	long maxIterations = DistributedNucleatingAssemblerOptions::getMaxIterations();

	Options::FileListType inputFiles = Options::getInputFiles();
	std::string contigFile = ContigExtenderOptions::getContigFile();

	ReadSet reads;
	LOG_VERBOSE(1, "Reading Input Files" );
	reads.appendAllFiles(inputFiles, world.rank(), world.size());

	LOG_VERBOSE(2, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount() << " Bases ");
	setGlobalReadSetOffsets(world, reads);

	std::string tmpDir = DistributedNucleatingAssemblerOptions::getVmatchIndexPath() + "/";
	if (world.rank() == 0)
		mkdir(tmpDir.c_str(), 0777);
	world.barrier();

	std::string rankOutputDir = getRankSubdir(world, tmpDir);

	Vmatch myVmatch = Vmatch(rankOutputDir + "/myReads", reads);

	SequenceLengthType minKmerSize, maxKmerSize, kmerStep;
	ContigExtender<KS>::getMinMaxKmerSize(reads, minKmerSize, maxKmerSize, kmerStep);
	maxKmerSize = boost::mpi::all_reduce(world, maxKmerSize, mpi::minimum<SequenceLengthType>());

	ReadSet finalContigs;
	short iteration = 0;
	while(++iteration <= maxIterations) {
		ReadSet contigs;
		contigs.appendFastaFile(contigFile, world.rank(), world.size());

		setGlobalReadSetOffsets(world, contigs);
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Iteration: " << iteration << ". Reading Contig File: " << contigFile << ". loaded " << contigs.getGlobalSize() << " Reads");
		if (contigs.getGlobalSize() == 0) {
			LOG_VERBOSE_OPTIONAL(1, true, "There are no contigs to extend in " << contigFile);
			break;
		}

		ReadSet::ReadSetVector contigReadSet(contigs.getGlobalSize(), ReadSet());
		Vmatch::MatchResults matches = myVmatch.match(contigFile, DistributedNucleatingAssemblerOptions::getVmatchOptions());

		std::vector< std::set<long> > contigReadHits(contigs.getGlobalSize());
		for(Vmatch::MatchResults::iterator match = matches.begin(); match != matches.end(); match++) {
			ReadSet::ReadSetSizeType globalContigIdx = match->queryNumber, readIdx = match->subjectNumber;
			contigReadHits[globalContigIdx].insert(readIdx);
			// include pairs
			if (readIdx % 2 == 0)
				contigReadHits[globalContigIdx].insert( readIdx + 1);
			else
				contigReadHits[globalContigIdx].insert( readIdx - 1 );
		}
		for(long globalContigIdx = 0; globalContigIdx < (long)contigReadHits.size(); globalContigIdx++) {
			LOG_DEBUG_OPTIONAL(3, true, "GlobalContig: " << globalContigIdx << " has " << contigReadHits[globalContigIdx].size() << " contigReadHits");
			for(std::set<long>::iterator it2 = contigReadHits[globalContigIdx].begin(); it2 != contigReadHits[globalContigIdx].end(); it2++) {
				long readIdx = *it2;
				contigReadSet[ globalContigIdx ].append( reads.getRead( readIdx ) );
			}
		}
		contigReadHits.clear();

		int sendBytes[world.size()], recvBytes[world.size()], sendDisp[world.size()], recvDisp[world.size()];
		int totalSend = 0, totalRecv = 0;

		for(int rank = 0; rank < world.size(); rank++) {
			sendBytes[rank] = 0;
			recvBytes[rank] = 0;
		}
		for(ReadSet::ReadSetSizeType globalContigIdx = 0; globalContigIdx < contigs.getGlobalSize(); globalContigIdx++) {
			int rank;
			ReadSet::ReadSetSizeType rankReadIdx;
			contigs.getRankReadForGlobalReadIdx(globalContigIdx, rank, rankReadIdx);
			int sendByteCount = contigReadSet[ globalContigIdx ].getStoreSize();
            sendBytes[ rank ] += sendByteCount;
            LOG_DEBUG_OPTIONAL(3, true, "GlobalContig: " << globalContigIdx << " sending " << contigReadSet[ globalContigIdx ].getSize() << " reads / "<< sendByteCount << " bytes to " << rank << " total: " << sendBytes[rank]);
		}
		for(int rank = 0; rank < world.size(); rank++) {
			sendDisp[rank] = totalSend;
			totalSend += sendBytes[rank];
			LOG_DEBUG_OPTIONAL(3, true, "Sending to rank " << rank << " sendBytes " << sendBytes[rank] << " sendDisp[] " << sendDisp[rank] << ". 0 == recvBytes[] "<< recvBytes[rank]);
		}

		MPI_Alltoall(sendBytes, 1, MPI_INT, recvBytes, 1, MPI_INT, world);

		for(int rank = 0; rank < world.size(); rank++) {
			recvDisp[rank] = totalRecv;
			totalRecv += recvBytes[rank];
			LOG_DEBUG_OPTIONAL(3, true, "to/from rank " << rank << ": sendDisp[] " << sendDisp[rank] << " sendBytes[] " << sendBytes[rank] << " recvDisp[] " << recvDisp[rank] << " recvBytes[] " << recvBytes[rank] );
		}
		LOG_DEBUG_OPTIONAL(3, true, "Sending " << totalSend << " Receiving " << totalRecv);
		world.barrier();
		char *sendBuf = (char*) malloc(totalSend == 0 ? 8 : totalSend);
		if (sendBuf == NULL)
			LOG_THROW("Could not allocate totalSend bytes! " << totalSend);
		char *tmp = sendBuf;
		for(ReadSet::ReadSetSizeType contigGlobalIndex = 0; contigGlobalIndex < contigs.getGlobalSize(); contigGlobalIndex++) {
			tmp += contigReadSet[contigGlobalIndex].store(tmp);
		}
		contigReadSet.clear();
		char *recvBuf = (char*) malloc(totalRecv == 0 ? 8 : totalRecv);
		if (recvBuf == NULL)
			LOG_THROW("Could not allocate totalRecv bytes! " << totalRecv);

		MPI_Alltoallv(sendBuf, sendBytes, sendDisp, MPI_BYTE, recvBuf, recvBytes, recvDisp, MPI_BYTE, world);
		free(sendBuf);
		// set contigReadSet to hold read sets for local set of contigs
		contigReadSet.resize(contigs.getSize(), ReadSet());

		for(int rank = 0; rank < world.size(); rank++) {
			ReadSet::ReadSetSizeType contigIdx = 0;
			tmp = recvBuf + recvDisp[rank];
			while (tmp != recvBuf + recvDisp[rank] + recvBytes[rank]) {
				ReadSet rs;
				tmp = (char*) rs.restore(tmp);
				contigReadSet[contigIdx].append(rs);
				LOG_DEBUG_OPTIONAL(3, true, "Restored from rank " << rank << " ReadSet " << rs.getSize() << " totaling " << contigReadSet[contigIdx].getSize());
				contigIdx++;
			}
		}
		free(recvBuf);

		ReadSet changedContigs;

        //#pragma omp parallel for
		for(ReadSet::ReadSetSizeType i = 0; i < contigs.getSize(); i++) {
			const Read &oldRead = contigs.getRead(i);
			LOG_VERBOSE_OPTIONAL(2, true, "Extending " << oldRead.getName() << " with " << contigReadSet[i].getSize() << " pool of reads");
			ReadSet myContig;
			myContig.append(oldRead);
			ReadSet newContig;
			SequenceLengthType myKmerSize = minKmerSize;
			SequenceLengthType newLen = 0;
			while (newLen <= oldRead.getLength() && myKmerSize <= maxKmerSize) {
				newContig = ContigExtender<KS>::extendContigs(myContig, contigReadSet[i], myKmerSize, myKmerSize);
				newLen = newContig.getRead(0).getLength();
				myKmerSize += kmerStep;
			}
			const Read &newRead = newContig.getRead(0);

			//#pragma omp critical
			{
				long deltaLen = newRead.getLength() - oldRead.getLength();
				LOG_VERBOSE_OPTIONAL(1, true, "Extended " << oldRead.getName() << " " << deltaLen << " bases to " << newRead.getLength() << ": " << newRead.getName() << " with " << contigReadSet[i].getSize() << " reads in the pool");
				if (deltaLen > 0) {
					changedContigs.append(newRead);
				} else {
					finalContigs.append(oldRead);
				}
			}
		}

		LOG_DEBUG_OPTIONAL(1, true, "Changed contigs: " << changedContigs.getSize() << " finalContigs: " << finalContigs.getSize());
		setGlobalReadSetOffsets(world, changedContigs);
		setGlobalReadSetOffsets(world, finalContigs);
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Changed contigs: " << changedContigs.getGlobalSize() << " finalContigs: " << finalContigs.getGlobalSize());

		{
			// write out the final contigs (so far) so we do not loose them
			DistributedOfstreamMap om(world, Options::getOutputFile(), "");
			om.setBuildInMemory();
			finalContigs.writeAll(om.getOfstream(""), FormatOutput::FASTA);
		}
		if (!Log::isDebug(1)) {
			// remove most recent contig file (if not debugging)
			if (ContigExtenderOptions::getContigFile().compare(contigFile) != 0) {
				if (world.rank() == 0) {
					LOG_VERBOSE_OPTIONAL(2, true, "Removing " << contigFile);
					unlink(contigFile.c_str());
				}
			}
		}

		if (changedContigs.getGlobalSize() == 0) {
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 1, "No more contigs to extend " << changedContigs.getSize());
			break;
		}

		{
			std::string filekey = "contig-" + boost::lexical_cast<std::string>(iteration);
			DistributedOfstreamMap om(world, Options::getOutputFile(), FormatOutput::getSuffix(FormatOutput::FASTA));
			om.setBuildInMemory();
			changedContigs.writeAll(om.getOfstream(filekey), FormatOutput::FASTA);
			std::string newContigFile = om.getRealFilePath(filekey);
			contigFile = newContigFile;
		}
	}

	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");
	MPI_Finalize();

	return 0;
}

