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
	static bool parseOpts(int argc, char *argv[]) {

		bool ret = ContigExtenderOptions::parseOpts(argc, argv);

		return ret;
	}
};

#include <stdlib.h>
class Vmatch {
public:
	class FieldsType {
	public:

		FieldsType(std::string &_line) : line(_line) {
			static const char *format = "%d %d %d %c %d %d %d %d %f %d %f";
			sscanf(line.c_str(), format, subjectLength, subjectNumber, subjectPosition, type, queryLength, queryNumber, queryPosition, distance, scoreValue, percentIdentity);
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
		std::string inputFile = indexName + ".tmp.input";
		OfstreamMap ofm(inputFile, "");
		inputs.writeAll(ofm.getOfstream(""), FormatOutput::FASTA);
		ofm.clear();
		std::string cmd("mkvtree -dna -allout -pl -indexname " + indexName + " -db " + inputFile);
		LOG_VERBOSE(1, "Building vmatch index " << indexName << " : " << cmd);
		int ret = system(cmd.c_str());
		if (ret != 0)
			throw;
		LOG_VERBOSE(1, "Removing temporary inputFile" << inputFile);
		unlink(inputFile.c_str());
	}
	MatchResults &match(std::string queryFile, std::string options = "-d -p -seedlength 10 -l 50 -e 3") {
		_results.clear();
		IPipestream vmatchOutput("vmatch " + options + " -q " + queryFile + " " + _indexName);
		std::string line;
		while (!vmatchOutput.eof()) {
			getline(vmatchOutput, line);
			if (line.length() == 0)
				break;
			if (line[0] == '#')
				continue;
			_results.push_back(FieldsType(line));
		}
		return _results;
	}
};

int main(int argc, char *argv[]) {
	// do not apply artifact filtering by default
	Options::getSkipArtifactFilter() = 1;
	// override the default output format!
	Options::getFormatOutput() = 3;
	Options::getMmapInput() = 0;
	Options::getVerbosity() = 2;

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

	Options::FileListType inputFiles = Options::getInputFiles();
	std::string contigFile = ContigExtenderOptions::getContigFile();

	ReadSet reads;
	LOG_VERBOSE(1, "Reading Input Files" );
	reads.appendAllFiles(inputFiles, world.rank(), world.size());

	LOG_VERBOSE(2, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount() << " Bases ");
	setGlobalReadSetOffsets(world, reads);

	std::string tmpDir = Options::getOutputFile() + "-tmp." + boost::lexical_cast<std::string>(world.size());
	if (world.rank() == 0)
		mkdir(tmpDir.c_str(), 0x777);
	world.barrier();

	std::string rankOutputDir = getRankSubdir(world, tmpDir);

	Vmatch myVmatch = Vmatch(rankOutputDir + "/myReads", reads);

	ReadSet finalContigs;
	int iteration = 0;
	while(++iteration) {
		ReadSet contigs;
		LOG_VERBOSE(1, "Reading Contig File" );
		contigs.appendFastaFile(contigFile, world.rank(), world.size());

		LOG_VERBOSE(2, "loaded " << contigs.getSize() << " Reads, " << contigs.getBaseCount() << " Bases ");
		setGlobalReadSetOffsets(world, contigs);
		if (contigs.getGlobalSize() == 0)
			break;

		ReadSet::ReadSetVector contigReadSet(contigs.getGlobalSize(), ReadSet());
		Vmatch::MatchResults matches = myVmatch.match(contigFile);
		std::vector< std::set<long> > contigReadHits;
		for(Vmatch::MatchResults::iterator match = matches.begin(); match != matches.end(); match++) {
			ReadSet::ReadSetSizeType contigIdx = match->queryNumber, readIdx = match->subjectNumber;
			contigReadHits[contigIdx].insert(readIdx);
		}
		for(long contigIdx = 0; contigIdx < (long)contigReadHits.size(); contigIdx++) {
			for(std::set<long>::iterator it2 = contigReadHits[contigIdx].begin(); it2 != contigReadHits[contigIdx].end(); it2++) {
				long readIdx = *it2;
				contigReadSet[ contigIdx ].append( reads.getRead( readIdx ) );
			}
		}

		int sendBytes[world.size()], recvBytes[world.size()], sendDisp[world.size()], recvDisp[world.size()], totalSend = 0, totalRecv = 0;
		for(int i = 0; i < world.size(); i++) {
			sendBytes[i] = 0;
			recvBytes[i] = 0;
		}
		for(ReadSet::ReadSetSizeType i = 0; i < contigs.getSize(); i++) {
			int rank;
			ReadSet::ReadSetSizeType rankReadIdx;
			contigs.getRankReadForGlobalReadIdx(i, rank, rankReadIdx);
			sendBytes[ rank ] += contigReadSet[ i ].getStoreSize();
		}
		MPI_Alltoall(sendBytes, world.size(), MPI_INT, recvBytes, world.size(), MPI_LONG, world);
		for(int i = 0; i < world.size(); i++) {
			sendDisp[i] = totalSend;
			totalSend += sendBytes[i];
			recvDisp[i] = totalRecv;
			totalRecv += recvBytes[i];
		}
		char *sendBuf = (char*) malloc(totalSend);
		if (sendBuf == NULL)
			throw;
		char *tmp = sendBuf;
		for(ReadSet::ReadSetSizeType i = 0; i < contigs.getSize(); i++) {
			tmp += contigReadSet[i].store(tmp);
		}
		contigReadSet.clear();
		char *recvBuf = (char*) malloc(totalRecv);
		if (recvBuf == NULL)
			throw;

		MPI_Alltoallv(sendBuf, sendBytes, sendDisp, MPI_BYTE, recvBuf, recvBytes, recvDisp, MPI_BYTE, world);
		free(sendBuf);
		contigReadSet.resize(contigs.getSize(), ReadSet());

		for(int i = 0; i < world.size(); i++) {
			ReadSet::ReadSetSizeType contigIdx = 0;
			tmp = recvBuf + recvDisp[i];
			while (tmp != recvBuf + recvDisp[i] + recvBytes[i]) {
				ReadSet rs;
				tmp = (char*) rs.restore(tmp);
				contigReadSet[contigIdx++].append(rs);
			}
		}
		free(recvBuf);

		ReadSet changedContigs;

		for(ReadSet::ReadSetSizeType i = 0; i < contigs.getSize(); i++) {
			const Read &oldRead = contigs.getRead(i);
			ReadSet myContig;
			myContig.append(oldRead);
			ReadSet newContig = ContigExtender<KS>::extendContigs(myContig, contigReadSet[i]);
			const Read &newRead = newContig.getRead(0);
			if (newRead.getLength() > oldRead.getLength())
				changedContigs.append(newRead);
			else
				finalContigs.append(oldRead);
		}

		LOG_VERBOSE(1, "Changed contigs: " << changedContigs.getSize());

		std::string filekey = "contig-" + boost::lexical_cast<std::string>(iteration);
		DistributedOfstreamMap om(world, tmpDir, FormatOutput::getSuffix(FormatOutput::FASTA));
		om.setBuildInMemory();
		changedContigs.writeAll(om.getOfstream(filekey), FormatOutput::FASTA);
		std::string newContigFile = om.getFilePath(filekey);
		om.clear();

		contigFile = newContigFile;
	}

	DistributedOfstreamMap om(world, Options::getOutputFile(), "");
	om.setBuildInMemory();
	finalContigs.writeAll(om.getOfstream(""), FormatOutput::FASTA);

}

