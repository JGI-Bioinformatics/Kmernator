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
#include "MemoryUtils.h"
#include "Utils.h"
#include "Log.h"

using namespace std;
typedef TrackingDataMinimal4f DataType;
typedef KmerSpectrum<DataType, DataType> KS;


class _DistributedNucleatingAssemblerOptions : public _ContigExtenderOptions  {
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
	static bool getVmatchPreload() {
		return getVarMap()["vmatch-preload"].as<bool>();
	}
	bool _parseOpts(po::options_description &desc, po::positional_options_description &p, po::variables_map &vm, int argc, char *argv[]) {

		desc.add_options()

		("max-iterations", po::value<int>()->default_value(1000),
				"the maximum number of rounds to extend the set of contigs")

		("vmatch-options", po::value<std::string>()->default_value("-d -p -seedlength 10 -l 50 -e 3"),
				"options with which to call vmatch")

		("vmatch-index-path", po::value<std::string>()->default_value("."),
				"top level directory under which to create the vmatch index directories for each rank")

		("vmatch-preload", po::value<bool>()->default_value(false),
				"pre-load the necessary vmatch files in a mmap")

				;

		bool ret = ContigExtenderOptions::parseOpts(argc, argv);

		if (Options::getOptions().getOutputFile().empty()) {
			LOG_WARN(1, "You must specify an --output");
			ret = -1;
		}

		return ret;
	}
};
typedef OptionsBaseTemplate< _DistributedNucleatingAssemblerOptions > DistributedNucleatingAssemblerOptions;

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
	Kmernator::MmapSourceVector _mmaps;

public:
	Vmatch(std::string indexName, ReadSet &inputs, bool cacheIndexes = DistributedNucleatingAssemblerOptions::getOptions().getVmatchPreload()) : _indexName(indexName) {
		if (FileUtils::getFileSize(indexName + ".suf") > 0) {
			LOG_DEBUG(1, "Vmatch(" << indexName << ", reads(" << inputs.getSize() << ")): Index already exists: " << indexName);
		} else {
			buildVmatchIndex(indexName, inputs);
		}
		if (cacheIndexes)
			mapIndexes();
	}
	virtual ~Vmatch() {
		clearMaps();
	}
	static void buildVmatchIndex(std::string indexName, ReadSet &inputs) {
		std::string inputFile = indexName + ".tmp.input";
		OfstreamMap ofm(inputFile, "");
		inputs.writeAll(ofm.getOfstream(""), FormatOutput::FASTA);
		ofm.clear();
		buildVmatchIndex(indexName, inputFile);
		LOG_DEBUG(1, "Removing temporary inputFile" << inputFile);
		unlink(inputFile.c_str());
	}
	static void buildVmatchIndex(std::string indexName, std::string inputFasta) {
		std::string cmd("mkvtree -dna -allout -pl -indexname " + indexName + " -db " + inputFasta);
		LOG_DEBUG(1, "Building vmatch index " << indexName << " : " << cmd);
		int ret = system(cmd.c_str());
		if (ret != 0)
			LOG_THROW("mkvtree failed to build(" << ret << "): " << cmd);
	}
	MatchResults &match(std::string queryFile, std::string options = "") {
		double time = MPI_Wtime();
		_results.clear();
		std::string cmd = "vmatch " + options + " -q " + queryFile + " " + _indexName;
		if (Log::isDebug(2))
			cmd = "strace -tt -T " + cmd;

		LOG_DEBUG_OPTIONAL(1, true, "Executing vmatch: " << cmd);
		IPipestream vmatchOutput(cmd, Log::isDebug(2));
		std::string line;
		while (!vmatchOutput.eof()) {
			getline(vmatchOutput, line);
			if (line.length() == 0)
				continue;
			if (line[0] == '#')
				continue;
			_results.push_back(FieldsType(line));
		}
		vmatchOutput.close();
		LOG_VERBOSE_OPTIONAL(1, true, "Vmatch::match(,): Found " << _results.size() << " results in " << (MPI_Wtime() - time) << " sec");
		if (Log::isDebug(2)) {
			LOG_DEBUG(1, "Vmatch::match() strace:\n" << vmatchOutput.getStdErr());
		}
		return _results;
	}
	void mapIndexes(bool _flush = true) {
		LOG_DEBUG(2, "memory mapping vmatch index files for " << _indexName);
		long size = 0;
		const std::string suffixes[] = {"prj", "al1", "tis", "ssp", "suf", "lcp", "bck", "sti1"};
		double time = MPI_Wtime(), time2, firstTime = 0;
		for(int i = 0 ; i < 8 ; i++) {
			std::string filePath = _indexName + "." + suffixes[i];
			LOG_DEBUG(3, "mmaping " << filePath);
			Kmernator::MmapSource mmap(filePath, FileUtils::getFileSize(filePath));
			madvise(const_cast<char*>(mmap.data()), mmap.size(), MADV_WILLNEED);
			size += mmap.size();

			if (_flush) {
				time2 = MPI_Wtime();
				flush(mmap.data(), mmap.size());
				firstTime += MPI_Wtime() - time2;
			}

			_mmaps.push_back(mmap);
		}
		LOG_DEBUG(2, MemoryUtils::getMmapUsage());
		LOG_VERBOSE(1, "Vmatch(): memory mapped " << size << " bytes for " << _indexName << " in " << (MPI_Wtime()-time) << " sec.  touched in " << firstTime << " sec");
	}
	void clearMaps() {
		_mmaps.clear();
	}
	static long flush(const char *data, long size) {
		const long *p =(const long*) data;
		size /= sizeof(long);
		long c = 0;
		for (long j = 0 ; j < size; j++)
			c += *(p++);
		return c;
	}
};


int main(int argc, char *argv[]) {
	// do not apply artifact filtering by default
	Options::getOptions().getSkipArtifactFilter() = 1;
	// override the default output format!
	Options::getOptions().getFormatOutput() = 3;
	Options::getOptions().getMmapInput() = 0;
	Options::getOptions().getVerbose() = 1;
	Options::getOptions().getMaxThreads() = 1;
	double timing1, timing2;

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

		if (Options::getOptions().getGatheredLogs())
			Logger::setWorld(&world, Options::getOptions().getDebug() >= 3);

	} catch (...) {
		std::cerr << OptionsBaseInterface::getDesc() << std::endl;
		std::cerr << std::endl << "Please fix the options and/or MPI environment" << std::endl;
		exit(1);
	}
	timing1 = MPI_Wtime();

	OptionsBaseInterface::FileListType inputFiles = Options::getOptions().getInputFiles();
	std::string contigFile = ContigExtenderOptions::getOptions().getContigFile();
	std::string finalContigFile;
	double minimumCoverage = ContigExtenderOptions::getOptions().getMinimumCoverage();
	long maxIterations = DistributedNucleatingAssemblerOptions::getOptions().getMaxIterations();
	std::string tmpDir = DistributedNucleatingAssemblerOptions::getOptions().getVmatchIndexPath() + "/";
	if (world.rank() == 0)
		mkdir(tmpDir.c_str(), 0777);

	world.barrier();

	ReadSet reads;
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files" );
	reads.appendAllFiles(inputFiles, world.rank(), world.size());

	setGlobalReadSetOffsets(world, reads);
	timing2 = MPI_Wtime();

	LOG_VERBOSE_OPTIONAL(1, world.rank() == 1, "loaded " << reads.getGlobalSize() << " Reads in " << (timing2-timing1) << " seconds" );

	std::string rankOutputDir = getRankSubdir(world, tmpDir);
	Vmatch myVmatch = Vmatch(rankOutputDir + "/myReads", reads);

	SequenceLengthType minKmerSize, maxKmerSize, kmerStep, maxExtend;
	ContigExtender<KS>::getMinMaxKmerSize(reads, minKmerSize, maxKmerSize, kmerStep);
	maxKmerSize = boost::mpi::all_reduce(world, maxKmerSize, mpi::minimum<SequenceLengthType>());
	LOG_VERBOSE(1, "Kmer size ranges: " << minKmerSize << "\t" << maxKmerSize << "\t" << kmerStep);
	maxExtend = maxKmerSize;

	timing1 = timing2; timing2 = MPI_Wtime();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Prepared vmatch indexes in " << (timing2-timing1) << " seconds");

	ReadSet finalContigs;
	ReadSet contigs;
	contigs.appendFastaFile(contigFile, world.rank(), world.size());

	short iteration = 0;
	while(++iteration <= maxIterations) {
		double startIterationTime = MPI_Wtime();

		setGlobalReadSetOffsets(world, contigs);

		double setOffsetsTime = MPI_Wtime() - startIterationTime;

		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Iteration: " << iteration << ". Contig File: " << contigFile << ". contains " << contigs.getGlobalSize() << " Reads");
		if (contigs.getGlobalSize() == 0) {
			LOG_VERBOSE_OPTIONAL(1, true, "There are no contigs to extend in " << contigFile);
			break;
		}

		ReadSet::ReadSetVector contigReadSet(contigs.getGlobalSize(), ReadSet());
		Vmatch::MatchResults matches = myVmatch.match(contigFile, DistributedNucleatingAssemblerOptions::getOptions().getVmatchOptions());

		double vmatchMatchTime = MPI_Wtime() - startIterationTime;

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

		double prepareSizesTime = MPI_Wtime() - startIterationTime;

		MPI_Alltoall(sendBytes, 1, MPI_INT, recvBytes, 1, MPI_INT, world);

		double exchangeSizesTime = MPI_Wtime() - startIterationTime;

		for(int rank = 0; rank < world.size(); rank++) {
			recvDisp[rank] = totalRecv;
			totalRecv += recvBytes[rank];
			LOG_DEBUG_OPTIONAL(3, true, "to/from rank " << rank << ": sendDisp[] " << sendDisp[rank] << " sendBytes[] " << sendBytes[rank] << " recvDisp[] " << recvDisp[rank] << " recvBytes[] " << recvBytes[rank] );
		}
		LOG_DEBUG_OPTIONAL(3, true, "Sending " << totalSend << " Receiving " << totalRecv);

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

		double prepareSendTime = MPI_Wtime() - startIterationTime;

		MPI_Alltoallv(sendBuf, sendBytes, sendDisp, MPI_BYTE, recvBuf, recvBytes, recvDisp, MPI_BYTE, world);

		double exchangeReadsTime = MPI_Wtime() - startIterationTime;

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

		std::stringstream extendLog;
        //#pragma omp parallel for
		for(ReadSet::ReadSetSizeType i = 0; i < contigs.getSize(); i++) {
			const Read &oldRead = contigs.getRead(i);
			Read newRead;
			SequenceLengthType oldLen = oldRead.getLength(), newLen = 0;
			ReadSet::ReadSetSizeType poolSize = contigReadSet[i].getSize();
			SequenceLengthType myKmerSize = minKmerSize;
			if (poolSize > minimumCoverage) {
				LOG_VERBOSE_OPTIONAL(2, true, "Extending " << oldRead.getName() << " with " << poolSize << " pool of reads");
				ReadSet myContig;
				myContig.append(oldRead);
				ReadSet newContig;

				while (newLen <= oldLen && myKmerSize <= maxKmerSize) {
					newContig = ContigExtender<KS>::extendContigs(myContig, contigReadSet[i], maxExtend, myKmerSize, myKmerSize);
					newLen = newContig.getRead(0).getLength();
					myKmerSize += kmerStep;
				}
				newRead = newContig.getRead(0);
			} else {
				newRead = oldRead;
			}
			long deltaLen = newLen - oldLen;
			if (deltaLen > 0) {
				extendLog << std::endl << "Extended " << oldRead.getName() << " " << deltaLen << " bases to " << newRead.getLength() << ": " << newRead.getName() << " with " << poolSize << " reads in the pool K " << (myKmerSize-kmerStep);
				//#pragma omp critical
				changedContigs.append(newRead);
			} else {
				extendLog << std::endl << "Did not extend " << oldRead.getName() << " with " << poolSize << " reads in the pool";
				//#pragma omp critical
				finalContigs.append(oldRead);

			}

		}

		double extendContigsTime = MPI_Wtime() - startIterationTime;
		LOG_VERBOSE(1, (extendLog.str()));

		LOG_DEBUG(1, "Changed contigs: " << changedContigs.getSize() << " finalContigs: " << finalContigs.getSize());
		setGlobalReadSetOffsets(world, changedContigs);
		setGlobalReadSetOffsets(world, finalContigs);
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Changed contigs: " << changedContigs.getGlobalSize() << " finalContigs: " << finalContigs.getGlobalSize());

		double calcChangedSizesTime = MPI_Wtime() - startIterationTime;

		std::string oldFinalContigFile = finalContigFile;
		std::string oldContigFile = contigFile;
		{
			// write out the state of the contig files (so far) so we do not loose them
			DistributedOfstreamMap om(world, Options::getOptions().getOutputFile(), "");
			om.setBuildInMemory();
			if (finalContigs.getGlobalSize() > 0) {
				std::string fileKey = "final-" + boost::lexical_cast<std::string>(iteration);
				finalContigs.writeAll(om.getOfstream(fileKey), FormatOutput::FASTA);
				finalContigFile = om.getRealFilePath(fileKey);
			}
			if (changedContigs.getGlobalSize() > 0) {
				std::string filekey = "-inputcontigs-" + boost::lexical_cast<std::string>(iteration) + ".fasta";
				changedContigs.writeAll(om.getOfstream(filekey), FormatOutput::FASTA);
				contigFile = om.getRealFilePath(filekey);
			}
			contigs = changedContigs;
		}

		double writeFinalTime = MPI_Wtime() - startIterationTime;

		if (!Log::isDebug(1) && world.rank() == 0) {
			// remove most recent contig files (if not debugging)
			if (!oldFinalContigFile.empty()) {
				LOG_VERBOSE_OPTIONAL(1, true, "Removing " << oldFinalContigFile);
				unlink(oldFinalContigFile.c_str());
			}

			if (ContigExtenderOptions::getOptions().getContigFile().compare(oldContigFile) != 0) {
				LOG_VERBOSE_OPTIONAL(1, true, "Removing " << oldContigFile);
				unlink(oldContigFile.c_str());
			}
		}

		if (changedContigs.getGlobalSize() == 0) {
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 1, "No more contigs to extend " << changedContigs.getSize());
			break;
		}

		double finishIterationTime = MPI_Wtime() - startIterationTime;
		if (Log::isVerbose(1)) {
			char buf[1024];
			sprintf(buf, "Time SO:%0.2f vm:%0.2f ps:%0.2f EST:%0.2f ps:%0.2f ER:%0.2f ec:%0.2f CCS:%0.2f WF:%0.2f FI:%0.2f total: %0.2f Mem: ",
					setOffsetsTime,
					vmatchMatchTime - setOffsetsTime,
					prepareSizesTime - vmatchMatchTime,
					exchangeSizesTime - prepareSizesTime,
					prepareSendTime - exchangeSizesTime,
					exchangeReadsTime - prepareSendTime,
					extendContigsTime - exchangeReadsTime,
					calcChangedSizesTime - extendContigsTime,
					writeFinalTime - calcChangedSizesTime,
					finishIterationTime - writeFinalTime,
					finishIterationTime);
			Log::Verbose(std::string(buf) + MemoryUtils::getMemoryUsage(), false);
		}
	}

	if (world.rank() == 0 && !Log::isDebug(1)) {
		if (ContigExtenderOptions::getOptions().getContigFile().compare(contigFile) != 0) {
			LOG_DEBUG_OPTIONAL(1, true, "Removing " << contigFile);
			unlink(contigFile.c_str());
		}
	}

	{
		// write final contigs (and any unfinished contigs still remaining)
		std::string tmpFinalFile = finalContigFile;
		DistributedOfstreamMap om(world, Options::getOptions().getOutputFile(), "");
		om.setBuildInMemory();
		finalContigs.append(contigs);
		std::string fileKey = "";
		finalContigs.writeAll(om.getOfstream(fileKey), FormatOutput::FASTA);
		om.clear();

		if (world.rank() == 0 && !finalContigFile.empty()) {
			LOG_DEBUG_OPTIONAL(1, true, "Removing " << finalContigFile);
			unlink(finalContigFile.c_str());
		}
		finalContigFile = om.getRealFilePath(fileKey);
	}
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Final contigs are in: " << finalContigFile);

	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");
	MPI_Finalize();

	return 0;
}

