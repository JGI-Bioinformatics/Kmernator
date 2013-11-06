/*
 * Vmatch.h
 *
 *  Created on: Sep 7, 2011
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

#ifndef VMATCH_H_
#define VMATCH_H_
#include <stdlib.h>
#include "Options.h"
#include "Utils.h"
#include "ReadSet.h"
#include "DistributedFunctions.h"
#include "MatcherInterface.h"

class _VmatchOptions : public OptionsBaseInterface {
public:
	static std::string getVmatchPath() {
		return getVarMap()["vmatch-path"].as<std::string>();
	}
	static std::string getVmatchOptions() {
		return getVarMap()["vmatch-options"].as<std::string> ();
	}
	static std::string getVmatchIndexPath() {
		return getVarMap()["vmatch-index-path"].as<std::string> ();
	}
	static bool getVmatchPreload() {
		return getVarMap()["vmatch-preload"].as<bool> ();
	}
	void _resetDefaults() {
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {

		po::options_description opts("Vmatch Options");
		opts.add_options()

				("vmatch-path", po::value<std::string>()->default_value(""), "if specified the path to the directory containing vmatch and mkvtree")

				("vmatch-options", po::value<std::string>()->default_value( "-d -p -seedlength 10 -l 50 -e 3"), "options with which to call vmatch")

				("vmatch-index-path", po::value<std::string>()->default_value("."), "top level directory under which to create the vmatch index directories for each rank")

				("vmatch-preload", po::value<bool>()->default_value(false), "pre-load the necessary vmatch files in a mmap")

				;
		desc.add(opts);

	}
};
typedef OptionsBaseTemplate< _VmatchOptions > VmatchOptions;

class Vmatch : public MatcherInterface {
public:
	typedef MatcherInterface::MatchResults MatchResults;

	class FieldsType {
	public:

		FieldsType(std::string &_line) :
			line(_line) {
			static const char *format = "%d %d %d %c %d %d %d %d %e %d %f";
			sscanf(line.c_str(), format, &subjectLength, &subjectNumber,
					&subjectPosition, &type, &queryLength, &queryNumber,
					&queryPosition, &distance, &eValue, &scoreValue,
					&percentIdentity);
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
	typedef std::vector<FieldsType> VmatchMatchResults;

private:
	std::string _indexName;
	VmatchMatchResults _results;
	Kmernator::MmapSourceVector _mmaps;
	std::string _binaryPath;

public:
	Vmatch(	mpi::communicator &world,
			std::string indexName,
			const ReadSet &target,
			bool cacheIndexes =
			VmatchOptions::getOptions().getVmatchPreload()
		) :  MatcherInterface(world, target), _indexName(indexName) {

		std::string	tmpDir = VmatchOptions::getOptions().getVmatchIndexPath() + "/";
		mkdir(tmpDir.c_str(), 0777); // all ranks need to mkdir if writing to local disks...
		tmpDir += _indexName + "/";
		mkdir(tmpDir.c_str(), 0777); // all ranks need to mkdir if writing to local disks...
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Created temporary directory for ranked vmatch files: " << tmpDir);

		std::string rankOutputDir = DistributedDirectoryManagement::getRankSubDir(world, tmpDir);
		_indexName = rankOutputDir + "/myReads";

		_binaryPath = VmatchOptions::getOptions().getVmatchPath();
		if (!_binaryPath.empty())
			_binaryPath += "/";

		if (FileUtils::getFileSize(_indexName + ".suf") > 0) {
			LOG_DEBUG(1, "Vmatch(" << _indexName << ", reads(" << target.getSize() << ")): Index already exists: " << _indexName);
		} else {
			buildVmatchIndex(target);
		}
		if (cacheIndexes)
			mapIndexes();
	}
	virtual ~Vmatch() {
		clearMaps();
	}
	void buildVmatchIndex(const ReadSet &target) {
		std::string inputFile = _indexName + ".tmp.input";
		OfstreamMap ofm(inputFile, "");
		target.writeAll(ofm.getOfstream(""), FormatOutput::FastaUnmasked());
		ofm.clear();
		buildVmatchIndex(inputFile);
		LOG_DEBUG(1, "Removing temporary inputFile" << inputFile);
		unlink(inputFile.c_str());
	}
	void buildVmatchIndex(std::string inputFasta) {
		std::string logFile = _indexName + "-mkvtree.log";
		std::string cmd(_binaryPath + "mkvtree -dna -allout -pl -indexname " + _indexName
				+ " -db " + inputFasta + " >" + logFile + " 2>&1");
		LOG_DEBUG(1, "Building vmatch index " << _indexName << " : " << cmd);
		int ret = ForkDaemon::system(cmd);
		if (ret != 0)
			LOG_THROW("mkvtree failed to build(" << ret << "): " << cmd << " Log:\n" << FileUtils::dumpFile(logFile));
	}

	MatchResults matchLocalImpl(std::string queryFile) {
		VmatchMatchResults matches = _match(queryFile);
		MatchResults contigReadHits;
		bool returnPairs = MatcherInterfaceOptions::getOptions().getIncludeMate();

		int myRank = getWorld().rank();
		for (VmatchMatchResults::iterator match = matches.begin(); match
		!= matches.end(); match++) {
			ReadSet::ReadSetSizeType globalContigIdx = match->queryNumber;
			ReadSet::ReadSetSizeType globalReadIdx = this->getTarget().getGlobalReadIdx(myRank, match->subjectNumber);

			if (globalContigIdx >= contigReadHits.size())
				contigReadHits.resize(globalContigIdx+1);

			contigReadHits[globalContigIdx].insert(globalReadIdx);
			// include pairs
			if (returnPairs) {
				if (globalReadIdx % 2 == 0)
					contigReadHits[globalContigIdx].insert(globalReadIdx + 1);
				else
					contigReadHits[globalContigIdx].insert(globalReadIdx - 1);
			}
		}
		return contigReadHits;
	}

	VmatchMatchResults &_match(std::string queryFile, std::string options = VmatchOptions::getOptions().getVmatchOptions()) {
		double time = MPI_Wtime();
		_results.clear();
		std::string cmd = _binaryPath + "vmatch " + options + " -q " + queryFile + " "
				+ _indexName;
		if (Log::isDebug(2))
			cmd = "strace -tt -T " + cmd;

		LOG_DEBUG_OPTIONAL(1, true, "Executing vmatch: " << cmd);
		IPipestream vmatchOutput(cmd, true);
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
		if (vmatchOutput.getExitStatus() != 0)
			LOG_THROW(cmd << " failed.");

		LOG_VERBOSE_OPTIONAL(1, true, "Vmatch::match(,): Found " << _results.size() << " results in " << (MPI_Wtime() - time) << " sec");
		if (Log::isDebug(2)) {
			LOG_DEBUG(1, "Vmatch::match() strace:\n" << vmatchOutput.getStdErr());
		}
		return _results;
	}
	void mapIndexes(bool _flush = true) {
		LOG_DEBUG(2, "memory mapping vmatch index files for " << _indexName);
		long size = 0;
		const std::string suffixes[] = { "prj", "al1", "tis", "ssp", "suf",
				"lcp", "bck", "sti1" };
		double time = MPI_Wtime(), time2, firstTime = 0;
		for (int i = 0; i < 8; i++) {
			std::string filePath = _indexName + "." + suffixes[i];
			LOG_DEBUG(3, "mmaping " << filePath);
			Kmernator::MmapSource mmap(filePath, FileUtils::getFileSize(
					filePath));
			madvise(const_cast<char*> (mmap.data()), mmap.size(), MADV_WILLNEED);
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
		const long *p = (const long*) data;
		size /= sizeof(long);
		long c = 0;
		for (long j = 0; j < size; j++)
			c += *(p++);
		return c;
	}
};


#endif /* VMATCH_H_ */
