/*
 * Vmatch.h
 *
 *  Created on: Sep 7, 2011
 *      Author: regan
 */

#ifndef VMATCH_H_
#define VMATCH_H_
#include <stdlib.h>
#include "Options.h"
#include "Utils.h"
#include "ReadSet.h"

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

		("vmatch-path", po::value<std::string>()->default_value(""),
				"if specified the path to the directory containing vmatch and mkvtree")

		("vmatch-options", po::value<std::string>()->default_value(
				"-d -p -seedlength 10 -l 50 -e 3"),
				"options with which to call vmatch")

		("vmatch-index-path",
				po::value<std::string>()->default_value("."),
				"top level directory under which to create the vmatch index directories for each rank")

		("vmatch-preload", po::value<bool>()->default_value(false),
				"pre-load the necessary vmatch files in a mmap")

		;
		desc.add(opts);

	}
};
typedef OptionsBaseTemplate< _VmatchOptions > VmatchOptions;

class Vmatch {
public:
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
	typedef std::vector<FieldsType> MatchResults;

private:
	std::string _indexName;
	MatchResults _results;
	Kmernator::MmapSourceVector _mmaps;
	std::string _binaryPath;

public:
	Vmatch(
			std::string indexName,
			ReadSet &inputs,
			bool cacheIndexes =
					VmatchOptions::getOptions().getVmatchPreload()) :
		_indexName(indexName) {
		_binaryPath = VmatchOptions::getOptions().getVmatchPath();
		if (!_binaryPath.empty())
			_binaryPath += "/";

		if (FileUtils::getFileSize(indexName + ".suf") > 0) {
			LOG_DEBUG(1, "Vmatch(" << indexName << ", reads(" << inputs.getSize() << ")): Index already exists: " << indexName);
		} else {
			buildVmatchIndex(inputs);
		}
		if (cacheIndexes)
			mapIndexes();
	}
	virtual ~Vmatch() {
		clearMaps();
	}
	void buildVmatchIndex(ReadSet &inputs) {
		std::string inputFile = _indexName + ".tmp.input";
		OfstreamMap ofm(inputFile, "");
		inputs.writeAll(ofm.getOfstream(""), FormatOutput::FastaUnmasked());
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
		int ret = system(cmd.c_str());
		if (ret != 0)
			LOG_THROW("mkvtree failed to build(" << ret << "): " << cmd << " Log:\n" << FileUtils::dumpFile(logFile));
	}
	MatchResults &match(std::string queryFile, std::string options = "") {
		double time = MPI_Wtime();
		_results.clear();
		std::string logFile = _indexName + "-vmatch.log";
		std::string cmd = _binaryPath + "vmatch " + options + " -q " + queryFile + " "
				+ _indexName + " 2>" + logFile;
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
