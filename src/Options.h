//
// Kmernator/src/Options.h
//
// Author: Rob Egan, Craig Furman
//
// Copyright 2010 The Regents of the University of California.
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
// binaries display the following acknowledgement:  "This product
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

#ifndef _OPTIONS_H
#define _OPTIONS_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>


#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "config.h"
#include "Log.h"

#define PASSES_LENGTH(length, readLength, minimumLength) ((minimumLength == MAX_SEQUENCE_LENGTH) ? (length == readLength) : (length >= minimumLength))

class GlobalOptions {
public:
	static inline bool &isCommentStored() {
		static bool _isCommentStored = false;
		return _isCommentStored;
	}
};

// The base class that returns the single, static set of options
class OptionsInstance {
public:
	static po::options_description &getDesc() {
		static po::options_description _desc;
		return _desc;
	}
	static po::positional_options_description &getPosDesc() {
		static po::positional_options_description _p;
		return _p;
	}
	static po::variables_map &getVarMap() {
		static po::variables_map _vm;
		return _vm;
	}
	static std::string &getOptionsErrorMsg() {
		static std::string _message;
		return _message;
	}
	static std::string setOptionsErrorMsg(std::string msg) {
		LOG_ERROR(1, "Invalid option: " << msg);
		getOptionsErrorMsg() += msg + "\n";
		return getOptionsErrorMsg();
	}
private:
	OptionsInstance();
	~OptionsInstance();
};

class OptionsBaseInterface {
public:
	typedef std::vector<std::string> StringListType;
	typedef boost::shared_ptr< std::ofstream > OStreamPtr;
	typedef StringListType FileListType;
	static std::string toString(const FileListType &fl) {
		std::stringstream ss;
		bool isNotFirst = false;
		for (FileListType::const_iterator it = fl.begin(); it != fl.end(); it++) {
			if (isNotFirst)
				ss << ", ";

			isNotFirst = true;
			ss << *it;
		}
		return ss.str();
	}

	static po::options_description &getDesc() {
		return OptionsInstance::getDesc();
	}
	static po::positional_options_description &getPosDesc() {
		return OptionsInstance::getPosDesc();
	}
	static po::variables_map &getVarMap() {
		return OptionsInstance::getVarMap();
	}
	static unsigned int &getVerbose() {
		return Log::Verbose().setLevel();
	}
	static unsigned int &getDebug() {
		return Log::Debug().setLevel();
	}

	// helper method to set external (or instance) values and log/print the option
	// expected to be called in _parseOptions(&vm)
	template<typename U>
	static void setOpt(std::string key, U &val, bool print = Log::printOptions()) {
		po::variables_map &vm = getVarMap();
		if (vm.count(key.c_str())) {
			val = vm[key.c_str()].as<U>();
			if (print) {
				LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), key << " is: " << val );
			}
		} else if (print) {
			LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), key << " was not specified.");
		}
	};
	void setOpt2(std::string key, StringListType &val, bool print = Log::printOptions()) {
		po::variables_map &vm = getVarMap();
		if (vm.count(key.c_str())) {
			val = vm[key.c_str()].as<StringListType>();
			if (print) {
				LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), key << " is: " << toString(val) );
			}
		} else if (print) {
			LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), key << " was not specified.");
		}
	}

	static void setOptionsErrorMsg(std::string msg) {
		OptionsInstance::setOptionsErrorMsg(msg);
	}
	static std::string getOptionsErrorMsg() {
		return OptionsInstance::getOptionsErrorMsg();
	}
	static bool hasOptionsErrorMsg() {
		return ! OptionsInstance::getOptionsErrorMsg().empty();
	}

	// TODO this could be moved somewhere else
	static std::string getHostname() {
		char hostname[256];
		gethostname(hostname, 256);
		return std::string(hostname);
	}

	// use to set/overrided any defaults on options that are stored persistently
	virtual void _resetDefaults() {}
	// use to set the description of all options
	virtual void _setOptions(po::options_description &desc, po::positional_options_description &p) {}
	// use to post-process options, returning true if everything is okay
	virtual bool _parseOptions(po::variables_map &vm) {
		return true;
	}
};

// Singleton template class
template <class T>
class OptionsBaseTemplate {
public:
	static T& getOptions() {
		static T _instance;
		return _instance;
	}
	static po::options_description &getDesc() {
		return OptionsInstance::getDesc();
	}
	static po::positional_options_description &getPosDesc() {
		return OptionsInstance::getPosDesc();
	}
	static po::variables_map &getVarMap() {
		return OptionsInstance::getVarMap();
	}

	// helper methods to pass through static methods to base class
	static void _resetDefaults() {
		getOptions()._resetDefaults();
	}
	static void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		getOptions()._setOptions(desc, p);
	}
	static bool _parseOptions(po::variables_map &vm) {
		return getOptions()._parseOptions(vm);
	}
	static bool hasOptionsErrorMsg() {
		return getOptions().hasOptionsErrorMsg();
	}
	static void setOptionsErrorMsg(std::string msg) {
		getOptions().setOptionsErrorMsg(msg);
	}
	static std::string getOptionsErrorMsg() {
		return getOptions().getOptionsErrorMsg();
	}
	static void parseOpts(int argc, char *argv[]) {
		// set any defaults
		bool ret = true;
		try {
			_resetDefaults();

			// load all the descriptions
			_setOptions(getDesc(), getPosDesc());

			// process the command line
			po::store(po::command_line_parser(argc, argv).options(getDesc()).positional(getPosDesc()).run(), getVarMap());

			// update the varmap
			po::notify( getVarMap() );

			// post-process

			ret = _parseOptions( getVarMap() );

			if (hasOptionsErrorMsg())
				LOG_THROW(getOptionsErrorMsg());
			if (!Log::getErrorMessages().empty())
				LOG_THROW("");
			if (!ret)
				LOG_THROW("one or more options are invalid");

		} catch (std::exception &e) {
			std::cerr << getDesc() << std::endl << std::endl;
			std::cerr << "Please fix the command line arguments:" << std::endl << std::endl << e.what();
			std::cerr << std::endl << std::endl;
			exit(1);
		} catch (...) {
			std::cerr << getDesc() << std::endl << std::endl;
			std::cerr << "Please fix the command line arguments: something is wrong... " << std::endl << std::endl;
			std::cerr << std::endl << std::endl;
			exit(1);
		}

	}


protected:

	OptionsBaseTemplate() {} // hidden constructor
	virtual ~OptionsBaseTemplate() {}  // hidden destructor
	OptionsBaseTemplate(OptionsBaseTemplate const&); // disabled
	OptionsBaseTemplate& operator=(OptionsBaseTemplate const&); // disabled
};

/* Example copy-paste template */

class _MySpecificOptions : public OptionsBaseInterface {
public:
	_MySpecificOptions() : myOption(1234) {
	}
	virtual ~_MySpecificOptions() {}
	void _resetDefaults() {
		// Other *::_resetDefaults();
	}
	int &getMyOption() {
		return myOption;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("My Specific Options");
		opts.add_options()("my-option", po::value<int>()->default_value(myOption), "my option");
		desc.add(opts);
		// Other *::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		setOpt("my-option", myOption);
		// Other ret &= *::_parseOptions(vm);
		return ret;
	}
private:
	int myOption;
};
typedef OptionsBaseTemplate< _MySpecificOptions > MySpecificOptions;


class _GeneralOptions : public OptionsBaseInterface {
public:
	_GeneralOptions() : maxThreads(OMP_MAX_THREADS_DEFAULT), tmpDir("/tmp"),
	formatOutput(0), keepReadComment(GlobalOptions::isCommentStored()), buildOutputInMemory(false),
	minQuality(3),  ignoreQual(false), mmapInput(true), gatheredLogs(true),
	batchSize(100000), separateOutputs(true)
	{
		char *tmpPath;
		tmpPath = getenv ("TMPDIR");
		if (tmpPath != NULL) {
			tmpDir = std::string(tmpPath);
		}
	}

	// make this final, so preserving the singleton state
private:
	~_GeneralOptions() {}
	friend class OptionsBaseTemplate< _GeneralOptions >;

private:
	int          maxThreads;
	FileListType inputFiles;
	FileListType inputFilePrefixes;
	std::string  outputFile;
	std::string  logFile;
	OStreamPtr   logFileStream;
	std::string  tmpDir;
	unsigned int formatOutput;
	bool         keepReadComment;
	bool         buildOutputInMemory;
	unsigned int minQuality;
	bool ignoreQual;
	bool mmapInput;
	bool gatheredLogs;
	unsigned int batchSize;
	bool separateOutputs;

public:
	void _resetOptions() {
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {

		po::options_description general("General Options");

		po::options_description logOpts("Logging and Messages");
		logOpts.add_options()
				("help", "produce help message")

		// logging & messages
				("verbose", po::value<unsigned int>()->default_value(getVerbose()), "level of verbosity (0+)")

				("debug", po::value<unsigned int>()->default_value(getDebug()), "level of debug verbosity (0+)")

				("log-file", po::value<std::string>()->default_value(logFile), "If set all INFO and DEBUG messages will be logged here (default stderr)")

				("gathered-logs", po::value<bool>()->default_value(gatheredLogs), "If set and MPI is enabled, VERBOSE1, VERBOSE2 and DEBUG1 logs will be gathered to the master before being output.")

				;

		general.add(logOpts);

		po::options_description miscOpts("Misc General");
		miscOpts.add_options()
#ifdef _USE_OPENMP
				("threads", po::value<int>()->default_value(maxThreads), "maximum number of threads")
#endif
				// Misc General Options
				("batch-size", po::value<unsigned int>()->default_value(batchSize), "default size of batches (reads, kmers, MPI, etc)")
				;
		general.add(miscOpts);

		po::options_description readOpts("Reading and Input Files");
		readOpts.add_options()
		// reading / input files
				("input-file", po::value<FileListType>(), "input file(s)")

				("mmap-input", po::value<bool>()->default_value(mmapInput), "If false, prevents input files from being mmaped, instead import reads into memory (somewhat faster if memory is abundant)")

				("ignore-quality", po::value<bool>()->default_value(ignoreQual), "ignore the quality score, to save memory or if they are untrusted")

				("min-quality-score", po::value<unsigned int>()->default_value(minQuality), "minimum quality score below which will be evaluated as Q=0 (i.e. 'N', prob = 0.0)")
				;
		general.add(readOpts);

		po::options_description writeOpts("Writing and Output Files");
		writeOpts.add_options()
		// writing / output files
				("output-file", po::value<std::string>(), "output file (sometimes a pattern prefix)")

				("format-output", po::value<unsigned int>()->default_value(formatOutput), "0: fastq, 1: fasta, 2: fastq unmasked, 3: fasta unmasked")

				("keep-read-comment", po::value<bool>()->default_value(keepReadComment), "If set, per-read comment will be preserved")

				("separate-outputs", po::value<bool>()->default_value(separateOutputs), "If set, each input (plus consensus) will generate a new outputfile.  If set false, all input files will be merged into one output file.")

				("build-output-in-memory", po::value<bool>()->default_value(buildOutputInMemory), "if set, all temporary output files will first be stored in memory (faster for MPI applications)")

				("temp-dir", po::value<std::string>()->default_value(tmpDir), "temporary directory to utilize")
				;
		general.add(writeOpts);

		desc.add(general);
	}
	void _setVerbosityAndLogs(bool print) {
		setOpt("verbose", getVerbose(), print);
		setOpt("debug", getDebug(), print);
		setOpt("log-file", getLogFile(), print);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		try {

			if (vm.count("help")) {
				std::cerr << getDesc() << std::endl;
				return false;
			}

			_setVerbosityAndLogs(false);
			if ( ! getLogFile().empty() ) {
				logFileStream.reset( new std::ofstream(getLogFile().c_str(), std::ios_base::out | std::ios_base::ate) );
				Log::Verbose().setOstream( *logFileStream );
				Log::Debug().setOstream( *logFileStream );
				LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), "log-file is: " << getLogFile().c_str());
			}

			if (getVerbose() > 0)
				LOG_VERBOSE(1, "Starting on " << OptionsBaseInterface::getHostname() << " version: " << VERSION);

			bool print = Log::printOptions();
			if (print) {
				_setVerbosityAndLogs(true);
			}

#ifdef _USE_OPENMP

			setOpt("threads", getMaxThreads(), print);

#endif
			validateOMPThreads();

			setOpt2("input-file", getInputFiles());

			setOpt("output-file", getOutputFile(), print);

			setOpt("temp-dir", getTmpDir(), print);

			setOpt("format-output", getFormatOutput(), print);

			setOpt("keep-read-comment", getKeepReadComment(), print);
			GlobalOptions::isCommentStored() = getKeepReadComment();

			setOpt("build-output-in-memory", getBuildOutputInMemory(), print);

			// set minimum quality score
			setOpt("min-quality-score", getMinQuality(), print);
			// TODO move to BimodalOptions
			//			if (getBimodalSigmas() >= 0 && (getMinReadLength() < getKmerSize() + 2)) {
			//				if(Logger::isMaster())
			//					LOG_WARN(1, "Bimodal Read Detection does not work unless min-read-length >= 2 + kmer-size");
			//			}

			// set the ignore quality value
			setOpt("ignore-quality", getIgnoreQual(), print);

			// set mmapInput
			setOpt("mmap-input", getMmapInput() , print);

			setOpt("gathered-logs", getGatheredLogs(), print);

			setOpt("batch-size", getBatchSize(), print);

			setOpt("separate-outputs", getSeparateOutputs(), print);

		} catch (std::exception& e) {
			LOG_ERROR(1,"Exception processing options" << std::endl << getDesc() << std::endl << e.what() << std::endl << "Exception processing options!" );
			return false;
		} catch (...) {
			LOG_ERROR(1, "Exception of unknown type!" << std::endl << getDesc() << std::endl );
			return false;
		}

		return ret;
	}

public:
	std::string &getInputFileSubstring(unsigned int fileIdx) {
		if (inputFilePrefixes.empty()) {
			// populate all the file indexes
			for(FileListType::iterator it = inputFiles.begin(); it != inputFiles.end(); it++) {
				size_t start = it->find_last_of('/');
				if (start == std::string::npos)
					start = 0;
				else
					start++;
				size_t end = it->find_last_of('.');
				if (end == std::string::npos)
					end = it->length() - 1;

				std::string fileprefix = it->substr( start, end-start );
				LOG_DEBUG_OPTIONAL(1, Logger::isMaster(), "InputFilePrefix: " << fileprefix );

				inputFilePrefixes.push_back( fileprefix );
			}
		}
		return inputFilePrefixes[fileIdx];
	}

	void validateOMPThreads() {
#ifdef _USE_OPENMP
		int maxThreads = omp_get_max_threads();
		LOG_DEBUG(2, "validating OpenMP threads: " << maxThreads);

		if (getMaxThreads() > maxThreads) {
			LOG_DEBUG(2, "Reducing the number of threads from " << getMaxThreads() << " to " << maxThreads);
			getMaxThreads() = maxThreads;
		}
		omp_set_num_threads(getMaxThreads());
#endif
	}


	unsigned int &getBatchSize()
	{
		return batchSize;
	}

	bool &getBuildOutputInMemory()
	{
		return buildOutputInMemory;
	}

	unsigned int &getFormatOutput()
	{
		return formatOutput;
	}

	bool &getKeepReadComment()
	{
		return keepReadComment;
	}

	bool &getGatheredLogs()
	{
		return gatheredLogs;
	}

	bool &getIgnoreQual()
	{
		return ignoreQual;
	}

	FileListType &getInputFilePrefixes()
	{
		return inputFilePrefixes;
	}

	FileListType &getInputFiles()
	{
		return inputFiles;
	}

	std::string &getLogFile()
	{
		return logFile;
	}

	OStreamPtr &getLogFileStream()
	{
		return logFileStream;
	}

	int &getMaxThreads()
	{
		return maxThreads;
	}


	unsigned int &getMinQuality()
	{
		return minQuality;
	}

	bool &getMmapInput()
	{
		return mmapInput;
	}

	std::string &getOutputFile()
	{
		return outputFile;
	}


	bool &getSeparateOutputs()
	{
		return separateOutputs;
	}


	std::string &getTmpDir()
	{
		return tmpDir;
	}

};

typedef OptionsBaseTemplate< _GeneralOptions > GeneralOptions;
typedef GeneralOptions Options;

#endif

