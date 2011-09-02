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
private:
	OptionsInstance();
	~OptionsInstance();
};

class OptionsBaseInterface {
public:
	typedef std::vector<std::string> StringListType;
	typedef boost::shared_ptr< std::ofstream > OStreamPtr;
	typedef StringListType FileListType;

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
	template<typename U>
	static void setOpt(std::string key, U &val, bool print = false) {
		po::variables_map &vm = getVarMap();
		if (vm.count(key.c_str())) {
			val = vm[key.c_str()].as<U>();
			if (print) {
				LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), key << " is: " << val );
			}
		} else if (print) {
			LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), key << " was not specified.");
		}
	}
	static std::string getHostname() {
		char hostname[256];
		gethostname(hostname, 256);
		return std::string(hostname);
	}

	virtual void _setOptions(po::options_description &desc, po::positional_options_description &p) {}
	virtual bool _parseOpts(po::options_description &desc, po::positional_options_description &p, po::variables_map &vm, int argc, char *argv[]) {
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
	static bool parseOpts(int argc, char *argv[]) {
		getOptions()._setOptions(getDesc(), getPosDesc());
		return getOptions()._parseOpts(getDesc(), getPosDesc(), getVarMap(),
				argc, argv);
	}


protected:

	OptionsBaseTemplate() {} // hidden constructor
	virtual ~OptionsBaseTemplate() {}  // hidden destructor
	OptionsBaseTemplate(OptionsBaseTemplate const&); // disabled
	OptionsBaseTemplate& operator=(OptionsBaseTemplate const&); // disabled
};

class _GeneralOptions : public OptionsBaseInterface {
public:
	_GeneralOptions() : maxThreads(OMP_MAX_THREADS_DEFAULT), tmpDir("/tmp"), formatOutput(0), buildOutputInMemory(false), kmerSize(21), minKmerQuality(0.10),
	minQuality(5), minDepth(2), depthRange(2), minReadLength(25), bimodalSigmas(-1.0), variantSigmas(-1.0), ignoreQual(0),
	periodicSingletonPurge(0), skipArtifactFilter(0), artifactFilterMatchLength(24), artifactFilterEditDistance(2), buildArtifactEditsInFilter(2),
	maskSimpleRepeats(1), phiXOutput(0), filterOutput(0),
	deDupMode(1), deDupSingle(0), deDupEditDistance(0), deDupStartOffset(0), deDupLength(16),
	mmapInput(1), saveKmerMmap(0), buildPartitions(0), gcHeatMap(1), gatheredLogs(1), batchSize(100000), separateOutputs(1)
	{
		 char *tmpPath;
		 tmpPath = getenv ("TMPDIR");
		 if (tmpPath != NULL) {
			 tmpDir = std::string(tmpPath);
		 }
	}
	virtual ~_GeneralOptions() {}

private:
    int          maxThreads;
	FileListType referenceFiles;
	FileListType inputFiles;
	FileListType inputFilePrefixes;
	std::string  outputFile;
	std::string  logFile;
	OStreamPtr   logFileStream;
	std::string  tmpDir;
	unsigned int formatOutput;
	bool         buildOutputInMemory;
	unsigned int kmerSize;
	double       minKmerQuality;
	unsigned int minQuality;
	unsigned int minDepth;
	unsigned int depthRange;
	unsigned int minReadLength;
	double       bimodalSigmas;
	double       variantSigmas;
	unsigned int ignoreQual;
	unsigned int periodicSingletonPurge;
	unsigned int skipArtifactFilter;
	unsigned int artifactFilterMatchLength;
	unsigned int artifactFilterEditDistance;
	unsigned int buildArtifactEditsInFilter;
	unsigned int maskSimpleRepeats;
	unsigned int phiXOutput;
	FileListType artifactReferenceFiles;
	unsigned int filterOutput;
	unsigned int deDupMode;
	unsigned int deDupSingle;
	unsigned int deDupEditDistance;
	unsigned int deDupStartOffset;
	unsigned int deDupLength;
	unsigned int mmapInput;
	unsigned int saveKmerMmap;
	std::string  loadKmerMmap;
	unsigned int buildPartitions;
	unsigned int gcHeatMap;
	unsigned int gatheredLogs;
	unsigned int batchSize;
	unsigned int separateOutputs;

public:
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {

		desc.add_options()("help", "produce help message")

		("verbose", po::value<unsigned int>()->default_value(getVerbose()),
				"level of verbosity (0+)")

		("debug", po::value<unsigned int>()->default_value(getDebug()),
				"level of debug verbosity (0+)")

#ifdef _USE_OPENMP
		("threads", po::value<int>()->default_value(maxThreads),
				"maximum number of threads")
#endif
		("reference-file", po::value<FileListType>(), "set reference file(s)")

		("kmer-size", po::value<unsigned int>()->default_value(kmerSize), "kmer size.  A size of 0 will skip k-mer calculations")

		("input-file", po::value<FileListType>(), "input file(s)")

		("output-file", po::value<std::string>(), "output file pattern")

		("format-output", po::value<unsigned int>()->default_value(formatOutput),
				"0: fastq, 1: fasta, 2: fastq unmasked, 3: fasta unmasked")

		("build-output-in-memory", po::value<bool>()->default_value(buildOutputInMemory),
				"if set, all temporary output files will first be stored in memory (faster for MPI applications)")

		("phix-output", po::value<unsigned int>()->default_value(phiXOutput),
		        "if set, artifact filter also screens for PhiX174, and any matching reads will be output into a separate file (requires --output-file set)")

		("filter-output", po::value<unsigned int>()->default_value(filterOutput),
				"if set, artifact filter reads will be output into a separate file. If not set, then affected reads will be trimmed and then output normally.  (requires --output-file set)")

		("log-file", po::value<std::string>()->default_value(logFile),
				"If set all INFO and DEBUG messages will be logged here (default stderr)")

		("temp-dir", po::value<std::string>()->default_value(tmpDir), "temporary directory to deposit mmap file")

		("min-read-length",
				po::value<unsigned int>()->default_value(minReadLength),
				"minimum (trimmed) read length of selected reads.  0: no minimum, 1: full read length")

		("min-kmer-quality", po::value<double>()->default_value(minKmerQuality),
				"minimum quality-adjusted kmer probability (0-1)")

		("min-quality-score", po::value<unsigned int>()->default_value(minQuality),
				"minimum quality score over entire kmer")

		("min-depth", po::value<unsigned int>()->default_value(minDepth),
				"minimum depth for a solid kmer")

		("depth-range", po::value<unsigned int>()->default_value(depthRange),
				"if > min-depth, then output will be created in cycles of files ranging from min-depth to depth-range")

		("bimodal-sigmas", po::value<double>()->default_value(bimodalSigmas),
				"Detect bimodal kmer-signatures across reads and trim at transition point if the two means are separated by bimodal-sigmas * stdDev (2.0 to 3.0 suggested).  disabled if < 0.0")

		("variant-sigmas", po::value<double>()->default_value(variantSigmas),
				"Detect and purge kmer-variants if >= variant-sigmas * Poisson-stdDev (2.0-3.0 suggested).  disabled if < 0.0")

		("ignore-quality", po::value<unsigned int>()->default_value(ignoreQual),
				"ignore the quality score, to save memory or if they are untrusted")

		("periodic-singleton-purge", po::value<unsigned int>()->default_value(periodicSingletonPurge),
				"Purge singleton memory structure every # of reads")

		("skip-artifact-filter", po::value<unsigned int>()->default_value(skipArtifactFilter),
				"Skip homo-polymer, primer-dimer and duplicated fragment pair filtering")

		("artifact-match-length", po::value<unsigned int>()->default_value(artifactFilterMatchLength),
				"Kmer match length to known artifact sequences")

		("artifact-edit-distance", po::value<unsigned int>()->default_value(artifactFilterEditDistance),
				"edit-distance to apply to artifact-match-length matches to know artifacts")

		("build-artifact-edits-in-filter", po::value<unsigned int>()->default_value(buildArtifactEditsInFilter),
				"0 - edits will be applied to reads on the fly, 1 - edits will be pre-build in the filter (needs more memory, less overall CPU), 2 - automatic based on size")

		("mask-simple-repeats", po::value<unsigned int>()->default_value(maskSimpleRepeats),
				"if filtering artifacts, also mask simple repeats")

		("artifact-reference-file", po::value<FileListType>(), "additional artifact reference file(s)")

		("dedup-mode", po::value<unsigned int>()->default_value(deDupMode),
				"if 0, no fragment de-duplication will occur.  if 1, single orientation (AB and BA are separated) will collapse to consensus. if 2, both orientations (AB and BA are the same) will collapse")

        ("dedup-single", po::value<unsigned int>()->default_value(deDupSingle),
		    "if 0, no single read de-duplication will occur.  if 1, then single read deduplication will occur")

		("dedup-edit-distance", po::value<unsigned int>()->default_value(deDupEditDistance),
				"if -1, no fragment de-duplication will occur, if 0, only exact match, ...")

		("dedup-start-offset", po::value<unsigned int>()->default_value(deDupStartOffset),
				"de-duplication start offset to find unique fragments, must be multiple of 4")

		("dedup-length", po::value<unsigned int>()->default_value(deDupLength),
				"de-duplication length to find unique fragments, must be multiple of 4 (length on each read of a paired fragment or doubled when in single-end mode)")

		("mmap-input", po::value<unsigned int>()->default_value(mmapInput),
				"If set to 0, prevents input files from being mmaped, instead import reads into memory (somewhat faster if memory is abundant)")

		("save-kmer-mmap", po::value<unsigned int>()->default_value(saveKmerMmap),
				"If set to 1, creates a memory map of the kmer spectrum for later use")

		("load-kmer-mmap", po::value<std::string>(), "Instead of generating kmer spectrum, load an existing one (read-only) named by this option")

		("build-partitions", po::value<unsigned int>()->default_value(buildPartitions),
				"If set, kmer spectrum will be computed in stages and then combined in mmaped files on disk.")

		("gc-heat-map", po::value<unsigned int>()->default_value(gcHeatMap),
				"If set, a GC Heat map will be output (requires --output)")

		("gathered-logs", po::value<unsigned int>()->default_value(gatheredLogs),
				"If set and MPI is enabled, VERBOSE1, VERBOSE2 and DEBUG1 logs will be gathered to the master before being output.")

		("batch-size", po::value<unsigned int>()->default_value(batchSize),
				"default size of batches (reads, kmers, MPI, etc)")

		("separate-outputs", po::value<unsigned int>()->default_value(separateOutputs),
				"If set, each input (plus consensus) will generate a new outputfile.  If set to 0, all input files will be merged into one output file.")

		;

	}
	bool _parseOpts(po::options_description &desc, po::positional_options_description &p, po::variables_map &vm, int argc, char *argv[]) {
		try {
			po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
			po::notify(vm);

			if (vm.count("help")) {
				std::cerr << desc << std::endl;
				return false;
			}


			setOpt<unsigned int>("verbose", getVerbose());
			setOpt<unsigned int>("debug", getDebug());
			setOpt<std::string>("log-file", getLogFile(), false);
			if ( ! getLogFile().empty() ) {
				logFileStream.reset( new std::ofstream(getLogFile().c_str(), std::ios_base::out | std::ios_base::ate) );
				Log::Verbose().setOstream( *logFileStream );
				Log::Debug().setOstream( *logFileStream );
				LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), "log-file is: " << getLogFile().c_str());
			}

			if (getVerbose() > 0)
				LOG_VERBOSE(1, "Starting on " << OptionsBaseInterface::getHostname());

			bool print = Logger::isMaster() && ((Log::isVerbose(1) || Log::isDebug(1)));

			std::ostream *output = NULL;
			if (print) {
				output = Log::Verbose("Options Set").getOstreamPtr();
			}

#ifdef _USE_OPENMP

			setOpt<int>("threads", getMaxThreads(), print);

#endif
			validateOMPThreads();

			if (vm.count("reference-file")) {
				if (print)
					Log::Verbose() << "Reference files are: ";
				FileListType & referenceFiles = getReferenceFiles()
						= vm["reference-file"].as<FileListType> ();
				if (print) {
					for (FileListType::iterator it = referenceFiles.begin(); it
						!= referenceFiles.end(); it++)
						*output << *it << ", ";
					*output << std::endl;
				}
			} else if (print) {
				Log::Verbose() << "reference was not set." << std::endl;
			}

			if (vm.count("kmer-size")) {
				setOpt<unsigned int>("kmer-size", getKmerSize(), print);
			} else {
				LOG_WARN(1, "There was no kmer size specified!");
			}
			if (vm.count("input-file")) {
				if (print) {
					Log::Verbose() << "Input files are: ";
				}
				FileListType inputs = getInputFiles() = vm["input-file"].as<FileListType> ();
				if (print) {
					for (FileListType::iterator it = inputs.begin(); it
							!= inputs.end(); it++)
						*output << *it << ", ";
					*output << std::endl;
				}
			} else  {
				LOG_WARN(1, "There were no input files specified!");
			}

			setOpt<std::string>("output-file", getOutputFile(), print);

			setOpt<std::string>("temp-dir", getTmpDir(), print);

			setOpt<unsigned int>("format-output", getFormatOutput(), print);

			setOpt<bool>("build-output-in-memory", getBuildOutputInMemory(), print);

			// set kmer quality
			setOpt<double>("min-kmer-quality", getMinKmerQuality(), print);

			// set minimum quality score
			setOpt<unsigned int>("min-quality-score", getMinQuality(), print);

			// set minimum depth
			setOpt<unsigned int>("min-depth", getMinDepth(), print);

			setOpt<unsigned int>("depth-range", getDepthRange(), print);

			// set read length
			setOpt<unsigned int>("min-read-length", getMinReadLength(), print);
			if (getMinReadLength() == 1) {
				getMinReadLength() = MAX_SEQUENCE_LENGTH;
			}

			setOpt<double>("bimodal-sigmas", getBimodalSigmas(), print);
			setOpt<double>("variant-sigmas", getVariantSigmas(), print);

			if (getBimodalSigmas() >= 0 && (getMinReadLength() < getKmerSize() + 2)) {
				if(Logger::isMaster())
					LOG_WARN(1, "Bimodal Read Detection does not work unless min-read-length >= 2 + kmer-size");
			}
			// set the ignore quality value
			setOpt<unsigned int>("ignore-quality", getIgnoreQual(), print);

			// set periodic singleton purge value
			setOpt<unsigned int>("periodic-singleton-purge", getPeriodicSingletonPurge(), print);

			// set skipArtifactFiltering
			setOpt<unsigned int>("skip-artifact-filter", getSkipArtifactFilter(), print);
			setOpt<unsigned int>("artifact-match-length", getArtifactFilterMatchLength(), print);
			setOpt<unsigned int>("artifact-edit-distance", getArtifactFilterEditDistance(), print);
			setOpt<unsigned int>("build-artifact-edits-in-filter", getBuildArtifactEditsInFilter(), print);

			// set simple repeat masking
			setOpt<unsigned int>("mask-simple-repeats", getMaskSimpleRepeats() , print);

			// set phix masking
			setOpt<unsigned int>("phix-output", getPhiXOutput() , print);

			if (vm.count("artifact-reference-file")) {
				if (print) {
					Log::Verbose() << "Artifact Reference files are: ";
				}
				FileListType artifacts = getArtifactReferenceFiles() = vm["artifact-reference-file"].as<FileListType> ();
				if (print) {
					for (FileListType::iterator it = artifacts.begin(); it
							!= artifacts.end(); it++)
						*output << *it << ", ";
					*output << std::endl;
				}
			}
			// set simple repeat masking
			setOpt<unsigned int>("filter-output", getFilterOutput() , print);

			// set dedup mode
			setOpt<unsigned int>("dedup-mode", getDeDupMode() , print);

			// set dedup single
			setOpt<unsigned int>("dedup-single", getDeDupSingle() , print);

			// set dedup edit distance
			setOpt<unsigned int>("dedup-edit-distance", getDeDupEditDistance() , print);
			if (getDeDupEditDistance() > 1) {
				LOG_ERROR(1, "Unsupported option dedup-edit-distance > 1!" << std::endl << desc << std::endl << "Unsupported option dedup-edit-distance > 1!");
			    return false;
			}
			setOpt<unsigned int>("dedup-start-offset", getDeDupStartOffset(), print);
			setOpt<unsigned int>("dedup-length", getDeDupLength(), print);
			if (getDeDupStartOffset() % 4 != 0 || getDeDupLength() % 4 != 0) {
				LOG_ERROR(1, "Unsupported option dedup-start-offset and dedup-length must both be mulitples of 4!" << std::endl << desc << std::endl << "Unsuppored option dedup-start-offset and dedup-length must both be mulitples of 4!");
				return false;
			}

			// set mmapInput
			setOpt<unsigned int>("mmap-input", getMmapInput() , print);

			setOpt<unsigned int>("save-kmer-mmap", getSaveKmerMmap(), print);

			setOpt<std::string>("load-kmer-mmap", getLoadKmerMmap(), print);

			// set buildPartitions
			setOpt<unsigned int>("build-partitions", getBuildPartitions(), print);

			setOpt<unsigned int>("gc-heat-map", getGCHeatMap(), print);

			setOpt<unsigned int>("gathered-logs", getGatheredLogs(), print);

			setOpt<unsigned int>("batch-size", getBatchSize(), print);

			setOpt<unsigned int>("separate-outputs", getSeparateOutputs(), print);

		} catch (std::exception& e) {
			LOG_ERROR(1,"Exception processing options" << std::endl << getDesc() << std::endl << e.what() << std::endl << "Exception processing options!" );
			return false;
		} catch (...) {
			LOG_ERROR(1, "Exception of unknown type!" << std::endl << getDesc() << std::endl );
			return false;
		}

		return true;
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

	unsigned int &getArtifactFilterEditDistance()
	{
	    return artifactFilterEditDistance;
	}

	unsigned int &getArtifactFilterMatchLength()
	{
	    return artifactFilterMatchLength;
	}

	FileListType &getArtifactReferenceFiles()
	{
	    return artifactReferenceFiles;
	}

	unsigned int &getBatchSize()
	{
	    return batchSize;
	}

	double &getBimodalSigmas()
	{
	    return bimodalSigmas;
	}

	unsigned int &getBuildArtifactEditsInFilter()
	{
	    return buildArtifactEditsInFilter;
	}

	bool &getBuildOutputInMemory()
	{
	    return buildOutputInMemory;
	}

	unsigned int &getBuildPartitions()
	{
	    return buildPartitions;
	}

	unsigned int &getDeDupEditDistance()
	{
	    return deDupEditDistance;
	}

	unsigned int &getDeDupLength()
	{
	    return deDupLength;
	}

	unsigned int &getDeDupMode()
	{
	    return deDupMode;
	}

	unsigned int &getDeDupSingle()
	{
	    return deDupSingle;
	}

	unsigned int &getDeDupStartOffset()
	{
	    return deDupStartOffset;
	}

	unsigned int &getDepthRange()
	{
	    return depthRange;
	}

	unsigned int &getFilterOutput()
	{
	    return filterOutput;
	}

	unsigned int &getFormatOutput()
	{
	    return formatOutput;
	}

	unsigned int &getGatheredLogs()
	{
	    return gatheredLogs;
	}

	unsigned int &getGCHeatMap()
	{
	    return gcHeatMap;
	}

	unsigned int &getIgnoreQual()
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

	unsigned int &getKmerSize()
	{
	    return kmerSize;
	}

	std::string &getLoadKmerMmap()
	{
	    return loadKmerMmap;
	}

	std::string &getLogFile()
	{
	    return logFile;
	}

	OStreamPtr &getLogFileStream()
	{
	    return logFileStream;
	}

	unsigned int &getMaskSimpleRepeats()
	{
	    return maskSimpleRepeats;
	}

	int &getMaxThreads()
	{
	    return maxThreads;
	}

	unsigned int &getMinDepth()
	{
	    return minDepth;
	}

	double &getMinKmerQuality()
	{
	    return minKmerQuality;
	}

	unsigned int &getMinQuality()
	{
	    return minQuality;
	}

	unsigned int &getMinReadLength()
	{
	    return minReadLength;
	}

	unsigned int &getMmapInput()
	{
	    return mmapInput;
	}

	std::string &getOutputFile()
	{
	    return outputFile;
	}

	unsigned int &getPeriodicSingletonPurge()
	{
	    return periodicSingletonPurge;
	}

	unsigned int &getPhiXOutput()
	{
	    return phiXOutput;
	}

	FileListType &getReferenceFiles()
	{
	    return referenceFiles;
	}

	unsigned int &getSaveKmerMmap()
	{
	    return saveKmerMmap;
	}

	unsigned int &getSeparateOutputs()
	{
	    return separateOutputs;
	}

	unsigned int &getSkipArtifactFilter()
	{
	    return skipArtifactFilter;
	}

	std::string &getTmpDir()
	{
	    return tmpDir;
	}

	double &getVariantSigmas()
	{
	    return variantSigmas;
	}


};

typedef OptionsBaseTemplate< _GeneralOptions > Options;

// extend the class for specific options for each application
// put common, universal options in this class
class _Options : public OptionsBaseInterface, public OptionsBaseTemplate<_Options> {
	friend class OptionsBaseTemplate<_Options>;
private:
	// cache of variables (for inline lookup and defaults)
    int          maxThreads;
	FileListType referenceFiles;
	FileListType inputFiles;
	FileListType inputFilePrefixes;
	std::string  outputFile;
	std::string  logFile;
	OStreamPtr   logFileStream;
	std::string  tmpDir;
	unsigned int formatOutput;
	bool         buildOutputInMemory;
	unsigned int kmerSize;
	double       minKmerQuality;
	unsigned int minQuality;
	unsigned int minDepth;
	unsigned int depthRange;
	unsigned int minReadLength;
	double       bimodalSigmas;
	double       variantSigmas;
	unsigned int ignoreQual;
	unsigned int periodicSingletonPurge;
	unsigned int skipArtifactFilter;
	unsigned int artifactFilterMatchLength;
	unsigned int artifactFilterEditDistance;
	unsigned int buildArtifactEditsInFilter;
	unsigned int maskSimpleRepeats;
	unsigned int phiXOutput;
	FileListType artifactReferenceFiles;
	unsigned int filterOutput;
	unsigned int deDupMode;
	unsigned int deDupSingle;
	unsigned int deDupEditDistance;
	unsigned int deDupStartOffset;
	unsigned int deDupLength;
	unsigned int mmapInput;
	unsigned int saveKmerMmap;
	std::string  loadKmerMmap;
	unsigned int buildPartitions;
	unsigned int gcHeatMap;
	unsigned int gatheredLogs;
	unsigned int batchSize;
	unsigned int separateOutputs;

	_Options() : maxThreads(OMP_MAX_THREADS_DEFAULT), tmpDir("/tmp"), formatOutput(0), buildOutputInMemory(false), kmerSize(21), minKmerQuality(0.10),
	minQuality(5), minDepth(2), depthRange(2), minReadLength(25), bimodalSigmas(-1.0), variantSigmas(-1.0), ignoreQual(0),
	periodicSingletonPurge(0), skipArtifactFilter(0), artifactFilterMatchLength(24), artifactFilterEditDistance(2), buildArtifactEditsInFilter(2),
	maskSimpleRepeats(1), phiXOutput(0), filterOutput(0),
	deDupMode(1), deDupSingle(0), deDupEditDistance(0), deDupStartOffset(0), deDupLength(16),
	mmapInput(1), saveKmerMmap(0), buildPartitions(0), gcHeatMap(1), gatheredLogs(1), batchSize(1000000), separateOutputs(1)
	{
		 char *tmpPath;
		 tmpPath = getenv ("TMPDIR");
		 if (tmpPath != NULL) {
			 tmpDir = std::string(tmpPath);
		 }
	}
	virtual ~_Options() {}

public:

	static inline int &getMaxThreads() {
		return getOptions().maxThreads;
	}
	static inline FileListType &getReferenceFiles() {
		return getOptions().referenceFiles;
	}
	static inline FileListType &getInputFiles() {
		return getOptions().inputFiles;
	}
	static inline std::string &getOutputFile() {
		return getOptions().outputFile;
	}
	static inline std::string &getLogFile() {
		return getOptions().logFile;
	}
	static inline std::string &getTmpDir() {
		return getOptions().tmpDir;
	}
	static inline unsigned int &getFormatOutput() {
		return getOptions().formatOutput;
	}
	static inline bool &getBuildOutputInMemory() {
		return getOptions().buildOutputInMemory;
	}
	static inline unsigned int &getKmerSize() {
		return getOptions().kmerSize;
	}
	static inline double &getMinKmerQuality() {
		return getOptions().minKmerQuality;
	}
	static inline unsigned int &getVerbose() {
		return Log::Verbose().setLevel();
	}
	static inline unsigned int &getDebug() {
		return Log::Debug().setLevel();
	}
	static inline unsigned int &getMinQuality() {
		return getOptions().minQuality;
	}
	static inline unsigned int &getMinDepth() {
		return getOptions().minDepth;
	}
	static inline unsigned int &getDepthRange() {
		return getOptions().depthRange;
	}
	static inline unsigned int &getMinReadLength() {
		return getOptions().minReadLength;
	}
	static inline double &getBimodalSigmas() {
		return getOptions().bimodalSigmas;
	}
	static inline double &getVariantSigmas() {
		return getOptions().variantSigmas;
	}
	static inline unsigned int &getIgnoreQual() {
		return getOptions().ignoreQual;
	}
	static inline unsigned int &getPeriodicSingletonPurge() {
		return getOptions().periodicSingletonPurge;
	}
	static inline unsigned int &getSkipArtifactFilter() {
		return getOptions().skipArtifactFilter;
	}
	static inline unsigned int &getArtifactFilterMatchLength() {
		return getOptions().artifactFilterMatchLength;
	}
	static inline unsigned int &getArtifactFilterEditDistance() {
		return getOptions().artifactFilterEditDistance;
	}
	static inline unsigned int &getBuildArtifactEditsInFilter() {
		return getOptions().buildArtifactEditsInFilter;
	}
	static inline unsigned int &getMaskSimpleRepeats() {
		return getOptions().maskSimpleRepeats;
	}
	static inline unsigned int &getPhiXOutput(){
		return getOptions().phiXOutput;
	}
	static inline FileListType &getArtifactReferenceFiles() {
		return getOptions().artifactReferenceFiles;
	}
	static inline unsigned int &getFilterOutput(){
		return getOptions().filterOutput;
	}
	static inline unsigned int &getDeDupMode() {
		return getOptions().deDupMode;
	}
	static inline unsigned int &getDeDupSingle() {
		return getOptions().deDupSingle;
	}
	static inline unsigned int &getDeDupEditDistance() {
		return getOptions().deDupEditDistance;
	}
	static inline unsigned int &getDeDupStartOffset() {
		return getOptions().deDupStartOffset;
	}
	static inline unsigned int &getDeDupLength() {
		return getOptions().deDupLength;
	}
	static inline unsigned int &getMmapInput() {
		return getOptions().mmapInput;
	}
	static inline unsigned int &getSaveKmerMmap() {
		return getOptions().saveKmerMmap;
	}
	static inline std::string  &getLoadKmerMmap() {
		return getOptions().loadKmerMmap;
	}
	static inline unsigned int &getBuildPartitions() {
		return getOptions().buildPartitions;
	}
	static inline unsigned int &getGCHeatMap() {
		return getOptions().gcHeatMap;
	}
	static inline unsigned int &getGatheredLogs() {
		return getOptions().gatheredLogs;
	}
	static inline unsigned int &getBatchSize() {
		return getOptions().batchSize;
	}
	static inline unsigned int &getSeparateOutputs() {
		return getOptions().separateOutputs;
	}


protected:
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {

	}

	bool _parseOpts(po::options_description &desc, po::positional_options_description &p, po::variables_map &vm, int argc, char *argv[]) {
		return false;
	}
};

#endif

