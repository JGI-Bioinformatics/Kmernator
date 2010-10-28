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

// put common, universal options in this class
// extend the class for specific options for each application
class Options {
public:
	typedef std::vector<std::string> FileListType;
	typedef boost::shared_ptr< std::ofstream > OStreamPtr;

protected:
	static inline Options &getOptions() {
		static Options singleton;
		return singleton;
	}

private:
	po::options_description desc;
	po::positional_options_description p;
	po::variables_map vm;

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
	unsigned int kmerSize;
	double       minKmerQuality;
	unsigned int minQuality;
	unsigned int minDepth;
	unsigned int depthRange;
	unsigned int minReadLength;
	unsigned int ignoreQual;
	unsigned int periodicSingletonPurge;
	unsigned int skipArtifactFilter;
	unsigned int artifactFilterMatchLength;
	unsigned int artifactFilterEditDistance;
	unsigned int maskSimpleRepeats;
	unsigned int phiXOutput;
	unsigned int filterOutput;
	unsigned int deDupMode;
	unsigned int deDupSingle;
	unsigned int deDupEditDistance;
	unsigned int deDupStartOffset;
	unsigned int deDupLength;
	unsigned int mmapInput;
	unsigned int buildPartitions;
	unsigned int gcHeatMap;

	Options() : maxThreads(OMP_MAX_THREADS_DEFAULT), tmpDir("/tmp"), formatOutput(0), kmerSize(21), minKmerQuality(0.10),
	minQuality(10), minDepth(2), depthRange(2), minReadLength(22), ignoreQual(0),
	periodicSingletonPurge(0), skipArtifactFilter(0), artifactFilterMatchLength(24), artifactFilterEditDistance(2),
	maskSimpleRepeats(1), phiXOutput(0), filterOutput(0),
	deDupMode(1), deDupSingle(0), deDupEditDistance(0), deDupStartOffset(0), deDupLength(16),
	mmapInput(1), buildPartitions(0), gcHeatMap(1) {
	}

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
	static inline unsigned int &getKmerSize() {
		return getOptions().kmerSize;
	}
	static inline double &getMinKmerQuality() {
		return getOptions().minKmerQuality;
	}
	static inline unsigned int &getVerbosity() {
		return Log::Verbose().getLevel();
	}
	static inline unsigned int &getDebug() {
		return Log::Debug().getLevel();
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
	static inline unsigned int &getMaskSimpleRepeats() {
		return getOptions().maskSimpleRepeats;
	}
	static inline unsigned int &getPhiXOutput(){
		return getOptions().phiXOutput;
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
	static inline unsigned int &getBuildPartitions() {
		return getOptions().buildPartitions;
	}
	static inline unsigned int &getGCHeatMap() {
		return getOptions().gcHeatMap;
	}
	const static unsigned int MAX_INT = (unsigned int) -1;

	static std::string &getInputFileSubstring(unsigned int fileIdx) {
		if (getOptions().inputFilePrefixes.empty()) {
			// populate all the file indexes
			for(FileListType::iterator it = getOptions().inputFiles.begin(); it != getOptions().inputFiles.end(); it++) {
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
				
				getOptions().inputFilePrefixes.push_back( fileprefix );
			}
		}
		return getOptions().inputFilePrefixes[fileIdx];
	}

	static void validateOMPThreads() {
#ifdef _USE_OPENMP
		int maxThreads = omp_get_max_threads();

		if (getMaxThreads() > maxThreads) {
			LOG_DEBUG(2, "Reducing the number of threads from " << getMaxThreads() << " to " << maxThreads);
			getMaxThreads() = maxThreads;
		}
		if ((getMaxThreads() & (getMaxThreads()-1)) != 0) {
			int start = getMaxThreads();
			int t = getMaxThreads();
			if (t > 32) {
				t=32;
			} else if (t > 16) {
				t=16;
			} else if (t > 8) {
				t=8;
			} else if (t > 4) {
				t=4;
			} else if (t > 2) {
				t=2;
			} else {
				t=1;
			}
			getMaxThreads() = t;
			LOG_DEBUG(2, "Reducing the number of threads from " << start << " to " << t );
		}
		omp_set_num_threads(getMaxThreads());
#endif
	}

protected:
	void setOptions() {

		desc.add_options()("help", "produce help message")

		("verbose", po::value<unsigned int>()->default_value(getVerbosity()),
				"level of verbosity (0+)")

		("debug", po::value<unsigned int>()->default_value(getDebug()),
				"level of debug verbosity (0+)")

#ifdef _USE_OPENMP
		("threads", po::value<int>()->default_value(maxThreads),
				"maximum number of threads. This must be a power of 2")
#endif
		("reference-file", po::value<FileListType>(), "set reference file(s)")

		("kmer-size", po::value<unsigned int>()->default_value(kmerSize), "kmer size.  A size of 0 will skip k-mer calculations")

		("input-file", po::value<FileListType>(), "input file(s)")

		("output-file", po::value<std::string>(), "output file pattern")

		("format-output", po::value<unsigned int>()->default_value(formatOutput),
				"0: fastq, 1: fasta, 2: fastq unmasked, 3: fasta unmasked")

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

		("mask-simple-repeats", po::value<unsigned int>()->default_value(maskSimpleRepeats),
				"if filtering artifacts, also mask simple repeats")

		("dedup-mode", po::value<unsigned int>()->default_value(deDupMode),
				"if 0, no fragment de-duplication will occur.  if 1, single orientation (AB and BA are separated) will collapse to consensus. if 2, both orientations (AB and BA are the same) will collapse")

        ("dedup-single", po::value<unsigned int>()->default_value(deDupSingle),
		    "if 0, no single read de-duplication will occur.  if 1, then single read deduplication will occur")

		("dedup-edit-distance", po::value<unsigned int>()->default_value(deDupEditDistance),
				"if -1, no fragment de-duplication will occur, if 0, only exact match, ...")

		("dedup-start-offset", po::value<unsigned int>()->default_value(deDupStartOffset),
				"de-duplication start offset to find unique fragments, must be multiple of 4")

		("dedup-length", po::value<unsigned int>()->default_value(deDupLength),
				"de-duplication length to find unique fragments, must be multiple of 4 (doubled when in single-end mode)")

		("mmap-input", po::value<unsigned int>()->default_value(mmapInput),
				"If set to 0, prevents input files from being mmaped, instead import reads into memory (somewhat faster if memory is abundant)")

		("build-partitions", po::value<unsigned int>()->default_value(buildPartitions),
				"If set, kmer spectrum will be computed in stages and then combined in mmaped files on disk.  Must be a power of 2")

		("gc-heat-map", po::value<unsigned int>()->default_value(gcHeatMap),
				"If set, a GC Heat map will be output (requires --output)")					;

	}

public:


	template<typename T>
	static void setOpt(std::string key, T &val, bool print = false) {
		po::variables_map &vm = getOptions().vm;
		if (vm.count(key.c_str())) {
			val = vm[key.c_str()].as<T>();
			if (print) {
				LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), key << " is: " << val );
			}
		} else if (print) {
			LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), key << " was not specified.");
		}

	}

	static po::options_description &getDesc() {
		return getOptions().desc;
	}
	static po::positional_options_description &getPosDesc() {
		return getOptions().p;
	}
	static po::variables_map &getVarMap() {
		return getOptions().vm;
	}

	static bool parseOpts(int argc, char *argv[]) {
		try {
			getOptions().setOptions();

			po::options_description & desc = getDesc();
			po::positional_options_description & p = getPosDesc();

			po::variables_map & vm = getVarMap();
			po::store(
					po::command_line_parser(argc, argv). options(desc).positional(
							p).run(), vm);
			po::notify(vm);

			if (vm.count("help")) {
				std::cerr << desc << std::endl;
				return false;
			}


			setOpt<unsigned int>("verbose", getVerbosity());
			{
				char hostname[128];
				gethostname(hostname, 128);
				LOG_VERBOSE(1, "Starting on " << hostname);
			}

			setOpt<unsigned int>("debug", getDebug());

			setOpt<std::string>("log-file", getLogFile(), false);
			if ( ! getLogFile().empty() ) {
				getOptions().logFileStream.reset( new std::ofstream(getLogFile().c_str(), std::ios_base::out | std::ios_base::ate) );
				Log::Verbose().setOstream( *getOptions().logFileStream );
				Log::Debug().setOstream( *getOptions().logFileStream );
				LOG_VERBOSE_OPTIONAL(1, Logger::isMaster(), "log-file is: " << getLogFile().c_str());
			}


			bool print = Logger::isMaster() && ((Log::isVerbose(1) || Log::isDebug(1)));
			std::ostream *output = NULL;
			if (print) {
				output = Log::Verbose("Options Set").getOstreamPtr();
			}

#ifdef _USE_OPENMP

			setOpt<int>("threads", getMaxThreads(), print);
			validateOMPThreads();

#endif

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
				LOG_ERROR(1, "There was no kmer size specified!");
				return false;
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
				LOG_ERROR(1, "There were no input files specified!" << std::endl << desc << std::endl << "There were no input files specified!");
				return false;
			}
			setOpt<std::string>("output-file", getOutputFile(), print);

			setOpt<std::string>("temp-dir", getTmpDir(), print);

			setOpt<unsigned int>("format-output", getFormatOutput(), print);

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
				getMinReadLength() = MAX_INT;
			}

			// set the ignore quality value
			setOpt<unsigned int>("ignore-quality", getIgnoreQual(), print);

			// set periodic singleton purge value
			setOpt<unsigned int>("periodic-singleton-purge", getPeriodicSingletonPurge(), print);

			// set skipArtifactFiltering
			setOpt<unsigned int>("skip-artifact-filter", getSkipArtifactFilter(), print);
			setOpt<unsigned int>("artifact-match-length", getArtifactFilterMatchLength(), print);
			setOpt<unsigned int>("artifact-edit-distance", getArtifactFilterEditDistance(), print);

			// set simple repeat masking
			setOpt<unsigned int>("mask-simple-repeats", getMaskSimpleRepeats() , print);

			// set phix masking
			setOpt<unsigned int>("phix-output", getPhiXOutput() , print);

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
				LOG_ERROR(1, "Unsuppored option dedup-start-offset and dedup-length must both be mulitples of 4!" << std::endl << desc << std::endl << "Unsuppored option dedup-start-offset and dedup-length must both be mulitples of 4!");
				return false;
			}

			// set mmapInput
			setOpt<unsigned int>("mmap-input", getMmapInput() , print);

			// set buildPartitions
			setOpt<unsigned int>("build-partitions", getBuildPartitions(), print);

			setOpt<unsigned int>("gc-heat-map", getGCHeatMap(), print);

		} catch (std::exception& e) {
			LOG_ERROR(1,"Exception processing options" << std::endl << getDesc() << std::endl << e.what() << std::endl << "Exception processing options!" );
			return false;
		} catch (...) {
			LOG_ERROR(1, "Exception of unknown type!" << std::endl << getDesc() << std::endl );
			return false;
		}

		return true;
	}
};

#endif

//
// $Log: Options.h,v $
// Revision 1.19  2010-08-18 17:50:39  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.18.4.1  2010-07-21 18:06:11  regan
// added options to change unique mask offset and length for de-duplication
//
// Revision 1.18  2010-06-23 22:15:15  regan
// added --threads option
//
// Revision 1.17  2010-06-23 20:58:02  regan
// fixed minimum read length logic
//
// Revision 1.16  2010-05-24 21:48:46  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.15.2.2  2010-05-19 22:43:49  regan
// bugfixes
//
// Revision 1.15.2.1  2010-05-19 00:20:46  regan
// refactored fomat output options
// added options to fastq2fasta
//
// Revision 1.15  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.14.6.4  2010-05-18 16:43:47  regan
// added GC heatmap output .. still refining
//
// Revision 1.14.6.3  2010-05-12 20:47:06  regan
// minor refactor.  adjusted defaults
//
// Revision 1.14.6.2  2010-05-12 19:52:10  regan
// refactored options
// removed obsolete parameters
// added new ones
//
// Revision 1.14.6.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.14  2010-05-04 23:49:21  regan
// merged changes for FixConsensusOutput-201000504
//
// Revision 1.13.6.2  2010-05-04 23:41:30  regan
// checkpoint
//
// Revision 1.13.6.1  2010-05-04 21:51:47  regan
// checkpoint
//
// Revision 1.13  2010-05-01 21:57:53  regan
// merged head with serial threaded build partitioning
//
// Revision 1.12.4.14  2010-05-01 05:57:40  regan
// made edit distance for dedup fragment pairs optional
//
// Revision 1.12.4.13  2010-04-30 16:33:21  regan
// better option descriptions
//
// Revision 1.12.4.12  2010-04-29 20:33:29  regan
// bugfix
//
// Revision 1.12.4.11  2010-04-29 16:54:07  regan
// bugfixes
//
// Revision 1.12.4.10  2010-04-29 06:58:55  regan
// modified default optoins
//
// Revision 1.12.4.9  2010-04-29 04:26:44  regan
// phix changes
//
// Revision 1.12.4.8  2010-04-28 22:28:10  regan
// refactored writing routines
//
// Revision 1.12.4.7  2010-04-28 16:57:00  regan
// bugfix in output filenames
//
// Revision 1.12.4.6  2010-04-27 23:17:50  regan
// fixed naming of output files
//
// Revision 1.12.4.5  2010-04-27 22:53:20  regan
// reworked defaults
//
// Revision 1.12.4.4  2010-04-27 18:25:20  regan
// bugfix in temp directory usage
//
// Revision 1.12.4.3  2010-04-27 05:48:46  regan
// added bulild in parts to main code and options
//
// Revision 1.12.4.2  2010-04-26 22:56:26  regan
// bugfix
//
// Revision 1.12.4.1  2010-04-26 22:53:08  regan
// added more global options
//
// Revision 1.12  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.11.2.2  2010-04-15 17:59:52  regan
// made mmap optional
//
// Revision 1.11.2.1  2010-04-15 17:47:19  regan
// added option to skip artifact filter
//
// Revision 1.11  2010-03-15 07:42:39  regan
// added kmer of 0 to skip kmer calculations
//
// Revision 1.10  2010-03-10 13:18:19  regan
// added quality ignore and singleton purge to save memory
//
// Revision 1.9  2010-03-03 17:08:52  regan
// fixed options to support -1 in min read length
//
// Revision 1.8  2010-03-02 15:03:33  regan
// added debug option and reformatted
//
// Revision 1.7  2010-02-26 13:01:16  regan
// reformatted
//
// Revision 1.6  2010-02-22 14:40:05  regan
// added optoin
//
// Revision 1.5  2010-01-16 01:06:28  regan
// refactored
//
// Revision 1.4  2010-01-14 01:17:48  regan
// working out options
//
// Revision 1.3  2010-01-14 00:47:44  regan
// added two common options
//
// Revision 1.2  2009-11-09 19:37:17  regan
// enhanced some debugging / analysis output
//
// Revision 1.1  2009-11-06 04:10:21  regan
// refactor of cmd line option handling
// added methods to evaluate spectrums
//
//
