//
// Kmernator/apps/FixPair.cpp
//
// Author: Rob Egan
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

#include <iostream>
#include <cstdlib>
#include <cstring>

#include <boost/shared_ptr.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "Options.h"
#include "Log.h"
#include "Utils.h"
#ifdef ENABLE_MPI
#include "DistributedFunctions.h"
#endif
 #include <unistd.h>

using namespace std;

class _SSOptions : public OptionsBaseInterface {
public:
	typedef OptionsBaseInterface::StringListType StringListType;

	static int &getDefaultNumFiles() {
		static int dummy = 1;
		return dummy;
	}
	static int &getDefaultFileNum() {
		static int dummy2 = 0;
		return dummy2;
	}
	int &getNumFiles() {
		return numFiles;
	}
	int &getFileNum() {
		return fileNum;
	}
	std::string &getSplitFile() {
		return splitFile;
	}
	int &getEvenChunks() {
		return evenChunks;
	}
	int &getFifoFile() {
		return fifoFile;
	}
	StringListType &getForkCommand() {
		return forkCommand;
	}
	std::string &getPipeCommand() {
		return pipeCommand;
	}
	StringListType &getMergeList() {
		return mergeList;
	}
	StringListType &getExtraFifo() {
		return extraFifo;
	}

	_SSOptions() : splitFile(), pipeCommand(), mergeList(), forkCommand(), extraFifo(), numFiles(getDefaultNumFiles()), fileNum(getDefaultFileNum()), evenChunks(1), fifoFile(0) {}

	std::string _replaceWithKeys(std::string input) {
		size_t pos;
		char buf[16];
		while ((pos = input.find("{Uniq}")) != std::string::npos) {
			sprintf(buf, "%06d", getFileNum());
			std::string buf2 = std::string(buf) + "of";
			sprintf(buf, "%06d", getNumFiles());
			buf2 += std::string(buf);
			LOG_DEBUG_OPTIONAL(2, true, "Replacing {Uniq} (" << buf2 << ") in " << input);
			input.replace(pos, 6, buf2);
		}

		while ((pos = input.find("{FileNum}")) != std::string::npos) {
			sprintf(buf, "%06d", getFileNum());
			std::string buf2(buf);
			LOG_DEBUG_OPTIONAL(2, true, "Replacing {FileNum} (" << buf2 << ") in " << input);
			input.replace(pos, 9, buf2);
		}
		while ((pos = input.find("{NumFiles}")) != std::string::npos) {
			sprintf(buf, "%06d", getNumFiles());
			std::string buf2(buf);
			LOG_DEBUG_OPTIONAL(2, true, "Replacing {NumFiles}  (" << buf2 << ") in " << input);
			input.replace(pos, 10, buf2);
		}
		LOG_DEBUG_OPTIONAL(2, true, "final string " << input);
		return input;
	}
	StringListType _replaceWithKeys(StringListType inputs) {
		StringListType output;
		for(unsigned int i = 0; i < inputs.size(); i++) {
			output.push_back( _replaceWithKeys(inputs[i]));
		}
		return output;
	}

	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		Options::getOptions().getVerbose() = 0;
		Options::getOptions().getDebug() = 0;
		Options::getOptions().getMmapInput() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("Split Sequence Options");
		opts.add_options()
						("num-files", po::value<int>()->default_value(getDefaultNumFiles()), "The number of files to split into N")
						("file-num",  po::value<int>()->default_value(getDefaultFileNum()), "The number of the file to output (0-(N-1))")
						("even-chunks", po::value<int>()->default_value(evenChunks), "if > 1 then the output of each partition will be spread out across the file (recommend 10 if the ordering is less important than predictable runtime of forked commands)")
						("pipe-command", po::value<std::string>(), "a command to pipe the portion of the file(s) into.  Use the keyword variables '{FileNum}' and '{NumFiles}' to replace with MPI derived values")
						("merge", po::value<StringListType>(), "two arguments.  First is per-mpi file (use keywords) second is final file; can be specified multiple times")
						("split-file", po::value<std::string>()->default_value(splitFile), "if set, paired (second) reads will be sent to this file.  first reads will go to --output-file")
						("output-fifo", po::value<int>()->default_value(fifoFile), "if set, --output-file (and --split-file) will be fifo files which will be created on start and deleted on exit")
						("extra-fifo", po::value<StringListType>(), "optional additional fifo file(s) (to potentially use with --fork-command)")
						("fork-command", po::value<StringListType>(), "optional command(s) to fork between opening output files and wait to finish before starting any merges")
						;

		desc.add(opts);
		GeneralOptions::_setOptions(desc, p);
		p.add("input-file", -1);

	}
	bool _parseOptions(po::variables_map &vm) {
		// set options specific to this program
		bool ret = GeneralOptions::_parseOptions(vm);
		setOpt<int>("num-files", getNumFiles());
		setOpt<int>("file-num", getFileNum());
		setOpt<int>("even-chunks", getEvenChunks());
		setOpt<string>("split-file", getSplitFile());
		setOpt<string>("pipe-command", getPipeCommand());
		setOpt<int>("output-fifo", getFifoFile());
		setOpt2("merge", getMergeList());
		setOpt2("fork-command", getForkCommand());
		setOpt2("extra-fifo", getExtraFifo());


		getPipeCommand() = _replaceWithKeys(getPipeCommand());
		getForkCommand() = _replaceWithKeys(getForkCommand());
		getExtraFifo() = _replaceWithKeys(getExtraFifo());
		getMergeList() = _replaceWithKeys(getMergeList());
		getSplitFile() = _replaceWithKeys(getSplitFile());
		GeneralOptions::getOptions().getOutputFile() = _replaceWithKeys(GeneralOptions::getOptions().getOutputFile());
		GeneralOptions::getOptions().getInputFiles() = _replaceWithKeys(GeneralOptions::getOptions().getInputFiles());

		if (Options::getOptions().getInputFiles().empty() || getNumFiles() == 0 || getFileNum() >= getNumFiles()) {
			ret = false;
			std::stringstream ss;
			ss <<"Please specify num-files, file-num and at least one input file.\nnum-files=" << getNumFiles() << "\nfile-num=" << getFileNum();
			setOptionsErrorMsg(ss.str());
		}
		if ((!getPipeCommand().empty()) && (!getSplitFile().empty())) {
			ret = false;
			setOptionsErrorMsg("You can not specify both --pipe-command and --split-file!");
		}
		if (Options::getOptions().getOutputFile().empty() && !getSplitFile().empty()) {
			ret = false;
			setOptionsErrorMsg("If you specify --split-file, you must also specify --output-file.");
		}
		return ret;
	}
private:
	std::string splitFile, pipeCommand;
	StringListType mergeList, forkCommand, extraFifo;
	int numFiles, fileNum, evenChunks, fifoFile;
};
typedef OptionsBaseTemplate< _SSOptions > SSOptions;
typedef OptionsBaseInterface::StringListType StringListType;

void outputRegularFiles(OptionsBaseInterface::FileListType inputs, std::ostream &output1) {

	FormatOutput outputFormat = FormatOutput::getDefault();
	int chunks = SSOptions::getOptions().getEvenChunks();
	for (unsigned int i = 0 ; i < inputs.size(); i++) {
		ReadFileReader reader(inputs[i], "");
		for(int j = 0 ; j < chunks; j++) {
			reader.seekToPartition( SSOptions::getOptions().getFileNum()+j*chunks, SSOptions::getOptions().getNumFiles()*chunks );
			std::string name, bases, quals, comment;
			while (reader.nextRead(name, bases, quals, comment)) {
				output1 << outputFormat.toString(name, bases, quals, comment);
			}
		}
	}

}
typedef boost::shared_ptr< Kmernator::MmapSource > MmapPtr;
typedef  boost::shared_ptr< ReadFileReader > RFRPtr;

void outputSplitFiles(OptionsBaseInterface::FileListType inputs, string outputFile1, string outputFile2) {

	int numInputs = inputs.size();
	std::vector< MmapPtr > mmaps;
	mmaps.reserve(numInputs);
	std::vector< RFRPtr > readers1, readers2;
	readers1.reserve(numInputs);
	readers2.reserve(numInputs);
	int chunks = SSOptions::getOptions().getEvenChunks();

	for(int i=0; i < numInputs; i++) {
		mmaps.push_back( MmapPtr(new Kmernator::MmapSource(inputs[i], FileUtils::getFileSize(inputs[i]))) );
		readers1.push_back(RFRPtr( new ReadFileReader(*mmaps[i])) );
		readers2.push_back(RFRPtr( new ReadFileReader(*mmaps[i])) );
	}


	FormatOutput outputFormat = FormatOutput::getDefault();
	omp_set_dynamic(0);
	omp_set_num_threads(2);
	unsigned long numReads1 = 0, numReads2 = 0;
	#pragma omp single
	for(int i=0; i < numInputs; i++)
		madvise(const_cast<char*>(mmaps[i]->data()), mmaps[i]->size(), MADV_DONTNEED);

    #pragma omp parallel num_threads(2)
	{
		if (omp_get_num_threads() != 2)
			LOG_THROW("Invalid number of threads! " << omp_get_num_threads());

		if (omp_get_thread_num() == 0) {
			std::ofstream output1(outputFile1.c_str());
			std::string name, bases, quals, comment;
			for(int i=0; i < numInputs; i++) {
				for(int j = 0 ; j < chunks; j++) {
					readers1[i]->seekToPartition( SSOptions::getOptions().getFileNum()+j*chunks, SSOptions::getOptions().getNumFiles()*chunks );
					//readers2[i]->seekToPartition( SSOptions::getOptions().getFileNum()+j*chunks, SSOptions::getOptions().getNumFiles()*chunks );

					madvise(const_cast<char*>(mmaps[i]->data()+readers1[i]->getPos()), readers1[i]->getLastPos() - readers1[i]->getPos(), MADV_SEQUENTIAL);

					while (readers1[i]->nextRead(name, bases, quals, comment)) {
						if ((numReads1 & 0x01) == 0) {
							output1 << outputFormat.toString(name, bases, quals, comment);

						}
						numReads1++;
						LOG_DEBUG_OPTIONAL(2, (numReads1 & 0xffff) == 0, "reading and writing split file " << i << " " << inputs[i] << " for read1 " << numReads1/2);
					}
				}
				LOG_DEBUG_OPTIONAL(2, true, "Finished reading and writing split file " << i << " " << inputs[i] << " for read1. " << numReads1 << " " << readers1[i]->getPos());
			}

			output1.close();
			LOG_DEBUG_OPTIONAL(1, true, "Finished reading and writing split files for read1. " << numReads1);

		} else {
			std::ofstream output2(outputFile2.c_str());
			std::string name, bases, quals, comment;
			for(int i=0; i < numInputs; i++) {
				for(int j = 0 ; j < chunks; j++) {
					//readers1[i]->seekToPartition( SSOptions::getOptions().getFileNum()+j*chunks, SSOptions::getOptions().getNumFiles()*chunks );
					readers2[i]->seekToPartition( SSOptions::getOptions().getFileNum()+j*chunks, SSOptions::getOptions().getNumFiles()*chunks );

					madvise(const_cast<char*>(mmaps[i]->data()+readers2[i]->getPos()), readers2[i]->getLastPos() - readers2[i]->getPos(), MADV_SEQUENTIAL);

					while (readers2[i]->nextRead(name, bases, quals, comment)) {
						if ((numReads2 & 0x01) == 1) {
							output2 << outputFormat.toString(name, bases, quals, comment);
						}
						numReads2++;
						LOG_DEBUG_OPTIONAL(2, (numReads2 & 0xffff) == 0, "reading and writing split file " << i << " " << inputs[i] << " for read2 " << numReads2/2);
					}
				}
				LOG_DEBUG_OPTIONAL(2, true, "Finished reading and writing split file " << i << " " << inputs[i] << " for read2. " << numReads2 << " " << readers2[i]->getPos());
			}

			output2.close();
			LOG_DEBUG_OPTIONAL(1, true, "Finished reading and writing split files for read2. " << numReads2);

		}

		LOG_DEBUG_OPTIONAL(1, true, "Finished splitFiles parallel " << numReads1 << " " << numReads2);
		#pragma omp barrier
	}



	if (numReads1 != numReads2) {
		LOG_THROW("Invalid number of paired reads! " << numReads1/2 << " vs " << numReads2/2);
	}
	LOG_DEBUG_OPTIONAL(1, true, "Finished outputSplitFiles");

}

std::vector< int > forkCommand() {
	std::vector< int > forks;
	StringListType forkCommand = SSOptions::getOptions().getForkCommand();
	if (forkCommand.empty())
		return forks;

	for(int i = 0; i < (int) forkCommand.size(); i++) {
		LOG_DEBUG_OPTIONAL(1, true, "Starting " << forkCommand[i]);

		int child = Fork::forkCommand(forkCommand[i]);
		forks.push_back(child);

		if (0) {
		int child = fork();
		if (child < 0)
			LOG_THROW("Could not fork child process");

		if (child == 0) {
			// child
			if (0) {
			std::vector<std::string> forkCommandArgs;
			boost::split(forkCommandArgs, forkCommand[i], boost::is_any_of("\t "));
			char *forkArgs[forkCommandArgs.size()];
			forkArgs[forkCommandArgs.size()-1] = NULL;
			for(int j = 0 ; j < (int) forkCommandArgs.size(); j++)
				forkArgs[j] = strdup(forkCommandArgs[j].c_str());

			std::stringstream ss;
			ss << "Executing ( " << getpid() << ") : " << forkArgs[0] << " ";
			for (int j=1; j < (int) forkCommandArgs.size(); j++)
				ss << forkArgs[j] << "\t";
			ss << std::endl;
			std::string msg = ss.str();
			std::cerr << msg;
			execvp(forkArgs[0], forkArgs);
			std::cerr << "Should not get here!" << std::endl;
			exit(1);
			} else {
			int ret = system(forkCommand[i].c_str());
			if (ret != 0)
				LOG_ERROR(1, "Failed to execute: " << forkCommand[i]);
			exit(ret);
			}
		} else {
			// parent
			LOG_DEBUG_OPTIONAL(1, true, "forkPid: " << child);
			forks.push_back( child );
			Cleanup::trackChild( child );
		}
		}
	}
	return forks;
}

int main(int argc, char *argv[]) {

	int exitStatus = 0;

#ifdef ENABLE_MPI
	MPI_Init(&argc, &argv);
	mpi::environment env(argc, argv);
	mpi::communicator world;
	SSOptions::getOptions().getDefaultNumFiles() = world.size();
	SSOptions::getOptions().getDefaultFileNum() = world.rank();
	Logger::setWorld(&world);
#endif
	SSOptions::parseOpts(argc, argv);
	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();

	std::string outputFile = Options::getOptions().getOutputFile();
	std::string splitFile = SSOptions::getOptions().getSplitFile();

	OPipestream *ops = NULL;
	ostream *out1 = &std::cout;
	ofstream *out1f = NULL;
	std::vector< int > forkPids;
	StringListType extraFifos = SSOptions::getOptions().getExtraFifo();

	try {

		if (!extraFifos.empty()) {
			for(int i = 0; i < (int) extraFifos.size(); i++) {
				LOG_DEBUG_OPTIONAL(1, true, "making extra fifo: " << extraFifos[i]);
				unlink(extraFifos[i].c_str());
				mkfifo(extraFifos[i].c_str(), 0700);
				Cleanup::addTemp(extraFifos[i]);
			}
		}

		if (SSOptions::getOptions().getFifoFile()) {

			if (!outputFile.empty()) {
				LOG_DEBUG_OPTIONAL(1, true, "making output fifo: " << outputFile);
				unlink(outputFile.c_str());
				mkfifo(outputFile.c_str(), 0700);
				Cleanup::addTemp(outputFile);
			}
			if (!splitFile.empty()) {
				LOG_DEBUG_OPTIONAL(1, true, "making split fifo: " << splitFile);
				unlink(splitFile.c_str());
				mkfifo(splitFile.c_str(), 0700);
				Cleanup::addTemp(splitFile);
			}

			forkPids = forkCommand();

		}

		if (!splitFile.empty()) {

			outputSplitFiles(inputs, outputFile, splitFile);

		} else {
			std::string pipeCommand = SSOptions::getOptions().getPipeCommand();
			if (!pipeCommand.empty()) {
				ops = new OPipestream(pipeCommand);
			}

			if (!outputFile.empty()) {
				out1 = out1f = new std::ofstream(outputFile.c_str());
			}

			assert(out1 != NULL || ops != NULL);
			outputRegularFiles(inputs, (ops == NULL ? *out1 : *ops));

			if (ops != NULL) {
				delete ops;
				ops = NULL;
			}
			if (out1f != NULL) {
				out1f->close();
				out1f = NULL;
			}

		}

		if (! SSOptions::getOptions().getFifoFile()) {
			forkPids = forkCommand();
		}

		if (!forkPids.empty()) {
			for(int i = 0; i < (int) forkPids.size(); i++) {
				LOG_DEBUG_OPTIONAL(1, true, "Waiting for forkCommand to finish: " << forkPids[i]);
				int status = Fork::wait(forkPids[i]);
				if (status != 0) {
					exitStatus += status;
					LOG_THROW("Child process " << forkPids[i] << " errored.  Aborting");
				}
			}
		}

	} catch (...) {
		LOG_WARN(1, "Cleaning up after exception");
		if (ops != NULL) {
			delete ops;
		}
		if (out1f != NULL) {
			out1f->close();
		}
		Cleanup::cleanup();
		sleep(1); // attempt to let other processes get the signal
#ifdef ENABLE_MPI
		MPI_Abort(world, 1);
#endif
		LOG_THROW("Found a problem...");
	}

	Cleanup::cleanup();

#ifdef ENABLE_MPI
	niceBarrier(world);
	OptionsBaseInterface::StringListType merges = SSOptions::getOptions().getMergeList();
	for(unsigned int i = 0; i < merges.size(); i++) {
		std::string &merge = merges[i];
		size_t pos = merge.find(' ');
		if (pos == std::string::npos) {
			LOG_WARN(1, "Could not parse merge string: " << merge);
		} else {
			LOG_DEBUG(2, "'" << merge << "' pos=" << pos);
			std::string myFile = merge.substr(0, pos);
			std::string mergeFile = merge.substr(pos+1);
			LOG_DEBUG(2, "'" << merge << "' pos=" << pos << " '" << myFile.c_str() << "' '" << mergeFile.c_str() << "'");
			LOG_DEBUG(2, "Merging '" << myFile << "' into '" << mergeFile << "'");
			DistributedOfstreamMap::mergeFiles(world, myFile, mergeFile, true);
		}
	}

	MPI_Finalize();
#endif

	return exitStatus;
}
