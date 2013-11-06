//
// Kmernator/apps/FixPair.cpp
//
// Author: Rob Egan
//
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

#undef _USE_OPENMP

#include <boost/thread.hpp>

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/ref.hpp>
#include <boost/container/vector.hpp>

#include "config.h"
#include "mpi.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "Options.h"
#include "Log.h"
#include "Utils.h"
#include "MPIUtils.h"
#include "BroadcastOstream.h"
#include "DistributedOfstreamMap.h"
#include <unistd.h>

using namespace std;

class _SSOptions : public OptionsBaseInterface {
public:
	typedef OptionsBaseInterface::StringListType StringListType;

	int &getNumFiles() {
		return numFiles;
	}
	int &getFileNum() {
		return fileNum;
	}
	int getFirstDim() {
		if (secondDim <= 1)
			return numFiles;
		else
			return numFiles / secondDim;
	}
	int getSecondDim() {
		if (secondDim <= 1)
			return 1;
		else
			return secondDim;
	}
	int getFirst() {
		if (secondDim <= 0)
			return fileNum;
		else
			return fileNum / secondDim;
	}
	int getSecond() {
		if (secondDim <= 0)
			return 1;
		else
			return fileNum % secondDim;
	}
	std::string &getSplitFile() {
		return splitFile;
	}
	int &getEvenChunks() {
		return evenChunks;
	}
	int &getMinReadLength() {
		return minReadLength;
	}
	bool &getFifoFile() {
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
	bool &getTrimPairInName() {
		return trimPairInName;
	}

	_SSOptions() : splitFile(), pipeCommand(), mergeList(), forkCommand(), extraFifo(), numFiles(-1), fileNum(-1), secondDim(0), evenChunks(1), minReadLength(0), fifoFile(false), trimPairInName(false) {}

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

		if (secondDim > 0) {
			while ((pos = input.find("{FirstNum}")) != std::string::npos) {
				sprintf(buf, "%06d", getFirst());
				std::string buf2(buf);
				LOG_DEBUG_OPTIONAL(2, true, "Replacing {FirstNum} (" << buf2 << ") in " << input);
				input.replace(pos, 10, buf2);
			}
			while ((pos = input.find("{SecondNum}")) != std::string::npos) {
				sprintf(buf, "%06d", getSecond());
				std::string buf2(buf);
				LOG_DEBUG_OPTIONAL(2, true, "Replacing {SecondNum} (" << buf2 << ") in " << input);
				input.replace(pos, 11, buf2);
			}
			while ((pos = input.find("{UniqFirst}")) != std::string::npos) {
				sprintf(buf, "%06d", getFirst());
				std::string buf2 = std::string(buf) + "of";
				sprintf(buf, "%06d", getFirstDim());
				buf2 += std::string(buf);
				LOG_DEBUG_OPTIONAL(2, true, "Replacing {UniqFirst} (" << buf2 << ") in " << input);
				input.replace(pos, 11, buf2);
			}
			while ((pos = input.find("{UniqSecond}")) != std::string::npos) {
				sprintf(buf, "%06d", getSecond());
				std::string buf2 = std::string(buf) + "of";
				sprintf(buf, "%06d", getSecondDim());
				buf2 += std::string(buf);
				LOG_DEBUG_OPTIONAL(2, true, "Replacing {UniqSecond} (" << buf2 << ") in " << input);
				input.replace(pos, 12, buf2);
			}
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

	void setFileDimensions(const MPI_Comm &comm) {
		if (numFiles == -1)
			MPI_Comm_size(comm, &numFiles);
		if (fileNum == -1)
			MPI_Comm_rank(comm, &fileNum);

		if (fileNum >= numFiles)
			LOG_THROW("You can not specify --file-num ge --num-files (" << fileNum << " " << numFiles);

		getPipeCommand() = _replaceWithKeys(getPipeCommand());
		getForkCommand() = _replaceWithKeys(getForkCommand());
		getExtraFifo() = _replaceWithKeys(getExtraFifo());
		getMergeList() = _replaceWithKeys(getMergeList());
		getSplitFile() = _replaceWithKeys(getSplitFile());
		GeneralOptions::getOptions().getOutputFile() = _replaceWithKeys(GeneralOptions::getOptions().getOutputFile());
		GeneralOptions::getOptions().getInputFiles() = _replaceWithKeys(GeneralOptions::getOptions().getInputFiles());

		if (secondDim > 0 && numFiles % secondDim != 0) {
			LOG_THROW("if you use --second-dim, it must be a factor of the number of files: " << numFiles << " which " << secondDim << " is not.");
		}

	}

	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		GeneralOptions::getOptions().getKeepReadComment() = false;

		Options::getOptions().getVerbose() = 0;
		Options::getOptions().getDebug() = 0;
		Options::getOptions().getMmapInput() = false;
		Options::getOptions().getMaxThreads() = 1;

	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("SplitSequenceOnTheFly <options> [input-file ...]\n\tNote --input-file can either be specified as a positional argument or within the options\n\nSplit Sequence Options");
		opts.add_options()
								("num-files", po::value<int>()->default_value(numFiles), "The number of files to split into N (automatically set under MPI)")
								("file-num",  po::value<int>()->default_value(fileNum), "The number of the file to output (0-(N-1)) (automatically set under MPI)")
								("second-dim", po::value<int>()->default_value(secondDim), "if > 0, then ranks will be divided by secondDim, and secondDim iterations of splitting will be performed, broadcasting secondDim chunks (of NumFiles) of the split data to each subrank group.  Use the keywords '{FirstNum}' '{SecondNum}' and '{UniqFirst}' and '{UniqSecond}'")
								("even-chunks", po::value<int>()->default_value(evenChunks), "if > 1 then the output of each partition will be spread out across the file (recommend 10 if the ordering is less important than predictable runtime of forked commands)")
								("pipe-command", po::value<std::string>(), "a command to pipe the portion of the file(s) into.  Use the keyword variables '{Uniq}', '{FileNum}' and '{NumFiles}' to replace with MPI derived values")
								("merge", po::value<StringListType>(), "two arguments.  First is per-mpi file (use keywords) second is final file; can be specified multiple times")
								("split-file", po::value<std::string>()->default_value(splitFile), "if set, paired (second) reads will be sent to this file.  first reads will go to --output-file")
								("output-fifo", po::value<bool>()->default_value(fifoFile), "if set, --output-file (and --split-file) will be fifo files which will be created on start and deleted on exit")
								("extra-fifo", po::value<StringListType>(), "optional additional fifo file(s) (to potentially use with --fork-command)")
								("fork-command", po::value<StringListType>(), "optional command(s) to fork between opening output files and wait to finish before starting any merges")
								("min-read-length", po::value<int>()->default_value(minReadLength), "minimum read length to output (applies to unpaired only)")
								("trim-pair-in-name", po::value<bool>()->default_value(trimPairInName), "remove /1 or /2 from the name.  Useful for some aligners like bowtie2")
								;

		desc.add(opts);
		GeneralOptions::_setOptions(desc, p);
		p.add("input-file", -1);

	}
	bool _parseOptions(po::variables_map &vm) {
		// set options specific to this program
		bool ret = GeneralOptions::_parseOptions(vm);
		setOpt("num-files", getNumFiles());
		setOpt("file-num", getFileNum());
		setOpt("second-dim", secondDim);
		setOpt("even-chunks", getEvenChunks());
		setOpt("min-read-length", getMinReadLength());
		setOpt("split-file", getSplitFile());
		setOpt("pipe-command", getPipeCommand());
		setOpt("output-fifo", getFifoFile());
		setOpt2("merge", getMergeList());
		setOpt2("fork-command", getForkCommand());
		setOpt2("extra-fifo", getExtraFifo());
		setOpt("trim-pair-in-name", getTrimPairInName());


		if (Options::getOptions().getInputFiles().empty() || getNumFiles() == 0) {
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
	int numFiles, fileNum, secondDim, evenChunks, minReadLength;
	bool fifoFile, trimPairInName;
};
typedef OptionsBaseTemplate< _SSOptions > SSOptions;
typedef OptionsBaseInterface::StringListType StringListType;

void trimName(std::string &name) {
	size_t pos = name.find_last_of("/");
	LOG_DEBUG(5, "trimName: " << name << " " << pos);
	if (pos != std::string::npos)
		name = name.substr(0, pos);
}

template<typename T>
void outputRegularFiles(OptionsBaseInterface::FileListType inputs, T &output1, int fileNum, int numFiles) {

	FormatOutput outputFormat = FormatOutput::getDefault();
	int chunks = SSOptions::getOptions().getEvenChunks();
	int minReadLength = SSOptions::getOptions().getMinReadLength();
	bool trimPairInName = SSOptions::getOptions().getTrimPairInName();

	for (unsigned int i = 0 ; i < inputs.size(); i++) {
		ReadFileReader reader(inputs[i], "");
		for(int j = 0 ; j < chunks; j++) {
			reader.seekToPartition( fileNum+j*chunks, numFiles*chunks );
			std::string name, bases, quals, comment;
			while (reader.nextRead(name, bases, quals, comment)) {
				if ((int) bases.length() >= minReadLength) {
					if (trimPairInName)
						trimName(name);
					output1 << outputFormat.toString(name, bases, quals, comment);
				}
			}
		}
	}
}

template<typename T>
class OutputSplitFiles {
private:
	const OptionsBaseInterface::FileListType &inputs;
	int readNum;
	T &output;
	int first, secondDim, numFiles;
	unsigned long &numReads;
	OutputSplitFiles() { throw; }
	OutputSplitFiles(const OutputSplitFiles& copy) { throw; }
	OutputSplitFiles &operator=(const OutputSplitFiles &copy) { throw; }
public:
	typedef boost::shared_ptr< Kmernator::MmapSource > MmapPtr;
	typedef boost::shared_ptr< ReadFileReader > RFRPtr;
	OutputSplitFiles(const OptionsBaseInterface::FileListType &_inputs, int _readNum, T &_output, int _first, int _secondDim, int _numFiles, unsigned long &_numReads)
	: inputs(_inputs), readNum(_readNum), output(_output), first(_first), secondDim(_secondDim), numFiles(_numFiles), numReads(_numReads) {
		if (readNum > 1)
			LOG_THROW("Invalid readNum: " << readNum);
		LOG_DEBUG_OPTIONAL(2, true, "OutputSplitFiles(): " << _readNum);
	}
	void operator()() {
		numReads = 0;
		if (secondDim > 1) {
			for(int secondIt = 0 ; secondIt < secondDim ; secondIt++) {
				LOG_DEBUG_OPTIONAL(1, true, "OutputSplitFiles()(): readNum: " << readNum << " second: " << secondIt << " mydim:" << first << ", " << secondIt << " : " << secondDim << " part " << first*secondDim + secondIt);
				numReads += outputSplitFiles(inputs, readNum, output, first*secondDim + secondIt, numFiles);
			}
		} else {
			numReads = outputSplitFiles(inputs, readNum, output, first, numFiles);
		}
		output.close();
	}
	unsigned long getNumReads() {
		return numReads;
	}
protected:

	unsigned long outputSplitFiles(const OptionsBaseInterface::FileListType &inputs, int readNum, T &output, int fileNum, int numFiles) {

		int numInputs = inputs.size();
		std::vector< RFRPtr > readers;

		readers.reserve(numInputs);
		int chunks = SSOptions::getOptions().getEvenChunks();

		for(int i=0; i < numInputs; i++) {
			readers.push_back(RFRPtr( new ReadFileReader(inputs[i], true)) );
		}

		assert(numInputs == (int)readers.size());

		FormatOutput outputFormat = FormatOutput::getDefault();
		unsigned long numReads = 0;

		_writeSplitReads(readers, readNum, chunks, fileNum, numFiles, numReads, output, outputFormat, inputs);

		LOG_DEBUG_OPTIONAL(1, true, "Finished splitFiles parallel " << numReads);

		return numReads;

	}
	void _writeSplitReads(std::vector<RFRPtr> & readers, unsigned long readNumOfPair, int chunks, int fileNum, int numFiles, unsigned long  & numReads, T & output, FormatOutput & outputFormat, const OptionsBaseInterface::FileListType & inputs)
	{
		bool trimPairInName = SSOptions::getOptions().getTrimPairInName();
		std::string name, bases, quals, comment;
		for(int i=0; i < (int)readers.size(); i++) {
			for(int j = 0 ; j < chunks; j++) {
				readers[i]->seekToPartition( fileNum+j*chunks, numFiles*chunks );

				while (readers[i]->nextRead(name, bases, quals, comment)) {
					if ((numReads & 0x01) == readNumOfPair) {
						if (trimPairInName)
							trimName(name);
						output << outputFormat.toString(name, bases, quals, comment);
					}
					numReads++;
					LOG_DEBUG_OPTIONAL(2, ((numReads & 0xfffff) == 0), "reading and writing split file " << i << " " << inputs[i] << " for read " << numReads/2);
				}
			}
			LOG_DEBUG_OPTIONAL(2, true, "Finished reading and writing split file " << i << " " << inputs[i] << " for read. " << numReads << " " << readers[i]->getPos());
		}
		LOG_DEBUG_OPTIONAL(1, true, "Finished reading and writing split files for read. " << numReads);
	}
};

class ForkCommandThread {
public:
	typedef boost::shared_ptr< ForkCommandThread > Ptr;
	typedef std::vector< Ptr > Vector;
	ForkCommandThread(std::string command) : _command(command), _status(0) {
		_thread = boost::thread(*this);
	}
	~ForkCommandThread() {
	}
	ForkCommandThread(ForkCommandThread &move) : _command(move._command), _status(move._status), _thread(boost::move(move._thread)) {}
	void operator()() {
		LOG_DEBUG_OPTIONAL(1, true, "ForkCommandThread(" << _command << ")::() Starting");
		_status = ForkDaemon::system(_command);
		LOG_DEBUG_OPTIONAL(1, true, "ForkCommandThread(" << _command << ")::() Finished: " << _status);
	}
	std::string toString() {
		return "ForkCommandThread(" + _command + "): " + boost::lexical_cast<std::string>(_status);
	}
	int getStatus() {
		return _status;
	}
	boost::thread &getThread() {
		return _thread;
	}
	int join() {
		_thread.join();
		return _status;
	}
	static int join(Vector v) {
		int status = 0;
		LOG_DEBUG_OPTIONAL(1, true, "Waiting for all " << v.size() << " threads to join");
		for(Vector::iterator it = v.begin(); it != v.end(); it++) {
			status += (*it)->join();
		}
		return status;
	}
private:
	std::string _command;
	int _status;
	boost::thread _thread;
};

ForkCommandThread::Vector forkCommand() {
	LOG_DEBUG_OPTIONAL(1, true, "forkCommand()");
	ForkCommandThread::Vector forks;
	StringListType forkCommand = SSOptions::getOptions().getForkCommand();

	forks.reserve(forkCommand.size());
	for(int i = 0; i < (int) forkCommand.size(); i++) {
		LOG_DEBUG_OPTIONAL(2, true, "forkCommand(): Starting " << forkCommand[i]);

		ForkCommandThread::Ptr ptr(new ForkCommandThread(forkCommand[i])); 
		forks.push_back(ptr);

	}
	return forks;
}

/*
std::vector< int > _forkCommand() {
	LOG_DEBUG_OPTIONAL(1, true, "forkCommand()");
	std::vector< int > forks;
	StringListType forkCommand = SSOptions::getOptions().getForkCommand();
	if (forkCommand.empty())
		return forks;

	for(int i = 0; i < (int) forkCommand.size(); i++) {
		LOG_DEBUG_OPTIONAL(1, true, "Starting " << forkCommand[i]);

		int child = Fork::forkCommand(forkCommand[i]);
		forks.push_back(child);

	}
	return forks;
}
*/

int main(int argc, char *argv[]) {

	ForkDaemon::initialize();

	int exitStatus = 0;

	ScopedMPIComm< SSOptions > world(argc, argv);

	Cleanup::prepare();

	SSOptions::getOptions().setFileDimensions(world);
	OptionsBaseInterface::FileListType inputs = Options::getOptions().getInputFiles();

	std::string outputFile = Options::getOptions().getOutputFile();
	std::string splitFile = SSOptions::getOptions().getSplitFile();

	OPipestream *ops = NULL;
	ostream *out1 = &std::cout;
	ofstream *out1f = NULL;
	ofstream *out2f = NULL;
	ForkCommandThread::Vector forkThreads;
	StringListType extraFifos = SSOptions::getOptions().getExtraFifo();

	try {
		int fileNum = SSOptions::getOptions().getFileNum();
		int numFiles = SSOptions::getOptions().getNumFiles();

		int second = SSOptions::getOptions().getSecond();
		int first = SSOptions::getOptions().getFirst();
		int secondDim = SSOptions::getOptions().getSecondDim();
		//int firstDim  = SSOptions::getOptions().getFirstDim();


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

			forkThreads = forkCommand();

		}

		if (!outputFile.empty()) {
			LOG_VERBOSE_OPTIONAL(1, true, "writing to output file: " << outputFile);
			out1 = out1f = new std::ofstream(outputFile.c_str());
		}

		if (!splitFile.empty()) {
			assert(out1f != NULL);
			LOG_VERBOSE_OPTIONAL(1, true, "writing to second output (split) file: " << splitFile);
			out2f = new std::ofstream(splitFile.c_str());
			unsigned long numReads1 = 0, numReads2 = 0;
			OutputSplitFiles<std::ofstream> osf1(inputs, 0, *out1f, first, secondDim, numFiles, numReads1);
			OutputSplitFiles<std::ofstream> osf2(inputs, 1, *out2f, first, secondDim, numFiles, numReads2);
			boost::thread t1(boost::ref(osf1)), t2(boost::ref(osf2));
			t1.join();
			t2.join();
			LOG_DEBUG_OPTIONAL(1, true, "Split " << numReads1 << " and " << numReads2 << " reads");
		} else {
			std::string pipeCommand = SSOptions::getOptions().getPipeCommand();
			if (!pipeCommand.empty()) {
				LOG_VERBOSE_OPTIONAL(1, true, "Outputting to pipeCommand: " << pipeCommand);
				ops = new OPipestream(pipeCommand);
			}

			assert(out1 != NULL || ops != NULL);
			if (secondDim > 0) {
				mpi::communicator local = world.split(first);
				int root = 0;
				BroadcastOstream bcastos(local, root, (ops == NULL ? *out1 : *ops));
				for(int secondIt = 0 ; secondIt < secondDim ; secondIt++) {
					LOG_DEBUG_OPTIONAL(1, true, "second: " << secondIt << " mydim:" << first << ", " << second << " : " << secondDim << " part " << first*secondDim + second);
					if (local.rank() == root)
						outputRegularFiles(inputs, bcastos, first*secondDim + secondIt, numFiles);
				}
			} else {
				outputRegularFiles(inputs, (ops == NULL) ? *out1 : *ops, fileNum, numFiles);
			}

		}

		LOG_DEBUG_OPTIONAL(1, true, "Closing output(s)");
		if (ops != NULL) {
			delete ops;
			ops = NULL;
		}

		if (out1f != NULL) {
			if (out1f->good())
				out1f->close();
			out1f = NULL;
		}
		if (out2f != NULL) {
			if (out2f->good())
				out2f->close();
			out2f = NULL;
		}

		if (! SSOptions::getOptions().getFifoFile()) {
			forkThreads = forkCommand();
		}

		if (!forkThreads.empty()) {
			for(int i = 0; i < (int) forkThreads.size(); i++) {
				LOG_DEBUG_OPTIONAL(1, true, "Waiting for forkCommand to finish: " << forkThreads[i]->toString());
				int status = forkThreads[i]->join();
				if (status != 0) {
					exitStatus += status;
					LOG_THROW("Child process " << forkThreads[i]->toString() << " errored.  Aborting");
				}
			}
		}

		ForkDaemon::finalize();

	} catch (...) {
		LOG_WARN(1, "Cleaning up after exception");
		if (ops != NULL) {
			delete ops;
		}
		if (out1f != NULL) {
			out1f->close();
		}
		Cleanup::cleanup(-1);
		sleep(1); // attempt to let other processes get the signal
		if (MPI::Is_thread_main())
			world.abort(1);
		LOG_THROW("Found a problem... should not get here");
	}

	Cleanup::cleanup();

	MPIUtils::niceBarrier(world);
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

	return exitStatus;
}
