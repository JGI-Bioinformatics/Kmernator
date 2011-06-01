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

#include "config.h"
#include "Sequence.h"
#include "ReadSet.h"
#include "Options.h"
#include "Log.h"
#include "Utils.h"
#ifdef ENABLE_MPI
#include "DistributedFunctions.h"
#endif

using namespace std;

class SSOptions : public Options {
public:
	typedef Options::StringListType StringListType;

	static int &getDefaultNumFiles() {
		static int dummy = 0;
		return dummy;
	}
	static int &getDefaultFileNum() {
		static int dummy2 = 0;
		return dummy2;
	}
	static int getNumFiles() {
		return getVarMap()["num-files"].as<int> ();
	}
	static int getFileNum() {
		return getVarMap()["file-num"].as<int> ();
	}
	static std::string _replaceWithKeys(std::string input) {
		size_t pos;
		while ((pos = input.find("{FileNum}")) != std::string::npos) {
			LOG_DEBUG_OPTIONAL(1, true, "Replacing {FileNum} in " << input);
			input.replace(pos, 9, boost::lexical_cast<std::string>(getFileNum()));
		}
		while ((pos = input.find("{NumFiles}")) != std::string::npos) {
			LOG_DEBUG_OPTIONAL(1, true, "Replacing {NumFiles} in " << input);
			input.replace(pos, 10, boost::lexical_cast<std::string>(getNumFiles()));
		}
		LOG_DEBUG_OPTIONAL(1, true, "final string " << input);
		return input;
	}
	static std::string getPipeCommand() {
		std::string pipeCmd = getVarMap().count("pipe-command") ? getVarMap()["pipe-command"].as<std::string>() : std::string();
		return _replaceWithKeys(pipeCmd);
	}
	static StringListType getMergeList() {
		StringListType input = getVarMap().count("merge") ? getVarMap()["merge"].as<StringListType>() : StringListType();
		StringListType output;
		for(unsigned int i = 0; i < input.size(); i++) {
			output.push_back( _replaceWithKeys(input[i]));
		}
		return output;
	}
	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("input-file", -1);
		getDesc().add_options()("help", "produce help message")
				("num-files", po::value<int>()->default_value(getDefaultNumFiles()), "The number of files to split into N")
				("file-num",  po::value<int>()->default_value(getDefaultFileNum()), "The number of the file to ouput (0-(N-1))")
				("pipe-command", po::value<std::string>(), "a command to pipe the portion of the file(s) into.  Use the keyword variables '{FileNum}' and '{NumFiles}' to replace with MPI derived values")
				("merge", po::value<StringListType>(), "two arguments.  First is per-mpi file (use keywords) second is final file; can be specified multiple times")
				;

		bool ret = Options::parseOpts(argc, argv);
		if (getInputFiles().empty() || getNumFiles() == 0 || getFileNum() >= getNumFiles()) {
			ret = false;
			LOG_ERROR(1, "Please specify num-files, file-num and at least one input file.\nnum-files=" << getNumFiles() << "\nfile-num=" << getFileNum());
		}
		return ret;
	}
};


int main(int argc, char *argv[]) {
	Options::getVerbosity() = 0;
	Options::getMmapInput() = 0;

#ifdef ENABLE_MPI
	MPI_Init(&argc, &argv);
	mpi::environment env(argc, argv);
	mpi::communicator world;
	SSOptions::getDefaultNumFiles() = world.size();
	SSOptions::getDefaultFileNum() = world.rank();
	Logger::setWorld(&world);
#endif
	if (!SSOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	OPipestream *ops = NULL;
	std::string pipeCommand = SSOptions::getPipeCommand();
	if (!pipeCommand.empty()) {
		ops = new OPipestream(pipeCommand);
	}

	Options::FileListType inputs = Options::getInputFiles();

	ostream &output = (ops == NULL ? std::cout : *ops);
	for (unsigned int i = 0 ; i < inputs.size(); i++) {
		ReadFileReader reader(inputs[i], "");
		unsigned long lastPos = reader.seekToPartition( SSOptions::getFileNum(), SSOptions::getNumFiles() );
		std::string name, bases, quals;
	    while (reader.nextRead(name, bases, quals)) {
	        Read read(name, bases, quals);
	        read.write(output);
            if (reader.getPos() >= lastPos)
            	break;
	    }
	}

	if (ops != NULL) {
		delete ops;
	}

#ifdef ENABLE_MPI
	niceBarrier(world);
	Options::StringListType merges = SSOptions::getMergeList();
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

}

// $Log: FixPair.cpp,v $
// Revision 1.4  2010-05-18 20:50:18  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.3.2.1  2010-05-07 22:59:29  regan
// refactored base type declarations
//
// Revision 1.3  2010-05-06 21:46:57  regan
// merged changes from PerformanceTuning-20100501
//
//
