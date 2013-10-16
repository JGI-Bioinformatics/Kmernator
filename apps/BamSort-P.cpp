/*
 * BamSort-P.cpp
 *
 *  Created on: Oct 11, 2012
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

#include <mpi.h>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

#include <string>

#include "SamUtils.h"
#include "Log.h"
#include "Options.h"
#include "MPIUtils.h"

class _BamSortOptions : public OptionsBaseInterface {
public:
	_BamSortOptions() : unmappedReads(), unmappedReadPairs(), keepUnmappedPairedRead(true), numPartitions(0) {
	}
	virtual ~_BamSortOptions() {}
	std::string &getUnmappedReads() {
		return unmappedReads;
	}
	std::string &getUnmappedReadPairs() {
		return unmappedReadPairs;
	}
	bool &getKeepUnmappedPairedRead() {
		return keepUnmappedPairedRead;
	}
	std::string &getOutputBam() {
		return outputBam;
	}
	FileListType &getInputBams() {
		return inputBams;
	}
	int &getNumPartitions() {
		return numPartitions;
	}
	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		GeneralOptions::getOptions().getDebug() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		p.add("output-bam", 1);
		p.add("input-bams", -1);

		po::options_description opts("BamSort-P <options> output.bam input.bam [...]\n\nOptions");
		opts.add_options()
				("output-bam", po::value<std::string>())
				("input-bams", po::value<FileListType>())
				("unmapped-read-pairs", po::value<std::string>()->default_value(unmappedReadPairs), "gzipped file to place unmapped read Pairs Fastqs (can be same as --unmapped-reads)")
				("unmapped-reads", po::value<std::string>()->default_value(unmappedReads), "gzipped file to place unmapped reads Fastqs (can be same as --unmapped-read-pairs)")
				("keep-unmapped-paired-read", po::value<bool>()->default_value(keepUnmappedPairedRead), "if unmapped-read-pairs file is specified, keep an unmapped read in the bam if its pair is mapped")
				("num-partitions", po::value<int>()->default_value(numPartitions), "The number of alignment-index partitions to merge. Input bams expected to come ordered grouped by in batches of num-partitions where each group has the exact same read counts in the exact same order");
		desc.add(opts);

		GeneralOptions::_setOptions(desc, p);
		// Other *::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		ret |= GeneralOptions::_parseOptions(vm);

		setOpt("unmapped-read-pairs", unmappedReadPairs);
		setOpt("unmapped-reads", unmappedReads);
		setOpt("keep-unmapped-paired-read", keepUnmappedPairedRead);
		setOpt("output-bam", outputBam);
		setOpt2("input-bams", inputBams);
		setOpt("num-partitions", numPartitions);

		if (outputBam.empty())
			setOptionsErrorMsg("You must specify at least the outputBam");

		if (inputBams.empty())
			setOptionsErrorMsg("You must specify at least one input bam");

		if (unmappedReads.empty() && !unmappedReadPairs.empty()) {
			unmappedReads = unmappedReadPairs;
		}
		if (unmappedReadPairs.empty() && !unmappedReads.empty()) {
			unmappedReadPairs = unmappedReads;
		}

		// Other ret &= *::_parseOptions(vm);
		return ret;
	}
private:
	std::string outputBam;
	FileListType inputBams;
	std::string unmappedReads;
	std::string unmappedReadPairs;
	bool keepUnmappedPairedRead;
	int numPartitions;
};
typedef OptionsBaseTemplate< _BamSortOptions > BamSortOptions;

int main(int argc, char **argv)
{
	ScopedMPIComm< BamSortOptions > world(argc, argv);

	Cleanup::prepare();

	try {

		int rank = world.rank(), size = world.size();

		BamVector reads;
		BamHeaderPtr header;

		std::string outputBam = BamSortOptions::getOptions().getOutputBam();
		std::string outputBamTmp = outputBam + ".tmp";
		OptionsBaseInterface::FileListType inputBams = BamSortOptions::getOptions().getInputBams();
		int partitions = BamSortOptions::getOptions().getNumPartitions();

		if (partitions > 1) {
			if (size != (int) inputBams.size()) {
				if (rank == 0)
					std::cerr << "The number of files is a mismatch to the job size.  When merging partitions, this is necessary." << BamSortOptions::getDesc() << std::endl;
				exit(1);
			}
			if (size % partitions != 0) {
				if (rank == 0)
					std::cerr << "The partitions " << partitions << " is not a factor of the size " << size << BamSortOptions::getDesc() << std::endl;
				exit(1);
			}
		}

		if (world.rank() == 0)
			unlink(outputBamTmp.c_str());
		world.barrier();
	
		LOG_VERBOSE_GATHER(1, "Reading input files");
		if (partitions > 1) {
			std::string myInputFile = inputBams[rank];
			int color = rank / partitions;
			mpi::communicator partitionWorld = world.split(color);
			LOG_VERBOSE_GATHER(1, "Input: " << myInputFile << " split color: " << color << " rank: " << partitionWorld.rank() << " of " << partitionWorld.size());
			SamUtils::MPIMergeSam mergeSam(partitionWorld, myInputFile, reads);
	
			if (color == 0)
				mergeSam.outputMergedHeader(outputBamTmp);
	
			world.barrier();
	
		} else {
			header = BamStreamUtils::readBamFile(world, inputBams, reads);
		}
	
		bool needsCollapse = false;
		bool keepUnmappedPairedRead = BamSortOptions::getOptions().getKeepUnmappedPairedRead();
		std::string unmappedReadPairFile = BamSortOptions::getOptions().getUnmappedReadPairs();
		std::string unmappedReadsFile  = BamSortOptions::getOptions().getUnmappedReads();
		BamVector unmappedReadSingles, unmappedReadPairs, unmappedPairedReads;
		if (!unmappedReadPairFile.empty() || !unmappedReadsFile.empty()) {

			std::string tmpFile(unmappedReadPairFile + ".tmp");
			std::string tmpFile2(unmappedReadsFile + ".tmp");
			if (world.rank() == 0) {
				unlink(tmpFile.c_str());
				unlink(tmpFile2.c_str());
			}
			world.barrier();

			SamUtils::splitUnmapped(reads, unmappedReadSingles, unmappedReadPairs, unmappedPairedReads, keepUnmappedPairedRead);
			if (!unmappedReadSingles.empty() || !unmappedReadPairs.empty() || (!keepUnmappedPairedRead && !unmappedPairedReads.empty())) {
				needsCollapse = true;
			}
			LOG_VERBOSE(1, "Purging unmapped read pairs: " << unmappedReadPairs.size());
			if (!unmappedReadPairFile.empty() && unmappedReadPairFile.compare("/dev/null") != 0) {
				LOG_VERBOSE(1, "Writing unmappedReadPairs " << unmappedReadPairs.size() << " to " << unmappedReadPairFile << " (" << tmpFile << ")");

				ScopedMPIFile pairedReads(world, tmpFile);
				SamUtils::writeFastqGz(world, unmappedReadPairs, pairedReads, true);
				world.barrier();
			}
			assert(unmappedReadPairs.empty());

			LOG_VERBOSE(1, "Purging unmapped read singles: " << unmappedPairedReads.size() << " and " << unmappedReadSingles.size());

			if (!unmappedReadsFile.empty() && unmappedReadsFile.compare("/dev/null") != 0) {
				LOG_VERBOSE(1, "Writing unmappedPairedReads " << unmappedPairedReads.size() << " to " << unmappedReadsFile << " (" << tmpFile2 << ")");

				ScopedMPIFile singleReads(world, tmpFile2);
				SamUtils::writeFastqGz(world, unmappedPairedReads, singleReads, !keepUnmappedPairedRead);

				LOG_VERBOSE(1, "Writing unmappedReadSingles " << unmappedReadSingles.size() << " to " << unmappedReadsFile << " (" << tmpFile2 << ")");

				SamUtils::writeFastqGz(world, unmappedReadSingles, singleReads, true);
				world.barrier();
			}
			assert(unmappedReadSingles.empty());
			assert(keepUnmappedPairedRead || unmappedPairedReads.empty());

			if (world.rank() == 0) {
				unlink(unmappedReadsFile.c_str());
				unlink(unmappedReadPairFile.c_str());
				rename(tmpFile.c_str(), unmappedReadPairFile.c_str());
				if (tmpFile.compare(tmpFile2) != 0) {
					rename(tmpFile2.c_str(), unmappedReadsFile.c_str());
				}
			}

		}

		if (needsCollapse) {
			LOG_DEBUG(1, "Collapsing read vector");
			long removed = SamUtils::collapseVector(reads);
			LOG_DEBUG(1, "Collapsed: " << removed);			
		}
		BamManager::destroyOrRecycleBamVector(unmappedReadSingles);
		BamManager::destroyOrRecycleBamVector(unmappedReadPairs);
		if (keepUnmappedPairedRead) {
			unmappedPairedReads.clear();
		} else {
			BamManager::destroyOrRecycleBamVector(unmappedPairedReads);
		}
	
		{
			LOG_VERBOSE_GATHER(1, "Redistributing reads before the sort:" << reads.size());
			BamStreamUtils::distributeReadsFinal(world, reads);
	
			LOG_VERBOSE_GATHER(1, "Sorting myreads: " << reads.size());
			SamUtils::MPISortBam sortem(world, reads, outputBamTmp, header.get());
			world.barrier();
			if (world.rank() == 0) {
				unlink(outputBam.c_str());
				rename(outputBamTmp.c_str(), outputBam.c_str());
			}
		}
	
		header.reset();
	
		BamManager::clearRecycledReads();
		LOG_VERBOSE_GATHER(1, "Finished");
	
		return 0;
	} catch (std::exception &e) {
		LOG_ERROR(1, "BamSort-P threw an exception! Aborting...\n\t" << e.what());
		world.abort(1);
	} catch (...) {
		LOG_ERROR(1, "BamSort-P threw an error! Aborting...\n");
		world.abort(1);
	}
}
