/*
 * MeraculousCounter-P.cpp
 *
 *  Created on: Sep 20, 2011
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

#include "config.h"
#include "Options.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "KmerTrackingData.h"
#include "Utils.h"
#include "Log.h"
#include "DistributedFunctions.h"
#include "Meraculous.h"

class _MeraculousCounterOptions : public OptionsBaseInterface {
public:
	// use to set/overrided any defaults on options that are stored persistently
	void _resetDefaults() {
		MeraculousOptions::_resetDefaults();
		MPIOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
		KmerBaseOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();

		GeneralOptions::getOptions().getMmapInput() = false;
		GeneralOptions::getOptions().getVerbose() = 2;
		GeneralOptions::getOptions().getMinQuality() = 2;

		KmerSpectrumOptions::getOptions().getMinKmerQuality() = 0;
		KmerSpectrumOptions::getOptions().getSaveKmerMmap() = 0;
	}
	// use to set the description of all options
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {

		MeraculousOptions::_setOptions(desc, p);
		MPIOptions::_setOptions(desc,p);
		KmerBaseOptions::_setOptions(desc,p);
		KmerSpectrumOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc, p);

	}
	// use to post-process options, returning true if everything is okay
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= MeraculousOptions::_parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);
		ret &= KmerBaseOptions::_parseOptions(vm);
		ret &= KmerSpectrumOptions::_parseOptions(vm);
		if (KmerBaseOptions::getOptions().getKmerSize() == 0) {
			setOptionsErrorMsg("The Kmer size can not be 0");
		}
		if (GeneralOptions::getOptions().getInputFiles().empty()) {
			setOptionsErrorMsg("You must specify at least one input file");
		}
		return ret;
	}
};
typedef OptionsBaseTemplate< _MeraculousCounterOptions > MeraculousCounterOptions;
typedef MeraculousDistributedKmerSpectrum KS;

int main(int argc, char *argv[]) {

	ScopedMPIComm< MeraculousCounterOptions > world(argc, argv);

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	ReadSet reads;

	try {
		OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files");

		// TODO save memory! read file and build spectrum in 100MB chunks
		reads.appendAllFiles(inputs, world.rank(), world.size());

		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		setGlobalReadSetOffsets(world, reads);
		long numBuckets = 0;
		numBuckets = KS::estimateWeakKmerBucketSize(reads);

		numBuckets = all_reduce(world, numBuckets, mpi::maximum<int>());
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "targeting " << numBuckets << " buckets for reads");

		KS spectrum(world, numBuckets);

		spectrum.buildKmerSpectrum(reads);
		if (Log::isVerbose(1)) {
			std::string hist = spectrum.getHistogram(false);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Kmer Histogram\n" << hist);
		}
		std::string outputFilenameBase = outputFilename + ".mercount.m" + boost::lexical_cast<std::string>(KmerSizer::getSequenceLength());
		spectrum.dumpCounts(outputFilenameBase);
		outputFilenameBase = outputFilename + ".mergraph.m" + boost::lexical_cast<std::string>(KmerSizer::getSequenceLength()) + ".D" + boost::lexical_cast<std::string>(KmerSpectrumOptions::getOptions().getMinDepth());
		spectrum.dumpGraphs(outputFilenameBase);
	} catch (...) {
		LOG_ERROR(1, "caught an error!" << StackTrace::getStackTrace());
	}
	world.barrier();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");

	return 0;

}

