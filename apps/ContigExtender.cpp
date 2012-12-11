//
// Kmernator/apps/ContigExtender.cpp
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

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "ReadSet.h"
#include "Options.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "FilterKnownOddities.h"
#include "DuplicateFragmentFilter.h"
#include "ContigExtender.h"
#include "Log.h"

using namespace std;
typedef TrackingDataMinimal4f DataType;
typedef KmerSpectrum<DataType, DataType> KS;

class _ContigExtenderOptions : public OptionsBaseInterface {
public:
	virtual ~_ContigExtenderOptions() {}
	void _resetDefaults() {
		ContigExtenderBaseOptions::_resetDefaults();
		GeneralOptions::_resetDefaults();
		FilterKnownOdditiesOptions::_resetDefaults();
		DuplicateFragmentFilterOptions::_resetDefaults();

		// do not apply artifact filtering by default
		FilterKnownOdditiesOptions::getOptions().getSkipArtifactFilter() = 1;
		// override the default output format!
		GeneralOptions::getOptions().getFormatOutput() = 3;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		ContigExtenderBaseOptions::_setOptions(desc, p);
		GeneralOptions::_setOptions(desc,p);
		FilterKnownOdditiesOptions::_setOptions(desc,p);
		DuplicateFragmentFilterOptions::_setOptions(desc,p);
	}
	// use to post-process options, returning true if everything is okay
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= FilterKnownOdditiesOptions::_parseOptions(vm);
		ret &= DuplicateFragmentFilterOptions::_parseOptions(vm);
		ret &= ContigExtenderBaseOptions::_parseOptions(vm);
		return ret;
	}
};
typedef OptionsBaseTemplate< _ContigExtenderOptions > ContigExtenderOptions;

int main(int argc, char *argv[]) {

	ContigExtenderOptions::parseOpts(argc, argv);

	OptionsBaseInterface::FileListType &inputFiles = Options::getOptions().getInputFiles();
	OptionsBaseInterface::FileListType contigFiles;
	contigFiles.push_back(ContigExtenderBaseOptions::getOptions().getContigFile());

	ReadSet reads;
	LOG_VERBOSE(1, "Reading Input Files" );
	reads.appendAllFiles(inputFiles);

	LOG_VERBOSE(1, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount() << " Bases ");

	ReadSet contigs;
	LOG_VERBOSE(1, "Reading Contig File" );
	contigs.appendAllFiles(contigFiles);
	LOG_VERBOSE(1, "loaded " << contigs.getSize() << " Reads, " << contigs.getBaseCount() << " Bases ");

	if (DuplicateFragmentFilterOptions::getOptions().getDeDupMode() > 0 && DuplicateFragmentFilterOptions::getOptions().getDeDupEditDistance() >= 0) {
		LOG_VERBOSE(2, "Applying DuplicateFragmentPair Filter to Input Files");
		unsigned long duplicateFragments = DuplicateFragmentFilter::filterDuplicateFragments(reads);
		LOG_VERBOSE(1, "filter removed duplicate fragment pair reads: " << duplicateFragments);
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
	}

	ReadSet newContigs = ContigExtender<KS>::extendContigs(contigs, reads);

	string outputFilename = Options::getOptions().getOutputFile();
	OfstreamMap ofmap(outputFilename,"");
	Options::getOptions().getFormatOutput() = FormatOutput::FASTA_UNMASKED;
	if (!outputFilename.empty()) {
		for(unsigned long i = 0; i < newContigs.getSize(); i++)
			newContigs.getRead(i).write(ofmap.getOfstream(""));
	}

	LOG_VERBOSE(1, "Finished");
}

