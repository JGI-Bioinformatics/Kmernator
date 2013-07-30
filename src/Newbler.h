/*
 * Newbler.h
 *
 *      Author: regan
 */
/*****************

Kmernator Copyright (c) 2013, The Regents of the University of California,
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

#ifndef NEWBLER_H
#define NEWBLER_H

#include "Options.h"
#include "Utils.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "Utils.h"

class _NewblerOptions : public OptionsBaseInterface {
public:
	_NewblerOptions() : newblerPath(), newblerOpts(), minOverlapLength(51), minOverlapIdentity(94), minLargeContigLength(151) {}
	virtual ~_NewblerOptions() {}
	std::string &getNewblerPath() {
		return newblerPath;
	}
	int &getMinOverlapLength() {
		return minOverlapLength;
	}
	int &getMinOverlapIdentity() {
		return minOverlapIdentity;
	}
	int &getMinLargeContigLength() {
		return minLargeContigLength;
	}
	string &getNewblerOpts() {
		return newblerOpts;
	}
	void _resetDefaults() {
		Options::getOptions().getMmapInput() = false;
	}

	void _setOptions(po::options_description &desc,	po::positional_options_description &p) {
		po::options_description opts("Newbler Options");

		opts.add_options()
				("newbler-path", po::value<std::string>()->default_value(newblerPath), "if set, ${newbler-path}/runAssembly will be used to extend contigs")
				("newbler-ml", po::value<int>()->default_value(minOverlapLength), "minimum overlap length")
				("newbler-mi", po::value<int>()->default_value(minOverlapIdentity), "minimum overlap identity")
				("newbler-l", po::value<int>()->default_value(minLargeContigLength), "minimum large contig length")
				("newbler-opts", po::value<std::string>()->default_value(newblerOpts), "extra options to send to newbler (use 'quotes')")
				;

		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		setOpt("newbler-path", newblerPath);
		setOpt("newbler-opts", newblerOpts);
		setOpt("newbler-ml", minOverlapLength);
		setOpt("newbler-mi", minOverlapIdentity);
		setOpt("newbler-l", minLargeContigLength);

		return ret;
	}
protected:
	std::string newblerPath, newblerOpts;
	int minOverlapLength, minOverlapIdentity, minLargeContigLength;
};

typedef OptionsBaseTemplate< _NewblerOptions > NewblerOptions;

class Newbler : public ExternalAssembler {
public:
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	Newbler() : ExternalAssembler("newbler", 25) {}
	virtual ~Newbler() {
	}
	Read extendContig(const Read &oldContig, const ReadSet &inputReads) {

		std::string prefix = "/.newbler-assembly";
		std::string outputDir = Cleanup::makeTempDir(Options::getOptions().getTmpDir(), prefix);
		LOG_VERBOSE_OPTIONAL(1, !GeneralOptions::getOptions().getKeepTempDir().empty(), "Saving Newbler working directory for " << oldContig.getName()
				<< " to " << GeneralOptions::getOptions().getKeepTempDir() << outputDir.substr(outputDir.find(prefix)));
		std::string baseName = outputDir + "/input.fastq";
		writeReads(inputReads, oldContig, baseName, FormatOutput::FastqUnmasked());
		std::string log = baseName + ".log";
		// TODO support CDNA & isotigs
		// -cdna -isplit //isotig traversal when depth spikes
		// -icc "$ICC" -icl "$ICL" -it $MAX_NUM_ISOTIGS_PER_ISOGROUP // maximum number of contigs in an isotig, min length, man num
		std::string cmd = NewblerOptions::getOptions().getNewblerPath() + "/runAssembly ";
		cmd += NewblerOptions::getOptions().getNewblerOpts();
		cmd += " -force -m -urt -nofe -noace -nobig -cpu " + boost::lexical_cast<string>(omp_get_num_threads());
		cmd += " -o " + outputDir + "/assembly";
		cmd += " -ml " + boost::lexical_cast<string>( NewblerOptions::getOptions().getMinOverlapLength() );
		cmd += " -mi " + boost::lexical_cast<string>( NewblerOptions::getOptions().getMinOverlapIdentity() );
		cmd += " -l " + boost::lexical_cast<string>( NewblerOptions::getOptions().getMinLargeContigLength() );
		cmd += baseName + " > " + log + " 2>&1";
		LOG_DEBUG_OPTIONAL(1, true, "Executing newbler for " << oldContig.getName() << "(" << inputReads.getSize() << " read pool): " << cmd);

		std::string newContigFile = outputDir + "/assembly/454Scaffolds.fna";
		Read bestRead = executeAssembly(cmd, newContigFile, oldContig);

		if (bestRead.empty())
			LOG_WARN(1, "Could not assemble " << oldContig.getName() << " with pool of " << inputReads.getSize() << " reads: " << FileUtils::dumpFile(log));

		Cleanup::removeTempDir(outputDir);
		return bestRead;
	}
};

#endif /* CAP3_H_ */
