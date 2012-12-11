/*
 * _VelvetOptimizer.h
 *
 *  Created on: Apr 4, 2012
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

#ifndef VELVET_OPTIMIZER_H_
#define VELVET_OPTIMIZER_H_

#include "Options.h"
#include "Utils.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"

class _VelvetOptimizerOptions : public OptionsBaseInterface {
public:
	_VelvetOptimizerOptions() : velvetOptPath() {}
	virtual ~_VelvetOptimizerOptions() {}

	void _resetDefaults() {
		Options::getOptions().getMmapInput() = 0;
	}
	void _setOptions(po::options_description &desc,
		po::positional_options_description &p) {
		po::options_description opts("Cap3 Options");

		opts.add_options()
				("velvet-optimizer-path", po::value<std::string>()->default_value(velvetOptPath), "The full path to the executable VelvetOptimizer.pl.  Set this to activate this module")
				;
		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		setOpt("velvet-optimizer-path", velvetOptPath);
		return ret;
	}
protected:
	std::string velvetOptPath;
};
typedef OptionsBaseTemplate< _VelvetOptimizerOptions > VelvetOptimizerOptions;

class VelvetOptimizer {
public:
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	static const int repeatContig = 4;
	static const int repeatContigCenter = 4;
	static const double minOverlapFraction = 0.9;
	static std::string getNewName(std::string oldName, int deltaLen) {
		std::string newName;
		size_t pos = oldName.find_last_of("+");
		if (pos != oldName.npos) {
			newName = oldName.substr(0, pos );
			deltaLen += atoi(oldName.substr(pos).c_str());
		} else {
			newName = oldName;
		}

		return newName + "+" + boost::lexical_cast<std::string>(deltaLen);
	}
	static Read selectBestContig(const ReadSet &candidateContigs, const Read &targetContig) {
		Read bestRead = targetContig;
		SequenceLengthType oldSize = KmerSizer::getSequenceLength();
		KmerSizer::set(16);
		KmerSpectrum<> tgtSpectrum(16);
		tgtSpectrum.setSolidOnly();
		tgtSpectrum.buildKmerSpectrum(targetContig);
		long tgtKmers = tgtSpectrum.solid.size();
		LOG_DEBUG(4, "Cap3::selectBestContig(): tgtContig: " << targetContig.getLength() << ", " << tgtKmers);

		for(ReadSetSizeType i = 0; i < candidateContigs.getSize(); i++) {
			const Read &tgtRead = candidateContigs.getRead(i);
			KmerWeights readKmers(tgtRead.getTwoBitSequence(), tgtRead.getLength(), true);
			tgtSpectrum.getCounts(readKmers, false);
			long readOverlap = readKmers.sumAll();
			if (readOverlap >= minOverlapFraction * tgtKmers && tgtRead.getLength() >= targetContig.getLength() && bestRead.getLength() < tgtRead.getLength()) {
				bestRead = tgtRead;
				bestRead.setName(getNewName(targetContig.getName(), tgtRead.getLength() - targetContig.getLength()));
			}
			LOG_DEBUG(4, "Cap3::selectBestContig(): tgtRead: " << tgtRead.getLength() << ", " << readOverlap << ": " << bestRead.getLength());
		}

		if (bestRead.getLength() > targetContig.getLength()) {
			LOG_DEBUG_OPTIONAL(2, true, "Cap3 new contig: " << bestRead.getName() << " from " << targetContig.getLength() << " to " << bestRead.getLength());
		} else {
			LOG_DEBUG_OPTIONAL(2, true, "Cap3 failed to extend: " << bestRead.getName());
		}
		KmerSizer::set(oldSize);
		return bestRead;
	}
	static Read extendContig(const Read &oldContig, const ReadSet &_inputReads) {
		ReadSet inputReads = _inputReads;

		for(int i = 0; i < repeatContig; i++)
			inputReads.append(oldContig);

		FormatOutput format = FormatOutput::FastaUnmasked();
		std::string outputDir = Options::getOptions().getTmpDir() + UniqueName::generateUniqueName("/.cap3-assembly");
		int status = mkdir(outputDir.c_str(), 0777);

		if (status != 0) {
			LOG_WARN(1, "Could not mkdir: " << outputDir << " bailing...");
			return Read();
		}
		std::string outputName = outputDir + "/input" + format.getSuffix();
		{
			OfstreamMap ofm(outputName, "");
			inputReads.writeAll(ofm.getOfstream(""), format);
		}
		std::string log = outputName + ".log";
		std::string cmd = Cap3Options::getOptions().getCap3Path() + " " + outputName + " > " + log + " 2>&1";
		LOG_DEBUG_OPTIONAL(1, true, "Executing: " << cmd);
		status = system(cmd.c_str());
		if (status == 0) {
			std::string newContigFile = outputName + ".cap.contigs";
			long fileSize = FileUtils::getFileSize(newContigFile);
			if (fileSize > 0) {
				ReadSet newContig;
				newContig.appendFastaFile(newContigFile);
				Read bestRead = selectBestContig(newContig, oldContig);
				clean(outputDir);
				return bestRead;
			}
		}
		LOG_WARN(1, "Could not assemble " << oldContig.getName() << " with pool of " << _inputReads.getSize() << " reads: " << FileUtils::dumpFile(log));

		clean(outputDir);
		return Read();
	}
	static void clean(std::string fileDir) {
		std::string cmd;
		cmd = "/bin/rm -rf " + fileDir;
		int status = system(cmd.c_str());
		if (status != 0)
			LOG_WARN(1, "Could not clean up directory: " + fileDir);
	}
};

#endif /* VELVET_OPTIMIZER_H_ */
