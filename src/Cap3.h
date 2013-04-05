/*
 * Cap3.h
 *
 *  Created on: Sep 7, 2011
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

#ifndef CAP3_H_
#define CAP3_H_

#include "Options.h"
#include "Utils.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "Utils.h"

class _Cap3Options : public OptionsBaseInterface {
public:
	_Cap3Options() : cap3Path() {}
	virtual ~_Cap3Options() {}
	std::string &getCap3Path() {
		return cap3Path;
	}
	void _resetDefaults() {
		Options::getOptions().getMmapInput() = false;
	}
	void _setOptions(po::options_description &desc,	po::positional_options_description &p) {
		po::options_description opts("Cap3 Options");

		opts.add_options()
				("cap3-path", po::value<std::string>()->default_value(cap3Path), "if set, cap3 will be used to extend contigs")

				;

		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		setOpt("cap3-path", cap3Path);

		return ret;
	}
protected:
	std::string cap3Path;
};
typedef OptionsBaseTemplate< _Cap3Options > Cap3Options;

class Cap3 {
public:
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	static const int repeatContig = 1;
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

	Read selectBestContig(const ReadSet &candidateContigs, const Read &targetContig) {
		LOG_DEBUG(2, "selectBestRead(): for " << targetContig.getName() << " " << targetContig.getLength());
		Read bestRead = targetContig;
		kmerReadUtils.buildReferenceMap(targetContig);
		SequenceLengthType contigLength = targetContig.getLength();
		double minExtensionFactor = ContigExtenderBaseOptions::getOptions().getMinimumExtensionFactor();

		SequenceLengthType minOverlap = contigLength * minExtensionFactor;
		int bestOverlap = minOverlap - 1;
		for(ReadSetSizeType i = 0; i < candidateContigs.getSize(); i++) {
			const Read &candidateRead = candidateContigs.getRead(i);
			if (candidateRead.getLength() <= contigLength) {
				LOG_DEBUG(2, "selectBestContig(): candidateRead is too short: " << candidateRead.getLength());
				continue;
			}
			long readOverlap = kmerReadUtils.getReferenceMapOverlap(candidateRead);
			if (readOverlap > bestOverlap) {
				bestOverlap = readOverlap;
				bestRead = candidateRead;
				bestRead.setName(getNewName(targetContig.getName(), candidateRead.getLength() - targetContig.getLength()));
				LOG_DEBUG(2, "Cap3::selectBestContig(): candidateRead: " << candidateRead.getLength() << ", " << readOverlap << ": " << bestRead.getLength() << " " << bestRead.getName());
			}
		}

		if (bestRead.getLength() > targetContig.getLength()) {
			LOG_DEBUG(2, "Cap3 new contig: " << bestRead.getName() << " from " << targetContig.getLength() << " to " << bestRead.getLength());
		} else {
			LOG_DEBUG(2, "Cap3 failed to extend: " << bestRead.getName());
		}
		return bestRead;
	}
	Read extendContig(const Read &oldContig, const ReadSet &_inputReads) {
		ReadSet inputReads = _inputReads;
		int status;

		for(int i = 0; i < repeatContig; i++)
			inputReads.append(oldContig);

		FormatOutput format = FormatOutput::FastaUnmasked();
		std::string prefix = "/.cap3-assembly";
		std::string outputDir = Cleanup::makeTempDir(Options::getOptions().getTmpDir(), prefix);
		LOG_VERBOSE_OPTIONAL(1, !GeneralOptions::getOptions().getKeepTempDir().empty(), "Saving Cap3 working directory for " << oldContig.getName()
				<< " to " << GeneralOptions::getOptions().getKeepTempDir() << outputDir.substr(outputDir.find(prefix)));
		std::string outputName = outputDir + "/input" + format.getSuffix();
		{
			OfstreamMap ofm(outputName, "");
			inputReads.writeAll(ofm.getOfstream(""), format);
		}
		std::string log = outputName + ".log";
		std::string cmd = Cap3Options::getOptions().getCap3Path() + " " + outputName + " > " + log + " 2>&1";
		LOG_DEBUG_OPTIONAL(1, true, "Executing cap3 for " << oldContig.getName() << "(" << inputReads.getSize() << " read pool): " << cmd);
		status = system(cmd.c_str());
		if (status == 0) {
			std::string newContigFile = outputName + ".cap.contigs";
			long fileSize = FileUtils::getFileSize(newContigFile);
			if (fileSize > 0) {
				ReadSet newContig;
				newContig.appendFastaFile(newContigFile);
				Read bestRead = selectBestContig(newContig, oldContig);
				Cleanup::removeTempDir(outputDir);
				return bestRead;
			}
		}
		LOG_WARN(1, "Could not assemble " << oldContig.getName() << " with pool of " << _inputReads.getSize() << " reads: " << FileUtils::dumpFile(log));

		Cleanup::removeTempDir(outputDir);
		return Read();
	}
private:
	KmerReadUtils kmerReadUtils;
};

#endif /* CAP3_H_ */
