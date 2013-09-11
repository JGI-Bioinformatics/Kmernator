/*
 * ExternalAssembler.h
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

#ifndef EXTERNAL_ASSEMBLER_H
#define EXTERNAL_ASSEMBLER_H

#include <stdlib.h>
#include <boost/date_time.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

#include "ReadSet.h"
#include "Log.h"
#include "Utils.h"
#include "KmerReadUtils.h"

class ExternalAssembler {
public:
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;

	int &getRepeatContig() {
		return repeatContig;
	}
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
	ExternalAssembler(std::string _name, int _repeatContig = 1, int _maxSeedLength = 1999, int _shredStep = 250):
		myName(_name), kmerReadUtils(), repeatContig(_repeatContig), maxSeedLength(_maxSeedLength), shredStep(_shredStep) {}
	virtual ~ExternalAssembler() {}

	void writeReads(const ReadSet &reads, const Read &oldContig, std::string outputName, FormatOutput format = FormatOutput::FastaUnmasked()) {
		OfstreamMap ofm(outputName, "");
		reads.writeAll(ofm.getOfstream(""), format);
		if (oldContig.getLength() > maxSeedLength) {
			LOG_DEBUG(1, "Shredding " << oldContig.getName() << " " << oldContig.getLength());
			for(int i = 0; i < repeatContig; i++) {
				ReadSet shreds = ReadSet::shred(oldContig, maxSeedLength, shredStep);
				for(int j = 0 ; j < (int) shreds.getSize(); j++) {
					shreds.getRead(j).write(ofm.getOfstream(""), format);
				}
			}
		} else {
			for(int i = 0; i < repeatContig; i++) {
				oldContig.write(ofm.getOfstream(""), format);
			}
		}
	}


	Read selectBestContig(const ReadSet &candidateContigs, const Read &targetContig) {
		LOG_DEBUG(2, myName << "::selectBestRead(): for " << targetContig.getName() << " " << targetContig.getLength());
		Read bestRead = targetContig;
		kmerReadUtils.buildReferenceMap(targetContig);
		SequenceLengthType contigLength = targetContig.getLength();
		double minExtensionFactor = ContigExtenderBaseOptions::getOptions().getMinimumExtensionFactor();

		SequenceLengthType minOverlap = contigLength * minExtensionFactor;
		int bestOverlap = minOverlap - 1;
		for(ReadSetSizeType i = 0; i < candidateContigs.getSize(); i++) {
			const Read &candidateRead = candidateContigs.getRead(i);
			if (candidateRead.getLength() <= contigLength) {
				LOG_DEBUG(2, myName << "::selectBestContig(): candidateRead is too short: " << candidateRead.getLength());
				continue;
			}
			long readOverlap = kmerReadUtils.getReferenceMapOverlap(candidateRead);
			if (readOverlap > bestOverlap) {
				bestOverlap = readOverlap;
				bestRead = candidateRead;
				bestRead.setName(getNewName(targetContig.getName(), candidateRead.getLength() - targetContig.getLength()));
				LOG_DEBUG(2, myName << "::selectBestContig(): candidateRead: " << candidateRead.getLength() << ", " << readOverlap << ": " << bestRead.getLength() << " " << bestRead.getName());
			}
		}

		if (bestRead.getLength() > targetContig.getLength()) {
			LOG_DEBUG(2, myName << " new contig: " << bestRead.getName() << " from " << targetContig.getLength() << " to " << bestRead.getLength());
		} else {
			LOG_DEBUG(2, myName << " failed to extend: " << bestRead.getName());
		}
		return bestRead;
	}

	Read executeAssembly(std::string cmd, std::string newContigFile, const Read &oldContig) {
		int status;
		status = system(cmd.c_str());
		if (status == 0) {

			long fileSize = FileUtils::getFileSize(newContigFile);
			if (fileSize > 0) {
				ReadSet newContig;
				newContig.appendFastaFile(newContigFile);
				Read bestRead = selectBestContig(newContig, oldContig);
				return bestRead;
			} else {
				LOG_WARN(1, myName << "::executeAssembly(): returned no file or zero byte file: " << newContigFile << " cmd: " << cmd);
			}
		} else {
			LOG_WARN(1, myName << "::executeAssembly(): cmd failed (" << status << "): " << cmd);
		}
		return Read();
	}

	std::string extendContigs(const ReadSet & contigs,
			ReadSet::ReadSetVector &contigReadSet, ReadSet & changedContigs,
			ReadSet & finalContigs, ReadSet::ReadSetSizeType minimumCoverage, int threadNum, int numThreads) {
		std::stringstream extendLog;
		int poolsWithoutMinimumCoverage = 0;

		long processedContigs = 0;
		for (long i = threadNum; i < (long) contigs.getSize(); i+=numThreads) {
			processedContigs++;
			const Read &oldRead = contigs.getRead(i);
			Read newRead = oldRead;
			SequenceLengthType oldLen = oldRead.getLength(), newLen = 0;

			ReadSet::ReadSetSizeType poolSize = contigReadSet[i].getSize();

			boost::posix_time::ptime extTime = boost::posix_time::microsec_clock::local_time();
			if (poolSize > minimumCoverage) {
				LOG_VERBOSE_OPTIONAL(2, true, myName << "::extendContigs(): Extending " << oldRead.getName() << " with " << poolSize << " pool of reads");
				newRead = extendContig(oldRead, contigReadSet[i]);
				newLen = newRead.getLength();
			} else {
				poolsWithoutMinimumCoverage++;
			}

			long microseconds =  (boost::posix_time::microsec_clock::local_time() - extTime).total_microseconds();
			long deltaLen = (long)newLen - (long)oldLen;
			if (deltaLen > 0) {
				extendLog << std::endl << myName << "::extendContigs(): Extended " << oldRead.getName() << " "
						<< deltaLen << " bases to " << newRead.getLength() << ": "
						<< newRead.getName() << " with " << poolSize
						<< " reads in the pool, in " << microseconds/1000  << " msec";
				//#pragma omp critical
				changedContigs.append(newRead);
			} else {
				extendLog << std::endl << "Did not extend " << oldRead.getName() << " with " << poolSize << " reads in the pool, in " << extTime << " sec";
				//#pragma omp critical
				finalContigs.append(oldRead);
			}
		}

		LOG_VERBOSE_OPTIONAL(2, true, myName << "::::extendContigs(): Extended " << contigs.getSize() - poolsWithoutMinimumCoverage << " contigs out of " << processedContigs << " total " << contigs.getSize());

		return extendLog.str();
	}

	virtual Read extendContig(const Read &oldContig, const ReadSet &_inputReads) = 0;

private:
	std::string myName;
	KmerReadUtils kmerReadUtils;
	int repeatContig;
	SequenceLengthType maxSeedLength, shredStep;

};

#endif
