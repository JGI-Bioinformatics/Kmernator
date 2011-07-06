//
// Kmernator/apps/ContigExtender.cpp
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
#include "ReadSet.h"
#include "Options.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "DuplicateFragmentFilter.h"
#include "Log.h"

using namespace std;
typedef TrackingDataWithDirection DataType;
typedef KmerSpectrum<DataType, DataType> KS;

class ContigExtenderOptions : public Options {
public:
	static std::string getContigFile () {
		return getVarMap()["contig-file"].as<std::string> ();
	}
	static double getMinimumConsensus() {
		return getVarMap()["minimum-consensus"].as<double>() / 100;
	}
	static double getMinimumCoverage() {
		return getVarMap()["minimum-coverage"].as<double>();
	}
	static bool parseOpts(int argc, char *argv[]) {
		// override the default output format!
		Options::getFormatOutput() = 3;

		// set options specific to this program
		getPosDesc().add("input-file", -1);

		getDesc().add_options()

		("minimum-consensus", po::value<double>()->default_value(95),
				"minimum percent consensus to call the next base")

		("minimum-coverage", po::value<double>()->default_value(9.9),
				"minimum (probability-weighted) coverage to continue calling the next base")

		("contig-file", po::value<std::string>(),
				"filename of input contigs.fa");

		bool ret = Options::parseOpts(argc, argv);

		if (getContigFile().empty() || getInputFiles().empty()) {
			LOG_ERROR(1, "you must specify the --contig-file and one or more input files");
			ret = false;
		}
		return ret;
	}
};

bool extendContig(std::string &fasta, bool toRight, double minimumCoverage, double minimumConsensus, KS &readSpectrum, KS &contigSpectrum) {

	SequenceLengthType kmerSize = KmerSizer::getSequenceLength();
	TEMP_KMER(tmp);
	std::string dir = toRight ? "Right" : "Left";

	std::string kmerString = toRight ? fasta.substr(fasta.length()-kmerSize, kmerSize) : fasta.substr(0, kmerSize);
	TwoBitSequence::compressSequence(kmerString, tmp.getTwoBitSequence());
	LOG_DEBUG(2, dir << " extending " << tmp.toFasta());

	KmerWeights ext = KmerWeights::extendKmer(tmp, toRight, true);
	double total = 0.0;
	for(unsigned int i = 0; i < ext.size(); i++) {
		KS::WeakElementType elem = readSpectrum.getIfExistsWeak(ext[i]);
		if (elem.isValid()) {
			total += ext.valueAt(i) = elem.value().getWeightedCount();
		}
	}
	bool wasExtended = false;

	double best = 0;
	if (total >= minimumCoverage) {
		for(unsigned int i = 0; i < ext.size(); i++) {
			double consensus = ext.valueAt(i) / total;
			best = std::max(consensus, best);
			if (consensus >= minimumConsensus) {
				// do not allow repeats
				if (contigSpectrum.getIfExistsSolid(ext[i]).isValid()) {
					LOG_VERBOSE(1, dir << " detected repeat: " << ext[i].toFasta() << " with kmer " << KmerSizer::getSequenceLength());
					break;
				}
				char base = TwoBitSequence::uncompressBase(i);
				fasta.insert(toRight ? fasta.length() : 0, 1, base);
				LOG_DEBUG(2, dir << " extended " << base << "\t" << fasta);
				wasExtended = true;
				break;
			}
		}
		if (! wasExtended ) {
			LOG_VERBOSE(1, "Ambiguous " << dir << " extension " << best << " " << ext.toString() << " with kmer " << KmerSizer::getSequenceLength())
		}
	} else {
		LOG_VERBOSE(1, "Not enough coverage to extend " << dir << " " << total << " " << ext.toString() << " with kmer " << KmerSizer::getSequenceLength());
	}

	return wasExtended;
}

void recordKmer(std::vector<KS> &contigSpectrums, bool toRight, std::string fasta, SequenceLengthType minKmerSize, SequenceLengthType maxKmerSize) {
	SequenceLengthType kmerSize;
	KmerSizer::set(maxKmerSize);
	TEMP_KMER(tmp);

	for(kmerSize = minKmerSize ; kmerSize <= maxKmerSize; kmerSize+=2) {
		if (fasta.length() < kmerSize)
			break;
		KmerSizer::set(kmerSize);
		SequenceLengthType pos = toRight ? fasta.length() - kmerSize : 0;
		tmp.set(fasta.substr(pos, kmerSize), true);
		contigSpectrums[kmerSize].appendSolid( tmp );
	}
}

ReadSet extendContigs(const ReadSet &contigs, const ReadSet &reads, SequenceLengthType minKmerSize, SequenceLengthType maxKmerSize) {
	SequenceLengthType kmerSize = minKmerSize;
	KmerSizer::set(kmerSize);

	// set scoping of spectrum arrays
	std::vector<KS> readSpectrums(maxKmerSize+1, KS()), contigSpectrums(maxKmerSize+1, KS());
	for(kmerSize = minKmerSize ; kmerSize <= maxKmerSize; kmerSize+=2) {
		KmerSizer::set(kmerSize);
		LOG_VERBOSE(2, "Building kmer spectrum for kmer sized: " << kmerSize);
		KS readSpectrum(KS::estimateWeakKmerBucketSize(reads, 64));
		readSpectrum.buildKmerSpectrum( reads, false );
		KS contigSpectrum(128);

		readSpectrums[kmerSize].swap(readSpectrum);
		contigSpectrums[kmerSize].swap(contigSpectrum);

	}

	double minimumConsensus = ContigExtenderOptions::getMinimumConsensus();
	double minimumCoverage = ContigExtenderOptions::getMinimumCoverage();

	ReadSet newContigs;

	for(unsigned int i = 0; i < contigs.getSize(); i++) {
		bool extendLeft = true;
		bool extendRight = true;
		Read read = contigs.getRead(i);

		ReadSet thisReadOnlySet;
		thisReadOnlySet.append(read);
		for(kmerSize = minKmerSize ; kmerSize <= maxKmerSize; kmerSize+=2) {
			contigSpectrums[kmerSize].reset(false);
			contigSpectrums[kmerSize].buildKmerSpectrum(thisReadOnlySet,true);
		}

		std::string fasta = read.getFasta();
		LOG_VERBOSE(1, "Extending " << read.getName() << " of length " << fasta.length());
		int leftTotal = 0, rightTotal = 0;
		while (extendLeft | extendRight) {
			SequenceLengthType len = fasta.length();
			if (len < minKmerSize)
				break;

			if (extendLeft) {
				bool toRight = false;
				for(kmerSize = minKmerSize; kmerSize <= maxKmerSize; kmerSize += 2) {
					LOG_DEBUG(1, "Extending for kmer size " << kmerSize);
					KmerSizer::set(kmerSize);
					extendLeft = extendContig(fasta, toRight, minimumCoverage, minimumConsensus, readSpectrums[kmerSize], contigSpectrums[kmerSize]);
					if (extendLeft) {
						recordKmer(contigSpectrums, toRight, fasta, minKmerSize, maxKmerSize);
						leftTotal++;
						break;
					}
				}
			}

			if (extendRight) {
				bool toRight = true;
				for(kmerSize = minKmerSize; kmerSize <= maxKmerSize; kmerSize += 2) {
					KmerSizer::set(kmerSize);
					extendRight =  extendContig(fasta, toRight, minimumCoverage, minimumConsensus, readSpectrums[kmerSize], contigSpectrums[kmerSize]);
					if (extendRight) {
						recordKmer(contigSpectrums, toRight, fasta, minKmerSize, maxKmerSize);
						rightTotal++;
						break;
					}
				}
			}
		}

		std::string newName = read.getName() + "-l" + boost::lexical_cast<std::string>(leftTotal) + "r" + boost::lexical_cast<std::string>(rightTotal);
		LOG_VERBOSE(1, "Extended " << newName << " : " << read.getName() << " left +" << leftTotal << " right +" << rightTotal << " to " << fasta.length());
		Read newContig(newName, fasta, std::string(fasta.length(), Read::REF_QUAL));
		newContigs.append(newContig);
	}

	for(kmerSize = minKmerSize ; kmerSize < maxKmerSize; kmerSize+=2) {
		KmerSizer::set(kmerSize);
		// must explicitly release memory as KmerSizer could cause memory leaks if set differently from when created
		readSpectrums[kmerSize-minKmerSize].reset(true);
		contigSpectrums[kmerSize-minKmerSize].reset(true);
	}
	return newContigs;
}

int main(int argc, char *argv[]) {
	Options::getSkipArtifactFilter() = 1;
	if (!ContigExtenderOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	Options::FileListType inputFiles = ContigExtenderOptions::getInputFiles();
	Options::FileListType contigFiles;
	contigFiles.push_back(ContigExtenderOptions::getContigFile());

	ReadSet reads;
	LOG_VERBOSE(1, "Reading Input Files" );
	reads.appendAllFiles(inputFiles);

	LOG_VERBOSE(1, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount() << " Bases ");
	reads.identifyPairs();

	ReadSet contigs;
	LOG_VERBOSE(1, "Reading Contig File" );
	contigs.appendAllFiles(contigFiles);
	LOG_VERBOSE(1, "loaded " << contigs.getSize() << " Reads, " << contigs.getBaseCount() << " Bases ");

	if (Options::getDeDupMode() > 0 && Options::getDeDupEditDistance() >= 0) {
	  LOG_VERBOSE(2, "Applying DuplicateFragmentPair Filter to Input Files");
	  unsigned long duplicateFragments = DuplicateFragmentFilter::filterDuplicateFragments(reads);
	  LOG_VERBOSE(1, "filter removed duplicate fragment pair reads: " << duplicateFragments);
	  LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
	}

	SequenceLengthType minKmerSize = Options::getKmerSize();
	SequenceLengthType maxKmerSize = std::min(reads.getMaxSequenceLength(), (SequenceLengthType) (reads.getBaseCount() / reads.getSize())) - 1;
	maxKmerSize = std::min(maxKmerSize, contigs.getMaxSequenceLength());
	maxKmerSize = maxKmerSize > 3*minKmerSize ? 3*minKmerSize : maxKmerSize;
	minKmerSize = maxKmerSize > minKmerSize ? minKmerSize : maxKmerSize - 2;
	if (maxKmerSize < 4) {
		LOG_ERROR(1, "There is not enough data to run a kmer walk in input files.");
		exit(1);
	}

	//maxKmerSize = minKmerSize;
	ReadSet newContigs = extendContigs(contigs, reads, minKmerSize, maxKmerSize);

	OfstreamMap ofmap;
	string outputFilename = Options::getOutputFile();
	Options::getFormatOutput() = FormatOutput::FASTA_UNMASKED;
	if (!outputFilename.empty()) {
		for(unsigned long i = 0; i < newContigs.getSize(); i++)
			newContigs.getRead(i).write(ofmap.getOfstream("-contigs"));
	}

	LOG_VERBOSE(1, "Finished");
}

