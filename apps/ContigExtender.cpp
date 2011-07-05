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
#include "Log.h"

using namespace std;
typedef TrackingDataWithDirection DataType;
typedef KmerSpectrum<DataType, DataType> KS;

class ContigExtenderOptions : public Options {
public:
	static std::string getContigFile () {
		return getVarMap()["contig-file"].as<std::string> ();
	}
	static bool parseOpts(int argc, char *argv[]) {
		// override the default output format!
		Options::getFormatOutput() = 3;

		// set options specific to this program
		getPosDesc().add("input-file", -1);

		getDesc().add_options()

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

int main(int argc, char *argv[]) {
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

	SequenceLengthType kmerSize = Options::getKmerSize();
	KmerSizer::set(kmerSize);

	KS readSpectrum(KS::estimateWeakKmerBucketSize(reads, 64));
	readSpectrum.buildKmerSpectrum( reads, false );
	LOG_VERBOSE(1, "Built read kmer spectrum");

	KS contigSpectrum(KS::estimateWeakKmerBucketSize(contigs, 64));
	contigSpectrum.buildKmerSpectrum(contigs, true);
	LOG_VERBOSE(1, "Built contig kmer spectrum");

	ReadSet newContigs;
	double minimumConsensus = 0.90;
	double minimumCoverage = 9.9;

	for(unsigned int i = 0; i < contigs.getSize(); i++) {
		SequenceLengthType len = 0;
		Read read = contigs.getRead(i);
		std::string fasta = read.getFasta();
		while (fasta.length() != len) {
			len = fasta.length();
			if (len < kmerSize)
				continue;

			LOG_VERBOSE(1, "Processing " << fasta);

			TEMP_KMER(tmp);
			TwoBitSequence::compressSequence(fasta.substr(0, kmerSize), tmp.getTwoBitSequence());
			LOG_VERBOSE(1, "Left extending " << tmp.toFasta());

			KmerWeights leftExt = KmerWeights::extendKmer(tmp, false, true);
			double leftTotal = 0.0;
			for(unsigned int i = 0; i < leftExt.size(); i++) {
				KS::WeakElementType elem = readSpectrum.getIfExistsWeak(leftExt[i]);
				if (elem.isValid()) {
					leftTotal += leftExt.valueAt(i) = elem.value().getWeightedCount();
					// do not allow repeats
					if (contigSpectrum.getIfExistsSolid(leftExt[i]).isValid()) {
						LOG_VERBOSE(1, "Left detected repeat: " << leftExt[i].toFasta());
						leftTotal = 0;
						break;
					}
				}
			}
			if (leftTotal >= minimumCoverage) {
				for(unsigned int i = 0; i < leftExt.size(); i++) {
					if (leftExt.valueAt(i) / leftTotal >= minimumConsensus) {
						char base = TwoBitSequence::uncompressBase(i);
						fasta.insert(0, 1, base);
						LOG_VERBOSE(1, "Left extended " << base);
						break;
					}
				}
			} else {
				LOG_VERBOSE(1, "Not enough coverage to extend left: " << leftTotal);
			}

			TwoBitSequence::compressSequence(fasta.substr(fasta.length()-kmerSize, kmerSize), tmp.getTwoBitSequence());
			LOG_VERBOSE(1, "Right extending " << tmp.toFasta());

			KmerWeights rightExt = KmerWeights::extendKmer(tmp, true, true);
			double rightTotal = 0.0;
			for(unsigned int i = 0; i < rightExt.size(); i++) {
				KS::WeakElementType elem = readSpectrum.getIfExistsWeak(rightExt[i]);
				if (elem.isValid()) {
					rightTotal += rightExt.valueAt(i) = elem.value().getWeightedCount();
					// do not allow repeats
					if (contigSpectrum.getIfExistsSolid(rightExt[i]).isValid()) {
						LOG_VERBOSE(1, "Right detected repeat: " << rightExt[i].toFasta());
						rightTotal = 0;
						break;
					}
				}
			}
			if (rightTotal >= minimumCoverage) {
				for(unsigned int i = 0; i < rightExt.size(); i++) {
					if (rightExt.valueAt(i) / rightTotal >= minimumConsensus) {
						char base = TwoBitSequence::uncompressBase(i);
						fasta.insert(fasta.length(), 1, base);
						LOG_VERBOSE(1, "Right extended " << base);
						break;
					}
				}
			} else {
				LOG_VERBOSE(1, "Not enough coverage to extend right: " << rightTotal);
			}

			LOG_VERBOSE(1, "Extended " << fasta << " to " << leftTotal << " " << leftExt.toString() << " and " << rightTotal << " " << rightExt.toString());
		}
		Read newContig(read.getName(), fasta, std::string(fasta.length(), Read::REF_QUAL));
		newContigs.append(newContig);
	}


	OfstreamMap ofmap;
	string outputFilename = Options::getOutputFile();
	Options::getFormatOutput() = FormatOutput::FASTA_UNMASKED;
	if (!outputFilename.empty()) {
		for(unsigned long i = 0; i < newContigs.getSize(); i++)
			newContigs.getRead(i).write(ofmap.getOfstream("-contigs"));
	}

	LOG_VERBOSE(1, "Finished");
}

// $Log: Fastq2Fasta.cpp,v $
// Revision 1.7  2010-05-24 21:48:50  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.6.2.1  2010-05-19 00:20:49  regan
// refactored fomat output options
// added options to fastq2fasta
//
// Revision 1.6  2010-05-18 20:50:18  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.5.2.1  2010-05-07 22:59:29  regan
// refactored base type declarations
//
// Revision 1.5  2010-05-06 21:46:57  regan
// merged changes from PerformanceTuning-20100501
//
//
