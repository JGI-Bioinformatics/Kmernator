/*
 * ContigExtender.h
 *
 *  Created on: Jul 7, 2011
 *      Author: regan
 */

#ifndef CONTIGEXTENDER_H_
#define CONTIGEXTENDER_H_

#include "config.h"
#include "ReadSet.h"
#include "Options.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "Log.h"

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



template <typename KS>
class ContigExtender {
protected:
	static void recordKmer(std::vector<KS> &contigSpectrums, bool toRight, std::string fasta, SequenceLengthType minKmerSize, SequenceLengthType maxKmerSize) {
		SequenceLengthType enteringKmerSize = KmerSizer::getSequenceLength();
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
		KmerSizer::set(enteringKmerSize);
	}

public:
	static ReadSet extendContigs(const ReadSet &contigs, const ReadSet &reads, SequenceLengthType minKmerSize, SequenceLengthType maxKmerSize) {
		SequenceLengthType enteringKmerSize = KmerSizer::getSequenceLength();
		SequenceLengthType kmerSize = minKmerSize;
		KmerSizer::set(kmerSize);

		double minimumConsensus = ContigExtenderOptions::getMinimumConsensus();
		double minimumCoverage = ContigExtenderOptions::getMinimumCoverage();

		LOG_VERBOSE(1, "Starting extendContigs with consensus fraction " << minimumConsensus << " and coverage " << minimumCoverage << " using kmers " << minKmerSize << " to " << maxKmerSize);

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
						extendLeft = readSpectrums[kmerSize].extendContig(fasta, toRight, minimumCoverage, minimumConsensus, &contigSpectrums[kmerSize]);
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
						extendRight =  readSpectrums[kmerSize].extendContig(fasta, toRight, minimumCoverage, minimumConsensus, &contigSpectrums[kmerSize]);
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
		KmerSizer::set(enteringKmerSize);
		return newContigs;
	}
};

#endif /* CONTIGEXTENDER_H_ */
