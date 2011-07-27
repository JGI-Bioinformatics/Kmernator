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
	static ReadSet extendContigs(const ReadSet &contigs, const ReadSet &reads) {
		SequenceLengthType minKmerSize, maxKmerSize;
		getMinMaxKmerSize(reads, minKmerSize, maxKmerSize);
		return extendContigs(contigs, reads, minKmerSize, maxKmerSize);
	}
	static ReadSet extendContigs(const ReadSet &contigs, const ReadSet &reads, SequenceLengthType minKmerSize, SequenceLengthType maxKmerSize, SequenceLengthType maxSteps = 5) {
		SequenceLengthType enteringKmerSize = KmerSizer::getSequenceLength();
		SequenceLengthType kmerSize = minKmerSize;
		KmerSizer::set(kmerSize);

		double minimumConsensus = ContigExtenderOptions::getMinimumConsensus();
		double minimumCoverage = ContigExtenderOptions::getMinimumCoverage();

		SequenceLengthType kmerStep = maxKmerSize - minKmerSize / maxSteps;
		// ensure is even
		kmerStep = (kmerStep & 1) == 1 ? kmerStep + 1 : kmerStep;
		kmerStep = std::min(2, kmerStep);

		LOG_DEBUG_OPTIONAL(1, true, "Starting extendContigs with consensus fraction " << minimumConsensus << " and coverage " << minimumCoverage << " using kmers " << minKmerSize << " to " << maxKmerSize << " step " << kmerStep << " and " << reads.getSize() << " reads");

		// set scoping of spectrum arrays
		std::vector<KS> readSpectrums(maxKmerSize+1, KS()), contigSpectrums(maxKmerSize+1, KS());
		for(kmerSize = minKmerSize ; kmerSize <= maxKmerSize; kmerSize+=kmerStep) {
			KmerSizer::set(kmerSize);
			LOG_DEBUG_OPTIONAL(1, true, "Building kmer spectrum for kmer sized: " << kmerSize << " over reads sized " << reads.getSize());
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

			std::string fasta = read.getFasta();
			LOG_VERBOSE_OPTIONAL(2, true, "Extending " << read.getName() << " of length " << fasta.length());

			ReadSet thisReadOnlySet;
			thisReadOnlySet.append(read);
			for(kmerSize = minKmerSize ; kmerSize <= maxKmerSize; kmerSize+=2) {
				contigSpectrums[kmerSize].reset(false);
				contigSpectrums[kmerSize].buildKmerSpectrum(thisReadOnlySet,true);
			}

			std::string::size_type leftTotal = 0, rightTotal = 0;
			while (extendLeft | extendRight) {
				SequenceLengthType len = fasta.length();
				if (len < minKmerSize)
					break;

				if (extendLeft) {
					bool toRight = false;
					for(kmerSize = minKmerSize; kmerSize <= maxKmerSize; kmerSize += 2) {
						LOG_DEBUG_OPTIONAL(2, true, "Extending for kmer size " << kmerSize);
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

			std::string newName = getNewName(read.getName(), leftTotal, rightTotal);
			LOG_DEBUG_OPTIONAL(1, true, "Extended " << newName << " : " << read.getName() << " left +" << leftTotal << " right +" << rightTotal << " to " << fasta.length());
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

	static std::string getNewName(std::string oldName, std::string::size_type leftTotal, std::string::size_type rightTotal) {
		std::string::size_type preLeftTotal = 0, preRightTotal = 0, pos, pos2;
		std::string newName = oldName;
		pos = oldName.find_last_of("-l");
		if (pos != std::string::npos) {
			pos2 = oldName.find_first_of("r", pos);
			if (pos2 != std::string::npos) {
				LOG_DEBUG_OPTIONAL(5, true, "ContigExtender::getNewName: " << oldName << " " << pos << " " << pos2);
				LOG_DEBUG_OPTIONAL(5, true, oldName.substr(pos+1, pos2-pos-1) << " " << oldName.substr(pos2+1));
				preLeftTotal = atoi(oldName.substr(pos+1, pos2-pos-1).c_str());
				preRightTotal = atoi(oldName.substr(pos2+1).c_str());
				newName = oldName.substr(0, pos-1);
			}
		}
		newName = newName + "-l" + boost::lexical_cast<std::string>(leftTotal+preLeftTotal) + "r" + boost::lexical_cast<std::string>(rightTotal+preRightTotal);
		LOG_DEBUG_OPTIONAL(1, true, "ContigExtender::getNewName(" << oldName << ", " << leftTotal << ", " << rightTotal << "): " << newName);
		return newName;
	}
	static void getMinMaxKmerSize(const ReadSet &reads, SequenceLengthType &minKmerSize, SequenceLengthType &maxKmerSize) {
		minKmerSize = Options::getKmerSize();
		SequenceLengthType maxLen = std::min(reads.getMaxSequenceLength(), reads.getBaseCount() / reads.getSize());
		maxKmerSize = std::min((SequenceLengthType) maxLen*0.80, maxLen - 1);
		maxKmerSize = std::max(minKmerSize, maxKmerSize);
	}
};

#endif /* CONTIGEXTENDER_H_ */
