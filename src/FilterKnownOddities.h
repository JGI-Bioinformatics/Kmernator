//
// Kmernator/src/FilterKnownOddities.h
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

#ifndef _FILTER_H
#define _FILTER_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>

#include "config.h"
#include "Kmer.h"
#include "ReadSet.h"
#include "KmerReadUtils.h"
#include "KmerSpectrum.h"
#include "Options.h"
#include "Log.h"
#include "ReadSelector.h"

class _FilterKnownOdditiesOptions : public OptionsBaseInterface {
public:
	_FilterKnownOdditiesOptions() : skipArtifactFilter(false), artifactFilterMatchLength(24),
	artifactFilterEditDistance(2), buildArtifactEditsInFilter(2),
	maskSimpleRepeats(false), phiXOutput(false), filterOutput(false) {}
	~_FilterKnownOdditiesOptions() {}

	unsigned int &getArtifactFilterEditDistance()
	{
		return artifactFilterEditDistance;
	}
	unsigned int &getArtifactFilterMatchLength()
	{
		return artifactFilterMatchLength;
	}
	FileListType &getArtifactReferenceFiles()
	{
		return artifactReferenceFiles;
	}
	unsigned int &getBuildArtifactEditsInFilter()
	{
		return buildArtifactEditsInFilter;
	}
	bool &getFilterOutput()
	{
		return filterOutput;
	}
	bool &getMaskSimpleRepeats()
	{
		return maskSimpleRepeats;
	}
	bool &getPhiXOutput()
	{
		return phiXOutput;
	}
	bool &getSkipArtifactFilter()
	{
		return skipArtifactFilter;
	}


	void _resetDefaults() {
		// *::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		// *::_setOptions(desc,p);
		po::options_description opts("Filter Known Artifacts and Oddities Options");
		opts.add_options()

	    				("phix-output", po::value<bool>()->default_value(phiXOutput), "if set, artifact filter also screens for PhiX174, and any matching reads will be output into a separate file (requires --output-file set)")

	    				("filter-output", po::value<bool>()->default_value(filterOutput), "if set, artifact filter reads will be output into a separate file. If not set, then affected reads will be trimmed and then output normally.  (requires --output-file set)")

	    				("skip-artifact-filter", po::value<bool>()->default_value(skipArtifactFilter), "if set, Skip homo-polymer, primer-dimer and duplicated fragment pair filtering")

	    				("artifact-match-length", po::value<unsigned int>()->default_value(artifactFilterMatchLength), "Kmer match length to known artifact sequences")

	    				("artifact-edit-distance", po::value<unsigned int>()->default_value(artifactFilterEditDistance), "edit-distance to apply to artifact-match-length matches to know artifacts")

	    				("build-artifact-edits-in-filter", po::value<unsigned int>()->default_value(buildArtifactEditsInFilter), "0 - edits will be applied to reads on the fly, 1 - edits will be pre-build in the filter (needs more memory, less overall CPU), 2 - automatic based on size")

	    				("mask-simple-repeats", po::value<bool>()->default_value(maskSimpleRepeats), "if set filtering artifacts will also mask simple repeats")

	    				("artifact-reference-file", po::value<FileListType>(), "additional artifact reference file(s)");

		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		// ret &= *::_parseOptions(vm);
		// set skipArtifactFiltering
		setOpt("skip-artifact-filter", skipArtifactFilter);
		setOpt("artifact-match-length", artifactFilterMatchLength);
		setOpt("artifact-edit-distance", artifactFilterEditDistance);
		setOpt("build-artifact-edits-in-filter", buildArtifactEditsInFilter);

		// set simple repeat masking
		setOpt("mask-simple-repeats", maskSimpleRepeats);

		// set phix masking
		setOpt("phix-output", phiXOutput);

		setOpt2("artifact-reference-file", artifactReferenceFiles);

		// set simple repeat masking
		setOpt("filter-output", filterOutput);

		return ret;
	}

protected:
	bool skipArtifactFilter;
	unsigned int artifactFilterMatchLength,
	artifactFilterEditDistance, buildArtifactEditsInFilter;
	bool maskSimpleRepeats, phiXOutput, filterOutput;
	FileListType artifactReferenceFiles;
};
typedef OptionsBaseTemplate< _FilterKnownOdditiesOptions > FilterKnownOdditiesOptions;


class FilterKnownOddities {
public:
	typedef Kmer::NumberType NumberType;
	typedef KmerMap<unsigned int> KM;
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	typedef std::vector< ReadSetSizeType > SequenceCounts;
	typedef std::vector< long > BaseCounts;
	typedef std::vector< ReadSet > ReadSets;

private:
	ReadSet sequences;
	unsigned short length;
	unsigned short twoBitLength;
	KM filter;
	SequenceCounts counts;
	int numErrors;
	ReadSets remnantReads;

public:
	FilterKnownOddities(int _length = FilterKnownOdditiesOptions::getOptions().getArtifactFilterMatchLength(), int _numErrors = FilterKnownOdditiesOptions::getOptions().getArtifactFilterEditDistance()) :
		length(_length), filter(512*1024), numErrors(_numErrors) {
		if (length > 28) {
			LOG_THROW("Invalid: FilterKnownOddities must use 7 bytes or less (<= 28 bases)");
		}
		twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(length);
		// T is 11, A is 00, so mask is all T's surrounded by A's

		// readIdx 0 is signal for no match!
		Read empty;
		sequences.append(empty);

		std::string fasta = getArtifactFasta();
		sequences.appendFastaData(fasta);
		if (FilterKnownOdditiesOptions::getOptions().getMaskSimpleRepeats()) {
			fasta = getSimpleRepeatFasta();
			getSimpleRepeatBegin() = sequences.getSize();
			sequences.appendFastaData(fasta);
			getSimpleRepeatEnd() = sequences.getSize();
		}
		if (FilterKnownOdditiesOptions::getOptions().getPhiXOutput()) {
			fasta = getPhiX();
			getPhiXReadIdx() = sequences.getSize();
			sequences.appendFastaData(fasta);
		}
		sequences.circularize(length);

		OptionsBaseInterface::FileListType artifacts = FilterKnownOdditiesOptions::getOptions().getArtifactReferenceFiles();
		if (!artifacts.empty()) {
			getReferenceReadIdx() = sequences.getSize();
			for(unsigned int i = 0 ; i < artifacts.size(); i++)
				sequences.appendAnyFile(artifacts[i]);
		}
		counts.resize( sequences.getSize()+1 );

		prepareMaps();

		remnantReads.resize( omp_get_max_threads(), ReadSet());

	}
	~FilterKnownOddities() {
		clear();
	}
	void clear() {
		unsigned long oldKmerLength = KmerSizer::getSequenceLength();
		KmerSizer::set(length);
		filter.clear();
		counts.clear();
		numErrors = 0;
		KmerSizer::set(oldKmerLength);
	}

	void prepareMaps() {
		unsigned long oldKmerLength = KmerSizer::getSequenceLength();
		assert(length % 4 == 0);
		KmerSizer::set(length);

		LOG_DEBUG(2, "Preparing exact match artifacts");
		KmerReadUtils kru;
		for (unsigned int i = 0; i < sequences.getSize(); i++) {
			const Read &read = sequences.getRead(i);
			KmerWeightedExtensions &kmers = kru.buildWeightedKmers(read, true);
			for (Kmer::IndexType j = 0; j < kmers.size(); j++) {
				filter.getOrSetElement( kmers[j] , i );
			}
			if (Log::isDebug(2) && i % 10000 == 0)
				LOG_DEBUG(2, "Processed" << i << " artifact reads." << filter.size() << "" << MemoryUtils::getMemoryUsage())
		}
		LOG_DEBUG(2, "Processed" << sequences.getSize() << " artifact reads." << filter.size() << "" << MemoryUtils::getMemoryUsage())

		int maxErrors = numErrors;
		int buildEdits = FilterKnownOdditiesOptions::getOptions().getBuildArtifactEditsInFilter();
		for (int error = 0; error < maxErrors; error++) {
			if (buildEdits == 1 || (buildEdits == 2 && (filter.size() < 750000))) {
				numErrors--;

				LOG_DEBUG(2,  "Preparing edit distance" << (error+1) << " match:" << filter.size() << "" << MemoryUtils::getMemoryUsage() );

				std::vector< KM::BucketType > tmpKmers;
				tmpKmers.reserve(filter.size());
				for(KM::Iterator it = filter.begin(); it != filter.end(); it++) {
					tmpKmers.push_back( KM::BucketType::permuteBases(it->key(), it->value(), true) );
				}
				for(std::vector< KM::BucketType >::iterator it = tmpKmers.begin(); it != tmpKmers.end(); it++) {
					KM::BucketType &kmers = *it;
					for (Kmer::IndexType j = 0; j < kmers.size(); j++) {
						filter.getOrSetElement( kmers[j] , kmers.valueAt(j) );
					}
				}
				tmpKmers.clear();
				LOG_DEBUG(2, "Prepared order" << (error+1) << ":" << filter.size() << "" << MemoryUtils::getMemoryUsage() );
			}
		}
		LOG_DEBUG(2, "filter is " << filter.size() << ".  Remaining edits is:" << numErrors);

		KmerSizer::set(oldKmerLength);
	}

	class FilterResults {
	public:
		KM::ValueType value; // the type of the artifact filter hit
		SequenceLengthType minPass, maxPass; // the range of the read to keep
		FilterResults() : value(0), minPass(0), maxPass(MAX_SEQUENCE_LENGTH) {}
		FilterResults(KM::ValueType &_v, SequenceLengthType &_min, SequenceLengthType &_max) :
			value(_v), minPass(_min), maxPass(_max) {}
	};
	class Recorder {
	public:
		OfstreamMap *omPhiX;
		OfstreamMap *omArtifact;
		SequenceCounts readCounts;
		SequenceCounts discardedCounts;
		BaseCounts baseCounts;
		float minimumReadLength;
		Recorder(const FilterKnownOddities &filter) : omPhiX(NULL), omArtifact(NULL), minimumReadLength(0) {
			const ReadSet &reads = filter.getSequences();
			readCounts.resize(reads.getSize() + 1);
			discardedCounts.resize(reads.getSize() + 1);
			baseCounts.resize(reads.getSize() + 1);
		}
		void recordDiscard(int value, Read &read, std::ostream *os, std::string label = "") {
			read.discard();
#pragma omp atomic
			discardedCounts[value]++;
			SequenceLengthType len = read.getLength();
#pragma omp atomic
			baseCounts[value] += len;
			if (os != NULL) {
				FilterKnownOddities::_writeFilterRead(*os, read, 0, len, label);
			}
		}
		void recordTrim(int value, Read &read, SequenceLengthType minPass, SequenceLengthType maxPass) {
			assert(maxPass > minPass);
			assert(!read.isDiscarded());
			SequenceLengthType seqLen = read.getLength();
			assert(maxPass - minPass < seqLen);

			LOG_DEBUG(3, "Trim read " << read.getName() << " " << minPass << "-" << maxPass);
			read = read.getTrimRead(minPass, maxPass - minPass, "AFTrim:" + boost::lexical_cast<string>(minPass) + "+" + boost::lexical_cast<string>(maxPass-minPass));

			#pragma omp atomic
			readCounts[value]++;
#pragma omp atomic
			baseCounts[value] += (seqLen - (maxPass - minPass));
		}
		long getDiscardedReads() const {
			long reads = 0;
			for(unsigned int i = 0; i < discardedCounts.size(); i++)
				reads += discardedCounts[i];
			return reads;
		}
		long getTrimmedReads() const {
			long reads = 0;
			for(unsigned int i = 0; i < readCounts.size(); i++)
				reads += readCounts[i];
			return reads;
		}
		long getTrimmedBases() const {
			long reads = 0;
			for(unsigned int i = 0; i < baseCounts.size(); i++)
				reads += baseCounts[i];
			return reads;
		}
	};

	FilterResults applyFilterToPair(ReadSet &reads, long pairIdx, Recorder &recorder) {
		ReadSet::Pair &pair = reads.getPair(pairIdx);
		FilterResults results1, results2;

		bool isRead1 = reads.isValidRead(pair.read1);
		bool isRead2 = reads.isValidRead(pair.read2);
		if (isRead1)
			results1 = applyFilterToRead(reads, pair.read1, recorder);
		if (isRead2)
			results2 = applyFilterToRead(reads, pair.read2, recorder);

		FilterResults results;
		if (isRead1 && isRead2) {
			recordAffectedRead(reads, recorder, results1, pair.read1, results2, pair.read2);
			int len1 = results1.maxPass - results1.minPass, len2 = results2.maxPass - results2.minPass;
			if (len1 > len2) {
				results = results1;
			} else {
				results = results2;
			}
			results.value = (isPhiX(results1.value) || isPhiX(results2.value)) ? getPhiXReadIdx() : std::max(results1.value, results2.value);
		} else {
			if (isRead1) {
				recordAffectedRead(reads, recorder, results1, pair.read1);
				results = results1;
			} else {
				recordAffectedRead(reads, recorder, results2, pair.read2);
				results = results2;
			}
		}
		return results;
	}


	FilterResults applyFilterToRead(ReadSet &reads, long readIdx, Recorder &recorder, bool recordEffects = false) {

		Read &read = reads.getRead(readIdx);
		FilterResults results;
		KM::ValueType &value = results.value;
		SequenceLengthType &minPass = results.minPass;
		SequenceLengthType &maxPass = results.maxPass;

		LOG_DEBUG(5, "Checking: " << read.getName());

		SequenceLengthType seqLen = read.getLength();
		long bytes = read.getTwoBitEncodingSequenceLength();
		value = 0;
		minPass = 0;
		maxPass = seqLen;

		TwoBitEncoding *ptr = read.getTwoBitSequence();

		// Validate quality scores.  Find the two best ranges of data.
		std::string quals = read.getQuals();
		SequenceLengthType qualsize = quals.size();
		assert(qualsize == seqLen);
		std::pair<long, long> best(0, 0), secondBest(0,0), test(0,0);

		char minQual = Read::FASTQ_START_CHAR + GeneralOptions::getOptions().getMinQuality();
		for(SequenceLengthType i = 0 ; i < qualsize ; i++) {
			test.second = i;
			if (quals[i] < minQual) {
				if (test.second - test.first > best.second - best.first) {
					std::swap(best, test);
				}
				if (test.second - test.first > secondBest.second - secondBest.first) {
					std::swap(secondBest, test);
				}
				test.first = test.second = i + 1;

				LOG_DEBUG(6, "QualityFilter(" << (int) minQual << ") hit baseIdx:" << i << "(" << quals[i] << ")\tmin: " << minPass <<  " max: " << maxPass << "\t" << read.toString());
			}
		}
		test.second = qualsize;
		if (test.second - test.first > best.second - best.first) {
			std::swap(best, test);
		}
		if (test.second - test.first > secondBest.second - secondBest.first) {
			std::swap(secondBest, test);
		}
		if (best.second > best.first) {
			minPass = best.first;
			maxPass = best.second;
		} else {
			minPass = 0;
			maxPass = 0;
		}
		LOG_DEBUG(3, "applyFilterToRead(): Trim: " << read.getName() << " best: " << best.first << "," << best.second << " second: " << secondBest.first << "," << secondBest.second);

		// Now find matches to known artifact reference

		long byteHops = ((maxPass+3)/4) - twoBitLength - ((seqLen & 0x03) == 0 ? 0 : 1);
		if (byteHops < 0 || byteHops > bytes)
			byteHops = 0;

		KM::ElementType elem;
		bool wasPhiX = false;
		KM::BucketType kmers;
		TEMP_KMER(leastKmer);

		LOG_DEBUG(3, "applyFilterToRad(): minPass: " << minPass << " maxPass: " << maxPass << " byteHops: " << byteHops << " twoBitLength: " << twoBitLength);
		SequenceLengthType minAffected = maxPass, maxAffected = minPass;
		for(long byteHop = minPass/4; byteHop <= byteHops; byteHop++) {

			const Kmer &fwd = (const Kmer&) *ptr;
			fwd.buildLeastComplement(leastKmer);

			elem = filter.getElementIfExists( leastKmer );
			if (elem.isValid()) {
				SequenceLengthType pos = byteHop*4;
				value = elem.value();
				wasPhiX |= isPhiX(value);
				if (minAffected > pos)
					minAffected = pos;
				if (maxAffected < pos+length)
					maxAffected = pos+length;
			}

			if (numErrors > 0) {
				KM::BucketType::permuteBases(leastKmer, kmers, numErrors, true);
				for(unsigned int i = 0; i < kmers.size(); i++) {
					elem = filter.getElementIfExists( kmers[i] );
					if (elem.isValid()) {
						SequenceLengthType pos = byteHop*4;
						value = elem.value();
						wasPhiX |= isPhiX(value);
						if (minAffected > pos)
							minAffected = pos;
						if (maxAffected < pos+length)
							maxAffected = pos+length;
					}
				}
			}

			ptr++;
		}
		LOG_DEBUG(3, "After artifact match: " << minAffected << " " << maxAffected << " value: " << value);

		if (wasPhiX) {
			value = getPhiXReadIdx();
		} else if (isSimpleRepeat(value)) {
			// allow simple repeats in the middle of a read with good edges
			bool isGoodMargin = true;
			if ((minAffected - minPass) < (long) 3*length/2)
				isGoodMargin = false;
			if ((maxPass - maxAffected) < (long) 3*length/2)
				isGoodMargin = false;

			if (isGoodMargin) {
				LOG_DEBUG(2, "Allowing simple repeat in middle:" << minAffected << "-" << maxAffected << "!  minPass: " << minPass << ", maxPass: " << maxPass << "\t" << read.toString());
				value = 0;
				minAffected = maxPass;
				maxAffected = minPass;
			}
		}

		if (value > 0 && minAffected <= maxAffected) {
			// trim out matches to artifacts

			if ((minAffected - minPass) >= (maxPass - maxAffected)) {
				// keep left side
				maxPass = minAffected;
			} else {
				// keep right side
				minPass = maxAffected;
			}
		}

		if (value == 0 && (maxPass - minPass) != seqLen) {
			value = sequences.getSize();
			if (ReadSelectorUtil::passesLength(secondBest.second - secondBest.first, seqLen, recorder.minimumReadLength)) {
				// there was only quality trim artifacts.
				// rescue the remaining read, if it is long enough
				int len = secondBest.second - secondBest.first;
				Read remnant = read.getTrimRead(secondBest.first, len, "AFTrim:" + boost::lexical_cast<string>(secondBest.first) + "+" + boost::lexical_cast<string>(len), "-qtrim");
				remnantReads[ omp_get_thread_num() ].append(remnant);
				LOG_DEBUG(1, "applyFilterToRead(): split read: " << read.getNameAndComment() << " to " << remnant.getNameAndComment());
			}
		}
		if (value > 0) {
			if (recordEffects) {
				recordAffectedRead(reads, recorder, results, readIdx);
				// trim the read now
			}
		}
		return results;
	}

	std::string getFilterName(KM::ValueType value) {
		assert(value > 0);
		if (value < sequences.getSize()) {
			return sequences.getRead(value).getName();
		} else {
			return std::string("MinQualityTrim" + boost::lexical_cast<std::string>(GeneralOptions::getOptions().getMinQuality()));
		}
	}
	void recordAffectedRead(ReadSet &reads, Recorder &recorder, FilterResults results1, ReadSet::ReadSetSizeType readIdx1, FilterResults results2 = FilterResults(), ReadSet::ReadSetSizeType readIdx2 = ReadSet::MAX_READ_IDX) {
		bool isRead1 = reads.isValidRead(readIdx1);
		bool isRead2 = reads.isValidRead(readIdx2);

		bool wasAffected = (results1.value != 0) | (results2.value != 0);
		bool wasPhiX = isPhiX(results1.value) | isPhiX(results2.value);
		bool wasReference = ((results1.value != sequences.getSize()) & isReference(results1.value))
						| ((results2.value != sequences.getSize()) & isReference(results2.value));

		if (wasAffected) {

			if (wasPhiX) {
				results1.value = results2.value = getPhiXReadIdx();
			}

			bool isRead1Affected = isRead1 && results1.value != 0;
			bool isRead2Affected = isRead2 && results2.value != 0;

			if (wasPhiX) {

				ostream *os = NULL;
				std::ostringstream *ss = NULL;
				if (recorder.omPhiX != NULL) {
					std::string fileSuffix = std::string("-") + reads.getReadFileNamePrefix(readIdx1);
					os = & (recorder.omPhiX->getOfstream( fileSuffix ) );
					ss = new std::ostringstream();
				}

				// always discard the read, as it contains some PhiX even if not output to a discard file
				if (isRead1) {
					Read &read = reads.getRead(readIdx1);
					recorder.recordDiscard(results1.value, read, ss);
				}
				if (isRead2) {
					Read &read = reads.getRead(readIdx2);
					recorder.recordDiscard(results2.value, read, ss);
				}

				if (os != NULL && ss->tellp() > 0)
				{
					std::string str = ss->str();
#pragma omp critical (writePhix)
					{ *os << str; }
					delete ss;
				}

			} else {

				ostream *os = NULL;
				std::ostringstream *ss = NULL;
				std::string label1, label2;
				if (recorder.omArtifact != NULL) {
					std::string fileSuffix = std::string("-") + reads.getReadFileNamePrefix( isRead1 ? readIdx1 : readIdx2);

					if (isRead1Affected)
						label1 = getFilterName(results1.value);
					if (isRead2Affected)
						label2 = getFilterName(results2.value);
					os = & (recorder.omArtifact->getOfstream( fileSuffix ) );
					ss = new std::ostringstream();
				}

				if (isRead1 && results1.value != 0) {
					Read &read = reads.getRead(readIdx1);
					SequenceLengthType len = read.getLength();
					int passLength = results1.maxPass - results1.minPass;
					if (wasReference || passLength <= 0 || !ReadSelectorUtil::passesLength(passLength, len, recorder.minimumReadLength)) {
						recorder.recordDiscard(results1.value, read, ss, label1);
					} else {
						recorder.recordTrim(results1.value, read, results1.minPass, results1.maxPass);
					}
				}
				if (isRead2 && results2.value != 0) {
					Read &read = reads.getRead(readIdx2);
					SequenceLengthType len = read.getLength();
					int passLength = results2.maxPass - results2.minPass;
					if (wasReference || passLength <= 0 || !ReadSelectorUtil::passesLength(passLength, len, recorder.minimumReadLength)) {
						recorder.recordDiscard(results2.value, read, ss, label2);
					} else {
						recorder.recordTrim(results2.value, read, results2.minPass, results2.maxPass);
					}
				}
				if (os != NULL && ss->tellp() > 0) {
					std::string str = ss->str();
#pragma omp critical (writeFilter)
					{ *os << str; }
					delete ss;
				}

			}

			if (Log::isDebug(5)){
				Read read;
				if (isRead1Affected) {
					read = reads.getRead(readIdx1);
					LOG_DEBUG(5, "FilterMatch1 to" << read.getName() << " "
							<< read.getFastaNoMarkup() << "" << read.getFasta() << " "
							<< wasPhiX << "" << results1.minPass << "-" << results1.maxPass
							<< ":" << getFilterName(results1.value));
				}
				if (isRead2Affected) {
					read = reads.getRead(readIdx2);
					LOG_DEBUG(5, "FilterMatch2 to" << read.getName() << " "
							<< read.getFastaNoMarkup() << "" << read.getFasta() << " "
							<< wasPhiX << " "  << results2.minPass << "-" << results2.maxPass
							<< ":"<< getFilterName(results2.value));

				}
			}
		}
	}

	unsigned long applyFilter(ReadSet &reads) {
		unsigned long oldKmerLength = KmerSizer::getSequenceLength();
		KmerSizer::set(length);

		Recorder recorder( *this );
		unsigned long affectedCount = 0;

		OfstreamMap _omPhiX(Options::getOptions().getOutputFile(), "-PhiX.fastq");
		if (FilterKnownOdditiesOptions::getOptions().getPhiXOutput()) {
			recorder.omPhiX = &_omPhiX;
		}
		OfstreamMap _omArtifact(Options::getOptions().getOutputFile(), "-Artifact.fastq");
		if (FilterKnownOdditiesOptions::getOptions().getFilterOutput()) {
			recorder.omArtifact = &_omArtifact;
		}

		float &minimumReadLength = recorder.minimumReadLength;
		minimumReadLength = ReadSelectorOptions::getOptions().getMinReadLength();

		bool byPair = reads.hasPairs();
		long size = reads.getSize();
		if (byPair)
			size = reads.getPairSize();

		LOG_VERBOSE(1, "Applying Artifact filter to " << size << " " << (byPair?"Pairs":"Reads"));
#pragma omp parallel for schedule(guided)
		for (long idx = 0; idx < size; idx++) {
			LOG_DEBUG(5, "filtering read" << (byPair ? "pair " : "" ) << idx);

			if (byPair)
				applyFilterToPair(reads, idx, recorder);
			else {
				applyFilterToRead(reads, idx, recorder, true);
			}
		}

		long remnants = 0;
		for(int i = 0; i < (int) remnantReads.size(); i++) {
			long s = remnantReads[i].getSize();
			remnants += s;
			if (s > 0) {
				reads.append(remnantReads[i]);
			}
		}
		LOG_VERBOSE_OPTIONAL(1, remnants > 0, "Rescued " << remnants << " reads from poor quality scores in the middle");

		for(ReadSetSizeType j = 0 ; j < recorder.readCounts.size() ; j++)
			affectedCount += recorder.readCounts[j];

		if (Log::isVerbose(1)) {
			std::stringstream ss ;
			ss << "Final Filter Matches to reads:" << affectedCount << std::endl;
			ss << "Discarded Reads:" << recorder.getDiscardedReads() << std::endl;
			ss << "Trimmed Reads:" << recorder.getTrimmedReads() << std::endl;
			ss << "Discarded/Trimmed Bases:" << recorder.getTrimmedBases() << std::endl;
			ss << std::endl;
			ss << "\tDiscarded\tAffected\tBasesRemoved\tMatch" << std::endl;
			// TODO sort
			for(unsigned long idx = 0 ; idx < recorder.readCounts.size(); idx++) {
				if (recorder.baseCounts[idx] > 0) {
					ss << "\t" <<  recorder.discardedCounts[idx] << "\t" << recorder.readCounts[idx] << "\t" << recorder.baseCounts[idx] << "\t" << getFilterName(idx) << std::endl;
				}
			}
			std::string s = ss.str();
			LOG_VERBOSE(1, s);
		}

		KmerSizer::set(oldKmerLength);
		return affectedCount;
	}

	static void _writeFilterRead(ostream &os, Read  &read, SequenceLengthType readOffset, SequenceLengthType readLength, std::string readLabel = "") {
		read.write(os, readOffset, readLength, readLabel, FormatOutput::FastqUnmasked());
	}

	const ReadSet &getSequences() const {
		return sequences;
	}

	static std::string getArtifactFasta() {
		std::stringstream ss;
		ss << ">PrimerDimer" << std::endl;
		ss << "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG" << std::endl;
		ss << ">RNA_Linker" << std::endl;
		ss << "ATCTCGTATGCCGTCTTCTGCTTGATCTCGTATGCCGTCTTCTGCTTG" << std::endl;
		ss << ">Homopolymer-A" << std::endl;
		ss << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" << std::endl;
		ss << ">Homopolymer-C" << std::endl;
		ss << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << std::endl;

		// from TagDust Lassmann T., et al. (2009) TagDust - A program to eliminate artifacts from next generation sequencing data. Bioinformatics.
		ss << ">Solexa_5_prime_adapter" << std::endl;
		ss << "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAG" << std::endl;
		ss << ">Solexa_3_prime_adapter" << std::endl;
		ss << "TTTTCGTATGCCGTCTTCTGCTTG" << std::endl;
		ss << ">Gex_Adapter_1" << std::endl;
		ss << "GATCGTCGGACTGTAGAACTCTGAAC" << std::endl;
		ss << ">Gex_Adapter_1_2" << std::endl;
		ss << "ACAGGTTCAGAGTTCTACAGTCCGAC" << std::endl;
		ss << ">Gex_Adapter_2" << std::endl;
		ss << "CAAGCAGAAGACGGCATACGANN" << std::endl;
		ss << ">Gex_Adapter_2_2" << std::endl;
		ss << "TCGTATGCCGTCTTCTGCTTG" << std::endl;
		ss << ">Gex_PCR_Primer_1" << std::endl;
		ss << "CAAGCAGAAGACGGCATACGA" << std::endl;
		ss << ">Gex_PCR_Primer_2" << std::endl;
		ss << "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA" << std::endl;
		ss << ">Gex_Sequencing_Primer" << std::endl;
		ss << "CGACAGGTTCAGAGTTCTACAGTCCGACGATC" << std::endl;
		ss << ">Adapters1" << std::endl;
		ss << "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG" << std::endl;
		ss << ">Adapters1_1" << std::endl;
		ss << "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" << std::endl;
		ss << ">PCR_Primers1" << std::endl;
		ss << "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" << std::endl;
		ss << ">PCR Primers1_1" << std::endl;
		ss << "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT" << std::endl;
		ss << ">Genomic_DNA_Sequencing_Primer" << std::endl;
		ss << "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" << std::endl;
		ss << ">PE_Adapters1" << std::endl;
		ss << "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG" << std::endl;
		ss << ">PE_Adapters1_" << std::endl;
		ss << "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" << std::endl;
		ss << ">PE_PCR_Primers1" << std::endl;
		ss << "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT" << std::endl;
		ss << ">PE_PCR_Primers1_1" << std::endl;
		ss << "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT" << std::endl;
		ss << ">PE_Sequencing_Primer" << std::endl;
		ss << "ACACTCTTTCCCTACACGACGCTCTTCCGATCT" << std::endl;
		ss << ">PE_Sequencing_Primer_1" << std::endl;
		ss << "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT" << std::endl;
		return ss.str();
	}


	static inline unsigned int &getSimpleRepeatBegin() {
		static unsigned int simpleRepeatBegin = -1;
		return simpleRepeatBegin;
	}
	static inline unsigned int &getSimpleRepeatEnd() {
		static unsigned int getSimpleRepeatEnd = -1;
		return getSimpleRepeatEnd;
	}

	static inline bool isSimpleRepeat(unsigned int readIdx) {
		return (readIdx >= getSimpleRepeatBegin() && readIdx < getSimpleRepeatEnd());
	}

	static std::string getSimpleRepeatFasta() {
		std::stringstream ss;
		ss << ">(CA)n#Simple_repeat" << std::endl;
		ss << "CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA" << std::endl;
		ss << "CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA" << std::endl;
		ss << ">(CAA)n#Simple_repeat" << std::endl;
		ss << "CAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAA" << std::endl;
		ss << "CAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAACAA" << std::endl;
		ss << ">(CAAA)n#Simple_repeat" << std::endl;
		ss << "CAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAA" << std::endl;
		ss << "CAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAACAAA" << std::endl;
		ss << ">(CAAAA)n#Simple_repeat" << std::endl;
		ss << "CAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAA" << std::endl;
		ss << "CAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAACAAAA" << std::endl;
		ss << ">(CAAAAA)n#Simple_repeat" << std::endl;
		ss << "CAAAAACAAAAACAAAAACAAAAACAAAAACAAAAACAAAAACAAAAACAAAAACAAAAA" << std::endl;
		ss << "CAAAAACAAAAACAAAAACAAAAACAAAAACAAAAACAAAAACAAAAACAAAAACAAAAA" << std::endl;
		ss << ">(CAAAC)n#Simple_repeat" << std::endl;
		ss << "CAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAAC" << std::endl;
		ss << "CAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAACCAAAC" << std::endl;
		ss << ">(CAAAG)n#Simple_repeat" << std::endl;
		ss << "CAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAG" << std::endl;
		ss << "CAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAGCAAAG" << std::endl;
		ss << ">(CAAAT)n#Simple_repeat" << std::endl;
		ss << "CAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAAT" << std::endl;
		ss << "CAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAATCAAAT" << std::endl;
		ss << ">(CAACC)n#Simple_repeat" << std::endl;
		ss << "CAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACC" << std::endl;
		ss << "CAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACCCAACC" << std::endl;
		ss << ">(CAACG)n#Simple_repeat" << std::endl;
		ss << "CAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACG" << std::endl;
		ss << "CAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACGCAACG" << std::endl;
		ss << ">(CAACT)n#Simple_repeat" << std::endl;
		ss << "CAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACT" << std::endl;
		ss << "CAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACTCAACT" << std::endl;
		ss << ">(CAAG)n#Simple_repeat" << std::endl;
		ss << "CAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAG" << std::endl;
		ss << "CAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAGCAAG" << std::endl;
		ss << ">(CAAGA)n#Simple_repeat" << std::endl;
		ss << "CAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGA" << std::endl;
		ss << "CAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGACAAGA" << std::endl;
		ss << ">(CAAGC)n#Simple_repeat" << std::endl;
		ss << "CAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGC" << std::endl;
		ss << "CAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGCCAAGC" << std::endl;
		ss << ">(CAAGG)n#Simple_repeat" << std::endl;
		ss << "CAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGG" << std::endl;
		ss << "CAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGGCAAGG" << std::endl;
		ss << ">(CAAGT)n#Simple_repeat" << std::endl;
		ss << "CAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGT" << std::endl;
		ss << "CAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGTCAAGT" << std::endl;
		ss << ">(CAAT)n#Simple_repeat" << std::endl;
		ss << "CAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAAT" << std::endl;
		ss << "CAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAATCAAT" << std::endl;
		ss << ">(CAATA)n#Simple_repeat" << std::endl;
		ss << "CAATACAATACAATACAATACAATACAATACAATACAATACAATACAATACAATACAATA" << std::endl;
		ss << "CAATACAATACAATACAATACAATACAATACAATACAATACAATACAATACAATACAATA" << std::endl;
		ss << ">(CAATC)n#Simple_repeat" << std::endl;
		ss << "CAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATC" << std::endl;
		ss << "CAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATCCAATC" << std::endl;
		ss << ">(CAATG)n#Simple_repeat" << std::endl;
		ss << "CAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATG" << std::endl;
		ss << "CAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATGCAATG" << std::endl;
		ss << ">(CAATT)n#Simple_repeat" << std::endl;
		ss << "CAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATT" << std::endl;
		ss << "CAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATTCAATT" << std::endl;
		ss << ">(CACAA)n#Simple_repeat" << std::endl;
		ss << "CACAACACAACACAACACAACACAACACAACACAACACAACACAACACAACACAACACAA" << std::endl;
		ss << "CACAACACAACACAACACAACACAACACAACACAACACAACACAACACAACACAACACAA" << std::endl;
		ss << ">(CACAC)n#Simple_repeat" << std::endl;
		ss << "CACACCACACCACACCACACCACACCACACCACACCACACCACACCACACCACACCACAC" << std::endl;
		ss << "CACACCACACCACACCACACCACACCACACCACACCACACCACACCACACCACACCACAC" << std::endl;
		ss << ">(CACAG)n#Simple_repeat" << std::endl;
		ss << "CACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAG" << std::endl;
		ss << "CACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAGCACAG" << std::endl;
		ss << ">(CACAT)n#Simple_repeat" << std::endl;
		ss << "CACATCACATCACATCACATCACATCACATCACATCACATCACATCACATCACATCACAT" << std::endl;
		ss << "CACATCACATCACATCACATCACATCACATCACATCACATCACATCACATCACATCACAT" << std::endl;
		ss << ">(CACCC)n#Simple_repeat" << std::endl;
		ss << "CACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCC" << std::endl;
		ss << "CACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCC" << std::endl;
		ss << ">(CACCG)n#Simple_repeat" << std::endl;
		ss << "CACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCG" << std::endl;
		ss << "CACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCGCACCG" << std::endl;
		ss << ">(CACCT)n#Simple_repeat" << std::endl;
		ss << "CACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCT" << std::endl;
		ss << "CACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCTCACCT" << std::endl;
		ss << ">(CACG)n#Simple_repeat" << std::endl;
		ss << "CACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACG" << std::endl;
		ss << "CACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACGCACG" << std::endl;
		ss << ">(CACGA)n#Simple_repeat" << std::endl;
		ss << "CACGACACGACACGACACGACACGACACGACACGACACGACACGACACGACACGACACGA" << std::endl;
		ss << "CACGACACGACACGACACGACACGACACGACACGACACGACACGACACGACACGACACGA" << std::endl;
		ss << ">(CACGC)n#Simple_repeat" << std::endl;
		ss << "CACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGC" << std::endl;
		ss << "CACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGC" << std::endl;
		ss << ">(CACGT)n#Simple_repeat" << std::endl;
		ss << "CACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGC" << std::endl;
		ss << "CACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGCCACGC" << std::endl;
		ss << ">(CACTA)n#Simple_repeat" << std::endl;
		ss << "CACTACACTACACTACACTACACTACACTACACTACACTACACTACACTACACTACACTA" << std::endl;
		ss << "CACTACACTACACTACACTACACTACACTACACTACACTACACTACACTACACTACACTA" << std::endl;
		ss << ">(CACTC)n#Simple_repeat" << std::endl;
		ss << "CACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTC" << std::endl;
		ss << "CACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTCCACTC" << std::endl;
		ss << ">(CACTG)n#Simple_repeat" << std::endl;
		ss << "CACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTG" << std::endl;
		ss << "CACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTGCACTG" << std::endl;
		ss << ">(CACTT)n#Simple_repeat" << std::endl;
		ss << "CACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTT" << std::endl;
		ss << "CACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTTCACTT" << std::endl;
		ss << ">(CAG)n#Simple_repeat" << std::endl;
		ss << "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG" << std::endl;
		ss << "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG" << std::endl;
		ss << ">(CAGA)n#Simple_repeat" << std::endl;
		ss << "CAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGA" << std::endl;
		ss << "CAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGACAGA" << std::endl;
		ss << ">(CAGAA)n#Simple_repeat" << std::endl;
		ss << "CAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAA" << std::endl;
		ss << "CAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAACAGAA" << std::endl;
		ss << ">(CAGAC)n#Simple_repeat" << std::endl;
		ss << "CAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGAC" << std::endl;
		ss << "CAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGACCAGAC" << std::endl;
		ss << ">(CAGAG)n#Simple_repeat" << std::endl;
		ss << "CAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAG" << std::endl;
		ss << "CAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAGCAGAG" << std::endl;
		ss << ">(CAGAT)n#Simple_repeat" << std::endl;
		ss << "CAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGAT" << std::endl;
		ss << "CAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGATCAGAT" << std::endl;
		ss << ">(CAGC)n#Simple_repeat" << std::endl;
		ss << "CAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGC" << std::endl;
		ss << "CAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGC" << std::endl;
		ss << ">(CAGCC)n#Simple_repeat" << std::endl;
		ss << "CAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCC" << std::endl;
		ss << "CAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCCCAGCC" << std::endl;
		ss << ">(CAGCG)n#Simple_repeat" << std::endl;
		ss << "CAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCG" << std::endl;
		ss << "CAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCGCAGCG" << std::endl;
		ss << ">(CAGCT)n#Simple_repeat" << std::endl;
		ss << "CAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCT" << std::endl;
		ss << "CAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCTCAGCT" << std::endl;
		ss << ">(CAGG)n#Simple_repeat" << std::endl;
		ss << "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGG" << std::endl;
		ss << "CAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGGCAGG" << std::endl;
		ss << ">(CAGGA)n#Simple_repeat" << std::endl;
		ss << "CAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGA" << std::endl;
		ss << "CAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGACAGGA" << std::endl;
		ss << ">(CAGGC)n#Simple_repeat" << std::endl;
		ss << "CAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGC" << std::endl;
		ss << "CAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGGC" << std::endl;
		ss << ">(CAGGG)n#Simple_repeat" << std::endl;
		ss << "CAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGG" << std::endl;
		ss << "CAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGGCAGGG" << std::endl;
		ss << ">(CAGGT)n#Simple_repeat" << std::endl;
		ss << "CAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGT" << std::endl;
		ss << "CAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGTCAGGT" << std::endl;
		ss << ">(CAGT)n#Simple_repeat" << std::endl;
		ss << "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT" << std::endl;
		ss << "CAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGT" << std::endl;
		ss << ">(CAGTA)n#Simple_repeat" << std::endl;
		ss << "CAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTA" << std::endl;
		ss << "CAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTACAGTA" << std::endl;
		ss << ">(CAGTC)n#Simple_repeat" << std::endl;
		ss << "CAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTC" << std::endl;
		ss << "CAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTCCAGTC" << std::endl;
		ss << ">(CAGTT)n#Simple_repeat" << std::endl;
		ss << "CAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTT" << std::endl;
		ss << "CAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTTCAGTT" << std::endl;
		ss << ">(CAT)n#Simple_repeat" << std::endl;
		ss << "CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT" << std::endl;
		ss << "CATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCAT" << std::endl;
		ss << ">(CATA)n#Simple_repeat" << std::endl;
		ss << "CATACATACATACATACATACATACATACATACATACATACATACATACATACATACATA" << std::endl;
		ss << "CATACATACATACATACATACATACATACATACATACATACATACATACATACATACATA" << std::endl;
		ss << ">(CATAA)n#Simple_repeat" << std::endl;
		ss << "CATAACATAACATAACATAACATAACATAACATAACATAACATAACATAACATAACATAA" << std::endl;
		ss << "CATAACATAACATAACATAACATAACATAACATAACATAACATAACATAACATAACATAA" << std::endl;
		ss << ">(CATAC)n#Simple_repeat" << std::endl;
		ss << "CATACCATACCATACCATACCATACCATACCATACCATACCATACCATACCATACCATAC" << std::endl;
		ss << "CATACCATACCATACCATACCATACCATACCATACCATACCATACCATACCATACCATAC" << std::endl;
		ss << ">(CATAG)n#Simple_repeat" << std::endl;
		ss << "CATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAG" << std::endl;
		ss << "CATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAGCATAG" << std::endl;
		ss << ">(CATAT)n#Simple_repeat" << std::endl;
		ss << "CATATCATATCATATCATATCATATCATATCATATCATATCATATCATATCATATCATAT" << std::endl;
		ss << "CATATCATATCATATCATATCATATCATATCATATCATATCATATCATATCATATCATAT" << std::endl;
		ss << ">(CATCC)n#Simple_repeat" << std::endl;
		ss << "CATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCC" << std::endl;
		ss << "CATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCCCATCC" << std::endl;
		ss << ">(CATCG)n#Simple_repeat" << std::endl;
		ss << "CATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCG" << std::endl;
		ss << "CATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCGCATCG" << std::endl;
		ss << ">(CATCT)n#Simple_repeat" << std::endl;
		ss << "CATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCT" << std::endl;
		ss << "CATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCTCATCT" << std::endl;
		ss << ">(CATG)n#Simple_repeat" << std::endl;
		ss << "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG" << std::endl;
		ss << "CATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG" << std::endl;
		ss << ">(CATGC)n#Simple_repeat" << std::endl;
		ss << "CATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGC" << std::endl;
		ss << "CATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGCCATGC" << std::endl;
		ss << ">(CATGT)n#Simple_repeat" << std::endl;
		ss << "CATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGT" << std::endl;
		ss << "CATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGTCATGT" << std::endl;
		ss << ">(CATTA)n#Simple_repeat" << std::endl;
		ss << "CATTACATTACATTACATTACATTACATTACATTACATTACATTACATTACATTACATTA" << std::endl;
		ss << "CATTACATTACATTACATTACATTACATTACATTACATTACATTACATTACATTACATTA" << std::endl;
		ss << ">(CATTC)n#Simple_repeat" << std::endl;
		ss << "CATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTC" << std::endl;
		ss << "CATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTC" << std::endl;
		ss << ">(CATTT)n#Simple_repeat" << std::endl;
		ss << "CATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTT" << std::endl;
		ss << "CATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTTCATTT" << std::endl;
		ss << ">(CCAA)n#Simple_repeat" << std::endl;
		ss << "CAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAAC" << std::endl;
		ss << "CAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAACCAAC" << std::endl;
		ss << ">(CCCA)n#Simple_repeat" << std::endl;
		ss << "CCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCA" << std::endl;
		ss << "CCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCACCCA" << std::endl;
		ss << ">(CCCCG)n#Simple_repeat" << std::endl;
		ss << "CCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCG" << std::endl;
		ss << "CCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCGCCCCG" << std::endl;
		ss << ">(CCCGA)n#Simple_repeat" << std::endl;
		ss << "CCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGA" << std::endl;
		ss << "CCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGACCCGA" << std::endl;
		ss << ">(CCGCG)n#Simple_repeat" << std::endl;
		ss << "CCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCG" << std::endl;
		ss << "CCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCGCCGCG" << std::endl;
		ss << ">(CCGCT)n#Simple_repeat" << std::endl;
		ss << "CCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCT" << std::endl;
		ss << "CCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCTCCGCT" << std::endl;
		ss << ">(CCGG)n#Simple_repeat" << std::endl;
		ss << "CGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGC" << std::endl;
		ss << "CGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGC" << std::endl;
		ss << ">(CCGAA)n#Simple_repeat" << std::endl;
		ss << "CCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAA" << std::endl;
		ss << "CCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAACCGAA" << std::endl;
		ss << ">(CCGAG)n#Simple_repeat" << std::endl;
		ss << "CCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAG" << std::endl;
		ss << "CCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAGCCGAG" << std::endl;
		ss << ">(CCGGA)n#Simple_repeat" << std::endl;
		ss << "CCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGA" << std::endl;
		ss << "CCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGACCGGA" << std::endl;
		ss << ">(CCGGG)n#Simple_repeat" << std::endl;
		ss << "CCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGG" << std::endl;
		ss << "CCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGGCCGGG" << std::endl;
		ss << ">(CCGTA)n#Simple_repeat" << std::endl;
		ss << "CCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTA" << std::endl;
		ss << "CCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTACCGTA" << std::endl;
		ss << ">(CCTAA)n#Simple_repeat" << std::endl;
		ss << "CCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAA" << std::endl;
		ss << "CCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAA" << std::endl;
		ss << ">(CCTAG)n#Simple_repeat" << std::endl;
		ss << "CCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAG" << std::endl;
		ss << "CCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAGCCTAG" << std::endl;
		ss << ">(CCTAT)n#Simple_repeat" << std::endl;
		ss << "CCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTAT" << std::endl;
		ss << "CCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTATCCTAT" << std::endl;
		ss << ">(CCTCG)n#Simple_repeat" << std::endl;
		ss << "CCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCG" << std::endl;
		ss << "CCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCGCCTCG" << std::endl;
		ss << ">(CG)n#Simple_repeat" << std::endl;
		ss << "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG" << std::endl;
		ss << "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG" << std::endl;
		ss << ">(CGA)n#Simple_repeat" << std::endl;
		ss << "CGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGA" << std::endl;
		ss << "CGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGA" << std::endl;
		ss << ">(CGAA)n#Simple_repeat" << std::endl;
		ss << "CGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAA" << std::endl;
		ss << "CGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAACGAA" << std::endl;
		ss << ">(CGAAA)n#Simple_repeat" << std::endl;
		ss << "CGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAA" << std::endl;
		ss << "CGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAACGAAA" << std::endl;
		ss << ">(CGAAG)n#Simple_repeat" << std::endl;
		ss << "CGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAG" << std::endl;
		ss << "CGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAGCGAAG" << std::endl;
		ss << ">(CGAAT)n#Simple_repeat" << std::endl;
		ss << "CGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAAT" << std::endl;
		ss << "CGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAATCGAAT" << std::endl;
		ss << ">(CGACG)n#Simple_repeat" << std::endl;
		ss << "CGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACG" << std::endl;
		ss << "CGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACGCGACG" << std::endl;
		ss << ">(CGAG)n#Simple_repeat" << std::endl;
		ss << "CGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAG" << std::endl;
		ss << "CGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAGCGAG" << std::endl;
		ss << ">(CGAGA)n#Simple_repeat" << std::endl;
		ss << "CGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGA" << std::endl;
		ss << "CGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGACGAGA" << std::endl;
		ss << ">(CGAGG)n#Simple_repeat" << std::endl;
		ss << "CGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGG" << std::endl;
		ss << "CGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGGCGAGG" << std::endl;
		ss << ">(CGAGT)n#Simple_repeat" << std::endl;
		ss << "CGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGT" << std::endl;
		ss << "CGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGTCGAGT" << std::endl;
		ss << ">(CGAT)n#Simple_repeat" << std::endl;
		ss << "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT" << std::endl;
		ss << "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT" << std::endl;
		ss << ">(CGATA)n#Simple_repeat" << std::endl;
		ss << "CGATACGATACGATACGATACGATACGATACGATACGATACGATACGATACGATACGATA" << std::endl;
		ss << "CGATACGATACGATACGATACGATACGATACGATACGATACGATACGATACGATACGATA" << std::endl;
		ss << ">(CGG)n#Simple_repeat" << std::endl;
		ss << "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG" << std::endl;
		ss << "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG" << std::endl;
		ss << ">(CGGA)n#Simple_repeat" << std::endl;
		ss << "CGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGA" << std::endl;
		ss << "CGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGACGGA" << std::endl;
		ss << ">(CGGG)n#Simple_repeat" << std::endl;
		ss << "CGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGG" << std::endl;
		ss << "CGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGGCGGG" << std::endl;
		ss << ">(CGGAA)n#Simple_repeat" << std::endl;
		ss << "CGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAA" << std::endl;
		ss << "CGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAACGGAA" << std::endl;
		ss << ">(CGGAG)n#Simple_repeat" << std::endl;
		ss << "CGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAG" << std::endl;
		ss << "CGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAGCGGAG" << std::endl;
		ss << ">(CGGGA)n#Simple_repeat" << std::endl;
		ss << "CGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGA" << std::endl;
		ss << "CGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGACGGGA" << std::endl;
		ss << ">(CGGT)n#Simple_repeat" << std::endl;
		ss << "CGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGT" << std::endl;
		ss << "CGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGTCGGT" << std::endl;
		ss << ">(CGTAA)n#Simple_repeat" << std::endl;
		ss << "CGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAA" << std::endl;
		ss << "CGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAACGTAA" << std::endl;
		ss << ">(CGTAG)n#Simple_repeat" << std::endl;
		ss << "CGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAG" << std::endl;
		ss << "CGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAGCGTAG" << std::endl;
		ss << ">(CTAA)n#Simple_repeat" << std::endl;
		ss << "CTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAA" << std::endl;
		ss << "CTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAA" << std::endl;
		ss << ">(CTAG)n#Simple_repeat" << std::endl;
		ss << "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG" << std::endl;
		ss << "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG" << std::endl;
		ss << ">(CTAAA)n#Simple_repeat" << std::endl;
		ss << "CTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAA" << std::endl;
		ss << "CTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAACTAAA" << std::endl;
		ss << ">(CTAAG)n#Simple_repeat" << std::endl;
		ss << "CTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAG" << std::endl;
		ss << "CTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAGCTAAG" << std::endl;
		ss << ">(CTAAT)n#Simple_repeat" << std::endl;
		ss << "CTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAAT" << std::endl;
		ss << "CTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAATCTAAT" << std::endl;
		ss << ">(CTACT)n#Simple_repeat" << std::endl;
		ss << "CTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACT" << std::endl;
		ss << "CTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACTCTACT" << std::endl;
		ss << ">(CTAGG)n#Simple_repeat" << std::endl;
		ss << "CTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGG" << std::endl;
		ss << "CTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGGCTAGG" << std::endl;
		ss << ">(CTAGT)n#Simple_repeat" << std::endl;
		ss << "CTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGT" << std::endl;
		ss << "CTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGTCTAGT" << std::endl;
		ss << ">(CTATA)n#Simple_repeat" << std::endl;
		ss << "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATA" << std::endl;
		ss << "CTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATACTATA" << std::endl;
		ss << ">(CTATT)n#Simple_repeat" << std::endl;
		ss << "CTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATT" << std::endl;
		ss << "CTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATTCTATT" << std::endl;
		ss << ">(CTTAT)n#Simple_repeat" << std::endl;
		ss << "CTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTAT" << std::endl;
		ss << "CTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTATCTTAT" << std::endl;
		ss << ">(CTTTA)n#Simple_repeat" << std::endl;
		ss << "CTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTA" << std::endl;
		ss << "CTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTACTTTA" << std::endl;
		ss << ">(GA)n#Simple_repeat" << std::endl;
		ss << "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA" << std::endl;
		ss << "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA" << std::endl;
		ss << ">(GAA)n#Simple_repeat" << std::endl;
		ss << "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA" << std::endl;
		ss << "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA" << std::endl;
		ss << ">(GAAA)n#Simple_repeat" << std::endl;
		ss << "GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA" << std::endl;
		ss << "GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA" << std::endl;
		ss << ">(GAAAA)n#Simple_repeat" << std::endl;
		ss << "GAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAA" << std::endl;
		ss << "GAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAA" << std::endl;
		ss << ">(GAGAA)n#Simple_repeat" << std::endl;
		ss << "GAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAA" << std::endl;
		ss << "GAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAAGAGAA" << std::endl;
		ss << ">(GGA)n#Simple_repeat" << std::endl;
		ss << "GGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGA" << std::endl;
		ss << "GGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGA" << std::endl;
		ss << ">(GGAA)n#Simple_repeat" << std::endl;
		ss << "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA" << std::endl;
		ss << "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA" << std::endl;
		ss << ">(GGAAA)n#Simple_repeat" << std::endl;
		ss << "GGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAA" << std::endl;
		ss << "GGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAAGGAAA" << std::endl;
		ss << ">(GGAGA)n#Simple_repeat" << std::endl;
		ss << "GGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGA" << std::endl;
		ss << "GGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGAGGAGA" << std::endl;
		ss << ">(GGGA)n#Simple_repeat" << std::endl;
		ss << "GGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGA" << std::endl;
		ss << "GGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGAGGGA" << std::endl;
		ss << ">(GGGAA)n#Simple_repeat" << std::endl;
		ss << "GGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAA" << std::endl;
		ss << "GGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGAA" << std::endl;
		ss << ">(GGGGA)n#Simple_repeat" << std::endl;
		ss << "GGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGA" << std::endl;
		ss << "GGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGAGGGGA" << std::endl;
		ss << ">(TA)n#Simple_repeat" << std::endl;
		ss << "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA" << std::endl;
		ss << "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA" << std::endl;
		ss << ">(TATAA)n#Simple_repeat" << std::endl;
		ss << "TATAATATAATATAATATAATATAATATAATATAATATAATATAATATAATATAATATAA" << std::endl;
		ss << "TATAATATAATATAATATAATATAATATAATATAATATAATATAATATAATATAATATAA" << std::endl;
		ss << ">(TAA)n#Simple_repeat" << std::endl;
		ss << "TAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAA" << std::endl;
		ss << "TAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAATAA" << std::endl;
		ss << ">(TAAA)n#Simple_repeat" << std::endl;
		ss << "TAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAA" << std::endl;
		ss << "TAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAATAAA" << std::endl;
		ss << ">(TAAAA)n#Simple_repeat" << std::endl;
		ss << "TAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAA" << std::endl;
		ss << "TAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAATAAAA" << std::endl;
		ss << ">(TAAG)n#Simple_repeat" << std::endl;
		ss << "TAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAG" << std::endl;
		ss << "TAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAGTAAG" << std::endl;
		ss << ">(TACG)n#Simple_repeat" << std::endl;
		ss << "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG" << std::endl;
		ss << "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG" << std::endl;
		ss << ">(TAG)n#Simple_repeat" << std::endl;
		ss << "TAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAG" << std::endl;
		ss << "TAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAGTAG" << std::endl;
		ss << ">(TAGA)n#Simple_repeat" << std::endl;
		ss << "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA" << std::endl;
		ss << "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA" << std::endl;
		ss << ">(TAGG)n#Simple_repeat" << std::endl;
		ss << "TAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGG" << std::endl;
		ss << "TAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGGTAGG" << std::endl;
		ss << ">(TAGGG)n#Simple_repeat" << std::endl;
		ss << "TAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGG" << std::endl;
		ss << "TAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGGTAGGG" << std::endl;
		ss << ">(TGAA)n#Simple_repeat" << std::endl;
		ss << "TGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAA" << std::endl;
		ss << "TGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAA" << std::endl;
		ss << ">(TGAG)n#Simple_repeat" << std::endl;
		ss << "TGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAG" << std::endl;
		ss << "TGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAGTGAG" << std::endl;
		ss << ">(TGG)n#Simple_repeat" << std::endl;
		ss << "TGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGG" << std::endl;
		ss << "TGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGG" << std::endl;
		ss << ">(TGGA)n#Simple_repeat" << std::endl;
		ss << "TGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGA" << std::endl;
		ss << "TGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGATGGA" << std::endl;
		ss << ">(TTAA)n#Simple_repeat" << std::endl;
		ss << "TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA" << std::endl;
		ss << "TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA" << std::endl;
		ss << ">(TTAAA)n#Simple_repeat" << std::endl;
		ss << "TTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAA" << std::endl;
		ss << "TTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAATTAAA" << std::endl;
		ss << ">(TTAAG)n#Simple_repeat" << std::endl;
		ss << "TTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAG" << std::endl;
		ss << "TTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAGTTAAG" << std::endl;
		ss << ">SUBTEL_sat#Satellite" << std::endl;
		ss << "GCGCCTCTCTGCGCCTGCGCCGGCGCSSCGCGCCTCTCTGCGCCTGCGCCGGCGCSSCGC" << std::endl;
		ss << "GCCTCTCTGCGCCTGCGCCGGCGCSSCGCGCCTCTCTGCGCCTGCGCCGGCGCSSC" << std::endl;
		ss << ">(CACCAT)n#Simple_repeat" << std::endl;
		ss << "CACCATCACCATCACCATCACCATCACCATCACCATCACCATCACCATCACCATCACCAT" << std::endl;
		ss << "CACCATCACCATCACCATCACCATCACCATCACCATCACCATCACCATCACCATCACCAT" << std::endl;
		ss << ">(CCCGAA)n#Simple_repeat" << std::endl;
		ss << "CCCGAACCCGAACCCGAACCCGAACCCGAACCCGAACCCGAACCCGAACCCGAACCCGAA" << std::endl;
		ss << "CCCGAACCCGAACCCGAACCCGAACCCGAACCCGAACCCGAACCCGAACCCGAACCCGAA" << std::endl;
		ss << ">(CCCCAA)n#Simple_repeat" << std::endl;
		ss << "CCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAA" << std::endl;
		ss << "CCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAA" << std::endl;
		ss << ">(CCCCAG)n#Simple_repeat" << std::endl;
		ss << "CCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAG" << std::endl;
		ss << "CCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAGCCCCAG" << std::endl;
		ss << ">(CCCTAA)n#Simple_repeat" << std::endl;
		ss << "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA" << std::endl;
		ss << "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA" << std::endl;
		ss << ">(GAAAAA)n#Simple_repeat" << std::endl;
		ss << "GAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAA" << std::endl;
		ss << "GAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAAGAAAAA" << std::endl;
		ss << ">(TAAAAA)n#Simple_repeat" << std::endl;
		ss << "TAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAA" << std::endl;
		ss << "TAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAATAAAAA" << std::endl;
		ss << ">(GGAGAA)n#Simple_repeat" << std::endl;
		ss << "GGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAA" << std::endl;
		ss << "GGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAAGGAGAA" << std::endl;
		ss << ">(GGGAGA)n#Simple_repeat" << std::endl;
		ss << "GGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGA" << std::endl;
		ss << "GGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGAGGGAGA" << std::endl;
		ss << ">(TGGGGG)n#Simple_repeat" << std::endl;
		ss << "TGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGG" << std::endl;
		ss << "TGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGGTGGGGG" << std::endl;
		ss << ">(CGGGGG)n#Simple_repeat" << std::endl;
		ss << "CGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGG" << std::endl;
		ss << "CGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGGCGGGGG" << std::endl;
		ss << ">(AGGGGG)n#Simple_repeat" << std::endl;
		ss << "AGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGG" << std::endl;
		ss << "AGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGGAGGGGG" << std::endl;
		ss << ">(CAGAGA)n#Simple_repeat" << std::endl;
		ss << "CAGAGACAGAGACAGAGACAGAGACAGAGACAGAGACAGAGACAGAGACAGAGACAGAGA" << std::endl;
		ss << "CAGAGACAGAGACAGAGACAGAGACAGAGACAGAGACAGAGACAGAGACAGAGACAGAGA" << std::endl;
		return ss.str();
	}

	static inline unsigned int &getPhiXReadIdx() {
		static unsigned int phiXReadIdx = -1;
		return phiXReadIdx;
	}
	static inline bool isPhiX(unsigned int readIdx) {
		return getPhiXReadIdx() == readIdx;
	}

	static std::string getPhiX() {
		std::stringstream ss;
		ss << ">gi|9626372|ref|NC_001422.1| Coliphage phiX174, complete genome" << std::endl;
		ss << "GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAAAAATTATCTT" << std::endl;
		ss << "GATAAAGCAGGAATTACTACTGCTTGTTTACGAATTAAATCGAAGTGGACTGCTGGCGGAAAATGAGAAA" << std::endl;
		ss << "ATTCGACCTATCCTTGCGCAGCTCGAGAAGCTCTTACTTTGCGACCTTTCGCCATCAACTAACGATTCTG" << std::endl;
		ss << "TCAAAAACTGACGCGTTGGATGAGGAGAAGTGGCTTAATATGCTTGGCACGTTCGTCAAGGACTGGTTTA" << std::endl;
		ss << "GATATGAGTCACATTTTGTTCATGGTAGAGATTCTCTTGTTGACATTTTAAAAGAGCGTGGATTACTATC" << std::endl;
		ss << "TGAGTCCGATGCTGTTCAACCACTAATAGGTAAGAAATCATGAGTCAAGTTACTGAACAATCCGTACGTT" << std::endl;
		ss << "TCCAGACCGCTTTGGCCTCTATTAAGCTCATTCAGGCTTCTGCCGTTTTGGATTTAACCGAAGATGATTT" << std::endl;
		ss << "CGATTTTCTGACGAGTAACAAAGTTTGGATTGCTACTGACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCT" << std::endl;
		ss << "TGCGTTTATGGTACGCTGGACTTTGTAGGATACCCTCGCTTTCCTGCTCCTGTTGAGTTTATTGCTGCCG" << std::endl;
		ss << "TCATTGCTTATTATGTTCATCCCGTCAACATTCAAACGGCCTGTCTCATCATGGAAGGCGCTGAATTTAC" << std::endl;
		ss << "GGAAAACATTATTAATGGCGTCGAGCGTCCGGTTAAAGCCGCTGAATTGTTCGCGTTTACCTTGCGTGTA" << std::endl;
		ss << "CGCGCAGGAAACACTGACGTTCTTACTGACGCAGAAGAAAACGTGCGTCAAAAATTACGTGCAGAAGGAG" << std::endl;
		ss << "TGATGTAATGTCTAAAGGTAAAAAACGTTCTGGCGCTCGCCCTGGTCGTCCGCAGCCGTTGCGAGGTACT" << std::endl;
		ss << "AAAGGCAAGCGTAAAGGCGCTCGTCTTTGGTATGTAGGTGGTCAACAATTTTAATTGCAGGGGCTTCGGC" << std::endl;
		ss << "CCCTTACTTGAGGATAAATTATGTCTAATATTCAAACTGGCGCCGAGCGTATGCCGCATGACCTTTCCCA" << std::endl;
		ss << "TCTTGGCTTCCTTGCTGGTCAGATTGGTCGTCTTATTACCATTTCAACTACTCCGGTTATCGCTGGCGAC" << std::endl;
		ss << "TCCTTCGAGATGGACGCCGTTGGCGCTCTCCGTCTTTCTCCATTGCGTCGTGGCCTTGCTATTGACTCTA" << std::endl;
		ss << "CTGTAGACATTTTTACTTTTTATGTCCCTCATCGTCACGTTTATGGTGAACAGTGGATTAAGTTCATGAA" << std::endl;
		ss << "GGATGGTGTTAATGCCACTCCTCTCCCGACTGTTAACACTACTGGTTATATTGACCATGCCGCTTTTCTT" << std::endl;
		ss << "GGCACGATTAACCCTGATACCAATAAAATCCCTAAGCATTTGTTTCAGGGTTATTTGAATATCTATAACA" << std::endl;
		ss << "ACTATTTTAAAGCGCCGTGGATGCCTGACCGTACCGAGGCTAACCCTAATGAGCTTAATCAAGATGATGC" << std::endl;
		ss << "TCGTTATGGTTTCCGTTGCTGCCATCTCAAAAACATTTGGACTGCTCCGCTTCCTCCTGAGACTGAGCTT" << std::endl;
		ss << "TCTCGCCAAATGACGACTTCTACCACATCTATTGACATTATGGGTCTGCAAGCTGCTTATGCTAATTTGC" << std::endl;
		ss << "ATACTGACCAAGAACGTGATTACTTCATGCAGCGTTACCATGATGTTATTTCTTCATTTGGAGGTAAAAC" << std::endl;
		ss << "CTCTTATGACGCTGACAACCGTCCTTTACTTGTCATGCGCTCTAATCTCTGGGCATCTGGCTATGATGTT" << std::endl;
		ss << "GATGGAACTGACCAAACGTCGTTAGGCCAGTTTTCTGGTCGTGTTCAACAGACCTATAAACATTCTGTGC" << std::endl;
		ss << "CGCGTTTCTTTGTTCCTGAGCATGGCACTATGTTTACTCTTGCGCTTGTTCGTTTTCCGCCTACTGCGAC" << std::endl;
		ss << "TAAAGAGATTCAGTACCTTAACGCTAAAGGTGCTTTGACTTATACCGATATTGCTGGCGACCCTGTTTTG" << std::endl;
		ss << "TATGGCAACTTGCCGCCGCGTGAAATTTCTATGAAGGATGTTTTCCGTTCTGGTGATTCGTCTAAGAAGT" << std::endl;
		ss << "TTAAGATTGCTGAGGGTCAGTGGTATCGTTATGCGCCTTCGTATGTTTCTCCTGCTTATCACCTTCTTGA" << std::endl;
		ss << "AGGCTTCCCATTCATTCAGGAACCGCCTTCTGGTGATTTGCAAGAACGCGTACTTATTCGCCACCATGAT" << std::endl;
		ss << "TATGACCAGTGTTTCCAGTCCGTTCAGTTGTTGCAGTGGAATAGTCAGGTTAAATTTAATGTGACCGTTT" << std::endl;
		ss << "ATCGCAATCTGCCGACCACTCGCGATTCAATCATGACTTCGTGATAAAAGATTGAGTGTGAGGTTATAAC" << std::endl;
		ss << "GCCGAAGCGGTAAAAATTTTAATTTTTGCCGCTGAGGGGTTGACCAAGCGAAGCGCGGTAGGTTTTCTGC" << std::endl;
		ss << "TTAGGAGTTTAATCATGTTTCAGACTTTTATTTCTCGCCATAATTCAAACTTTTTTTCTGATAAGCTGGT" << std::endl;
		ss << "TCTCACTTCTGTTACTCCAGCTTCTTCGGCACCTGTTTTACAGACACCTAAAGCTACATCGTCAACGTTA" << std::endl;
		ss << "TATTTTGATAGTTTGACGGTTAATGCTGGTAATGGTGGTTTTCTTCATTGCATTCAGATGGATACATCTG" << std::endl;
		ss << "TCAACGCCGCTAATCAGGTTGTTTCTGTTGGTGCTGATATTGCTTTTGATGCCGACCCTAAATTTTTTGC" << std::endl;
		ss << "CTGTTTGGTTCGCTTTGAGTCTTCTTCGGTTCCGACTACCCTCCCGACTGCCTATGATGTTTATCCTTTG" << std::endl;
		ss << "GATGGTCGCCATGATGGTGGTTATTATACCGTCAAGGACTGTGTGACTATTGACGTCCTTCCTCGTACGC" << std::endl;
		ss << "CGGGCAATAATGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCCGCGGATTGGT" << std::endl;
		ss << "TTCGCTGAATCAGGTTATTAAAGAGATTATTTGTCTCCAGCCACTTAAGTGAGGTGATTTATGTTTGGTG" << std::endl;
		ss << "CTATTGCTGGCGGTATTGCTTCTGCTCTTGCTGGTGGCGCCATGTCTAAATTGTTTGGAGGCGGTCAAAA" << std::endl;
		ss << "AGCCGCCTCCGGTGGCATTCAAGGTGATGTGCTTGCTACCGATAACAATACTGTAGGCATGGGTGATGCT" << std::endl;
		ss << "GGTATTAAATCTGCCATTCAAGGCTCTAATGTTCCTAACCCTGATGAGGCCGCCCCTAGTTTTGTTTCTG" << std::endl;
		ss << "GTGCTATGGCTAAAGCTGGTAAAGGACTTCTTGAAGGTACGTTGCAGGCTGGCACTTCTGCCGTTTCTGA" << std::endl;
		ss << "TAAGTTGCTTGATTTGGTTGGACTTGGTGGCAAGTCTGCCGCTGATAAAGGAAAGGATACTCGTGATTAT" << std::endl;
		ss << "CTTGCTGCTGCATTTCCTGAGCTTAATGCTTGGGAGCGTGCTGGTGCTGATGCTTCCTCTGCTGGTATGG" << std::endl;
		ss << "TTGACGCCGGATTTGAGAATCAAAAAGAGCTTACTAAAATGCAACTGGACAATCAGAAAGAGATTGCCGA" << std::endl;
		ss << "GATGCAAAATGAGACTCAAAAAGAGATTGCTGGCATTCAGTCGGCGACTTCACGCCAGAATACGAAAGAC" << std::endl;
		ss << "CAGGTATATGCACAAAATGAGATGCTTGCTTATCAACAGAAGGAGTCTACTGCTCGCGTTGCGTCTATTA" << std::endl;
		ss << "TGGAAAACACCAATCTTTCCAAGCAACAGCAGGTTTCCGAGATTATGCGCCAAATGCTTACTCAAGCTCA" << std::endl;
		ss << "AACGGCTGGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAGGTTAGTGCTGAGGTTGAC" << std::endl;
		ss << "TTAGTTCATCAGCAAACGCAGAATCAGCGGTATGGCTCTTCTCATATTGGCGCTACTGCAAAGGATATTT" << std::endl;
		ss << "CTAATGTCGTCACTGATGCTGCTTCTGGTGTGGTTGATATTTTTCATGGTATTGATAAAGCTGTTGCCGA" << std::endl;
		ss << "TACTTGGAACAATTTCTGGAAAGACGGTAAAGCTGATGGTATTGGCTCTAATTTGTCTAGGAAATAACCG" << std::endl;
		ss << "TCAGGATTGACACCCTCCCAATTGTATGTTTTCATGCCTCCAAATCTTGGAGGCTTTTTTATGGTTCGTT" << std::endl;
		ss << "CTTATTACCCTTCTGAATGTCACGCTGATTATTTTGACTTTGAGCGTATCGAGGCTCTTAAACCTGCTAT" << std::endl;
		ss << "TGAGGCTTGTGGCATTTCTACTCTTTCTCAATCCCCAATGCTTGGCTTCCATAAGCAGATGGATAACCGC" << std::endl;
		ss << "ATCAAGCTCTTGGAAGAGATTCTGTCTTTTCGTATGCAGGGCGTTGAGTTCGATAATGGTGATATGTATG" << std::endl;
		ss << "TTGACGGCCATAAGGCTGCTTCTGACGTTCGTGATGAGTTTGTATCTGTTACTGAGAAGTTAATGGATGA" << std::endl;
		ss << "ATTGGCACAATGCTACAATGTGCTCCCCCAACTTGATATTAATAACACTATAGACCACCGCCCCGAAGGG" << std::endl;
		ss << "GACGAAAAATGGTTTTTAGAGAACGAGAAGACGGTTACGCAGTTTTGCCGCAAGCTGGCTGCTGAACGCC" << std::endl;
		ss << "CTCTTAAGGATATTCGCGATGAGTATAATTACCCCAAAAAGAAAGGTATTAAGGATGAGTGTTCAAGATT" << std::endl;
		ss << "GCTGGAGGCCTCCACTATGAAATCGCGTAGAGGCTTTGCTATTCAGCGTTTGATGAATGCAATGCGACAG" << std::endl;
		ss << "GCTCATGCTGATGGTTGGTTTATCGTTTTTGACACTCTCACGTTGGCTGACGACCGATTAGAGGCGTTTT" << std::endl;
		ss << "ATGATAATCCCAATGCTTTGCGTGACTATTTTCGTGATATTGGTCGTATGGTTCTTGCTGCCGAGGGTCG" << std::endl;
		ss << "CAAGGCTAATGATTCACACGCCGACTGCTATCAGTATTTTTGTGTGCCTGAGTATGGTACAGCTAATGGC" << std::endl;
		ss << "CGTCTTCATTTCCATGCGGTGCACTTTATGCGGACACTTCCTACAGGTAGCGTTGACCCTAATTTTGGTC" << std::endl;
		ss << "GTCGGGTACGCAATCGCCGCCAGTTAAATAGCTTGCAAAATACGTGGCCTTATGGTTACAGTATGCCCAT" << std::endl;
		ss << "CGCAGTTCGCTACACGCAGGACGCTTTTTCACGTTCTGGTTGGTTGTGGCCTGTTGATGCTAAAGGTGAG" << std::endl;
		ss << "CCGCTTAAAGCTACCAGTTATATGGCTGTTGGTTTCTATGTGGCTAAATACGTTAACAAAAAGTCAGATA" << std::endl;
		ss << "TGGACCTTGCTGCTAAAGGTCTAGGAGCTAAAGAATGGAACAACTCACTAAAAACCAAGCTGTCGCTACT" << std::endl;
		ss << "TCCCAAGAAGCTGTTCAGAATCAGAATGAGCCGCAACTTCGGGATGAAAATGCTCACAATGACAAATCTG" << std::endl;
		ss << "TCCACGGAGTGCTTAATCCAACTTACCAAGCTGGGTTACGACGCGACGCCGTTCAACCAGATATTGAAGC" << std::endl;
		ss << "AGAACGCAAAAAGAGAGATGAGATTGAGGCTGGGAAAAGTTACTGTAGCCGACGTTTTGGCGGCGCAACC" << std::endl;
		ss << "TGTGACGACAAATCTGCTCAAATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA"     << std::endl;
		return ss.str();
	}

	static inline unsigned int &getReferenceReadIdx() {
		static unsigned int referenceReadId = -1;
		return referenceReadId;
	}
	static inline bool isReference(unsigned int readIdx) {
		return getReferenceReadIdx() <= readIdx;
	}
	static inline bool hasReference() {
		return getReferenceReadIdx() != (unsigned int) -1;
	}

};

#endif

