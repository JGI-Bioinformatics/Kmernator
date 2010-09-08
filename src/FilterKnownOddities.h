//
// Kmernator/src/FilterKnownOddities.h
//
// Author: Rob Egan, Craig Furman
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
#include "Log.h"

class FilterKnownOddities {
public:
	typedef Kmer::NumberType NumberType;
	typedef KmerMap<unsigned short> KM;
	typedef std::vector< ReadSet::ReadSetSizeType > SequenceCounts;

private:
	ReadSet sequences;
	unsigned short length;
	unsigned short twoBitLength;
	KM filter;
	SequenceCounts counts;

public:
	FilterKnownOddities(int _length = 24, int numErrors = 2) :
		length(_length), filter(32*1024*1024) {
		if (length > 28) {
			throw std::invalid_argument("FilterKnownOddities must use 7 bytes or less (<= 28 bases)");
		}
		twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(length);
		// T is 11, A is 00, so mask is all T's surrounded by A's

		// readIdx 0 is signal for no match!
		Read empty;
		sequences.append(empty);

		std::string fasta = getArtifactFasta();
		sequences.appendFastaFile(fasta);
		if (Options::getMaskSimpleRepeats()) {
			fasta = getSimpleRepeatFasta();
			getSimpleRepeatBegin() = sequences.getSize();
			sequences.appendFastaFile(fasta);
			getSimpleRepeatEnd() = sequences.getSize();
		}
		if (Options::getPhiXOutput()) {
			fasta = getPhiX();
			getPhiXReadIdx() = sequences.getSize();
			sequences.appendFastaFile(fasta);
		}
		counts.resize( sequences.getSize() );

		prepareMaps(numErrors);
	}
	~FilterKnownOddities() {
		clear();
	}
	void clear() {
		filter.clear();
		counts.clear();
	}

	void prepareMaps(int numErrors) {
		unsigned long oldKmerLength = KmerSizer::getSequenceLength();
		assert(length % 4 == 0);
		KmerSizer::set(length);

		sequences.circularize(length);

		for (unsigned short i = 0; i < sequences.getSize(); i++) {
			const Read read = sequences.getRead(i);
			KmerWeights kmers = KmerReadUtils::buildWeightedKmers(read, true);
			for (Kmer::IndexType j = 0; j < kmers.size(); j++) {
				filter.getOrSetElement( kmers[j] , i );
			}
		}
		LOG_VERBOSE(1,  "Prepared exact match: " << filter.size() << " " << MemoryUtils::getMemoryUsage() );

		for (int error = 0; error < numErrors; error++) {
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
			LOG_VERBOSE(1, "Prepared order " << (error+1) << ": " << filter.size() << " " << MemoryUtils::getMemoryUsage() );
		}

		KmerSizer::set(oldKmerLength);
	}

	class FilterResults {
	public:
		KM::ValueType value;
		SequenceLengthType minAffected;
		SequenceLengthType maxAffected;
		FilterResults() : value(0), minAffected(0), maxAffected(0) {}
		FilterResults(KM::ValueType &_v, SequenceLengthType &_m, SequenceLengthType &_x) : value(_v), minAffected(_m), maxAffected(_x){}
	};
	class Recorder {
	public:
		OfstreamMap *omPhiX;
		OfstreamMap *omArtifact;
		SequenceCounts *threadCounts;
		SequenceLengthType minReadPos;
		Recorder() : omPhiX(NULL), omArtifact(NULL), threadCounts(NULL), minReadPos(0) {}
	};

	FilterResults applyFilterToPair(ReadSet &reads, long pairIdx, SequenceCounts *threadCounts, Recorder &recorder) {
		ReadSet::Pair &pair = reads.getPair(pairIdx);
		FilterResults results1, results2;

		bool isRead1 = reads.isValidRead(pair.read1);
		bool isRead2 = reads.isValidRead(pair.read2);
		if (isRead1)
			results1 = applyFilterToRead(reads, pair.read1, threadCounts, recorder);
		if (isRead2)
			results2 = applyFilterToRead(reads, pair.read2, threadCounts, recorder);

		FilterResults results;
		if (isRead1 && isRead2) {
			recordAffectedRead(reads, recorder, results1, pair.read1, results2, pair.read2);
			results.minAffected = std::min(results1.minAffected, results2.minAffected);
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

	FilterResults applyFilterToRead(ReadSet &reads, long readIdx, SequenceCounts *threadCounts, Recorder &recorder, bool recordEffects = false) {

		Read &read = reads.getRead(readIdx);
		FilterResults results;
		KM::ValueType &value = results.value;
		SequenceLengthType &minAffected = results.minAffected;
		minAffected = MAX_SEQUENCE_LENGTH;
		SequenceLengthType &maxAffected = results.maxAffected;
		maxAffected = 0;

		LOG_DEBUG_MT(2, "Checking " << read.getName() << "\t" << read.getFasta() );

		SequenceLengthType seqLen = read.getLength();
		TwoBitEncoding *ptr = read.getTwoBitSequence();
		long bytes = read.getTwoBitEncodingSequenceLength();
		TwoBitEncoding revcomp[bytes+1];
		TwoBitEncoding *revPtr = revcomp;
		SequenceLengthType seqLenByteBoundary = seqLen & ~((SequenceLengthType) 0x03);
		TwoBitSequence::reverseComplement(ptr, revPtr, seqLenByteBoundary);
		long byteHops = bytes - twoBitLength - ((seqLen & 0x03) == 0 ? 0 : 1);
		if (byteHops < 0)
			return results;

		KM::ElementType elem;
		bool wasPhiX = false;


		for(long byteHop = 0; byteHop <= byteHops; byteHop++) {

			// need to test both forward and reverse paths since only one is stored in filter

			const Kmer &fwd = (const Kmer&) *ptr;

			elem = filter.getElementIfExists( fwd );
			if (elem.isValid()) {
				SequenceLengthType pos = byteHop*4;
				value = elem.value();
				wasPhiX |= isPhiX(value);
				if (minAffected > pos)
					minAffected = pos;
				if (maxAffected < pos+length)
					maxAffected = pos+length;
			}

			const Kmer &rev = (const Kmer&) *revPtr;

			elem = filter.getElementIfExists( rev );
			if (elem.isValid()) {
			    SequenceLengthType pos = 0;
			    SequenceLengthType posMinus = length + byteHop*4;
			    if (seqLen > posMinus) {
			        pos = seqLen - posMinus;
			    }
				value = elem.value();
				wasPhiX |= isPhiX(value);
				if (minAffected > pos)
					minAffected = pos;
				if (maxAffected < pos+length)
					maxAffected = pos+length;
			}

			ptr++;
			revPtr++;
		}

		if (wasPhiX) {
			value = getPhiXReadIdx();
		} else if (isSimpleRepeat(value)) {
			// allow simple repeats in the middle of a read with good edges
			if (minAffected > recorder.minReadPos && seqLen - maxAffected >= length) {
				value = 0;
				minAffected = 0;
				maxAffected = 0;
			}
		}
		if (value > 0) {
			read.markupBases(minAffected , maxAffected - minAffected, 'X');
			if (recordEffects)
				recordAffectedRead(reads, recorder, results, readIdx);
		}
		return results;
	}

	void recordAffectedRead(ReadSet &reads, Recorder &recorder, FilterResults results1, ReadSet::ReadSetSizeType readIdx1, FilterResults results2 = FilterResults(), ReadSet::ReadSetSizeType readIdx2 = ReadSet::MAX_READ_IDX) {
		bool isRead1 = reads.isValidRead(readIdx1);
		bool isRead2 = reads.isValidRead(readIdx2);

		bool wasAffected = (results1.value != 0) | (results2.value != 0);
		bool wasPhiX = isPhiX(results1.value) | isPhiX(results2.value);

		if (wasAffected) {

			if (wasPhiX) {
				results1.value = results2.value = getPhiXReadIdx();
			}

			bool isRead1Affected = isRead1 && results1.value != 0;
			bool isRead2Affected = isRead2 && results2.value != 0;

			int threadNum = omp_get_thread_num();

			if (isRead1Affected)
				recorder.threadCounts[threadNum][ results1.value ]++;
			if (isRead2Affected)
				recorder.threadCounts[threadNum][ results2.value ]++;

			if (wasPhiX && recorder.omPhiX != NULL) {

				std::string fileSuffix = std::string("-") + reads.getReadFileNamePrefix(readIdx1);
				#pragma omp critical (writePhix)
				{
					ostream &os = recorder.omPhiX->getOfstream( fileSuffix );
		            // always discard the read, as it contains some PhiX and was sorted
					if (isRead1) {
						Read &read = reads.getRead(readIdx1);
						_writeFilterRead(os, read, read.getLength());
						read.discard();
					}
					if (isRead2) {
						Read &read = reads.getRead(readIdx2);
						_writeFilterRead(os, read, read.getLength());
						read.discard();
					}
				}
			} else if ( (!wasPhiX) && recorder.omArtifact != NULL) {

				std::string fileSuffix = std::string("-") + reads.getReadFileNamePrefix( isRead1 ? readIdx1 : readIdx2);
				std::string label1, label2;
				if (isRead1Affected)
					label1 = sequences.getRead(results1.value).getName();
				if (isRead2Affected)
					label2 = sequences.getRead(results2.value).getName();

				#pragma omp critical (writeFilter)
				{
					ostream &os = recorder.omArtifact->getOfstream( fileSuffix );
					// if the read length is less than the desired minimum, discard it completely regardless if it was output or not
					if (isRead1 && results1.value != 0) {
						Read &read = reads.getRead(readIdx1);
						_writeFilterRead(os, read, read.getLength(), label1);
						if (results1.minAffected == 0 || results1.minAffected < recorder.minReadPos) {
							read.discard();
						}
					}
					if (isRead2 && results2.value != 0) {
						Read &read = reads.getRead(readIdx2);
						_writeFilterRead(os, read, read.getLength(), label2);
						if (results2.minAffected == 0 || results2.minAffected < recorder.minReadPos) {
							read.discard();
						}
					}
				}
			}

			if (Log::isDebug(1)) {

				#pragma omp critical (stderr)
				{
				  ostream &debug = Log::Debug();
				  Read read;
				  if (isRead1Affected) {
					  read = reads.getRead(readIdx1);
					  debug << "FilterMatch1 to " << read.getName() << " "
						  << read.getFastaNoMarkup() << " " << read.getFasta() << " "
							  << wasPhiX << " " << sequences.getRead(results1.value).getName() << std::endl;
				  }
				  if (isRead2Affected) {
					  read = reads.getRead(readIdx2);
					  debug << "FilterMatch2 to " << read.getName() << " "
						  << read.getFastaNoMarkup() << " " << read.getFasta() << " "
							  << wasPhiX << " " << sequences.getRead(results2.value).getName() << std::endl;

				  }
			    }
			}
		}
	}

	unsigned long applyFilter(ReadSet &reads) {
		unsigned long oldKmerLength = KmerSizer::getSequenceLength();
		KmerSizer::set(length);

		Recorder recorder;
		unsigned long affectedCount = 0;

		OfstreamMap _omPhiX(Options::getOutputFile(), "-PhiX.fastq");
		if (Options::getPhiXOutput()) {
			recorder.omPhiX = &_omPhiX;
		}
		OfstreamMap _omArtifact(Options::getOutputFile(), "-Artifact.fastq");
		if (Options::getFilterOutput()) {
			recorder.omArtifact = &_omArtifact;
		}

		SequenceLengthType &minReadPos = recorder.minReadPos;
		minReadPos = Options::getMinReadLength();
		if (minReadPos != 0 && minReadPos != MAX_SEQUENCE_LENGTH) {
			minReadPos -= 1;
		}

		int numThreads = omp_get_max_threads();
		SequenceCounts     threadCounts[numThreads];
		recorder.threadCounts = threadCounts;

        for (int i = 0; i < numThreads; i++) {
            threadCounts[i].resize( sequences.getSize() );
        }
		// start with any existing state
		counts.swap(threadCounts[0]);

		bool byPair = reads.hasPairs();
		long size = reads.getSize();
		if (byPair)
			size = reads.getPairSize();

		#pragma omp parallel for schedule(dynamic)
		for (long idx = 0; idx < size; idx++) {

			if (byPair)
				applyFilterToPair(reads, idx, threadCounts, recorder);
			else {
				applyFilterToRead(reads, idx, threadCounts, recorder, true);
			}
		}

		// consolidate threaded instances of global variables
		for(int i = 1 ; i < numThreads; i++) {
			for(ReadSetSizeType j = 0 ; j < threadCounts[i].size(); j++) {
				if (threadCounts[0].size() <= j)
					threadCounts[0].resize(j+1);
				threadCounts[0][j] += threadCounts[i][j];
			}
		}

		// restore state
		counts.swap(threadCounts[0]);
		for(ReadSetSizeType j = 0 ; j < counts.size() ; j++)
			affectedCount += counts[j];

		if (Log::isVerbose(1)) {
			ostream &log = Log::Verbose();
			log << "Final Filter Matches to reads: " << affectedCount << std::endl;
			// TODO sort
			for(unsigned long idx = 0 ; idx < counts.size(); idx++) {
				if (counts[idx] > 0) {
				    log << "\t" << counts[idx] << "\t" << sequences.getRead(idx).getName() << std::endl;
				}
			}
		}

		KmerSizer::set(oldKmerLength);
		return affectedCount;
	}

	static void _writeFilterRead(ostream &os, Read  &read, SequenceLengthType readLength, std::string readLabel = "") {
	    read.write(os, readLength, readLabel, FormatOutput::FASTQ_UNMASKED);
	}

	const ReadSet &getSequences() const {
		return sequences;
	}

	typedef TrackingDataWithAllReads TD;
	typedef KmerSpectrum<TD, TD, TD> KS;
	typedef TrackingData::ReadPositionWeightVector RPW;
	typedef KS::WeakMapType::ElementType KSElementType;
	typedef std::vector< KSElementType > KSElementVector;
	typedef RPW::iterator RPWIterator;
	typedef ReadSet::Pair Pair;
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;

	static void _buildDuplicateFragmentMap(KS::Vector &ksv, ReadSet &reads, unsigned char bytes, bool useReverseComplement, bool paired, unsigned int startOffset = Options::getDeDupStartOffset()) {
		// build one KS per thread, then merge, skipping singletons
        // no need to include quality scores

		// build the kmer spectrum with the concatenated prefixes
		// from each read in the pair
		int numThreads = ksv.size();
		KmerWeights::Vector tmpKmerv(numThreads);
		for(int i = 0; i < numThreads; i++) {
		  ksv[i] = KS(reads.getPairSize() / 64 / numThreads, false);
		  tmpKmerv[i].resize(1);
		  tmpKmerv[i].valueAt(0) = 1.0;
		}
		long pairSize =  reads.getPairSize();

		ReadSet::madviseMmapsSequential();
		if (!paired)
			bytes *= 2;
		SequenceLengthType sequenceLength = bytes * 4;

	    #pragma omp parallel for
		for(long pairIdx = 0; pairIdx < pairSize; pairIdx++) {
			Pair &pair = reads.getPair(pairIdx);
			int threadNum = omp_get_thread_num();
			KmerWeights &kmerWeights = tmpKmerv[threadNum];
			Kmer &kmer = kmerWeights[0];

			if (paired && pair.isPaired() ) {
				if(reads.isValidRead(pair.read1) && reads.isValidRead(pair.read2)) {
					const Read &read1 = reads.getRead(pair.read1);
					const Read &read2 = reads.getRead(pair.read2);
					if (read1.isDiscarded() || read2.isDiscarded())
						continue;

					// create read1 + the reverse complement of read2 (1:rev2)
					// when useReverseComplement, it is represented as a kmer, and the leastcomplement of 1:rev2 and 2:rev1 will be stored
					// and properly account for duplicate fragment pairs

					Sequence::BaseLocationVectorType markups = read1.getMarkups();
					if (TwoBitSequence::firstMarkupX(markups) + startOffset < sequenceLength) {
						memcpy(kmer.getTwoBitSequence()       , read1.getTwoBitSequence() + (startOffset/4), bytes);
					} else {
						continue;
					}
					markups = read2.getMarkups();
					if (TwoBitSequence::firstMarkupX(markups) + startOffset < sequenceLength) {
						TwoBitSequence::reverseComplement( read2.getTwoBitSequence() + (startOffset/4), kmer.getTwoBitSequence() + bytes, sequenceLength);
					} else {
						continue;
					}

					long myPairIdx = pairIdx;
					if (useReverseComplement) {
						// choose orientation and flag in pairIdx
						TEMP_KMER(tmpRevComp);
						if (! kmer.buildLeastComplement(tmpRevComp) ) {
							kmer = tmpRevComp;
							myPairIdx = pairIdx + pairSize;
						}
					}
					// store the pairIdx (not readIdx)
					ksv[threadNum].append(kmerWeights, myPairIdx);
				}
			} else if ( pair.isSingle() && (!paired) ) {
				ReadSetSizeType readIdx = pair.lesser();
				if (reads.isValidRead(readIdx)) {
					const Read &read1 = reads.getRead(readIdx);
					if (read1.isDiscarded())
						continue;
					Sequence::BaseLocationVectorType markups = read1.getMarkups();
					if (TwoBitSequence::firstMarkupX(markups) + startOffset < sequenceLength) {
						memcpy(kmer.getTwoBitSequence()        , read1.getTwoBitSequence() + (startOffset/4), bytes);
					} else {
						continue;
					}
					// store the readIdx (not the pairIdx)
					ksv[threadNum].append(kmerWeights, readIdx);
				}
			}
		}
		if (Log::isDebug(3)) {
		  for (int i = 0; i < numThreads; i++) {
			  Log::Debug() << "spectrum " << i << std::endl;
			  ksv[i].printHistograms();
		  }
		} else {
			LOG_VERBOSE(1, "merging duplicate fragment spectrums" );
		}

		KS::mergeVector(ksv, 1);
	}

	// TODO make useWeights an Option::
	static void _mergeNodesWithinEditDistance(KS &ks, unsigned int cutoffThreshold, unsigned int editDistance, bool useWeights = true) {
		// TODO honor edit distance > 1
		LOG_VERBOSE(1, "Merging kmers within edit-distance of " << editDistance << " " << MemoryUtils::getMemoryUsage());

		// create a sorted set of elements with > cutoffThreshold count
		// (singletons will not be included in this round)
		KSElementVector elems;

		#pragma omp parallel private(elems)
		{
			for(KS::WeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
				KSElementType &elem = *it;
				if (elem.isValid() && elem.value().getCount() >= cutoffThreshold)
					elems.push_back(elem);
			}

			#pragma omp single
			{
				LOG_VERBOSE(1, "Sorting all nodes >= " << cutoffThreshold << " count: " << elems.size() << " " << MemoryUtils::getMemoryUsage());
			}
			std::sort(elems.begin(), elems.end());


			#pragma omp single
			{
				LOG_VERBOSE(1, "Merging elements. " << MemoryUtils::getMemoryUsage());
			}
			for(KSElementVector::reverse_iterator it = elems.rbegin(); it != elems.rend(); it++) {
				KSElementType &elem = *it;
				if (elem.isValid() && elem.value().getCount() >= cutoffThreshold) {
					ks.consolidate(elem.key(), useWeights);
				}
			}

			#pragma omp single
			{
				LOG_VERBOSE(1, "Merging elements below cutoffThreshold. " << MemoryUtils::getMemoryUsage() );
			}

			// do not clear elements until all merging has completed
			elems.clear();

			// now merge those less than cutoff (singletons by default)
			for(KS::WeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
				KSElementType &elem = *it;
				if (elem.isValid()) {
					TD &value = elem.value();
					if (value > 0 && value.getCount() < cutoffThreshold) {
						ks.consolidate(elem.key(), useWeights);
					}
				}
			}
		} // omp parallel

		LOG_VERBOSE(1, "Merged histogram. " << MemoryUtils::getMemoryUsage() );
		ks.printHistograms();

	}
	static ReadSetSizeType _buildConsensusUnPairedReads(KS &ks, ReadSet &reads, ReadSet &newReads, unsigned int cutoffThreshold) {
		ReadSet::madviseMmapsRandom();
		LOG_VERBOSE(1, "Building consensus reads. " << MemoryUtils::getMemoryUsage());

		ReadSetSizeType affectedCount = 0;

		#pragma omp parallel reduction(+:affectedCount)
		for(KS::WeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
		    if (it->value().getCount() >= cutoffThreshold) {
		    	RPW rpw = it->value().getEachInstance();

		    	ReadSet tmpReadSet1;

		    	for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {

		    		ReadSetSizeType readIdx = rpwit->readId;

		    		const Read &read1 = reads.getRead(readIdx);
		    		tmpReadSet1.append( read1 );

		    	}

		    	Read consensus1 = tmpReadSet1.getConsensusRead();

		    	#pragma omp critical (BCUR_newReads)
		    	{
		    	   newReads.append(consensus1);
		    	}
		    	affectedCount += rpw.size();

		    	for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {

		    		ReadSetSizeType readIdx = rpwit->readId;

		    	    Read &read1 = reads.getRead(readIdx);
		    		read1.discard();
		    	}

		    }
		}
		LOG_VERBOSE(1, "Clearing duplicate pair map: " << MemoryUtils::getMemoryUsage() );
		ks.reset();
		LOG_VERBOSE(1, "Built " << newReads.getSize() << " new consensus reads: " <<  MemoryUtils::getMemoryUsage() );
		newReads.identifyPairs();
		return affectedCount;
	}

	static ReadSetSizeType _buildConsensusPairedReads(KS &ks, ReadSet &reads, ReadSet &newReads, unsigned int cutoffThreshold) {
		ReadSet::madviseMmapsRandom();
		LOG_VERBOSE(1, "Building consensus reads. " << MemoryUtils::getMemoryUsage());

		ReadSetSizeType affectedCount = 0;
		ReadSetSizeType pairSize = reads.getPairSize();
		int numThreads = omp_get_max_threads();
		ReadSet _threadNewReads[numThreads];

		#pragma omp parallel reduction(+:affectedCount)
		for(KS::WeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
		    if (it->value().getCount() >= cutoffThreshold) {
		    	RPW rpw = it->value().getEachInstance();

		    	ReadSet tmpReadSet1;
		    	ReadSet tmpReadSet2;
		    	ReadSet &threadNewReads = _threadNewReads[omp_get_thread_num()];

		    	for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {

		    		// iterator readId is actually the pairIdx built above
		    		ReadSetSizeType pairIdx = rpwit->readId;

		    		// correct orientation
		    		bool isCorrectOrientation = true;
		    		if (pairIdx >= pairSize) {
		    			isCorrectOrientation = false;
		    			pairIdx = pairIdx - pairSize;
		    		}
		    		Pair &pair = reads.getPair(pairIdx);
		    		const Read &read1 = reads.getRead(isCorrectOrientation ? pair.read1 : pair.read2);
		    		const Read &read2 = reads.getRead(isCorrectOrientation ? pair.read2 : pair.read1);

		    		tmpReadSet1.append( read1 );
		    		tmpReadSet2.append( read2 );

		    	}

		    	Read consensus1 = tmpReadSet1.getConsensusRead();
		    	Read consensus2 = tmpReadSet2.getConsensusRead();

		    	threadNewReads.append(consensus1);
		    	threadNewReads.append(consensus2);

		    	affectedCount += 2 * rpw.size();

		    	for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {

		    	    ReadSetSizeType pairIdx = rpwit->readId;

		    	    // orientation does not matter here, but correcting the index is important!
		    	    if (pairIdx >= pairSize) {
		    	    	pairIdx = pairIdx - pairSize;
		    	    }
		    	    Pair &pair = reads.getPair(pairIdx);
		    	    Read &read1 = reads.getRead(pair.read1);
		    	    Read &read2 = reads.getRead(pair.read2);
		    		read1.discard();
		    		read2.discard();
		    	}

		    }
		}
		for(int i = 0 ; i < numThreads; i++)
			newReads.append(_threadNewReads[i]);

		LOG_VERBOSE(1, "Clearing duplicate pair map: " << MemoryUtils::getMemoryUsage() );
		ks.reset();
		LOG_VERBOSE(1, "Built " << newReads.getSize() << " new consensus reads: " <<  MemoryUtils::getMemoryUsage() );
		newReads.identifyPairs();
		return affectedCount;
	}

	static ReadSetSizeType _filterDuplicateFragments(ReadSet &reads, unsigned char bytes, unsigned int cutoffThreshold, unsigned int editDistance, bool paired) {

		int numThreads = omp_get_max_threads();
		KS::Vector ksv(numThreads);

		LOG_VERBOSE(1, "Building " << (paired?"Paired":"Un-Paired") << " Duplicate Fragment Spectrum" );

		SequenceLengthType affectedCount = 0;

        bool useReverseComplement = (Options::getDeDupMode() == 2);

		// build the paired duplicate fragment map
		_buildDuplicateFragmentMap(ksv, reads, bytes, useReverseComplement, paired);

		KS &ks = ksv[0];
		// analyze the spectrum
		ks.printHistograms();

		if (editDistance > 0) {
			_mergeNodesWithinEditDistance(ks, cutoffThreshold, editDistance);
		}

		ReadSet newReads;
		if (paired) {
			affectedCount += _buildConsensusPairedReads(ks, reads, newReads, cutoffThreshold);
		} else {
			affectedCount += _buildConsensusUnPairedReads(ks, reads, newReads, cutoffThreshold);
		}
		reads.append(newReads);

		return affectedCount;
	}

	static ReadSetSizeType filterDuplicateFragments(ReadSet &reads, unsigned char sequenceLength = Options::getDeDupLength(), unsigned int cutoffThreshold = 2, unsigned int editDistance = Options::getDeDupEditDistance()) {

	  if ( Options::getDeDupMode() == 0 || editDistance == (unsigned int) -1) {
		  LOG_VERBOSE(1, "Skipping filter and merge of duplicate fragments");
		  return 0;
	  }
	  ReadSetSizeType affectedCount = 0;

	  // select the number of bytes from each pair to scan
	  unsigned char bytes = sequenceLength / 4;
	  if (bytes == 0) {
			bytes = 1;
	  }
	  SequenceLengthType oldKmerSize = KmerSizer::getSequenceLength();
	  KmerSizer::set(bytes * 4 * 2);

	  affectedCount += _filterDuplicateFragments(reads, bytes, cutoffThreshold, editDistance, true);

	  if (Options::getDeDupSingle() == 1)
		  affectedCount += _filterDuplicateFragments(reads, bytes, cutoffThreshold, editDistance, false);

      KmerSizer::set(oldKmerSize);
      ReadSet::madviseMmapsSequential();

	  return affectedCount;
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
	    ss << ">Gex_Adapter_1             " << std::endl;
	    ss << "GATCGTCGGACTGTAGAACTCTGAAC " << std::endl;
	    ss << ">Gex_Adapter_1_2" << std::endl;
	    ss << "ACAGGTTCAGAGTTCTACAGTCCGAC " << std::endl;
	    ss << ">Gex_Adapter_2" << std::endl;
	    ss << "CAAGCAGAAGACGGCATACGANN " << std::endl;
	    ss << ">Gex_Adapter_2_2" << std::endl;
	    ss << "TCGTATGCCGTCTTCTGCTTG " << std::endl;
	    ss << ">Gex_PCR_Primer_1" << std::endl;
	    ss << "CAAGCAGAAGACGGCATACGA " << std::endl;
	    ss << ">Gex_PCR_Primer_2" << std::endl;
	    ss << "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA " << std::endl;
	    ss << ">Gex_Sequencing_Primer" << std::endl;
	    ss << "CGACAGGTTCAGAGTTCTACAGTCCGACGATC " << std::endl;
	    ss << ">Adapters1" << std::endl;
	    ss << "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG  " << std::endl;
	    ss << ">Adapters1_1" << std::endl;
	    ss << "ACACTCTTTCCCTACACGACGCTCTTCCGATCT   " << std::endl;
	    ss << ">PCR_Primers1" << std::endl;
	    ss << "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT " << std::endl;
	    ss << ">PCR Primers1_1" << std::endl;
	    ss << "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT   " << std::endl;
	    ss << ">Genomic_DNA_Sequencing_Primer" << std::endl;
	    ss << "ACACTCTTTCCCTACACGACGCTCTTCCGATCT " << std::endl;
	    ss << ">PE_Adapters1" << std::endl;
	    ss << "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG " << std::endl;
	    ss << ">PE_Adapters1_" << std::endl;
	    ss << "ACACTCTTTCCCTACACGACGCTCTTCCGATCT " << std::endl;
	    ss << ">PE_PCR_Primers1" << std::endl;
	    ss << "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT " << std::endl;
	    ss << ">PE_PCR_Primers1_1" << std::endl;
	    ss << "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT " << std::endl;
	    ss << ">PE_Sequencing_Primer" << std::endl;
	    ss << "ACACTCTTTCCCTACACGACGCTCTTCCGATCT " << std::endl;
	    ss << ">PE_Sequencing_Primer_1" << std::endl;
	    ss << "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT" << std::endl;
		return ss.str();
	}


	static inline unsigned short &getSimpleRepeatBegin() {
		static unsigned short simpleRepeatBegin = -1;
		return simpleRepeatBegin;
	}
	static inline unsigned short &getSimpleRepeatEnd() {
		static unsigned short getSimpleRepeatEnd = -1;
		return getSimpleRepeatEnd;
	}

	static inline bool isSimpleRepeat(unsigned short readIdx) {
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

	static inline unsigned short &getPhiXReadIdx() {
		static unsigned short phiXReadIdx = -1;
		return phiXReadIdx;
	}
	static inline bool isPhiX(unsigned short readIdx) {
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

};

#endif

// $Log: FilterKnownOddities.h,v $
// Revision 1.28  2010-08-18 17:50:40  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.27.4.2  2010-07-21 18:06:11  regan
// added options to change unique mask offset and length for de-duplication
//
// Revision 1.27.4.1  2010-07-21 17:27:48  regan
// refactored
// filter by pairs, if available
// output PhiX by pairs
// allow simple repeats if edge sequences are unmasked
//
// Revision 1.27  2010-06-23 20:58:02  regan
// fixed minimum read length logic
//
// Revision 1.26  2010-06-22 23:06:31  regan
// merged changes in CorruptionBugfix-20100622 branch
//
// Revision 1.25.4.1  2010-06-22 23:02:51  regan
// named all critical sections
// made the code a bit clearer to follow
//
// Revision 1.25  2010-05-24 21:48:46  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.24.2.6  2010-05-24 21:44:44  regan
// save memory and scanning time by outputting reads as they are filtered (if output is requested).
//
// Revision 1.24.2.5  2010-05-20 18:26:43  regan
// attempt to fix a race condition when consolidating/merging edit-distance spectrums
//
// Revision 1.24.2.4  2010-05-20 03:42:24  regan
// added RNA_Linker to filter sequences
// optimized parallel performance in consensus generation
//
// Revision 1.24.2.3  2010-05-19 22:43:49  regan
// bugfixes
//
// Revision 1.24.2.2  2010-05-19 21:53:20  regan
// bugfixes
//
// Revision 1.24.2.1  2010-05-19 21:36:54  regan
// refactored duplicate fragment filter code
// added duplicate fragment on single ended reads
//
// Revision 1.24  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.23.2.4  2010-05-12 22:45:00  regan
// added readset circularize method
//
// Revision 1.23.2.3  2010-05-12 20:47:12  regan
// bugfix in names of output files
//
// Revision 1.23.2.2  2010-05-10 17:57:24  regan
// fixing types
//
// Revision 1.23.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.23  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.22  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
//
