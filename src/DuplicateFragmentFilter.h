/*
 * DuplicateFragmentFilter.h
 *
 *  Created on: Oct 26, 2010
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

#ifndef DUPLICATEFRAGMENTFILTER_H_
#define DUPLICATEFRAGMENTFILTER_H_

#include "FilterKnownOddities.h"
#include "Options.h"


class _DuplicateFragmentFilterOptions : public OptionsBaseInterface {
public:
	_DuplicateFragmentFilterOptions() : deDupMode(0), deDupSingle(false), deDupConsensus(true), deDupEditDistance(0), deDupStartOffset(0),
	deDupLength(24) {}
	virtual ~_DuplicateFragmentFilterOptions() {}

	unsigned int &getDeDupEditDistance()
	{
		return deDupEditDistance;
	}

	unsigned int &getDeDupLength()
	{
		return deDupLength;
	}

	unsigned int &getDeDupMode()
	{
		return deDupMode;
	}

	bool &getDeDupSingle()
	{
		return deDupSingle;
	}

	bool &getDeDupConsensus() {
		return deDupConsensus;
	}

	unsigned int &getDeDupStartOffset()
	{
		return deDupStartOffset;
	}

	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		FilterKnownOdditiesOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		// *::_setOptions(desc,p);
		// This class does not support non-default options for:
		// GeneralOptions or FilterKnownOdditiesOptions
		po::options_description opts("Duplicate Fragment Filter (PCR duplicates)");
		opts.add_options()
						("dedup-mode", po::value<unsigned int>()->default_value(deDupMode), "if 0, no fragment de-duplication will occur.  if 1, single orientation (AB and BA are separated) will collapse to consensus. if 2, both orientations (AB and BA are the same) will collapse")

						("dedup-single", po::value<bool>()->default_value(deDupSingle), "if not set, no single read de-duplication will occur.  if set, then single read deduplication will occur")

						("dedup-consensus", po::value<bool>()->default_value(deDupConsensus), "if set then a single consensus read will be calculated, if not set, a random read will be selected")

						("dedup-edit-distance", po::value<unsigned int>()->default_value(deDupEditDistance), "if -1, no fragment de-duplication will occur, if 0, only exact match, ...")

						("dedup-start-offset", po::value<unsigned int>()->default_value(deDupStartOffset), "de-duplication start offset to find unique fragments, must be multiple of 4")

						("dedup-length", po::value<unsigned int>()->default_value(deDupLength), "de-duplication length to find unique fragments, must be multiple of 4 (length on each read of a paired fragment or doubled when in single-end mode)");
		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		// This class does not support non-default options for:
		// GeneralOptions or FilterKnownOdditiesOptions
		// ret &= *::_parseOptions(vm);
		// set dedup mode
		setOpt("dedup-mode", getDeDupMode());

		// set dedup single
		setOpt("dedup-single", getDeDupSingle());

		setOpt("dedup-consensus", getDeDupConsensus());

		// set dedup edit distance
		setOpt("dedup-edit-distance", getDeDupEditDistance());
		if (getDeDupEditDistance() > 1) {
			LOG_ERROR(1, "Unsupported option dedup-edit-distance > 1!" << std::endl << getDesc() << std::endl << "Unsupported option dedup-edit-distance > 1!");
			ret = false;
		}
		setOpt("dedup-start-offset", getDeDupStartOffset());
		setOpt("dedup-length", getDeDupLength());
		if (getDeDupStartOffset() % 4 != 0 || getDeDupLength() % 4 != 0) {
			LOG_ERROR(1, "Unsupported option dedup-start-offset and dedup-length must both be mulitples of 4!" << std::endl << getDesc() << std::endl << "Unsuppored option dedup-start-offset and dedup-length must both be mulitples of 4!");
			ret = false;
		}

		return ret;
	}
protected:
	unsigned int deDupMode;
	bool deDupSingle, deDupConsensus;
	unsigned int deDupEditDistance, deDupStartOffset, deDupLength;

};
typedef OptionsBaseTemplate< _DuplicateFragmentFilterOptions > DuplicateFragmentFilterOptions;

// NOTE
//
//  _KmerSpectrum must be <TD, TD, TD>;

class DuplicateFragmentFilter
{
public:
	typedef TrackingData::ReadPositionWeightVector RPW;
	typedef TrackingDataWithAllReads TD;
	typedef KmerSpectrum<TD, TD, TD> KS;
	typedef KS::Vector KSV;
	typedef KS::WeakMapType::ElementType KSElementType;
	typedef std::vector< KSElementType > KSElementVector;
	typedef KS::WeakIterator KSWeakIterator;
	typedef KSElementVector::reverse_iterator KSElementVector__reverse_iterator;
	typedef RPW::iterator RPWIterator;
	typedef ReadSet::Pair Pair;
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;

	static void _buildDuplicateFragmentMap(KSV &ksv, ReadSet &reads, unsigned char bytes, bool useReverseComplement, bool paired, unsigned int startOffset = DuplicateFragmentFilterOptions::getOptions().getDeDupStartOffset()) {
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
		long skippedDiscard = 0, skippedInvalid = 0, skippedTooShort = 0, skippedUnpaired = 0;

#pragma omp parallel for reduction(+:skippedDiscard) reduction(+:skippedTooShort) reduction(+:skippedUnpaired) reduction(+:skippedInvalid)
		for(long pairIdx = 0; pairIdx < pairSize; pairIdx++) {
			Pair &pair = reads.getPair(pairIdx);
			int threadNum = omp_get_thread_num();
			KmerWeights &kmerWeights = tmpKmerv[threadNum];
			Kmer &kmer = kmerWeights[0];

			if (paired && pair.isPaired() ) {
				if(reads.isValidRead(pair.read1) && reads.isValidRead(pair.read2)) {
					const Read &read1 = reads.getRead(pair.read1);
					const Read &read2 = reads.getRead(pair.read2);
					if (read1.isDiscarded() || read2.isDiscarded()) {
						skippedDiscard++;
						LOG_DEBUG(6, "Skipped Discarded Reads: \n" << read1.toFastq() << read2.toFastq());
						continue;
					}

					// create read1 + the reverse complement of read2 (1:rev2)
					// when useReverseComplement, it is represented as a kmer, and the leastcomplement of 1:rev2 and 2:rev1 will be stored
					// and properly account for duplicate fragment pairs

					SequenceLengthType readLength;
					readLength = read1.getFirstMarkupXLength();
					if (readLength >= sequenceLength + startOffset) {
						memcpy(kmer.getTwoBitSequence()       , read1.getTwoBitSequence() + (startOffset/4), bytes);
					} else {
						skippedTooShort++;
						LOG_DEBUG(6, "Skipped Read1 TooShort: \n" << read1.toFastq() << read2.toFastq());
						continue;
					}

					readLength = read2.getFirstMarkupXLength();
					if (readLength >= sequenceLength + startOffset) {
						TwoBitSequence::reverseComplement( read2.getTwoBitSequence() + (startOffset/4), kmer.getTwoBitSequence() + bytes, sequenceLength);
					} else {
						skippedTooShort++;
						LOG_DEBUG(6, "Skipped Read2 TooShort: \n" << read1.toFastq() << read2.toFastq());
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
				} else {
					skippedInvalid++;
					LOG_DEBUG(6, "Skipped Read(s) invalid");
					continue;
				}
			} else if ( pair.isSingle() && (!paired) ) {
				ReadSetSizeType readIdx = pair.lesser();
				if (reads.isValidRead(readIdx)) {
					const Read &read1 = reads.getRead(readIdx);
					if (read1.isDiscarded()) {
						skippedDiscard++;
						LOG_DEBUG(6, "Skipped (single) Discarded : \n" << read1.toFastq());
						continue;
					}

					SequenceLengthType readLength = read1.getFirstMarkupXLength();
					if (readLength >= sequenceLength + startOffset) {
						memcpy(kmer.getTwoBitSequence()        , read1.getTwoBitSequence() + (startOffset/4), bytes);
					} else {
						skippedTooShort++;
						LOG_DEBUG(6, "Skipped (single) TooShort: \n" << read1.toFastq());
						continue;
					}
					// store the readIdx (not the pairIdx)
					ksv[threadNum].append(kmerWeights, readIdx);
				} else {
					skippedInvalid++;
					LOG_DEBUG(6, "Skipped Read(s) invalid");
					continue;
				}
			} else {
				LOG_DEBUG(6, "Skipped Unpaired: " << pairIdx);
				skippedUnpaired++;
				continue;
			}
		}
		if (Log::isDebug(3)) {
			for (int i = 0; i < numThreads; i++) {
				LOG_DEBUG(3, "spectrum " << i);
				ksv[i].printHistograms();
			}
		} else {
			LOG_VERBOSE(2, "merging duplicate fragment spectrums" );
		}

		LOG_VERBOSE(1, "Duplicate Detection skipped " << (paired?"pairs":"reads") << ": " << (skippedDiscard + skippedInvalid + skippedTooShort + skippedUnpaired) << "\n"
				<< "\tDiscarded " << skippedDiscard << " TooShort " << skippedTooShort << " UnPaired " << skippedUnpaired << " Invalid " << skippedInvalid);
		KS::mergeVector(ksv, 1);
	}

	// TODO make useWeights an Option::
	static void _mergeNodesWithinEditDistance(KS &ks, unsigned int cutoffThreshold, unsigned int editDistance) {
		// TODO honor edit distance > 1
		LOG_VERBOSE(1, "Merging kmers within edit-distance of " << editDistance << " " << MemoryUtils::getMemoryUsage());

		// create a sorted set of elements with > cutoffThreshold count
		// (singletons will not be included in this round)
		KSElementVector elems;

#pragma omp parallel private(elems)
		{
			int threadId = omp_get_thread_num();
			for(KSWeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
				KSElementType &elem = *it;
				if (elem.isValid() && elem.value().getCount() >= cutoffThreshold)
					elems.push_back(elem);
			}

			if (threadId == 0)
				LOG_VERBOSE(2, "Sorting all nodes >= " << cutoffThreshold << " count: " << elems.size() << " " << MemoryUtils::getMemoryUsage());

			std::sort(elems.begin(), elems.end());


			if (threadId == 0)
				LOG_VERBOSE(2, "Merging elements. " << MemoryUtils::getMemoryUsage());

			for(KSElementVector__reverse_iterator it = elems.rbegin(); it != elems.rend(); it++) {
				KSElementType &elem = *it;
				if (elem.isValid() && elem.value().getCount() >= cutoffThreshold) {
					ks.consolidate(elem.key());
				}
			}

			if (threadId == 0)
				LOG_VERBOSE(2, "Merging elements below cutoffThreshold. " << MemoryUtils::getMemoryUsage() );

			// do not clear elements until all merging has completed
			elems.clear();

			// now merge those less than cutoff (singletons by default)
			for(KSWeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
				KSElementType &elem = *it;
				if (elem.isValid()) {
					TD &value = elem.value();
					if (value > 0 && value.getCount() < cutoffThreshold) {
						ks.consolidate(elem.key());
					}
				}
			}
		} // omp parallel

		LOG_DEBUG(1, MemoryUtils::getMemoryUsage() );
		ks.printHistograms();

	}
	static ReadSetSizeType _buildConsensusUnPairedReads(KS &ks, ReadSet &reads, ReadSet &newReads, unsigned int cutoffThreshold) {
		ReadSet::madviseMmapsRandom();
		LOG_VERBOSE(1, "Building consensus reads. ");
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		unsigned char minQual = Options::getOptions().getMinQuality();

		ReadSetSizeType affectedCount = 0;

#pragma omp parallel reduction(+:affectedCount)
		for(KSWeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
			if (it->value().getCount() >= cutoffThreshold) {
				RPW rpw = it->value().getEachInstance();

				ReadSet tmpReadSet1;

				ReadSetSizeType randomIdx = rpw.size();
                                if (DuplicateFragmentFilterOptions::getOptions().getDeDupConsensus()) {
					for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {

						ReadSetSizeType readIdx = rpwit->readId;

						const Read &read1 = reads.getRead(readIdx);
						tmpReadSet1.append( read1 );
	
					}

					Read consensus1 = tmpReadSet1.getConsensusRead(minQual);

					#pragma omp critical (BCUR_newReads)
					{
						newReads.append(consensus1);
					}
				} else {
					randomIdx = LongRand::rand() % tmpReadSet1.getSize();
					LOG_DEBUG_OPTIONAL(2, true, "Selected " << randomIdx << " out of " << rpw.size() << reads.getRead( reads.getPair((rpw.begin() + randomIdx)->readId).read1 ).getName());
				}
				affectedCount += rpw.size();

				ReadSetSizeType count = 0;
				for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {
					if (count++ == randomIdx)
						continue;	

					ReadSetSizeType readIdx = rpwit->readId;

					Read &read1 = reads.getRead(readIdx);
					read1.discard();
				}

			}
		}
		LOG_DEBUG(1, "Clearing duplicate pair map: " << MemoryUtils::getMemoryUsage() );
		ks.reset();
		LOG_VERBOSE(1, "Built " << newReads.getSize() << " new consensus reads: " <<  MemoryUtils::getMemoryUsage() );
		newReads.identifyPairs();
		return affectedCount;
	}

	static ReadSetSizeType _buildConsensusPairedReads(KS &ks, ReadSet &reads, ReadSet &newReads, unsigned int cutoffThreshold) {
		ReadSet::madviseMmapsRandom();
		LOG_VERBOSE(1, "Building consensus reads. ");
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		unsigned char minQual = Options::getOptions().getMinQuality();

		ReadSetSizeType affectedCount = 0;
		ReadSetSizeType pairSize = reads.getPairSize();
		int numThreads = omp_get_max_threads();
		ReadSet _threadNewReads[numThreads];

#pragma omp parallel num_threads(numThreads) reduction(+:affectedCount)
		for(KSWeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
			if (it->value().getCount() >= cutoffThreshold) {
				RPW rpw = it->value().getEachInstance();

				ReadSetSizeType randomIdx = rpw.size();

				if (DuplicateFragmentFilterOptions::getOptions().getDeDupConsensus()) {
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

					Read consensus1 = tmpReadSet1.getConsensusRead(minQual);
					Read consensus2 = tmpReadSet2.getConsensusRead(minQual);

					threadNewReads.append(consensus1);
					threadNewReads.append(consensus2);
				} else {
					randomIdx = LongRand::rand() % rpw.size();
					LOG_DEBUG_OPTIONAL(2, true, "Selected " << randomIdx << " out of " << rpw.size() << reads.getRead( reads.getPair((rpw.begin() + randomIdx)->readId).read1 ).getName());
				}

				affectedCount += 2 * rpw.size();
				ReadSetSizeType count = 0;
				for(RPWIterator rpwit = rpw.begin(); rpwit != rpw.end(); rpwit++) {
					if (count++ == randomIdx)
						continue;

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

		LOG_DEBUG(1, "Clearing duplicate pair map: " << MemoryUtils::getMemoryUsage() );
		ks.reset();
		LOG_VERBOSE(1, "Built " << newReads.getSize() << " new consensus reads: " <<  MemoryUtils::getMemoryUsage() );
		newReads.identifyPairs();
		return affectedCount;
	}

	static ReadSetSizeType _filterDuplicateFragments(ReadSet &reads, unsigned char bytes, unsigned int cutoffThreshold, unsigned int editDistance, bool paired) {

		int numThreads = omp_get_max_threads();
		KSV ksv(numThreads);

		LOG_VERBOSE(1, "Building " << (paired?"Paired":"Un-Paired") << " Duplicate Fragment Spectrum" );

		SequenceLengthType affectedCount = 0;

		bool useReverseComplement = (DuplicateFragmentFilterOptions::getOptions().getDeDupMode() == 2);

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

		if (FilterKnownOdditiesOptions::getOptions().getSkipArtifactFilter() == 0) {
			LOG_VERBOSE(1, "Preparing artifact filter on new consensus reads: ");
			FilterKnownOddities filter;
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			// ignore user settings for outputing any filtered new consensus reads
			int oldFilterOutput = FilterKnownOdditiesOptions::getOptions().getFilterOutput();
			FilterKnownOdditiesOptions::getOptions().getFilterOutput() = 0;

			LOG_VERBOSE(2, "Applying sequence artifact filter to new consensus reads");
			unsigned long filtered = filter.applyFilter(newReads);
			LOG_VERBOSE(1, "filter affected (trimmed/removed) " << filtered << " Reads ");;
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			// reset user settings
			FilterKnownOdditiesOptions::getOptions().getFilterOutput() = oldFilterOutput;
		}

		reads.append(newReads);

		// force release all memory before KmerSizer is called
		for (int i = 0 ; i < numThreads ; i++)
			ksv[i].reset();

		return affectedCount;
	}

	static ReadSetSizeType filterDuplicateFragments(ReadSet &reads, unsigned char sequenceLength = DuplicateFragmentFilterOptions::getOptions().getDeDupLength(), unsigned int cutoffThreshold = 2, unsigned int editDistance = DuplicateFragmentFilterOptions::getOptions().getDeDupEditDistance()) {

		if ( DuplicateFragmentFilterOptions::getOptions().getDeDupMode() == 0 || editDistance == (unsigned int) -1) {
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

		if (DuplicateFragmentFilterOptions::getOptions().getDeDupSingle() == 1)
			affectedCount += _filterDuplicateFragments(reads, bytes, cutoffThreshold, editDistance, false);

		KmerSizer::set(oldKmerSize);
		ReadSet::madviseMmapsNormal();

		return affectedCount;
	}

};

#endif /* DUPLICATEFRAGMENTFILTER_H_ */
