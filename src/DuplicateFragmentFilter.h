/*
 * DuplicateFragmentFilter.h
 *
 *  Created on: Oct 26, 2010
 *      Author: regan
 */

#ifndef DUPLICATEFRAGMENTFILTER_H_
#define DUPLICATEFRAGMENTFILTER_H_


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

	static void _buildDuplicateFragmentMap(KSV &ksv, ReadSet &reads, unsigned char bytes, bool useReverseComplement, bool paired, unsigned int startOffset = Options::getDeDupStartOffset()) {
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
                long skipped = 0;

	    #pragma omp parallel for reduction(+:skipped)
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
						skipped++;
						continue;
					}

					// create read1 + the reverse complement of read2 (1:rev2)
					// when useReverseComplement, it is represented as a kmer, and the leastcomplement of 1:rev2 and 2:rev1 will be stored
					// and properly account for duplicate fragment pairs

					Sequence::BaseLocationVectorType markups = read1.getMarkups();
					if (TwoBitSequence::firstMarkupX(markups) + startOffset < sequenceLength) {
						memcpy(kmer.getTwoBitSequence()       , read1.getTwoBitSequence() + (startOffset/4), bytes);
					} else {
						skipped++;
						continue;
					}
					markups = read2.getMarkups();
					if (TwoBitSequence::firstMarkupX(markups) + startOffset < sequenceLength) {
						TwoBitSequence::reverseComplement( read2.getTwoBitSequence() + (startOffset/4), kmer.getTwoBitSequence() + bytes, sequenceLength);
					} else {
						skipped++;
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
					if (read1.isDiscarded()) {
						skipped++;
						continue;
					}
					Sequence::BaseLocationVectorType markups = read1.getMarkups();
					if (TwoBitSequence::firstMarkupX(markups) + startOffset < sequenceLength) {
						memcpy(kmer.getTwoBitSequence()        , read1.getTwoBitSequence() + (startOffset/4), bytes);
					} else {
						skipped++;
						continue;
					}
					// store the readIdx (not the pairIdx)
					ksv[threadNum].append(kmerWeights, readIdx);
				}
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

		LOG_VERBOSE(1, "Duplicate Detection skipped: " << skipped);
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
					ks.consolidate(elem.key(), useWeights);
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
						ks.consolidate(elem.key(), useWeights);
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

		ReadSetSizeType affectedCount = 0;

		#pragma omp parallel reduction(+:affectedCount)
		for(KSWeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
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

		ReadSetSizeType affectedCount = 0;
		ReadSetSizeType pairSize = reads.getPairSize();
		int numThreads = omp_get_max_threads();
		ReadSet _threadNewReads[numThreads];

		#pragma omp parallel num_threads(numThreads) reduction(+:affectedCount)
		for(KSWeakIterator it = ks.weak.beginThreaded(); it != ks.weak.endThreaded(); it++) {
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

		// force release all memory before KmerSizer is called
		for (int i = 0 ; i < numThreads ; i++)
			ksv[i].reset();

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

};

#endif /* DUPLICATEFRAGMENTFILTER_H_ */
