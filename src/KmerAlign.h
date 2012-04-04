/*
 * KmerAlign.h
 *
 *  Created on: Apr 3, 2012
 *      Author: regan
 */

#ifndef KMER_ALIGN_H_
#define KMER_ALIGN_H_

#include "Kmer.h"
#include "KmerTrackingData.h"
#include "KmerSpectrum.h"
#include "DistributedFunctions.h"
#include "ReadSet.h"
#include "Options.h"

class _KmerAlignOptions  : public OptionsBaseInterface {
public:
	_KmerAlignOptions() {}
	~_KmerAlignOptions() {}
	void _resetDefaults() {
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("Kmer-Align Options");
		opts.add_options()

				;
		desc.add(opts);

	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		return ret;
	}
protected:
};
typedef OptionsBaseTemplate< _KmerAlignOptions > KmerAlignOptions;

class AlignmentRecord {
public:
	AlignmentRecord() : startPos(0), endPos(0) {}
	AlignmentRecord(SequenceLengthType _s, SequenceLengthType _e) : startPos(_s), endPos(_e) {}
	bool isAligned() const {
		return startPos != endPos;
	}
	SequenceLengthType getOverlap() const {
		return isReversed() ? (startPos - endPos) : (endPos - startPos);
	}
	bool isReversed() const {
		return (startPos < endPos);
	}
	bool isAtEnd(const Read &read) const {
		if (!isAligned())
			return false;
		SequenceLengthType len = read.getLength();
		return isAtEnd(len);
	}
	bool isAtEnd(SequenceLengthType len) const {
		if (isReversed()) {
			if (endPos == 0 || startPos == len-1)
				return true;
		} else {
			if (startPos == 0 || endPos == len-1)
				return true;
		}
		return false;
	}
	SequenceLengthType startPos, endPos;
};
class Alignment {
public:
	Alignment() : targetAln(), queryAln(), mismatches(0) {}
	Alignment(AlignmentRecord target, AlignmentRecord query) : targetAln(target), queryAln(query), mismatches(0) {}
	bool isAligned() const {
		return (targetAln.isAligned() & queryAln.isAligned());
	}
	SequenceLengthType getOverlap() const {
		return std::min(targetAln.getOverlap(), queryAln.getOverlap());
	}
	SequenceLengthType getMisMatches() const {
		return mismatches;
	}
	float getIdentity() const {
		return 1.0 - (float) mismatches / (float) getOverlap();
	}
	bool isAtEnd(const Read &target, const Read &query) const {
		return (targetAln.isAtEnd(target) | queryAln.isAtEnd(query));
	}

	AlignmentRecord targetAln;
	AlignmentRecord queryAln;
	SequenceLengthType mismatches;
};

class KmerAlign {
public:
	typedef TrackingDataWithLastRead TD;
	typedef KmerSpectrum<TD, TD, TD> KS;

	KmerAlign(const Read &target) : _read(target), _spectrum(target.getLength()) {
		_spectrum.setSolidOnly();
		_spectrum.buildKmerSpectrum(target);
	}
	Alignment getAlignment(const Read &query) {
		return getAlignment(_spectrum, query, query);
	}
	const Read &getTarget() const {
		return _read;
	}
	static Alignment getBestAlignment(Alignment a, Alignment b) {
		if (a.getOverlap()*a.getIdentity() > b.getOverlap()*b.getIdentity())
			return a;
		else
			return b;
	}
	// targetKS needs to have been built with Read, so that readIDs match
	static Alignment getAlignment(const KS &targetKS, const Read &target, const Read &query) {
		assert(targetKS.hasSolids);
		Alignment bestAlignment;
		KmerWeights kmers(query.getTwoBitSequence(), query.getLength(), true);
		for(unsigned int j = 0; j < kmers.size(); j++) {
			KS::SolidElementType element = targetKS.getIfExistsSolid( kmers[j] );
			if (element.isValid()) {
				TrackingData::ReadPositionWeightVector rpwv = element.value().getEachInstance();
				for(TrackingData::ReadPositionWeightVector::iterator it = rpwv.begin(); it != rpwv.end(); it++) {
					ReadSet::ReadSetSizeType globalReadIdx = it->readId;
					assert(globalReadIdx == 0); // target is the only read
					bestAlignment = getBestAlignment(bestAlignment, getAlignment(target, it->position, query, j, KmerSizer::getSequenceLength()));
				}
			}
		}
		return bestAlignment;
	}

	// fast zipper alignment to the end of one of the reads
	static Alignment getAlignment(const Read &target, SequenceLengthType tpos, const Read &query, SequenceLengthType qpos, SequenceLengthType minLen) {
		Alignment alignment;

		std::string tseq = target.getFastaTrim();
		std::string qseq = query.getFastaTrim();
		SequenceLengthType targetLen = tseq.length();
		SequenceLengthType queryLen = qseq.length();

		if (tpos + minLen > targetLen || qpos + minLen > queryLen) {
			LOG_DEBUG_OPTIONAL(2, true, "no alignment returned because match is out of bounds");
			return alignment;
		}

		std::string tseqM = tseq.substr(tpos, minLen);
		bool revcomp = false;
		if (tseqM.compare( qseq.substr(qpos, minLen )) != 0 ) {
			revcomp = true;
			qseq = TwoBitSequence::getReverseComplementFasta(query.getTwoBitSequence(), queryLen);
			qpos = queryLen - qpos - minLen;
			if (tseqM.compare( qseq.substr(qpos, minLen )) != 0 ) {
				LOG_WARN(1, "did not find a match between " << target.toString() << " and " << query.toString());
				return alignment;
			}
		}

		AlignmentRecord &q = alignment.queryAln;
		AlignmentRecord &t = alignment.targetAln;
		q = AlignmentRecord(qpos, qpos+minLen);
		t = AlignmentRecord(tpos, tpos+minLen);

		// extend left
		while (q.startPos > 0 && t.startPos > 0) {
			q.startPos--; t.startPos--;
			if (tseq[t.startPos] != qseq[q.startPos])
				alignment.mismatches++;
		}

		// extend right
		while (q.endPos < (queryLen-1) && t.endPos < (targetLen-1)) {
			q.endPos++; t.endPos++;
			if (tseq[t.startPos] != qseq[q.startPos])
				alignment.mismatches++;
		}

		// fix direction
		if (revcomp) {
			q = AlignmentRecord( queryLen - q.endPos, queryLen - q.startPos);
		}
		return alignment;
	}
private:
	Read _read;
	KS _spectrum;
};
#endif /* KMER_ALIGN_H_ */
