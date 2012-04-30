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
	AlignmentRecord(const AlignmentRecord &copy) {
		*this = copy;
	}
	AlignmentRecord &operator=(const AlignmentRecord &copy) {
		startPos = copy.startPos;
		endPos = copy.endPos;
		return *this;
	}
	bool isAligned() const {
		return startPos != endPos;
	}
	SequenceLengthType getOverlap() const {
		if (!isAligned())
			return 0;
		return (isReversed() ? (startPos - endPos) : (endPos - startPos)) + 1;
	}
	bool isReversed() const {
		return (startPos > endPos);
	}
	bool isAtEnd(const Read &read, SequenceLengthType dist = 0) const {
		if (!isAligned())
			return false;
		SequenceLengthType len = read.getLength();
		return isAtEnd(len, dist);
	}
	bool isAtEnd(SequenceLengthType len, SequenceLengthType dist = 0) const {
		bool ret = false;
		if (dist > len)
			dist = len-1;
		if (isReversed()) {
			if (endPos <= dist || startPos >= len-1-dist)
				ret = true;
		} else {
			if (startPos <= dist || endPos >= len-1-dist)
				ret = true;
		}
		LOG_DEBUG(5, "isAtEnd(" << len << "," << dist <<"): " << ret << " from " << toString());
		return ret;
	}
	bool contains(SequenceLengthType pos) const {
		bool ret;
		if (!isAligned()) {
			ret = false;
		} else {
			if (isReversed())
				ret = ((endPos <= pos) & (pos <= startPos));
			else
				ret = ((startPos <= pos) & (pos <= endPos));
		}
		return ret;
	}
	std::string toString() const {
		return "{" + boost::lexical_cast<std::string>(startPos) + "," + boost::lexical_cast<std::string>(endPos) + "}";
	}
	SequenceLengthType startPos, endPos;
};
class Alignment {
public:
	Alignment() : targetAln(), queryAln(), mismatches(0) {}
	Alignment(AlignmentRecord target, AlignmentRecord query) : targetAln(target), queryAln(query), mismatches(0) {}
	Alignment(const Alignment &copy) {
		*this = copy;
	}
	Alignment &operator=(const Alignment &copy) {
		targetAln = copy.targetAln;
		queryAln = copy.queryAln;
		mismatches = copy.mismatches;
		return *this;
	}
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
		if (!isAligned())
			return 0.0;
		return 1.0 - (float) mismatches / (float) getOverlap();
	}
	bool isAtEnd(const Read &target, const Read &query) const {
		bool ret = (targetAln.isAtEnd(target) | queryAln.isAtEnd(query));
		return ret;
	}
	std::string toString() const {
		return "Align{" + targetAln.toString() + "," + queryAln.toString() + "mis" + boost::lexical_cast<std::string>(mismatches) + "}";
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
		return getAlignment(_spectrum, _read, query);
	}
	const Read &getTarget() const {
		return _read;
	}
	static Alignment getBestAlignment(Alignment a, Alignment b) {
		if (a.getOverlap()*a.getIdentity() >= b.getOverlap()*b.getIdentity())
			return a;
		else
			return b;
	}
	// targetKS needs to have been built with Read, so that readIDs match
	static Alignment getAlignment( KS &targetKS, const Read &target, const Read &query) {
		assert(targetKS.hasSolids);
		Alignment bestAlignment;
		LOG_DEBUG(4, "getAlignment(): tgt: " << target.getName() << " query: " << query.getName());
		KmerWeights kmers(query.getTwoBitSequence(), query.getLength(), true);
		for(unsigned int j = 0; j < kmers.size(); j++) {
			KS::SolidElementType element = targetKS.getIfExistsSolid( kmers[j] );
			if (element.isValid()) {
				TrackingData::ReadPositionWeightVector rpwv = element.value().getEachInstance();
				LOG_DEBUG(5, "getAlignment(): isvalid " << j << " rpwv: " << rpwv.size());
				for(TrackingData::ReadPositionWeightVector::iterator it = rpwv.begin(); it != rpwv.end(); it++) {
					assert(it->readId == 0); // target is the only read
					Alignment test;
					SequenceLengthType ksize = KmerSizer::getSequenceLength();
					if ( !( bestAlignment.targetAln.contains(it->position) & bestAlignment.queryAln.contains(j) ) ) {
						test = getAlignment(target, it->position, query, j, ksize);
						bestAlignment = getBestAlignment(bestAlignment, test);
					}
				}
			} else {
				LOG_DEBUG(5, "getAlignment(): " << kmers[j].toFasta() << " was not present");
			}
		}
		LOG_DEBUG(4, "" << (bestAlignment.isAligned()? "Aligned ":"Did not Align ") << target.toFasta() << " to " << query.toFasta() << " " << bestAlignment.toString());
		return bestAlignment;
	}

	// fast zipper alignment to the end of one of the reads
	static Alignment getAlignment(const Read &target, SequenceLengthType tpos, const Read &query, SequenceLengthType qpos, SequenceLengthType minLen) {
		Alignment alignment;
		LOG_DEBUG(5, "getAlignment(" << target.getName() << ", " << tpos << ", " << query.getName() << ", " << qpos << ", " << minLen << ")");

		std::string tseq = target.getFasta();
		std::string qseq = query.getFasta();
		SequenceLengthType targetLen = tseq.length();
		SequenceLengthType queryLen = qseq.length();

		if (tpos + minLen > targetLen || qpos + minLen > queryLen) {
			LOG_DEBUG(3, "no alignment returned because match is out of bounds: " << tpos << " / " << targetLen << ", " << qpos << "/ " << queryLen << " len: " << minLen << " " << target.toString() << " " << query.toString());
			return alignment;
		}

		std::string tseqM = tseq.substr(tpos, minLen);
		bool revcomp = false;
		if (tseqM.compare( qseq.substr(qpos, minLen )) != 0 ) {
			std::string qseq2 = TwoBitSequence::getReverseComplementFasta(query.getTwoBitSequence(), queryLen);
			SequenceLengthType qpos2 = queryLen - qpos - minLen;
			if (tseqM.compare( qseq2.substr(qpos2, minLen )) != 0 ) {
				LOG_DEBUG_OPTIONAL(2, qpos2 + minLen < query.getFirstMarkupLength(), "did not find a match at " << tpos << "," << qpos << "'" << tseqM << "' vs '" << qseq.substr(qpos,minLen) << "' or '" << qseq2.substr(qpos2,minLen) << "' between " << target.toFasta() << " and " << query.toFasta());
				return alignment;
			}
			revcomp = true;
			qseq = qseq2;
			qpos = qpos2;
		}

		AlignmentRecord &q = alignment.queryAln;
		AlignmentRecord &t = alignment.targetAln;
		q = AlignmentRecord(qpos, qpos+minLen-1);
		t = AlignmentRecord(tpos, tpos+minLen-1);

		// extend left
		while (q.startPos > 0 && t.startPos > 0) {
			q.startPos--; t.startPos--;
			if (tseq[t.startPos] != qseq[q.startPos])
				alignment.mismatches++;
		}

		// extend right
		while (q.endPos < (queryLen-1) && t.endPos < (targetLen-1)) {
			q.endPos++; t.endPos++;
			if (tseq[t.endPos] != qseq[q.endPos])
				alignment.mismatches++;
		}

		// fix direction
		if (revcomp) {
			q = AlignmentRecord( queryLen - 1 - q.startPos, queryLen - 1 - q.endPos);
		}
		if (!alignment.isAligned())
			LOG_WARN(1, "Could not find the alignment for " << tpos << ", " << qpos << " " << target.toString() << query.toString());
		return alignment;
	}
private:
	Read _read;
	KS _spectrum;
};
#endif /* KMER_ALIGN_H_ */
