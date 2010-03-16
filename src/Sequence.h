// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.h,v 1.23 2010-03-16 06:42:50 regan Exp $
//
#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <string>
#include <vector>

#include <tr1/memory>

#include "TwoBitSequence.h"

class Sequence {
public:
	static const char REF_QUAL = 0xff;
	static const char FASTQ_START_CHAR = 64;
	static const char PRINT_REF_QUAL = 126;

public:
	typedef TwoBitSequenceBase::SequenceLengthType SequenceLengthType;
	const static SequenceLengthType MAX_SEQUENCE_LENGTH =
			(SequenceLengthType) -1;
	typedef TwoBitSequenceBase::BaseLocationType BaseLocationType;
	typedef TwoBitSequenceBase::BaseLocationVectorType BaseLocationVectorType;

private:
	inline const Sequence &constThis() const {
		return *this;
	}

protected:
	SequenceLengthType _length;

	/*
	 _data contains a composite of 1 fields:
	 +0                 : the sequence as NCBI 2NA (2 bits per base ACGT)
	 += (length +3)/4   : non-ACGT bases: count followed by array of markups
	 */
	std::tr1::shared_ptr<unsigned char> _data;

	SequenceLengthType *_getMarkupBasesCount();
	const SequenceLengthType *_getMarkupBasesCount() const;
	BaseLocationType *_getMarkupBases();
	const BaseLocationType *_getMarkupBases() const;

	void reset();

	void setSequence(std::string fasta, unsigned int extraBytes);
public:

	Sequence();
	Sequence(std::string fasta);

	~Sequence();

	void setSequence(std::string fasta);

	SequenceLengthType getLength() const;
	std::string getFasta(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH) const;
	std::string getFastaNoMarkup() const;
	SequenceLengthType getMarkupBasesCount() const { return *_getMarkupBasesCount(); }
	BaseLocationVectorType getMarkups() const;

	SequenceLengthType getTwoBitEncodingSequenceLength() const;
	TwoBitEncoding *getTwoBitSequence();
	const TwoBitEncoding *getTwoBitSequence() const;

};

class BaseQual {
public:
	char base;
	char qual;
	BaseQual(char _base, char _qual) : base(_base), qual(_qual) {}
	BaseQual(char _base, double prob) : base(_base), qual() {
		qual = getQualChar(prob);
	}
	char getQualChar(double prob) {
		if (prob > 0.999) {
			return Sequence::FASTQ_START_CHAR + 30;
		} else if (prob > 0.99) {
			return Sequence::FASTQ_START_CHAR + 20;
		} else if (prob > 0.9) {
			return Sequence::FASTQ_START_CHAR + 10;
		} else {
			return Sequence::FASTQ_START_CHAR + 1;
		}
	}
};

class ProbabilityBase {
public:

	double a,c,g,t;
	double top;
	char best;
	ProbabilityBase() : a(0.0), c(0.0), g(0.0), t(0.0), top(0.0), best(' ') {}
	ProbabilityBase(const ProbabilityBase &copy) : a(copy.a), c(copy.c), g(copy.g), t(copy.t), top(copy.top), best(copy.best) {
		setTop(*this);
	}
	ProbabilityBase &operator=(const ProbabilityBase &copy) {
		a = copy.a;
		c = copy.c;
		g = copy.g;
		t = copy.t;
		setTop(*this);
		return *this;
	}
	void setTop(const ProbabilityBase &base) {
		setTop(base.a, base.c, base.g, base.t);
	}
	void setTop(double _a, double _c, double _g, double _t) {
		if (top < _a) {
			top = _a;
			best = 'A';
		}
		if (top < _c) {
			top = _c;
			best = 'C';
		}
		if (top < _g) {
			top = _g;
			best = 'G';
		}
		if (top < _t) {
			top = _t;
			best = 'T';
		}
	}
	ProbabilityBase operator+(const ProbabilityBase &other) {
		ProbabilityBase tmp(*this);
		tmp.a += other.a;
		tmp.c += other.c;
		tmp.g += other.g;
		tmp.t += other.t;
		tmp.setTop(other);
		return tmp;
	}
	ProbabilityBase &operator+=(const ProbabilityBase &other) {
		*this = *this + other;
		setTop(*this);
		return *this;
	}
	ProbabilityBase operator*(const ProbabilityBase &other) {
		ProbabilityBase tmp(*this);
		tmp.a *= other.a;
		tmp.c *= other.c;
		tmp.g *= other.g;
		tmp.t *= other.t;
		tmp.setTop(other);
		return tmp;
	}
	ProbabilityBase &operator*=(const ProbabilityBase &other) {
		*this = *this * other;
		setTop(*this);
		return *this;
	}
	ProbabilityBase &operator*=(double factor) {
		a *= factor;
		c *= factor;
		g *= factor;
		t *= factor;
		return *this;
	}
	double sum() const {
		return a+c+g+t;
	}
	double getA() const {
		if (best == 'A')
			return top;
		else
			return a;
	}
	double getC() const {
		if (best == 'C')
			return top;
		else
			return c;
	}
	double getG() const {
		if (best == 'G')
			return top;
		else
			return g;
	}
	double getT() const {
		if (best == 'T')
			return top;
		else
			return t;
	}

	BaseQual getBaseQual() const {
    	if (a > c) {
    		// A/G/T
    		if (a > g) {
    			// A/T
    			if (a > t) {
    				return BaseQual('A', getA());// A
    			} else {
    				return BaseQual('T', getT());// T
    			}
    		} else {
    			// G/T
    			if (g > t) {
    				return BaseQual('G', getG());// G
    			} else {
    				return BaseQual('T', getT());// T
    			}
    		}
    	} else {
    		// C/G/T
    		if (c > g) {
    			// C/T
    			if (c > t) {
    				return BaseQual('C', getC());// C
    			} else {
    				return BaseQual('T', getT());// T
    			}
    		} else {
    			// G/T
    			if (g > t) {
    				return BaseQual('G', getG());// G
    			} else {
    				return BaseQual('T', getT());// T
    			}
    		}
    	}
    	throw;
	}
};

class ProbabilityBases {
private:
	std::vector< ProbabilityBase > _bases;
public:
	ProbabilityBases(size_t _size) : _bases(_size) {}
	void resize(size_t _size) { _bases.resize(_size); }
	size_t size() const {
		return _bases.size();
	}
	double sum() const {
		double sum = 0.0;
		for(int i = 0; i < (int) _bases.size(); i++)
			sum += _bases[i].sum();
		return sum;
	}
	ProbabilityBase &operator[](size_t idx) {
		return _bases[idx];
	}
	const ProbabilityBase &operator[](size_t idx) const {
		return _bases[idx];
	}

	ProbabilityBases &operator+=(const ProbabilityBases &other) {
		if (_bases.size() < other._bases.size()) {
			resize(other._bases.size());
		}
		for(Sequence::SequenceLengthType i = 0 ; i < other._bases.size(); i++) {
			_bases[i] += other._bases[i];
		}
		return *this;
	}
	ProbabilityBases &operator*=(const ProbabilityBases &other) {
		if (_bases.size() < other._bases.size() ) {
			resize(other._bases.size());
		}
		for(int i = 0; i < (int) other._bases.size(); i++) {
			_bases[i] *= other._bases[i];
		}
		return *this;
	}
	ProbabilityBases &operator*=(double factor) {
		for(int i = 0; i < (int) _bases.size(); i++) {
		    _bases[i] *= factor;
		}
		return *this;
	}
};

class Read : public Sequence {


private:
	inline const Read &constThis() const {
		return *this;
	}

private:

	/*
	 _data is inherited and now contains a composite of 4 fields:
	 +0                    : the sequence as NCBI 2NA (2 bits per base ACGT)
	 += (length +3)/4      :  non-ACGT bases: count followed by array of markups
	 += getMarkupLength()  : qualities as 1 byte per base, 0 = N 1..255 Phred Quality Score.
	 += length             : null terminated name.
	 */

	char * _getQual();
	const char * _getQual() const;
	char * _getName();
	const char * _getName() const;
	SequenceLengthType _qualLength() const;

	static int qualityToProbabilityInitialized;
	static int
			initializeQualityToProbability(unsigned char minQualityScore = 0);

public:

	Read() :
		Sequence() {
	}
	;
	Read(std::string name, std::string fasta, std::string qualBytes);

	void setRead(std::string name, std::string fasta, std::string qualBytes);

	std::string getName() const;
	std::string getQuals(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH,
			bool forPrinting = false) const;
	void markupBases(SequenceLengthType offset, SequenceLengthType length, char mask = 'X');

	std::string toFastq(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH,
			std::string label = "") const;
	std::string getFormattedQuals(SequenceLengthType trimOffset =
			MAX_SEQUENCE_LENGTH) const;

	ProbabilityBases getProbabilityBases() const {
		std::string fasta = getFasta();
		std::string quals = getQuals();
		ProbabilityBases probs(fasta.length());
		for(int i = 0; i < (int) fasta.length(); i++) {
			double prob = qualityToProbability[ (unsigned char) quals[i] ];
			if (prob < 0.2501) {
				prob = 0.2501; // slightly better than random...
			}
			double otherProb = (1.0 - prob) / 3.0;
			ProbabilityBase &base = probs[i];
			switch (fasta[i]) {
			case 'A': base.a += prob; base.c += otherProb; base.g += otherProb; base.t += otherProb; break;
			case 'C': base.c += prob; base.a += otherProb; base.g += otherProb; base.t += otherProb; break;
			case 'G': base.g += prob; base.a += otherProb; base.c += otherProb; base.t += otherProb; break;
			case 'T': base.t += prob; base.a += otherProb; base.c += otherProb; base.g += otherProb; break;
			}
		}
		return probs;
	}
	double scoreProbabilityBases(const ProbabilityBases &probs) const {
		ProbabilityBases myprobs = getProbabilityBases();
		myprobs *= probs;
		return myprobs.sum();
	}

	static double qualityToProbability[256];
	static void setMinQualityScore(unsigned char minQualityScore);
};

#endif

//
// $Log: Sequence.h,v $
// Revision 1.23  2010-03-16 06:42:50  regan
// bugfixes
//
// Revision 1.22  2010-03-15 18:35:17  regan
// minor refactor and added consensus read
//
// Revision 1.21  2010-03-14 16:56:59  regan
// added probability base classes
//
// Revision 1.20  2010-03-08 22:14:38  regan
// replaced zero bases with markup bases to mask out reads that match the filter pattern
// bugfix in overrunning the mask
//
// Revision 1.19  2010-03-03 17:38:48  regan
// fixed quality scores
//
// Revision 1.18  2010-02-26 13:01:16  regan
// reformatted
//
// Revision 1.17  2010-02-26 11:09:43  regan
// bugfix
//
// Revision 1.16  2010-02-22 14:41:03  regan
// bugfix in printing
//
// Revision 1.15  2010-01-14 18:04:14  regan
// bugfixes
//
// Revision 1.14  2010-01-13 23:46:46  regan
// made const class modifications
//
// Revision 1.13  2010-01-13 21:16:00  cfurman
// added setMinQualityScore
//
// Revision 1.12  2010-01-13 00:24:30  regan
// use less memory for reference sequences and those without quality
//
// Revision 1.11  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.10  2010-01-05 06:44:39  regan
// fixed warnings
//
// Revision 1.9  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.8  2009-11-07 00:28:41  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
// Revision 1.7  2009-11-02 18:27:00  regan
// added getMarkups()
// added quality to probability lookup table
//
// Revision 1.6  2009-10-22 07:04:06  regan
// added a few unit tests
// minor refactor
//
// Revision 1.5  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
// Revision 1.4  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
// Revision 1.3  2009-10-21 00:00:58  cfurman
// working on kmers....
//
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//
