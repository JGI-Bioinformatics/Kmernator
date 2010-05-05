// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.h,v 1.27 2010-05-05 06:28:35 regan Exp $
//
#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <list>

#include "config.h"
#include "TwoBitSequence.h"
#include "Utils.h"
#include "LRUCache.h"

class Sequence {
public:
	static const char REF_QUAL = KoMer::REF_QUAL;
	static const char FASTQ_START_CHAR = KoMer::FASTQ_START_CHAR;
	static const char PRINT_REF_QUAL = KoMer::PRINT_REF_QUAL;

	static long constructCount;
	static long destructCount;

public:
	typedef TwoBitSequenceBase::SequenceLengthType SequenceLengthType;
	typedef TwoBitSequenceBase::SequenceLengthType1 SequenceLengthType1;
	typedef TwoBitSequenceBase::SequenceLengthType2 SequenceLengthType2;
	typedef TwoBitSequenceBase::BaseLocationType BaseLocationType;
	typedef TwoBitSequenceBase::BaseLocationType1 BaseLocationType1;
	typedef TwoBitSequenceBase::BaseLocationType2 BaseLocationType2;
	typedef TwoBitSequenceBase::BaseLocationVectorType BaseLocationVectorType;
	typedef TwoBitSequenceBase::MarkupElementSizeType MarkupElementSizeType;

	const static SequenceLengthType MAX_SEQUENCE_LENGTH = (SequenceLengthType) -1;
	typedef KoMer::RecordPtr RecordPtr;
	typedef TwoBitSequenceBase::TwoBitEncodingPtr DataPtr;

	typedef boost::shared_ptr< Sequence > SequencePtr;
	typedef std::pair< void*, SequencePtr > _CachedSequenceElement;
	typedef LRUCache< _CachedSequenceElement::first_type, _CachedSequenceElement::second_type > CachedSequences;
	typedef std::vector< CachedSequences > CachedSequencesVector;
	const static long maxCachePerThread = 31;

private:
	class DataPtrListVector
	{
		typedef std::list< DataPtr > DataPtrList;
		typedef std::vector< DataPtrList > _DataPtrListVector;
		_DataPtrListVector _vec;
		bool _isValid;
		DataPtrList &getList() {
			return _vec[omp_get_thread_num()];
		}
	public:
		static const long size = 511;
		DataPtrListVector() : _isValid(false) {
			_vec = _DataPtrListVector(OMP_MAX_THREADS, DataPtrList());
			_isValid = true;
		}
		~DataPtrListVector() {
			_isValid = false;
		}
		inline bool isValid() const { return _isValid; }
		void reset() {
			_isValid = false;
			_vec = _DataPtrListVector(OMP_MAX_THREADS, DataPtrList());
			_isValid = true;
		}
		DataPtr retrieveDataPtr();
		void returnDataPtr(DataPtr &ptr);
	};

	// dangling pointer
    // it is very important that ~Sequence::DataPtrListVector() NEVER gets called, as ~Sequence() inserts into it.
	static DataPtrListVector *preAllocatedDataPtrs;

public:
	static void clearCaches() {
		threadCacheSequences = CachedSequencesVector(OMP_MAX_THREADS, CachedSequences(Sequence::maxCachePerThread));
		preAllocatedDataPtrs->reset();
	}


private:
	inline const Sequence &constThis() const {
		return *this;
	}

protected:
	static CachedSequencesVector threadCacheSequences;

protected:
	char _flags;
	static const char MMAPED       = 0x80;
	static const char MARKUPS1     = 0x40;
	static const char MARKUPS2     = 0x20;
	static const char MARKUPS4     = MARKUPS1|MARKUPS2; // 0x60
	static const char HASQUALS     = 0x10;
	static const char PAIRED       = 0x08;
	static const char DISCARDED    = 0x04;
	static const char PREALLOCATED = 0x02;
	static const char INVALID      = 0x01;
	// 0x80 - mmaped data
	// 0x40 - markups in unsigned char
	// 0x20 - markups in unsigned short
	//      - 0x40|0x20 (0x60) markups in unsigned int
	// 0x10 - hasQuals
	// 0x08 - mmaped quals
	// 0x04 - isDiscarded
	// 0x01 - invalid state

	DataPtr _data;
	// if mmaped, _data consists of:
	//   +0 : char * = pointerToMmapedFileRecordStart
	//   if mmaped && hasquals
	//     + sizeof(char*) : char * = pointerToMmapedQualRecordStart
	//
	//   if markups4
	//     +sizeof(char*) : unsigned int = markupsCount
	//     +sizeof(unsigned int) : = BaseLocationType[markupsCount]
	//   else if markups1
	//     +sizeof(char*) : unsigned char = markupsCount
	//     +sizeof(unsigned char) : = BaseLocationType1[markupsCount]
	//   else if markups2
	//     +sizeof(char*) : unsigned short = markupsCount
	//     +sizeof(unsigned short) : = BaseLocationType2[markupsCount]
    //
	// else, _data consists of:
	//   +0 : SequenceLengthType = sequenceLength
	//   +sizeof(SequenceLengthType) : = TwoBitEncoding[ (sequenceLength+3)/4 ]
	//   if markups 1,2 or 4....
	//     +(sequenceLength+3)/4 : markupsCount 1,2 or 4 & BaseLocationType 1,2 or 4 [markupsCount]
	//

	void setFlag(char f);
	void unsetFlag(char f);

	RecordPtr *_getRecord();
	const RecordPtr *_getRecord() const;

	RecordPtr *_getQualRecord();
	const RecordPtr *_getQualRecord() const;

	SequenceLengthType *_getLength();
	const SequenceLengthType *_getLength() const;

	TwoBitEncoding *_getTwoBitSequence();
	const TwoBitEncoding *_getTwoBitSequence() const;

	SequenceLengthType _getStoredMarkupBasesCount() const;

	SequenceLengthType *_getMarkupBasesCount();
	const SequenceLengthType *_getMarkupBasesCount() const;

	SequenceLengthType1 *_getMarkupBasesCount1();
	const SequenceLengthType1 *_getMarkupBasesCount1() const;

	SequenceLengthType2 *_getMarkupBasesCount2();
	const SequenceLengthType2 *_getMarkupBasesCount2() const;

	BaseLocationType *_getMarkupBases();
	const BaseLocationType *_getMarkupBases() const;

	BaseLocationType1 *_getMarkupBases1();
	const BaseLocationType1 *_getMarkupBases1() const;

	BaseLocationType2 *_getMarkupBases2();
	const BaseLocationType2 *_getMarkupBases2() const;

	BaseLocationVectorType _getMarkups() const;

	void reset(char flags = 0);

	void setSequence(std::string fasta, long extraBytes, bool usePreAllocation = false);
	void setSequence(RecordPtr mmapRecordStart, const BaseLocationVectorType &markups, long extraBytes, RecordPtr mmapQualRecordStart = NULL);
	void setSequence(RecordPtr mmapRecordStart, long extraBytes, RecordPtr mmapQualRecordStart = NULL);

	void setMarkups(MarkupElementSizeType markupElementSize, const BaseLocationVectorType &markups);

	static CachedSequences &getCachedSequencesForThread() {
		int threadNum = omp_get_thread_num();
		return threadCacheSequences[threadNum];
	}
	SequencePtr getCache() const {
		assert(isMmaped());

		CachedSequences &cache = getCachedSequencesForThread();
		SequencePtr cachedSequence = cache.fetch( (void*)_data.get()->get() );
		if ( cachedSequence.get() == NULL ) {
			return setCache();
		} else {
			return cachedSequence;
		}
	}
	SequencePtr setCache() const {
		SequencePtr ptr;
		return setCache(ptr);
	}
	SequencePtr &setCache(SequencePtr &expandedSequence) const {
		assert(isMmaped());

		if (expandedSequence.get() == NULL)
			expandedSequence = readMmaped(true);
		CachedSequences &cache = getCachedSequencesForThread();
		cache.insert( (void*)_data.get()->get(), expandedSequence);
		return expandedSequence;
	}

public:
	static void setThreadCache(Sequence &mmapedSequence, SequencePtr &expandedSequence) {
		mmapedSequence.setCache(expandedSequence);
	}

public:

	Sequence();
	Sequence(const Sequence &copy) {
		*this = copy;
	}
	Sequence(std::string fasta, bool usePreAllocation = false);
	Sequence(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	Sequence &operator=(const Sequence &other) {
		if (this == &other)
			return *this;
		_flags = other._flags;
		_data = other._data;
		return *this;
	}
	Sequence clone() const {
		return Sequence(getFasta());
	}

	~Sequence();


public:
	inline bool isMmaped()       const { return (_flags & MMAPED)        == MMAPED; }
	inline bool isMarkups4()     const { return (_flags & MARKUPS4)      == MARKUPS4; }
	inline bool isMarkups2()     const { return (_flags & MARKUPS4)      == MARKUPS2; }
	inline bool isMarkups1()     const { return (_flags & MARKUPS4)      == MARKUPS1; }
	inline bool hasMarkups()     const { return (_flags & MARKUPS4)      != 0; }
	inline bool hasQuals()       const { return (_flags & HASQUALS)      == HASQUALS; }
	inline bool isPaired()       const { return (_flags & PAIRED)        == PAIRED; }
	inline bool isDiscarded()    const { return (_flags & DISCARDED)     == DISCARDED; }
	inline bool isPreAllocated() const { return (_flags & PREALLOCATED)  == PREALLOCATED; }
	inline bool isValid()        const { return (_flags & INVALID)       == 0; }

	void setSequence(std::string fasta, bool usePreAllocation = false);
	void setSequence(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	void setSequence(RecordPtr mmapRecordStart, const BaseLocationVectorType &markups, RecordPtr mmapQualRecordStart = NULL);

	void discard() {
		setFlag(DISCARDED);
	}
	void unDiscard() {
		unsetFlag(DISCARDED);
	}
	void markPaired() {
		setFlag(PAIRED);
	}

	SequenceLengthType getLength() const;

	const RecordPtr getRecord() const {
		assert(isMmaped());
		return *_getRecord();
	}
	RecordPtr getRecord() {
		return const_cast<RecordPtr> (constThis().getRecord());
	}
	const RecordPtr getQualRecord() const {
		assert(isMmaped() && hasQuals());
		return *_getQualRecord();
	}
	RecordPtr getQualRecord() {
		return const_cast<RecordPtr> (constThis().getQualRecord());
	}

	std::string getFasta(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH) const;
	std::string getFastaNoMarkup(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH) const;

	SequenceLengthType getMarkupBasesCount() const;
	BaseLocationVectorType getMarkups() const;

	SequenceLengthType getTwoBitEncodingSequenceLength() const;
	TwoBitEncoding *getTwoBitSequence() { return _getTwoBitSequence(); }
	const TwoBitEncoding *getTwoBitSequence() const { return _getTwoBitSequence(); }

	void readMmaped(std::string &name, std::string &bases, std::string &quals) const {
		assert(isMmaped());
		// TODO fix hack on NULL lastPtr.  Presently only works for single-lined fastas
		RecordPtr record(getRecord()), lastRecord(NULL), qualRecord(NULL), lastQualRecord(NULL);
		if (hasQuals()) {
		    qualRecord = getQualRecord();
			lastQualRecord = NULL;
		}
		SequenceRecordParser::parse(record, lastRecord, name, bases, quals, qualRecord, lastQualRecord);
	}
	SequencePtr readMmaped(bool usePreAllocation = false) const {
		std::string name, bases, quals;
		readMmaped(name, bases, quals);
		return SequencePtr(new Sequence(bases, usePreAllocation));
	}

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

public:
	typedef boost::shared_ptr< Sequence > ReadPtr;

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

	Read() : Sequence() {}
	Read(const Read &copy) {
		*this = copy;
	}
	Read(std::string name, std::string fasta, std::string qualBytes, bool usePreAllocation = false);
	Read(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	Read(RecordPtr mmapRecordStart, std::string markupFasta, RecordPtr mmapQualRecordStart = NULL);
	Read &operator=(const Read &other) {
		Sequence::operator=(other);
	    // there are no extra data members
	    return *this;
	}
	Read clone() const {
		return Read(getName(), getFasta(), getQuals());
	}

	void setRead(std::string name, std::string fasta, std::string qualBytes, bool usePreAllocation = false);
	void setRead(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	void setRead(RecordPtr mmapRecordStart, std::string markupFasta, RecordPtr mmapQualRecordStart);

	ReadPtr readMmaped(bool usePreAllocation = false) const {
		BaseLocationVectorType markups = _getMarkups();
		std::string name, bases, quals;
		Sequence::readMmaped(name, bases, quals);
		TwoBitSequence::applyMarkup(bases, markups);
		return ReadPtr(new Read(name, bases, quals, usePreAllocation));
	}

	bool recordHasQuals() const {
		assert(isMmaped());
		if (hasQuals())
			return true;
		else
			// TODO make this more general
		    return *getRecord() == '@'; // FASTQ
	}


	std::string getName() const;
	std::string getQuals(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH,
			bool forPrinting = false, bool unmasked = false) const;
	void markupBases(SequenceLengthType offset, SequenceLengthType length, char mask = 'X');

	std::string toFastq(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH,
			std::string label = "", bool unmasked = false) const;
	std::string toFasta(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH,
			std::string label = "", bool unmasked = false) const;
	std::string toQual(SequenceLengthType trimOffset, std::string label) const;

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
			case 'A':
			case 'a': base.a += prob; base.c += otherProb; base.g += otherProb; base.t += otherProb; break;
			case 'C':
			case 'c': base.c += prob; base.a += otherProb; base.g += otherProb; base.t += otherProb; break;
			case 'G':
			case 'g': base.g += prob; base.a += otherProb; base.c += otherProb; base.t += otherProb; break;
			case 'T':
			case 't': base.t += prob; base.a += otherProb; base.c += otherProb; base.g += otherProb; break;
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

	// format == 0 fastq
	// format == 1 fasta
	// format == 2 fastq unmasked
	// format == 3 fasta unmasked
	inline static std::ostream &write(std::ostream &os, const Read &read,
			SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH, std::string label = "", int format = 0) {
		switch (format) {
		case 0: os << read.toFastq(trimOffset, label); break;
		case 1: os << read.toFasta(trimOffset, label); break;
		case 2: os << read.toFastq(trimOffset, label, true) ; break;
		case 3: os << read.toFasta(trimOffset, label, true); break;
		}
		return os;
	}
	inline std::ostream &write(std::ostream &os,
			SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH, std::string label = "", int format = 0) const {
		return write(os, *this, trimOffset, label, format);
	}

};

#endif

//
// $Log: Sequence.h,v $
// Revision 1.27  2010-05-05 06:28:35  regan
// merged changes from FixPairOutput-20100504
//
// Revision 1.26.4.1  2010-05-05 05:57:53  regan
// fixed pairing
// fixed name to exclude labels and comments after whitespace
// applied some performance optimizations from other branch
// created FixPair application
//
// Revision 1.26.2.3  2010-05-03 21:34:07  regan
// fixed clone method
//
// Revision 1.26.2.2  2010-05-02 05:39:47  regan
// added method to use PAIRED flag
//
// Revision 1.26.2.1  2010-05-02 04:38:30  regan
// replaced mmap quals flag with paired
//
// Revision 1.26  2010-05-01 21:57:54  regan
// merged head with serial threaded build partitioning
//
// Revision 1.25.2.8  2010-05-01 21:28:50  regan
// bugfix
//
// Revision 1.25.2.7  2010-04-30 23:53:14  regan
// attempt to fix a bug.  clearing Sequence caches when it makes sense
//
// Revision 1.25.2.6  2010-04-30 23:33:37  regan
// replaced read cache with LRU 3rd party cache
//
// Revision 1.25.2.5  2010-04-30 22:29:58  regan
// added comments about dangling pointer
//
// Revision 1.25.2.4  2010-04-30 21:53:52  regan
// reuse memory efficiently for cache lookups
//
// Revision 1.25.2.3  2010-04-29 06:58:42  regan
// reworked output parameters to include option to print unmasked reads
//
// Revision 1.25.2.2  2010-04-29 04:26:32  regan
// bugfix in output filenames and content
//
// Revision 1.25.2.1  2010-04-28 22:28:11  regan
// refactored writing routines
//
// Revision 1.25  2010-04-21 00:33:19  regan
// merged with branch to detect duplicated fragment pairs with edit distance
//
// Revision 1.24.2.1  2010-04-16 23:46:06  regan
// checkpoint
//
// Revision 1.24  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.23.2.12.2.3  2010-04-16 21:38:18  regan
// removed commented code
//
// Revision 1.23.2.12.2.2  2010-04-16 17:49:42  regan
// fixed comments
//
// Revision 1.23.2.12.2.1  2010-04-16 05:29:59  regan
// checkpoint.. broke it
//
// Revision 1.23.2.12  2010-04-15 21:31:50  regan
// bugfix in markups and duplicate fragment filter
//
// Revision 1.23.2.11  2010-04-15 17:29:02  regan
// checkpoint, working with some optimizations
//
// Revision 1.23.2.10  2010-04-14 22:36:06  regan
// round of bugfixes
//
// Revision 1.23.2.9  2010-04-14 17:51:43  regan
// checkpoint
//
// Revision 1.23.2.8  2010-04-14 05:35:37  regan
// checkpoint. compiles but segfaults
//
// Revision 1.23.2.7  2010-04-14 03:51:19  regan
// checkpoint. compiles but segfaults
//
// Revision 1.23.2.6  2010-04-12 22:37:47  regan
// checkpoint
//
// Revision 1.23.2.5  2010-04-12 20:59:45  regan
// mmap checkpoint
//
// Revision 1.23.2.4  2010-04-05 05:42:53  regan
// checkpoint mmaping input files
//
// Revision 1.23.2.3  2010-04-05 02:56:08  regan
// bugfixes
//
// Revision 1.23.2.2  2010-04-04 15:58:02  regan
// fixed assertion code to obey debug rules
//
// Revision 1.23.2.1  2010-04-04 15:31:27  regan
// checkpoint - refactored underlying data structure t compress markups
//
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
