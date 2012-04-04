//
// Kmernator/src/Sequence.h
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

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <cstring>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <list>

#include "config.h"
#include "TwoBitSequence.h"
#include "Utils.h"
#include "LRUCache.h"

class Sequence {
public:
	static const char REF_QUAL = Kmernator::REF_QUAL;
	static const char PRINT_REF_QUAL = Kmernator::PRINT_REF_QUAL;
	static char FASTQ_START_CHAR;

public:
	typedef TwoBitSequenceBase::SequenceLengthType SequenceLengthType;
	typedef TwoBitSequenceBase::SequenceLengthType1 SequenceLengthType1;
	typedef TwoBitSequenceBase::SequenceLengthType2 SequenceLengthType2;
	typedef TwoBitSequenceBase::BaseLocationType BaseLocationType;
	typedef TwoBitSequenceBase::BaseLocationType1 BaseLocationType1;
	typedef TwoBitSequenceBase::BaseLocationType2 BaseLocationType2;
	typedef TwoBitSequenceBase::BaseLocationVectorType BaseLocationVectorType;
	typedef TwoBitSequenceBase::MarkupElementSizeType MarkupElementSizeType;

	typedef Kmernator::RecordPtr RecordPtr;
	typedef TwoBitSequenceBase::TwoBitEncodingPtr DataPtr;

	typedef boost::shared_ptr< Sequence > SequencePtr;
	typedef std::pair< const void*, SequencePtr > _CachedSequenceElement;
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
		inline DataPtrList &getList() {
			return _vec[omp_get_thread_num()];
		}
	public:
		static const unsigned long size = 511;
		DataPtrListVector() : _isValid(false) {
			_vec = _DataPtrListVector(omp_get_max_threads(), DataPtrList());
			_isValid = true;
		}
		~DataPtrListVector() {
			_isValid = false;
		}
		inline bool isValid() const { return _isValid; }
		void reset();
		DataPtr retrieveDataPtr();
		void returnDataPtr(DataPtr &ptr);
	};

	// dangling pointer
	// it is very important that ~Sequence::DataPtrListVector() NEVER gets called, as ~Sequence() inserts into it.
	static DataPtrListVector *preAllocatedDataPtrs;

public:
	static void clearCaches();


private:
	inline const Sequence &constThis() const {
		return *this;
	}

protected:
	static CachedSequencesVector threadCacheSequences;

protected:
	// TODO if mmaped w/o markups, can this be a regular pointer with no allocation?
	DataPtr _data;
	// TODO can _flags be embedded into _data -- save in memory alignment padding
	// -- move flags into ReadSet? as second parallel vector?
	char _flags;

	static const char MMAPED       = 0x80;
	static const char MARKUPS1     = 0x40;
	static const char MARKUPS2     = 0x20;
	static const char MARKUPS4     = MARKUPS1|MARKUPS2; // 0x60
	static const char HASQUALS     = 0x10;
	static const char PAIRED       = 0x08;
	static const char DISCARDED    = 0x04;
	static const char PREALLOCATED = 0x02;
	static const char HASFASTAQUAL = 0x01;
	// 0x80 - mmaped data
	// 0x40 - markups in unsigned char
	// 0x20 - markups in unsigned short
	//      - 0x40|0x20 (0x60) markups in unsigned int
	// 0x10 - hasQuals
	// 0x08 - mmaped quals
	// 0x04 - isDiscarded
	// 0x01 - hasFastaQual

	// if mmaped, _data consists of:
	//   +0 : char * = pointerToMmapedFileRecordStart
	//   if mmaped && hasFastaQual
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

	void *_getData();
	const void *_getData() const;

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

	virtual const void *_getEnd() const;

	BaseLocationVectorType _getMarkups() const;

	void reset(char flags = 0);

	void setSequence(std::string fasta, long extraBytes, bool usePreAllocation = false);
	void setSequence(RecordPtr mmapRecordStart, const BaseLocationVectorType &markups, long extraBytes, RecordPtr mmapQualRecordStart = NULL);
	void setSequence(RecordPtr mmapRecordStart, long extraBytes, RecordPtr mmapQualRecordStart = NULL);

	void setMarkups(MarkupElementSizeType markupElementSize, const BaseLocationVectorType &markups);

	static CachedSequences &getCachedSequencesForThread();

	SequencePtr getCache() const;
	SequencePtr setCache() const;
	SequencePtr &setCache(SequencePtr &expandedSequence) const;

public:
	static void setThreadCache(Sequence &mmapedSequence, SequencePtr &expandedSequence);

public:

	Sequence();
	Sequence(const Sequence &copy);
	Sequence(std::string fasta, bool usePreAllocation = false);
	Sequence(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	Sequence &operator=(const Sequence &other);
	Sequence clone(bool usePreAllocation = false) const;

	virtual ~Sequence();


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
	inline bool hasFastaQual()   const { return (_flags & HASFASTAQUAL)  == HASFASTAQUAL; }
	inline bool isValid()        const { return ( _getData() != NULL ); }

	virtual long getStoreSize() const;
	virtual long store(void *dst) const;
	virtual void *restore(void *src, long size);

	void setSequence(std::string fasta, bool usePreAllocation = false);
	void setSequence(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	void setSequence(RecordPtr mmapRecordStart, const BaseLocationVectorType &markups, RecordPtr mmapQualRecordStart = NULL);

	inline void discard()        {  setFlag(DISCARDED); }
	inline void unDiscard()      { unsetFlag(DISCARDED);}
	inline void markPaired()     { setFlag(PAIRED);     }

	SequenceLengthType getLength() const;

	const RecordPtr getRecord() const;
	inline RecordPtr getRecord() {
		return const_cast<RecordPtr> (constThis().getRecord());
	}

	const RecordPtr getQualRecord() const;
	inline RecordPtr getQualRecord() {
		return const_cast<RecordPtr> (constThis().getQualRecord());
	}

	std::string getFasta(SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH) const;
	std::string getFastaNoMarkup(SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH) const;
	std::string getFastaTrim() const {
		return getFasta(0, getFirstMarkupLength());
	}

	SequenceLengthType getMarkupBasesCount() const;
	BaseLocationVectorType getMarkups() const;

	SequenceLengthType getFirstMarkupLength() const;
	SequenceLengthType getFirstMarkupNLength() const;
	SequenceLengthType getFirstMarkupXLength() const;
	SequenceLengthType getFirstMarkupNorXLength() const;

	SequenceLengthType getTwoBitEncodingSequenceLength() const;
	inline TwoBitEncoding *getTwoBitSequence() { return _getTwoBitSequence(); }
	inline const TwoBitEncoding *getTwoBitSequence() const { return _getTwoBitSequence(); }

	void readMmaped(std::string &name, std::string &bases, std::string &quals) const;
	SequencePtr readMmaped(bool usePreAllocation = false) const;

};

class BaseQual {
public:
	char base;
	char qual;
	BaseQual(char _base, char _qual) : base(_base), qual(_qual) {}
	BaseQual(char _base, double prob) : base(_base), qual() {
		qual = getQualChar(prob);
	}
	static char getQualChar(double prob, bool ignoreLow = false);
};

class ProbabilityBase {
public:
	double a,c,g,t;
	double top;
	char best;
	short count;

public:
	ProbabilityBase() : a(0.0), c(0.0), g(0.0), t(0.0), top(0.0), best(' '), count(0) {}
	ProbabilityBase(const ProbabilityBase &copy) : a(copy.a), c(copy.c), g(copy.g), t(copy.t), top(copy.top), best(copy.best), count(copy.count) {
		setTop(*this);
	}
	ProbabilityBase &operator=(const ProbabilityBase &copy);

	void observe(char nuc, double prob);

	void setTop(double _a, double _c, double _g, double _t);
	void setTop(const ProbabilityBase &base) {
		setTop(base.a, base.c, base.g, base.t);
	}
	ProbabilityBase operator+(const ProbabilityBase &other);
	ProbabilityBase &operator+=(const ProbabilityBase &other);
	ProbabilityBase operator*(const ProbabilityBase &other);
	ProbabilityBase &operator*=(const ProbabilityBase &other);
	ProbabilityBase &operator*=(double factor);
	inline double sum() const {
		return a+c+g+t;
	}
	double getA() const;
	double getC() const;
	double getG() const;
	double getT() const;

	BaseQual getBaseQual() const;
};

class ProbabilityBases {
private:
	std::vector< ProbabilityBase > _bases;
public:
	ProbabilityBases(size_t _size) : _bases(_size) {}
	void resize(size_t _size) { _bases.resize(_size); }
	inline size_t size() const {
		return _bases.size();
	}
	double sum() const;
	ProbabilityBase &operator[](size_t idx);
	const ProbabilityBase &operator[](size_t idx) const;

	ProbabilityBases &operator+=(const ProbabilityBases &other);
	ProbabilityBases &operator*=(const ProbabilityBases &other);
	ProbabilityBases &operator*=(double factor);
	std::string toString() const;
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

	virtual const void *_getEnd() const;

	static int qualityToProbabilityInitialized;
	static int
	initializeQualityToProbability(unsigned char minQualityScore = 0, char startChar = Kmernator::FASTQ_START_CHAR_ILLUMINA);

public:

	Read() : Sequence() {}
	Read(const Read &copy) {
		*this = copy;
	}
	Read(std::string name, std::string fasta, std::string qualBytes, bool usePreAllocation = false);
	Read(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	Read(RecordPtr mmapRecordStart, std::string markupFasta, RecordPtr mmapQualRecordStart = NULL);
	virtual ~Read() {}

	Read &operator=(const Read &other);
	Read clone(bool usePreAllocation = false) const;

	void setRead(std::string name, std::string fasta, std::string qualBytes, bool usePreAllocation = false);
	void setRead(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	void setRead(RecordPtr mmapRecordStart, std::string markupFasta, RecordPtr mmapQualRecordStart);

	ReadPtr readMmaped(bool usePreAllocation = false) const;
	bool recordHasQuals() const ;

	std::string getName() const;
	void setName(const std::string name);

	std::string getQuals(SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH,
			bool forPrinting = false, bool unmasked = false) const;
	void markupBases(SequenceLengthType offset, SequenceLengthType length, char mask = 'X');

	std::string toFastq(SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH,
			std::string label = "", bool unmasked = false) const;
	std::string toFasta(SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH,
			std::string label = "", bool unmasked = false) const;
	std::string toQual(SequenceLengthType trimOffset, SequenceLengthType trimLength, std::string label) const;

	std::string getFormattedQuals(SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH) const;

	ProbabilityBases getProbabilityBases(unsigned char minQual = Options::getOptions().getMinQuality()) const;
	double scoreProbabilityBases(const ProbabilityBases &probs) const;

	static double qualityToProbability[256];
	static void setMinQualityScore(unsigned char minQualityScore, char startChar = Kmernator::FASTQ_START_CHAR_ILLUMINA);

	// format == 0 fastq
	// format == 1 fasta
	// format == 2 fastq unmasked
	// format == 3 fasta unmasked
	inline static std::ostream &write(std::ostream &os, const Read &read,
			SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH, std::string label = "", FormatOutput format = FormatOutput::getDefault()) {
		switch (format.getType()) {
		case FormatOutput::FASTQ: os << read.toFastq(trimOffset, trimLength, label); break;
		case FormatOutput::FASTA: os << read.toFasta(trimOffset, trimLength, label); break;
		case FormatOutput::FASTQ_UNMASKED: os << read.toFastq(trimOffset, trimLength, label, true) ; break;
		case FormatOutput::FASTA_UNMASKED: os << read.toFasta(trimOffset, trimLength, label, true); break;
		default: LOG_THROW("Invalid format for Sequence::write(): " << format.getType());
		}
		return os;
	}
	inline std::ostream &write(std::ostream &os,
			SequenceLengthType trimOffset = 0, SequenceLengthType trimLength = MAX_SEQUENCE_LENGTH, std::string label = "", FormatOutput format = FormatOutput::getDefault()) const {
		return write(os, *this, trimOffset, trimLength, label, format);
	}

	static long getIntendedWriteSize(const Read &read, SequenceLengthType sequenceLength, std::string label = "", FormatOutput format = FormatOutput::getDefault()) {
		long length = 0;
		switch (format.getType()) {
		case FormatOutput::FASTQ: length += 2 + 4 + sequenceLength*2 + read.getName().length() + (label.length() > 0 ? 1 + label.length() : 0); break;
		case FormatOutput::FASTA: length += 1 + 2 + sequenceLength   + read.getName().length() + (label.length() > 0 ? 1 + label.length() : 0); break;
		case FormatOutput::FASTQ_UNMASKED: length += 2 + 4 + read.getLength()*2 + read.getName().length() + (label.length() > 0 ? 1 + label.length() : 0); break;
		case FormatOutput::FASTA_UNMASKED: length += 1 + 2 + read.getLength()   + read.getName().length() + (label.length() > 0 ? 1 + label.length() : 0); break;
		default: LOG_THROW("Invalid format for Sequence::getIntendedWriteSize(): " << format.getType());
		}
		return length;
	}
	long getIntendedWriteSize(SequenceLengthType sequenceLength, std::string label = "", FormatOutput format = FormatOutput::getDefault()) const {
		return getIntendedWriteSize(*this, sequenceLength, label, format);
	}

	std::string toString() const;

};

#endif

