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
		static const long size = 511;
		DataPtrListVector() : _isValid(false) {
			_vec = _DataPtrListVector(OMP_MAX_THREADS, DataPtrList());
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
	Sequence clone() const;

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
	inline bool hasFastaQual()   const { return (_flags & HASFASTAQUAL)  == HASFASTAQUAL; }
	inline bool isValid()        const { return ( _getData() != NULL ); }

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

	std::string getFasta(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH) const;
	std::string getFastaNoMarkup(SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH) const;

	SequenceLengthType getMarkupBasesCount() const;
	BaseLocationVectorType getMarkups() const;

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
	Read &operator=(const Read &other);
	Read clone() const;

	void setRead(std::string name, std::string fasta, std::string qualBytes, bool usePreAllocation = false);
	void setRead(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart = NULL);
	void setRead(RecordPtr mmapRecordStart, std::string markupFasta, RecordPtr mmapQualRecordStart);

	ReadPtr readMmaped(bool usePreAllocation = false) const;
	bool recordHasQuals() const ;

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

	ProbabilityBases getProbabilityBases() const;
	double scoreProbabilityBases(const ProbabilityBases &probs) const;

	static double qualityToProbability[256];
	static void setMinQualityScore(unsigned char minQualityScore, char startChar = Kmernator::FASTQ_START_CHAR_ILLUMINA);

	// format == 0 fastq
	// format == 1 fasta
	// format == 2 fastq unmasked
	// format == 3 fasta unmasked
	inline static std::ostream &write(std::ostream &os, const Read &read,
			SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH, std::string label = "", FormatOutput format = FormatOutput::getDefault()) {
		switch (format.getType()) {
		case FormatOutput::FASTQ: os << read.toFastq(trimOffset, label); break;
		case FormatOutput::FASTA: os << read.toFasta(trimOffset, label); break;
		case FormatOutput::FASTQ_UNMASKED: os << read.toFastq(trimOffset, label, true) ; break;
		case FormatOutput::FASTA_UNMASKED: os << read.toFasta(trimOffset, label, true); break;
		default: throw std::invalid_argument("Invalid format");
		}
		return os;
	}
	inline std::ostream &write(std::ostream &os,
			SequenceLengthType trimOffset = MAX_SEQUENCE_LENGTH, std::string label = "", FormatOutput format = FormatOutput::getDefault()) const {
		return write(os, *this, trimOffset, label, format);
	}

	std::string toString() const;

};

#endif

//
// $Log: Sequence.h,v $
// Revision 1.34  2010-08-18 17:50:39  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.33.4.1  2010-07-20 20:02:56  regan
// autodetect fastq quality range
//
// Revision 1.33  2010-06-22 23:06:31  regan
// merged changes in CorruptionBugfix-20100622 branch
//
// Revision 1.32.4.1  2010-06-22 22:58:38  regan
// added has fasta qual flag to differentiate from has quals
//
// Revision 1.32  2010-05-24 21:48:46  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.31.2.1  2010-05-19 00:20:46  regan
// refactored fomat output options
// added options to fastq2fasta
//
// Revision 1.31  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.30.2.3  2010-05-10 22:40:20  regan
// minor refactor -- replaced invalid flag
//
// Revision 1.30.2.2  2010-05-10 20:44:35  regan
// minor refactor moved code into cpp
//
// Revision 1.30.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.30  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.29  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
// Revision 1.28.2.1  2010-05-06 18:45:35  regan
// broke it...
//
// Revision 1.28  2010-05-06 16:43:56  regan
// merged changes from ConsensusTesting-20100505
//
// Revision 1.27.2.1  2010-05-05 23:46:22  regan
// checkpoint... seems to compile
//
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
