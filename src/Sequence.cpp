//
// Kmernator/src/Sequence.cpp
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

#include <iostream>
#include <stdexcept>
#include <sstream>

#include "Sequence.h"
#include "Log.h"

#include <cstdlib>
#include <cmath>

using namespace std;

/*----------------------------- SEQUENCE -------------------------------------------*/
char Sequence::FASTQ_START_CHAR = Kmernator::FASTQ_START_CHAR_ILLUMINA;

// thread safe caches
Sequence::CachedSequencesVector Sequence::threadCacheSequences(OMP_MAX_THREADS_DEFAULT, Sequence::CachedSequences(Sequence::maxCachePerThread));

// dangling pointer!!
Sequence::DataPtrListVector *Sequence::preAllocatedDataPtrs = new Sequence::DataPtrListVector();

// static methods of Sequence
void Sequence::clearCaches() {
	threadCacheSequences = CachedSequencesVector(OMP_MAX_THREADS_DEFAULT, CachedSequences(Sequence::maxCachePerThread));
	preAllocatedDataPtrs->reset();
}

Sequence::CachedSequences &Sequence::getCachedSequencesForThread() {
	int threadNum = omp_get_thread_num();
	return threadCacheSequences[threadNum];
}

void Sequence::setThreadCache(Sequence &mmapedSequence, SequencePtr &expandedSequence) {
	mmapedSequence.setCache(expandedSequence);
}

// instance based cache methods
Sequence::SequencePtr Sequence::getCache() const {
	assert(isMmaped());

	CachedSequences &cache = getCachedSequencesForThread();
	SequencePtr cachedSequence = cache.fetch( getRecord() );
	if ( cachedSequence.get() == NULL ) {
		return setCache();
	} else {
		return cachedSequence;
	}
}
Sequence::SequencePtr Sequence::setCache() const {
	SequencePtr ptr;
	return setCache(ptr);
}
Sequence::SequencePtr &Sequence::setCache(Sequence::SequencePtr &expandedSequence) const {
	assert(isMmaped());

	if (expandedSequence.get() == NULL)
		expandedSequence = readMmaped(true);
	CachedSequences &cache = getCachedSequencesForThread();
	cache.insert( getRecord(), expandedSequence);
	return expandedSequence;
}


// Sequence::DataPtrListVector methods

Sequence::DataPtr Sequence::DataPtrListVector::retrieveDataPtr() {
	DataPtrList &list = getList();
	DataPtr ptr;
	if (!list.empty()) {
		ptr = list.front();
		list.pop_back();
	} else {
		ptr = DataPtr( TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(size) );
	}
	assert(ptr.get()->count() == 1);
	return ptr;
}
void Sequence::DataPtrListVector::returnDataPtr(DataPtr &dataPtr) {
	// claim ownership
	DataPtr newPtr;
	newPtr.swap(dataPtr);

	TwoBitSequenceBase::_TwoBitEncodingPtr *ptr = newPtr.get();
	if (ptr == NULL || ptr->count() != 1 || !isValid())
		return;

	DataPtrList &list = getList();
	list.push_back( newPtr );

	assert(newPtr.get()->count() == 2);
}

void Sequence::DataPtrListVector::reset() {
	_isValid = false;
	_vec = _DataPtrListVector(omp_get_max_threads(), DataPtrList());
	_isValid = true;
}


// constructors and operators

Sequence::Sequence() : _flags(0) {
	reset(0);
}

Sequence::Sequence(const Sequence &copy)  {
	*this = copy;
}
Sequence::Sequence(std::string fasta, bool usePreAllocation) :
			_flags(0) {
	setSequence(fasta, usePreAllocation);
}

Sequence::Sequence(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart) :
			_flags(0) {
	setSequence(mmapRecordStart, mmapQualRecordStart);
}

Sequence &Sequence::operator=(const Sequence &other) {
	if (this == &other)
		return *this;
	_flags = other._flags;
	_data = other._data;
	return *this;
}

Sequence Sequence::clone(bool usePreAllocation) const {
	return Sequence(getFasta(), usePreAllocation);
}

Sequence::~Sequence() {
	reset(0);
}

long Sequence::getStoreSize() const {
	if (isMmaped()) {
		// need to build it first!
		return clone(true).getStoreSize();
	}
	const char *start = (char*) this->_getData();
	const char *end = (const char*) this->_getEnd();
	return sizeof(char) + (end-start);
}

long Sequence::store(void *_dst) const {
	if (isMmaped()) {
		return clone(true).store(_dst);
	}
	char *dst = (char*) _dst;
	*(dst++) = _flags;
	const char *start = (const char*) this->_getData();
	const char *end = (const char*) this->_getEnd();
	memcpy(dst, start, (end-start));
	return sizeof(char) + (end-start);
}

void *Sequence::restore(void *_src, long size) {
	char *src = (char*) _src;
	reset(*(src++));
	size -= sizeof(char);
	try {
		_data = DataPtr( TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(size) );
	} catch (...) {
		LOG_THROW(
				"RuntimeError: Cannot allocate memory in Sequence::restore()");
	}
	memcpy(_data.get(), src, size);
	return src + size;
}

void Sequence::setSequence(std::string fasta, bool usePreAllocation) {
	setSequence(fasta, 0, usePreAllocation);
}

void Sequence::setSequence(std::string fasta, long extraBytes, bool usePreAllocation) {
	reset(0);
	SequenceLengthType length = fasta.length();
	unsigned long buffSize = TwoBitSequence::fastaLengthToTwoBitLength(length);

	bool needMalloc = buffSize > MAX_STACK_SIZE;
	TwoBitEncoding _stackBuffer[ needMalloc ? 0 : buffSize ];
	TwoBitEncoding *buffer = _stackBuffer;
	if (needMalloc) {
		buffer = new TwoBitEncoding[ buffSize ];
	}

	if (buffer == NULL)
		LOG_THROW(
				"RuntimeError: Could not allocate buffer memory in Sequence::setSequence");

	BaseLocationVectorType markupBases = TwoBitSequence::compressSequence(fasta, buffer);
	long totalMarkupSize = 0;
	MarkupElementSizeType markupSizes = TwoBitSequence::getMarkupElementSize(markupBases, totalMarkupSize);
	if (totalMarkupSize == 0 && markupBases.size() != 0) {
		// markups exist, but totalMarkupSize will not be stored.
		// signal for discarded record
		setFlag(DISCARDED);
	}
	try {
		unsigned int size =
				sizeof(SequenceLengthType)
				+ buffSize
				+ totalMarkupSize
				+ extraBytes;
		if (usePreAllocation && size <= DataPtrListVector::size) {
			_data = preAllocatedDataPtrs->retrieveDataPtr();
			setFlag(PREALLOCATED);
		} else {
			_data = DataPtr( TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(size) );
		}
	} catch (...) {
		LOG_THROW(
				"RuntimeError: Cannot allocate memory in Sequence::setSequence()");
	}

	*_getLength() = length;
	memcpy(_getTwoBitSequence(), buffer, buffSize);

	// free the buffer
	if (needMalloc)
		delete [] buffer;

	if (totalMarkupSize > 0)
		setMarkups(markupSizes, markupBases);

}

void Sequence::setSequence(Sequence::RecordPtr mmapRecordStart, Sequence::RecordPtr mmapQualRecordStart) {
	setSequence(mmapRecordStart, 0, mmapQualRecordStart);
}
void Sequence::setSequence(Sequence::RecordPtr mmapRecordStart, const BaseLocationVectorType &markups, long extraBytes, Sequence::RecordPtr mmapQualRecordStart) {
	long markupBytes = 0;
	MarkupElementSizeType markupElementSize = TwoBitSequence::getMarkupElementSize(markups, markupBytes);
	extraBytes += markupBytes;
	setSequence(mmapRecordStart, extraBytes, mmapQualRecordStart);
	if (markupBytes > 0)
		setMarkups(markupElementSize, markups);
}
void Sequence::setSequence(Sequence::RecordPtr mmapRecordStart, long extraBytes, Sequence::RecordPtr mmapQualRecordStart) {
	reset(MMAPED);
	long size = sizeof(Sequence::RecordPtr)+extraBytes;
	if (mmapQualRecordStart != NULL)
		setFlag(HASFASTAQUAL);
	setFlag(HASQUALS);
	size += sizeof(Sequence::RecordPtr);
	try {
		_data = DataPtr( TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(size) );
	} catch (...) {
		LOG_THROW("RuntimeError: Can not allocate memory for Sequence::setSequence(mmap)");
	}

	*_getRecord() = mmapRecordStart;
	if (hasFastaQual()) {
		*_getQualRecord() = mmapQualRecordStart;
	}

}

void Sequence::reset(char flags) {
	if (isPreAllocated()) {
		preAllocatedDataPtrs->returnDataPtr(_data);
	}
	_flags = flags;
	_data.reset();
}

void Sequence::setFlag(char f) {
	if ((_flags & f) != f) {
		_flags |= f;
	}
}
void Sequence::unsetFlag(char f) {
	if ((_flags & f) != 0) {
		_flags &= ~f;
	}
}
void Sequence::setMarkups(MarkupElementSizeType markupElementSize, const BaseLocationVectorType &markups) {
	// assume memory has been allocated already
	SequenceLengthType size = markups.size();
	if (size == 0) {
		unsetFlag(MARKUPS4);
		return;
	}
	// if the first base is masked, the entire sequence is masked
	// so discard the read
	if (markupElementSize.first == 0 || (size > 0 && markups[0].first == 'X' && markups[0].second == 0) ) {
		// markups exist, yet no markups will be stored.
		// signal for discarded record
		size = 0;
		unsetFlag(MARKUPS4);
		setFlag(DISCARDED);
		return;
	}
	if (markupElementSize.second == sizeof(BaseLocationType1)) {
		unsetFlag(MARKUPS4); setFlag(MARKUPS1);
		SequenceLengthType1 markupBasesSize = (SequenceLengthType1) size;
		memcpy(_getMarkupBasesCount1(), &markupBasesSize, sizeof(markupBasesSize));
		BaseLocationType1 *ptr = _getMarkupBases1();
		for (SequenceLengthType1 i = 0; i < markupBasesSize; i++) {
			BaseLocationType1 tmp(markups[i].first, markups[i].second);
			memcpy(ptr++, &tmp, sizeof(tmp));
		}
	} else if (markupElementSize.second == sizeof(BaseLocationType2)) {
		unsetFlag(MARKUPS4); setFlag(MARKUPS2);
		SequenceLengthType2 markupBasesSize = (SequenceLengthType2) size;
		memcpy(_getMarkupBasesCount2(), &markupBasesSize, sizeof(markupBasesSize));
		BaseLocationType2 *ptr = _getMarkupBases2();
		for (SequenceLengthType2 i = 0; i < markupBasesSize; i++) {
			BaseLocationType2 tmp(markups[i].first, markups[i].second);
			memcpy(ptr++, &tmp, sizeof(tmp));
		}
	} else {
		setFlag(MARKUPS4);
		SequenceLengthType markupBasesSize = (SequenceLengthType) size;
		memcpy(_getMarkupBasesCount(), &markupBasesSize, sizeof(markupBasesSize));
		BaseLocationType *ptr = _getMarkupBases();
		for (SequenceLengthType i = 0; i < markupBasesSize; i++) {
			BaseLocationType tmp(markups[i].first, markups[i].second);
			memcpy(ptr++, &tmp, sizeof(tmp));
		}
	}
}

string Sequence::getFastaNoMarkup(SequenceLengthType trimOffset, SequenceLengthType trimLength) const {
	if ( !isValid() )
		return string("");

	if (isMmaped()) {
		string name, bases, quals, comment;
		readMmaped(name, bases, quals, comment);
		SequenceLengthType len = bases.length();
		assert(trimOffset <= len);
		if (trimLength < len - trimOffset)
			return bases.substr(trimOffset, trimLength);
		else
			return bases.substr(trimOffset);
	} else {
		SequenceLengthType len = getLength();
		assert(trimOffset <= len);
		if (trimLength > len - trimOffset)
			trimLength = len - trimOffset;
		return TwoBitSequence::getFasta(getTwoBitSequence(), trimOffset, trimLength);
	}
}
string Sequence::getFasta(SequenceLengthType trimOffset, SequenceLengthType trimLength) const {
	if ( !isValid() )
		return string("");
	if (isDiscarded() || trimLength <= 1) {
		// to support printing paired reads where 1 read is trimmed to 0
		return string(1, 'N');
	}
	string fasta;
	if (isMmaped()) {
		SequencePtr sequencePtr = getCache();
		fasta = sequencePtr->getFasta(trimOffset, trimLength);
	} else {
		SequenceLengthType len = getLength();
		assert(trimOffset <= len);
		if (trimLength > len - trimOffset)
			trimLength = len - trimOffset;
		if (trimLength <= 1) {
			// to support printing paired reads where 1 read is trimmed to 0
			return string(1, 'N');
		}
		fasta = TwoBitSequence::getFasta(getTwoBitSequence(), trimOffset, trimLength);
	}
	BaseLocationVectorType markups = getMarkups();
	TwoBitSequence::applyMarkup(fasta, markups);

	return fasta;
}

const Sequence::RecordPtr Sequence::getRecord() const {
	assert(isMmaped());
	return *_getRecord();
}
const Sequence::RecordPtr Sequence::getQualRecord() const {
	assert(isMmaped() && hasFastaQual());
	return *_getQualRecord();
}

const void *Sequence::_getData() const {
	return _data.get();
}
void *Sequence::_getData() {
	return const_cast<void*>(constThis()._getData());
}

const Sequence::RecordPtr *Sequence::_getRecord() const {
	assert(isValid() && isMmaped());
	return (RecordPtr *) _getData();
}
Sequence::RecordPtr *Sequence::_getRecord() {
	return const_cast<RecordPtr *> (constThis()._getRecord());
}
const Sequence::RecordPtr *Sequence::_getQualRecord() const {
	assert(isMmaped() && hasFastaQual());
	return (_getRecord() + 1);
}
Sequence::RecordPtr *Sequence::_getQualRecord() {
	return const_cast<RecordPtr *> (constThis()._getQualRecord());
}

const SequenceLengthType *Sequence::_getLength() const {
	if (isMmaped()) {
		SequencePtr sequencePtr = getCache();
		return sequencePtr->_getLength();
	} else {
		return (SequenceLengthType *) _getData();
	}
}

SequenceLengthType *Sequence::_getLength() {
	return const_cast<SequenceLengthType*> (constThis()._getLength());
}

const TwoBitEncoding *Sequence::_getTwoBitSequence() const {
	if (isMmaped()) {
		SequencePtr sequencePtr = getCache();
		return sequencePtr->_getTwoBitSequence();
	} else {
		return (TwoBitEncoding *) (_getLength() + 1);
	}
}
TwoBitEncoding *Sequence::_getTwoBitSequence() {
	return const_cast<TwoBitEncoding*> (constThis()._getTwoBitSequence());
}

const SequenceLengthType *Sequence::_getMarkupBasesCount() const {
	if (isMmaped()) {
		if (hasFastaQual()) {
			return (SequenceLengthType *) (_getQualRecord() + 1);
		} else {
			return (SequenceLengthType *) (_getRecord() + 1);
		}
	} else {
		return (SequenceLengthType *) (_getTwoBitSequence()
				+ getTwoBitEncodingSequenceLength());
	}

}
SequenceLengthType *Sequence::_getMarkupBasesCount() {
	return const_cast<SequenceLengthType*> (constThis()._getMarkupBasesCount());
}
const SequenceLengthType1 *Sequence::_getMarkupBasesCount1() const {
	return (SequenceLengthType1 *) (_getMarkupBasesCount());
}
SequenceLengthType1  *Sequence::_getMarkupBasesCount1() {
	return const_cast<SequenceLengthType1 *> (constThis()._getMarkupBasesCount1());
}
const SequenceLengthType2 *Sequence::_getMarkupBasesCount2() const {
	return (SequenceLengthType2 *) (_getMarkupBasesCount());
}
SequenceLengthType2 *Sequence::_getMarkupBasesCount2() {
	return const_cast<SequenceLengthType2 *> (constThis()._getMarkupBasesCount2());
}

const BaseLocationType *Sequence::_getMarkupBases() const {
	return (BaseLocationType *) (_getMarkupBasesCount() + 1);
}
BaseLocationType *Sequence::_getMarkupBases() {
	return const_cast<BaseLocationType *> (constThis()._getMarkupBases());
}

const BaseLocationType1 *Sequence::_getMarkupBases1() const {
	return (BaseLocationType1 *) (_getMarkupBasesCount1() + 1);
}
BaseLocationType1 *Sequence::_getMarkupBases1() {
	return const_cast<BaseLocationType1 *> (constThis()._getMarkupBases1());
}

const BaseLocationType2 *Sequence::_getMarkupBases2() const {
	return (BaseLocationType2 *) (_getMarkupBasesCount2() + 1);
}
BaseLocationType2 *Sequence::_getMarkupBases2() {
	return const_cast<BaseLocationType2 *> (constThis()._getMarkupBases2());
}

const void *Sequence::_getEnd() const {
	const void *end;
	if (hasMarkups()) {
		if (isMarkups4())
			end =  _getMarkupBases() + *_getMarkupBasesCount();
		else if (isMarkups2())
			end =  _getMarkupBases2() + *_getMarkupBasesCount2();
		else
			end =  _getMarkupBases1() + *_getMarkupBasesCount1();
	} else {
		end = _getMarkupBasesCount();
	}
	return end;
}

SequenceLengthType Sequence::getLength() const {
	if ( !isValid() )
		return 0;
	if (isMmaped()) {
		assert(isValid());
		return getFastaNoMarkup().length();
	} else {
		return *_getLength();
	}
}
SequenceLengthType Sequence::getFirstMarkupLength() const {
	BaseLocationVectorType markups = getMarkups();
	SequenceLengthType markupLen = TwoBitSequence::firstMarkup(markups);
	if (markupLen == 0)
		return getLength();
	else
		return markupLen - 1;
}
SequenceLengthType Sequence::getFirstMarkupNLength() const {
	BaseLocationVectorType markups = getMarkups();
	SequenceLengthType markupLen = TwoBitSequence::firstMarkupN(markups);
	if (markupLen == 0)
		return getLength();
	else
		return markupLen - 1;
}

SequenceLengthType Sequence::getFirstMarkupXLength() const {
	BaseLocationVectorType markups = getMarkups();
	SequenceLengthType markupLen = TwoBitSequence::firstMarkupX(markups);
	if (markupLen == 0)
		return getLength();
	else
		return markupLen - 1;
}

SequenceLengthType Sequence::getFirstMarkupNorXLength() const {
	BaseLocationVectorType markups = getMarkups();
	SequenceLengthType markupLen = TwoBitSequence::firstMarkupNorX(markups);
	if (markupLen == 0)
		return getLength();
	else
		return markupLen - 1;
}


SequenceLengthType Sequence::getTwoBitEncodingSequenceLength() const {
	// does not need to be valid, yet, but at least _getLength() needs to be initialized
	if ( !isValid() )
		return 0;

	return TwoBitSequence::fastaLengthToTwoBitLength( getLength() );;
}

SequenceLengthType Sequence::_getStoredMarkupBasesCount() const {
	assert(isValid());
	if (isMarkups4()) {
		return *_getMarkupBasesCount();
	} else if (isMarkups1()) {
		return *_getMarkupBasesCount1();
	} else if(isMarkups2()) {
		return *_getMarkupBasesCount2();
	} else {
		return 0;
	}
}
SequenceLengthType Sequence::getMarkupBasesCount() const {
	assert(isValid());
	if (isMmaped())
		return getMarkups().size();
	else
		return _getStoredMarkupBasesCount();
}
BaseLocationVectorType Sequence::_getMarkups() const {
	assert(isValid());
	BaseLocationVectorType markups;
	SequenceLengthType size = _getStoredMarkupBasesCount();
	if (size > 0) {
		markups.reserve(size);

		if (isMarkups4()) {
			const BaseLocationType *ptr = _getMarkupBases();
			for (SequenceLengthType i = 0; i < size; i++) {
				markups.push_back(*(ptr++));
			}
		} else if (isMarkups1()) {
			const BaseLocationType1 *ptr = _getMarkupBases1();
			for (SequenceLengthType i = 0; i < size; i++) {
				markups.push_back(BaseLocationType( ptr->first, ptr->second) );
				ptr++;
			}
		} else {
			const BaseLocationType2 *ptr = _getMarkupBases2();
			for (SequenceLengthType i = 0; i < size; i++) {
				markups.push_back(BaseLocationType( ptr->first, ptr->second));
				ptr++;
			}
		}
	}
	// if the sequence is discarded
	// or the first base is masked, then the whole sequence is discarded and masked
	if (isDiscarded() || (size>0 && markups[0].first == 'X' && markups[0].second == 0)) {
		markups.clear();
		markups.reserve(getLength());
		for(SequenceLengthType i = 0 ; i < getLength(); i++)
			markups.push_back(BaseLocationType('X', i));
	}
	return markups;
}
BaseLocationVectorType Sequence::getMarkups() const {
	BaseLocationVectorType markups = _getMarkups();
	if (isMmaped()) {
		// we must find hidden markups (Ns) in the mmaped record
		// TODO optimize this
		std::string fasta = getFastaNoMarkup();
		TwoBitSequence::applyMarkup(fasta, markups);
		markups = TwoBitSequence::compressSequence(fasta, NULL);
	}
	return markups;
}

void Sequence::readMmaped(std::string &name, std::string &bases, std::string &quals, std::string &comment) const {
	assert(isMmaped());
	// TODO fix hack on NULL lastPtr.  Presently only works for single-lined fastas
	RecordPtr record(getRecord()), lastRecord(NULL), qualRecord(NULL), lastQualRecord(NULL);
	if (hasFastaQual()) {
		qualRecord = getQualRecord();
		lastQualRecord = NULL;
	}
	bool isGood = SequenceRecordParser::parse(record, lastRecord, name, bases, quals, comment, qualRecord, lastQualRecord, FASTQ_START_CHAR);
	if (!isGood)
		this->discard();
}
Sequence::SequencePtr Sequence::readMmaped(bool usePreAllocation) const {
	std::string name, bases, quals, comment;
	readMmaped(name, bases, quals, comment);
	return SequencePtr(new Sequence(bases, usePreAllocation));
}

/*------------------------------------ READ ----------------------------------------*/

double Read::qualityToProbability[256];

int Read::initializeQualityToProbability(unsigned char minQualityScore, char startChar) {
#pragma omp critical (FastqStartChar)
	{
		if (startChar != FASTQ_START_CHAR) {
			LOG_DEBUG_OPTIONAL(1, true, "Switching quality scale for FASTQ (std vs Illumina) to " << (int) startChar);
		}
		FASTQ_START_CHAR = startChar;
		for (int i = 0; i < 256; i++) {
			qualityToProbability[i] = 0.0;
		}
		int start = FASTQ_START_CHAR;
		for (int i = start + minQualityScore; i < PRINT_REF_QUAL; i++)
			qualityToProbability[i] = 1.0 - pow(10.0, ((start - i) / 10.0));
		for (int i = PRINT_REF_QUAL ; i < 256; i++)
			qualityToProbability[i] = 1.0; // for reads with no quality data
	}
	return 1;
}
int Read::qualityToProbabilityInitialized =
		Read::initializeQualityToProbability(0, Kmernator::FASTQ_START_CHAR_ILLUMINA);

void Read::setMinQualityScore(unsigned char minQualityScore, char startChar) {
	Read::initializeQualityToProbability(minQualityScore, startChar);
}

// Read constructors and operators

Read::Read(std::string name, std::string fasta, std::string qualBytes, std::string comment, bool usePreAllocation) {
	setRead(name, fasta, qualBytes, comment, usePreAllocation);
}
Read::Read(Sequence::RecordPtr mmapRecordStart, Sequence::RecordPtr mmapQualRecordStart) {
	setRead(mmapRecordStart, mmapQualRecordStart);
}

Read::Read(Sequence::RecordPtr mmapRecordStart, std::string markupFasta, Sequence::RecordPtr mmapQualRecordStart) {
	setRead(mmapRecordStart, markupFasta, mmapQualRecordStart);
}
Read &Read::operator=(const Read &other) {
	Sequence::operator=(other);
	// there are no extra data members
	return *this;
}
Read Read::clone(bool usePreAllocation) const {
	return Read(getName(), getFasta(), getQuals(), getComment(), usePreAllocation);
}


Read::ReadPtr Read::readMmaped(bool usePreAllocation) const {
	BaseLocationVectorType markups = _getMarkups();
	std::string name, bases, quals, comment;
	Sequence::readMmaped(name, bases, quals, comment);
	TwoBitSequence::applyMarkup(bases, markups);
	return ReadPtr(new Read(name, bases, quals, comment, usePreAllocation));
}

bool Read::recordHasQuals() const {
	assert(isMmaped());
	if (hasQuals() || hasFastaQual())
		return true;
	else
		// TODO make this more general
		return *getRecord() == '@'; // FASTQ
}

ProbabilityBases Read::getProbabilityBases(unsigned char minQuality) const {
	std::string fasta = getFasta();
	std::string quals = getQuals();
	ProbabilityBases probs(fasta.length());
	for(int i = 0; i < (int) fasta.length(); i++) {
		unsigned char q = quals[i];
		if (q < minQuality + Read::FASTQ_START_CHAR)
			break;
		double prob = qualityToProbability[ q ];
		if (prob < 0.2501) {
			prob = 0.2501; // slightly better than random...
		}
		ProbabilityBase &base = probs[i];
		base.observe(fasta[i], prob);
	}

	return probs;
}
double Read::scoreProbabilityBases(const ProbabilityBases &probs) const {
	ProbabilityBases myprobs = getProbabilityBases();
	myprobs *= probs;
	return myprobs.sum();
}

const char * Read::_getQual() const {
	assert(!isMmaped());
	if (isMarkups1()) {
		return (char *) (_getMarkupBases1() + *_getMarkupBasesCount1());
	} else if (isMarkups2()) {
		return (char *) (_getMarkupBases2() + *_getMarkupBasesCount2());
	} else if (isMarkups4()) {
		return (char *) (_getMarkupBases() + *_getMarkupBasesCount());
	} else {
		// there are no markups, so there is no markup count...
		return (char *) (_getMarkupBasesCount());
	}
}
char * Read::_getQual() {
	return const_cast<char*> (constThis()._getQual());
}

SequenceLengthType Read::_qualLength() const {
	assert(!isMmaped());
	if (hasQuals()) {
		SequenceLengthType len = getLength();
		if (len == 0)
			return 0;
		else if (*(_getQual()) == REF_QUAL)
			return 1;
		else
			return len;
	} else {
		return 0;
	}
}

const char * Read::_getName() const {
	assert(!isMmaped());
	return (char *) (_getQual() + _qualLength());
}
char * Read::_getName() {
	return const_cast<char*> (constThis()._getName());
}

const char * Read::_getComment() const {
	const char *name = _getName();
	int len = strlen(name);
	if (GlobalOptions::isCommentStored()) {
		return name + len+1;
	} else {
		return name + len; // null termination
	}
}
char * Read::_getComment() {
	return const_cast<char *> (constThis()._getComment());
}

const void *Read::_getEnd() const {
	assert(!isMmaped());
	if (GlobalOptions::isCommentStored()) {
		return (const void *) (_getComment() + strlen(_getComment())+1);
	} else {
		return (const void *) (_getName() + strlen(_getName())+1);
	}
}

void Read::setRead(std::string name, std::string fasta, std::string qualBytes, std::string comment, bool usePreAllocation) {
	assert(!isMmaped());
	if (fasta.length() != qualBytes.length())
		LOG_THROW(
				"InvalidFormat for setRead(): fasta length != qual length for name = " << name << " " << fasta << " " << qualBytes);

	// do not store quals if it is a reference
	if (qualBytes.length() > 1 && qualBytes[0] == REF_QUAL) {
		qualBytes.clear();
	}

	int extraLength = name.length() + 1;
	if (GlobalOptions::isCommentStored())
		extraLength += comment.length() + 1;
	Sequence::setSequence(fasta, qualBytes.length() + extraLength, usePreAllocation);

	if (qualBytes.empty()) {
		unsetFlag(HASQUALS);
	} else {
		setFlag(HASQUALS);
		memcpy(_getQual(), qualBytes.c_str(), qualBytes.length());
	}
	strcpy(_getName(), name.c_str());

	if (GlobalOptions::isCommentStored())
		strcpy(_getComment(), comment.c_str());

}
void Read::setRead(Sequence::RecordPtr mmapRecordStart, Sequence::RecordPtr mmapQualRecordStart) {
	setSequence(mmapRecordStart, 0, mmapQualRecordStart);
	if (recordHasQuals())
		setFlag(HASQUALS);
}
void Read::setRead(Sequence::RecordPtr mmapRecordStart, std::string markupFasta, Sequence::RecordPtr mmapQualRecordStart) {
	BaseLocationVectorType markups = TwoBitSequence::compressSequence(markupFasta.c_str(), NULL);
	setSequence(mmapRecordStart, markups, 0, mmapQualRecordStart);
	if (recordHasQuals())
		setFlag(HASQUALS);
}

void Read::markupBases(SequenceLengthType offset, SequenceLengthType length, char mask) {
	if (isDiscarded())
		return;

	string fasta = getFasta();
	string origFasta = fasta;
	SequenceLengthType len = getLength();
	if (offset >= len) {
		LOG_THROW("Read::markupBases(): offset can not be greater than length of sequence");
	}
	if (offset + length > len)
		length = len - offset;

	// save some memory if the first base is being masked, the whole sequence is discarded
	// so only need to mark the first base
	if ((fasta.length()>0 && fasta[0] == 'X')
			|| (offset == 0 && mask == 'X')) {
		setFlag(DISCARDED);
		return;
	}
	fasta.replace(offset, length, length, mask);


	if (isMmaped()) {
		RecordPtr record = getRecord();
		RecordPtr qualRecord = NULL;
		if (hasFastaQual())
			qualRecord = getQualRecord();

		reset();
		setRead(record, fasta, qualRecord);

	} else {
		string name  = getName();
		string qual  = getQuals();
		string comment = getComment();

		reset();
		setRead(name, fasta, qual, comment);
	}
}

string Read::getName() const {
	if ( !isValid() ) {
		return string("");
	} else if (isMmaped()) {
		string name, bases, quals, comment;
		Sequence::readMmaped(name, bases, quals, comment);
		return name;
	} else {
		return string(_getName());
	}
}

void Read::setName(const std::string name) { // inefficient!
	setRead(name, getFasta(), getQuals(), getComment());
}

string Read::getComment() const {
	if ( !isValid() ) {
		return string("");
	} else if (isMmaped()) {
		string name, bases, quals, comment;
		Sequence::readMmaped(name, bases, quals, comment);
		return comment;
	} else {
		return string(_getComment());
	}
}

void Read::setComment(const std::string comment) { // inefficient!
	setRead(getName(), getFasta(), getQuals(), comment);
}

string Read::getQuals(SequenceLengthType trimOffset, SequenceLengthType trimLength, bool forPrinting, bool unmasked) const {
	if ((isDiscarded() || trimLength <= 1) && !unmasked) {
		// to support printing paired reads where 1 read is trimmed to 0
		return string(1, FASTQ_START_CHAR+1);
	}
	if (isMmaped()) {
		string name, bases, quals, comment;
		Sequence::readMmaped(name, bases, quals, comment);
		SequenceLengthType len = bases.length();
		assert(trimOffset <= len);
		if (trimLength > len - trimOffset)
			trimLength = len - trimOffset;

		if (trimLength > 0) {
			if ( (!hasQuals()) || quals[0] == REF_QUAL) {
				if (forPrinting)
					return string(trimLength, PRINT_REF_QUAL);
				else
					return string(trimLength, REF_QUAL);
			} else {
				return quals.substr(trimOffset, trimLength);
			}
		} else {
			// to support printing paired reads where 1 read is trimmed to 0
			return string(1, FASTQ_START_CHAR+1);
		}
	} else {
		SequenceLengthType len = getLength();
		const char * qualPtr = NULL;
		assert(trimOffset <= len);
		if (trimLength > len - trimOffset)
			trimLength = len - trimOffset;

		if (trimLength > 0) {
			qualPtr = _getQual() + trimOffset;

			if ( (!hasQuals()) || *qualPtr == REF_QUAL) {
				if (forPrinting)
					return string(trimLength, PRINT_REF_QUAL);
				else
					return string(trimLength, REF_QUAL);
			} else {
				return string(qualPtr, trimLength);
			}
		} else if (_getData() == NULL) {
			// This is an invalid / empty sequence...
			return string("");
		} else {
			// to support printing paired reads where 1 read is trimmed to 0
			return string(1, FASTQ_START_CHAR+1);
		}
	}
}

string Read::toFastq(SequenceLengthType trimOffset, SequenceLengthType trimLength, std::string label, bool unmasked) const {
	std::string fasta;
	if (unmasked)
		fasta = getFastaNoMarkup(trimOffset, trimLength);
	else
		fasta = getFasta(trimOffset, trimLength);
	return string('@' + getNameAndComment() + (label.length() > 0 ? " " + label : "")
			+ "\n" + fasta + "\n+\n"
			+ getQuals(trimOffset, trimLength, true, unmasked) + "\n");
}
string Read::toFasta(SequenceLengthType trimOffset, SequenceLengthType trimLength, std::string label, bool unmasked) const {
	std::string fasta;
	if (unmasked)
		fasta = getFastaNoMarkup(trimOffset, trimLength);
	else
		fasta = getFasta(trimOffset, trimLength);
	return string('>' + getNameAndComment() + (label.length() > 0 ? " " + label : "")
			+ "\n" + fasta + "\n");
}
string Read::toQual(SequenceLengthType trimOffset, SequenceLengthType trimLength, std::string label) const {
	return string('>' + getNameAndComment() + (label.length() > 0 ? " " + label : "")
			+ "\n" + getFormattedQuals(trimOffset, trimLength) + "\n");
}

// TODO format with linebreaks
string Read::getFormattedQuals(SequenceLengthType trimOffset, SequenceLengthType trimLength) const {
	string quals = getQuals(trimOffset, trimLength, true);
	stringstream ss;
	SequenceLengthType len = quals.length();
	assert(trimOffset < len);
	if (trimLength > len - trimOffset)
		trimLength = len - trimOffset;

	for (SequenceLengthType i = trimOffset; i < trimOffset+trimLength; i++) {
		ss << (int) (quals[i] - FASTQ_START_CHAR) << ' ';
	}
	return ss.str();
}

std::string Read::toString() const {
	return getNameAndComment() + "\t" + getFasta() + "\t" + getQuals(0, MAX_SEQUENCE_LENGTH, true, true) + "\t" + (isDiscarded()? " discarded" : "");
}


/*------------------------------------ BaseQual ----------------------------------------*/

char probToQual(double prob) {
	return (char) (-10. * std::log10(1.0-prob));
}

char BaseQual::getQualChar(double prob, bool ignoreLow) {
	if (prob >= 0.9999) {
		return Sequence::FASTQ_START_CHAR + 40;
	} else if (ignoreLow && prob <= 0.25) {
		return ' ';
	} else {
		return Sequence::FASTQ_START_CHAR + probToQual(prob);
	}
}

/*------------------------------------ ProbabilityBase ----------------------------------------*/

ProbabilityBase &ProbabilityBase::operator=(const ProbabilityBase &copy) {
	a = copy.a;
	c = copy.c;
	g = copy.g;
	t = copy.t;
	top = copy.top;
	best = copy.best;
	count = copy.count;
	return *this;
}
ProbabilityBase ProbabilityBase::operator+(const ProbabilityBase &other) {
	ProbabilityBase tmp(*this);
	tmp.a += other.a;
	tmp.c += other.c;
	tmp.g += other.g;
	tmp.t += other.t;
	tmp.count += other.count;
	tmp.setTop(other);
	return tmp;
}
ProbabilityBase &ProbabilityBase::operator+=(const ProbabilityBase &other) {
	*this = *this + other;
	// setTop already called
	return *this;
}
ProbabilityBase ProbabilityBase::operator*(const ProbabilityBase &other) {
	ProbabilityBase tmp(*this);
	tmp.a *= other.a;
	tmp.c *= other.c;
	tmp.g *= other.g;
	tmp.t *= other.t;
	tmp.setTop(other);
	return tmp;
}
ProbabilityBase &ProbabilityBase::operator*=(const ProbabilityBase &other) {
	*this = *this * other;
	// setTop already called
	return *this;
}
ProbabilityBase &ProbabilityBase::operator*=(double factor) {
	a *= factor;
	c *= factor;
	g *= factor;
	t *= factor;
	return *this;
}

void ProbabilityBase::observe(char nuc, double prob) {
	double otherProb = (1.0 - prob) / 3.0;
	ProbabilityBase &base = *this;
	switch (nuc) {
	case 'A':
	case 'a': base.a += prob; base.c += otherProb; base.g += otherProb; base.t += otherProb; break;
	case 'C':
	case 'c': base.c += prob; base.a += otherProb; base.g += otherProb; base.t += otherProb; break;
	case 'G':
	case 'g': base.g += prob; base.a += otherProb; base.c += otherProb; base.t += otherProb; break;
	case 'T':
	case 't': base.t += prob; base.a += otherProb; base.c += otherProb; base.g += otherProb; break;
	}
	base.count++;
}

void ProbabilityBase::setTop(double _a, double _c, double _g, double _t) {
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

double ProbabilityBase::getA() const {
	if (best == 'A' && top < a * count)
		return top;
	else
		return a;
}
double ProbabilityBase::getC() const {
	if (best == 'C' && top < c * count)
		return top;
	else
		return c;
}
double ProbabilityBase::getG() const {
	if (best == 'G' && top < g * count)
		return top;
	else
		return g;
}
double ProbabilityBase::getT() const {
	if (best == 'T' && top < t * count)
		return top;
	else
		return t;
}

BaseQual ProbabilityBase::getBaseQual() const {
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
	LOG_THROW("Should not get here: ProbabilityBase::getBaseQual()");
}


/*------------------------------------ ProbabilityBases ----------------------------------------*/

double ProbabilityBases::sum() const {
	double sum = 0.0;
	for(int i = 0; i < (int) _bases.size(); i++)
		sum += _bases[i].sum();
	return sum;
}

ProbabilityBase &ProbabilityBases::operator[](size_t idx) {
	return _bases[idx];
}
const ProbabilityBase &ProbabilityBases::operator[](size_t idx) const {
	return _bases[idx];
}

ProbabilityBases &ProbabilityBases::operator+=(const ProbabilityBases &other) {
	if (_bases.size() < other._bases.size()) {
		resize(other._bases.size());
	}

	for(size_t i = 0 ; i < other._bases.size(); i++) {
		_bases[i] += other._bases[i];
	}
	return *this;
}
ProbabilityBases &ProbabilityBases::operator*=(const ProbabilityBases &other) {
	if (_bases.size() < other._bases.size() ) {
		resize(other._bases.size());
	}
	for(size_t i = 0; i < other._bases.size(); i++) {
		_bases[i] *= other._bases[i];
	}
	return *this;
}
ProbabilityBases &ProbabilityBases::operator*=(double factor) {
	for(size_t i = 0; i <  _bases.size(); i++) {
		_bases[i] *= factor;
	}
	return *this;
}
std::string ProbabilityBases::toString() const {
	std::stringstream ss;
	for(size_t i = 0; i < _bases.size(); i++) {
		ss << _bases[i].getBaseQual().base;
	}
	ss << std::endl;
	for(size_t i = 0; i < _bases.size(); i++) {
		ss << _bases[i].best;
	}
	ss << std::endl;
	for(size_t i = 0; i < _bases.size(); i++) {
		ss << _bases[i].getBaseQual().qual;
	}
	ss << std::endl;
	for(size_t i = 0; i < _bases.size(); i++) {
		ss << BaseQual::getQualChar( _bases[i].getA(), true );
	}
	ss << std::endl;
	for(size_t i = 0; i < _bases.size(); i++) {
		ss << BaseQual::getQualChar( _bases[i].getC(), true );
	}
	ss << std::endl;
	for(size_t i = 0; i < _bases.size(); i++) {
		ss << BaseQual::getQualChar( _bases[i].getG(), true );
	}
	ss << std::endl;
	for(size_t i = 0; i < _bases.size(); i++) {
		ss << BaseQual::getQualChar( _bases[i].getT(), true );
	}
	return ss.str();
}
