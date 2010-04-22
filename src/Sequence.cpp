// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.cpp,v 1.30 2010-04-22 23:41:32 regan Exp $
//

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include "Sequence.h"
#include <cstdlib>
#include <cmath>

using namespace std;

/*----------------------------- SEQUENCE -------------------------------------------*/

Sequence::CachedSequenceVector Sequence::threadCacheSequence(OMP_MAX_THREADS, CachedSequence());
long Sequence::constructCount = 0;
long Sequence::destructCount = 0;

Sequence::Sequence() {
#pragma omp atomic
	constructCount++;
	reset(INVALID);
}

Sequence::Sequence(std::string fasta) :
	_flags(0) {
#pragma omp atomic
	constructCount++;
	setSequence(fasta);
}

Sequence::Sequence(RecordPtr mmapRecordStart, RecordPtr mmapQualRecordStart) :
	_flags(0) {
#pragma omp atomic
	constructCount++;
	setSequence(mmapRecordStart, mmapQualRecordStart);
}

Sequence::~Sequence() {
#pragma omp atomic
	destructCount++;
	reset(INVALID);
}

void Sequence::setSequence(std::string fasta) {
	setSequence(fasta, 0);
}

void Sequence::setSequence(std::string fasta, long extraBytes) {
	reset(INVALID);
	SequenceLengthType length = fasta.length();
	bool freeBuffer = false;
	const unsigned long STACK_SIZE = 1024;
	TwoBitEncoding _stackBuffer[STACK_SIZE];

	unsigned long buffSize = TwoBitSequence::fastaLengthToTwoBitLength(length);
	TwoBitEncoding *buffer;
	if (buffSize > STACK_SIZE) {
	    buffer = (TwoBitEncoding*) malloc(buffSize);
	    freeBuffer = true;
	} else {
	    buffer = _stackBuffer;
	}

	if (buffer == NULL)
		throw std::runtime_error(
				"Could not allocate buffer memory in Sequence::setSequence");

	BaseLocationVectorType markupBases = TwoBitSequence::compressSequence(fasta, buffer);
	long totalMarkupSize = 0;
	MarkupElementSizeType markupSizes = TwoBitSequence::getMarkupElementSize(markupBases, totalMarkupSize);
	try {
		unsigned int size =
			    sizeof(SequenceLengthType)
		        + buffSize
				+ totalMarkupSize
				+ extraBytes;
		_data = DataPtr( TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(size) );
	} catch (...) {
		throw new std::runtime_error(
				"Cannot allocate memory in Sequence::setSequence()");
	}

    *_getLength() = length;
	memcpy(_getTwoBitSequence(), buffer, buffSize);

    // free the buffer
	if (freeBuffer)
	    free(buffer);

	setMarkups(markupSizes, markupBases);
	unsetFlag(INVALID);

}

void Sequence::setSequence(Sequence::RecordPtr mmapRecordStart, Sequence::RecordPtr mmapQualRecordStart) {
	setSequence(mmapRecordStart, 0, mmapQualRecordStart);
}
void Sequence::setSequence(Sequence::RecordPtr mmapRecordStart, const BaseLocationVectorType &markups, long extraBytes, Sequence::RecordPtr mmapQualRecordStart) {
	long markupBytes = 0;
	MarkupElementSizeType markupElementSize = TwoBitSequence::getMarkupElementSize(markups, markupBytes);
	extraBytes += markupBytes;
	setSequence(mmapRecordStart, extraBytes, mmapQualRecordStart);
	setMarkups(markupElementSize, markups);
}
void Sequence::setSequence(Sequence::RecordPtr mmapRecordStart, long extraBytes, Sequence::RecordPtr mmapQualRecordStart) {
	reset(MMAPED|INVALID);
	long size = sizeof(Sequence::RecordPtr)+extraBytes;
	if (mmapQualRecordStart != NULL)
		size += sizeof(Sequence::RecordPtr);
	try {
		_data = DataPtr( TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(size) );
	} catch (...) {
		throw new std::runtime_error("Can not allocate memory for Sequence::setSequence(mmap)");
	}
	unsetFlag(INVALID);
	*_getRecord() = mmapRecordStart;
	if (mmapQualRecordStart != NULL) {
		setFlag(MMAPED_QUALS);
		*_getQualRecord() = mmapQualRecordStart;
	}

}

void Sequence::reset(char flags) {
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
	if (size > 0 && markups[0].first == 'X' && markups[0].second == 0) {
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

string Sequence::getFastaNoMarkup() const {
	if (_data.get() == NULL)
	    return string("");
	if (isMmaped()) {
		string name, bases, quals;
		readMmaped(name, bases, quals);
		return bases;
	} else {
	    return TwoBitSequence::getFasta(getTwoBitSequence(), getLength());
	}
}
string Sequence::getFasta(SequenceLengthType trimOffset) const {
	if (_data.get() == NULL)
		return string("");
	string fasta;
	if (isMmaped()) {
		SequencePtr sequencePtr = getCache();
		fasta = sequencePtr->getFasta(trimOffset);
	} else {
	    SequenceLengthType len = getLength();
	    if (trimOffset < len)
		    len = trimOffset;
	    if (len == 0) {
		    // to support printing paired reads where 1 read is trimmed to 0
		    return string(1, 'N');
	    }
	    fasta = TwoBitSequence::getFasta(getTwoBitSequence(), len);
	}
	BaseLocationVectorType markups = getMarkups();
	TwoBitSequence::applyMarkup(fasta, markups);

	return fasta;
}

const Sequence::RecordPtr *Sequence::_getRecord() const {
	assert(isValid() && isMmaped());
	return (RecordPtr *) _data.get();
}
Sequence::RecordPtr *Sequence::_getRecord() {
	return const_cast<RecordPtr *> (constThis()._getRecord());
}
const Sequence::RecordPtr *Sequence::_getQualRecord() const {
	assert(isQualMmaped());
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
		return (SequenceLengthType *) _data.get();
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
		if (isQualMmaped()) {
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

SequenceLengthType Sequence::getLength() const {
	if (_data.get() == NULL)
		return 0;
	if (isMmaped()) {
		assert(isValid());
		return getFastaNoMarkup().length();
	} else {
		return *_getLength();
	}
}

SequenceLengthType Sequence::getTwoBitEncodingSequenceLength() const {
	// does not need to be valid, yet, but at least _getLength() needs to be initialized
	if (_data.get() == NULL)
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
		    for (unsigned int i = 0; i < size; i++) {
			    markups.push_back(*(ptr++));
		    }
		} else if (isMarkups1()) {
			const BaseLocationType1 *ptr = _getMarkupBases1();
			for (unsigned int i = 0; i < size; i++) {
			    markups.push_back(BaseLocationType( ptr->first, ptr->second) );
			    ptr++;
			}
		} else {
			const BaseLocationType2 *ptr = _getMarkupBases2();
		    for (unsigned int i = 0; i < size; i++) {
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
		for(unsigned int i = 0 ; i < getLength(); i++)
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

/*------------------------------------ READ ----------------------------------------*/

//Read::CachedReadVector Read::threadCacheRead(OMP_MAX_THREADS, CachedRead());
double Read::qualityToProbability[256];

int Read::initializeQualityToProbability(unsigned char minQualityScore) {
	for (int i = 0; i < 256; i++) {
		qualityToProbability[i] = 0;
	}
	int start = FASTQ_START_CHAR;
	for (int i = start + minQualityScore; i <= PRINT_REF_QUAL; i++)
		qualityToProbability[i] = 1.0 - pow(10.0, ((start - i) / 10.0));

	qualityToProbability[255] = 1.0; // for reads with no quality data
	return 1;
}
int Read::qualityToProbabilityInitialized =
		Read::initializeQualityToProbability(0);

void Read::setMinQualityScore(unsigned char minQualityScore) {
	Read::initializeQualityToProbability(minQualityScore);
}

Read::Read(std::string name, std::string fasta, std::string qualBytes) {
	setRead(name, fasta, qualBytes);
}
Read::Read(Sequence::RecordPtr mmapRecordStart, Sequence::RecordPtr mmapQualRecordStart) {
	setRead(mmapRecordStart, mmapQualRecordStart);
}

Read::Read(Sequence::RecordPtr mmapRecordStart, std::string markupFasta, Sequence::RecordPtr mmapQualRecordStart) {
	setRead(mmapRecordStart, markupFasta, mmapQualRecordStart);
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

void Read::setRead(std::string name, std::string fasta, std::string qualBytes) {
	assert(!isMmaped());
	if (fasta.length() != qualBytes.length())
		throw new std::invalid_argument(
				"fasta length != qual length for name = " + name);

	// do not store quals if it is a reference
	if (qualBytes.length() > 1 && qualBytes[0] == REF_QUAL) {
		qualBytes.clear();
	}

	Sequence::setSequence(fasta, qualBytes.length() + (name.length() + 1));

	if (qualBytes.empty()) {
		unsetFlag(HASQUALS);
	} else {
	    setFlag(HASQUALS);
		memcpy(_getQual(), qualBytes.c_str(), qualBytes.length());
	}
	strcpy(_getName(), name.c_str());
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
  string fasta = getFasta();
  SequenceLengthType len = getLength();
  if (offset + length > len)
    length = len - offset;

  // save some memory if the first base is being masked, the whole sequence is discarded
  // so only need to mark the first base
  if ((fasta.length()>0 && fasta[0] == 'X')
		|| (offset == 0 && mask == 'X')) {
		fasta = getFastaNoMarkup();
		offset = 0;
		length = 1;
		mask = 'X';
  }
  fasta.replace(offset, length, length, mask);

  if (isMmaped()) {
	RecordPtr record = getRecord();
	RecordPtr qualRecord = NULL;
	if (isQualMmaped())
		qualRecord = getQualRecord();

	reset();
	setRead(record, fasta, qualRecord);

  } else {
	string name  = getName();
	string qual  = getQuals();

	reset();
	setRead(name, fasta, qual);
  }
}

string Read::getName() const {
	if (_data.get() == NULL) {
		return string("");
	} else if (isMmaped()) {
		string name, bases, quals;
		Sequence::readMmaped(name, bases, quals);
		return name;
	} else {
		return string(_getName());
	}
}

string Read::getQuals(SequenceLengthType trimOffset, bool forPrinting) const {
  if (isMmaped()) {
	string name, bases, quals;
	Sequence::readMmaped(name, bases, quals);
	SequenceLengthType len = bases.length();
	len = (trimOffset <= len ? trimOffset : len);

	if (len > 0) {
		  if ( (!hasQuals()) || quals[0] == REF_QUAL) {
			if (forPrinting)
				return string(len, PRINT_REF_QUAL);
			else
				return string(len, REF_QUAL);
		  } else {
		      return quals.substr(0,len);
		  }
	} else {
		// to support printing paired reads where 1 read is trimmed to 0
		return string(1, FASTQ_START_CHAR+1);
	}
  } else {
	SequenceLengthType len = getLength();
	const char * qualPtr = NULL;
	len = (trimOffset <= len ? trimOffset : len);
	if (len > 0) {
	  qualPtr = _getQual();

	  if ( (!hasQuals()) || *qualPtr == REF_QUAL) {
		if (forPrinting)
			return string(len, PRINT_REF_QUAL);
		else
			return string(len, REF_QUAL);
	  } else {
		return string(qualPtr, len);
	  }
	} else if (_data.get() == NULL) {
		// This is an invalid / empty sequence...
	    return string("");
    } else {
		// to support printing paired reads where 1 read is trimmed to 0
		return string(1, FASTQ_START_CHAR+1);
	}
  }
}

string Read::toFastq(SequenceLengthType trimOffset, std::string label) const {
	return string('@' + getName() + (label.length() > 0 ? " " + label : "")
			+ "\n" + getFasta(trimOffset) + "\n+\n"
			+ getQuals(trimOffset, true) + "\n");
}
string Read::getFormattedQuals(SequenceLengthType trimOffset) const {
	string quals = getQuals(trimOffset, true);
	stringstream ss;
	SequenceLengthType len = quals.length();
	if (trimOffset < len)
		len = trimOffset;
	for (unsigned int i = 0; i < len; i++) {
		ss << (int) quals[i] - 64 << ' ';
	}
	return ss.str();
}
//
// $Log: Sequence.cpp,v $
// Revision 1.30  2010-04-22 23:41:32  regan
// fixed a few bugs
//
// Revision 1.29  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.28.2.9.2.2  2010-04-16 21:38:09  regan
// minor performance change
//
// Revision 1.28.2.9.2.1  2010-04-16 05:30:00  regan
// checkpoint.. broke it
//
// Revision 1.28.2.9  2010-04-15 21:31:50  regan
// bugfix in markups and duplicate fragment filter
//
// Revision 1.28.2.8  2010-04-15 17:29:02  regan
// checkpoint, working with some optimizations
//
// Revision 1.28.2.7  2010-04-14 22:36:06  regan
// round of bugfixes
//
// Revision 1.28.2.6  2010-04-14 20:53:49  regan
// checkpoint and passes unit tests!
//
// Revision 1.28.2.5  2010-04-14 03:51:20  regan
// checkpoint. compiles but segfaults
//
// Revision 1.28.2.4  2010-04-12 22:37:47  regan
// checkpoint
//
// Revision 1.28.2.3  2010-04-12 20:59:45  regan
// mmap checkpoint
//
// Revision 1.28.2.2  2010-04-05 02:56:08  regan
// bugfixes
//
// Revision 1.28.2.1  2010-04-04 15:31:27  regan
// checkpoint - refactored underlying data structure t compress markups
//
// Revision 1.28  2010-03-16 06:42:50  regan
// bugfixes
//
// Revision 1.27  2010-03-15 14:58:42  regan
// fixed major bug in markups
//
// Revision 1.26  2010-03-08 22:14:38  regan
// replaced zero bases with markup bases to mask out reads that match the filter pattern
// bugfix in overrunning the mask
//
// Revision 1.25  2010-03-03 17:38:48  regan
// fixed quality scores
//
// Revision 1.24  2010-03-03 17:10:49  regan
// let zero trimmed reads print out with a single N to support pairs staying together
//
// Revision 1.23  2010-02-26 13:01:16  regan
// reformatted
//
// Revision 1.22  2010-02-22 14:41:03  regan
// bugfix in printing
//
// Revision 1.21  2010-01-14 18:04:14  regan
// bugfixes
//
// Revision 1.20  2010-01-13 23:46:46  regan
// made const class modifications
//
// Revision 1.19  2010-01-13 21:16:00  cfurman
// added setMinQualityScore
//
// Revision 1.18  2010-01-13 00:24:30  regan
// use less memory for reference sequences and those without quality
//
// Revision 1.17  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.16  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.15  2009-11-28 01:00:07  regan
// fixed bugs and warnings
//
// Revision 1.14  2009-11-21 15:58:29  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.13  2009-11-10 07:04:59  regan
// bugfix in quality lookup table -- use float!
//
// Revision 1.12  2009-11-07 00:28:41  cfurman
// ReadSet now takes fasta, fastq or  fasta+qual files.
//
// Revision 1.11  2009-11-06 04:07:13  regan
// bugfix when stack size is limited
//
// Revision 1.10  2009-11-04 19:32:03  cfurman
// now reads in fasta (with optional qual) files
//
// Revision 1.9  2009-11-02 18:27:00  regan
// added getMarkups()
// added quality to probability lookup table
//
// Revision 1.8  2009-10-30 00:51:40  regan
// bug fix and working on executable
//
// Revision 1.7  2009-10-22 20:49:15  cfurman
// tests added
//
// Revision 1.6  2009-10-22 07:04:06  regan
// added a few unit tests
// minor refactor
//
// Revision 1.5  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
// Revision 1.4  2009-10-21 00:00:58  cfurman
// working on kmers....
//
// Revision 1.3  2009-10-20 20:56:27  cfurman
// Got it to compile!
//
// Revision 1.2  2009-10-20 17:25:50  regan
// added CVS tags
//
//


