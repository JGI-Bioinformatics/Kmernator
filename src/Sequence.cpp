// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.cpp,v 1.34 2010-05-18 20:50:24 regan Exp $
//

#include <iostream>
#include <stdexcept>
#include <sstream>

#include "Sequence.h"
#include <cstdlib>
#include <cmath>

using namespace std;

/*----------------------------- SEQUENCE -------------------------------------------*/

// thread safe caches
Sequence::CachedSequencesVector Sequence::threadCacheSequences(OMP_MAX_THREADS, Sequence::CachedSequences(Sequence::maxCachePerThread));

// dangling pointer!!
Sequence::DataPtrListVector *Sequence::preAllocatedDataPtrs = new Sequence::DataPtrListVector();

// static methods of Sequence
void Sequence::clearCaches() {
	threadCacheSequences = CachedSequencesVector(OMP_MAX_THREADS, CachedSequences(Sequence::maxCachePerThread));
	preAllocatedDataPtrs->reset();
}

Sequence::CachedSequences &Sequence::getCachedSequencesForThread() {
	int threadNum = omp_get_thread_num();
	return threadCacheSequences[threadNum];
}

void Sequence::setThreadCache(Sequence::Sequence &mmapedSequence, Sequence::SequencePtr &expandedSequence) {
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
	_vec = _DataPtrListVector(OMP_MAX_THREADS, DataPtrList());
	_isValid = true;
}


// constructors and operators

Sequence::Sequence() : _flags(0) {
	reset(0);
}

Sequence::Sequence(const Sequence::Sequence &copy)  {
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

Sequence::Sequence &Sequence::operator=(const Sequence::Sequence &other) {
	if (this == &other)
		return *this;
	_flags = other._flags;
	_data = other._data;
	return *this;
}

Sequence::Sequence Sequence::clone() const {
	return Sequence(getFasta());
}

Sequence::~Sequence() {
	reset(0);
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
		if (usePreAllocation && size <= DataPtrListVector::size) {
			_data = preAllocatedDataPtrs->retrieveDataPtr();
			setFlag(PREALLOCATED);
		} else {
		    _data = DataPtr( TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(size) );
	    }
	} catch (...) {
		throw std::runtime_error(
				"Cannot allocate memory in Sequence::setSequence()");
	}

    *_getLength() = length;
	memcpy(_getTwoBitSequence(), buffer, buffSize);

    // free the buffer
	if (needMalloc)
	    delete [] buffer;

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
	setMarkups(markupElementSize, markups);
}
void Sequence::setSequence(Sequence::RecordPtr mmapRecordStart, long extraBytes, Sequence::RecordPtr mmapQualRecordStart) {
	reset(MMAPED);
	long size = sizeof(Sequence::RecordPtr)+extraBytes;
	if (mmapQualRecordStart != NULL)
		size += sizeof(Sequence::RecordPtr);
	try {
		_data = DataPtr( TwoBitSequenceBase::_TwoBitEncodingPtr::allocate(size) );
	} catch (...) {
		throw std::runtime_error("Can not allocate memory for Sequence::setSequence(mmap)");
	}

	*_getRecord() = mmapRecordStart;
	if (mmapQualRecordStart != NULL) {
		setFlag(HASQUALS);
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

string Sequence::getFastaNoMarkup(SequenceLengthType trimOffset) const {
	if ( !isValid() )
	    return string("");

	if (isMmaped()) {
		string name, bases, quals;
		readMmaped(name, bases, quals);
		SequenceLengthType len = bases.length();
		if (trimOffset < len)
			return bases.substr(0,len);
		else
			return bases;
	} else {
		SequenceLengthType len = getLength();
		if (trimOffset < len)
			len = trimOffset;
	    return TwoBitSequence::getFasta(getTwoBitSequence(), len);
	}
}
string Sequence::getFasta(SequenceLengthType trimOffset) const {
	if ( !isValid() )
		return string("");
	if (isDiscarded() || trimOffset <= 1) {
		// to support printing paired reads where 1 read is trimmed to 0
		return string(1, 'N');
	}
	string fasta;
	if (isMmaped()) {
		SequencePtr sequencePtr = getCache();
		fasta = sequencePtr->getFasta(trimOffset);
	} else {
	    SequenceLengthType len = getLength();
	    if (trimOffset < len)
		    len = trimOffset;
	    if (len <= 1) {
		    // to support printing paired reads where 1 read is trimmed to 0
		    return string(1, 'N');
	    }
	    fasta = TwoBitSequence::getFasta(getTwoBitSequence(), len);
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
	assert(isMmaped() && hasQuals());
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
	assert(isMmaped() && hasQuals());
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
		if (hasQuals()) {
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
	if ( !isValid() )
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

void Sequence::readMmaped(std::string &name, std::string &bases, std::string &quals) const {
	assert(isMmaped());
	// TODO fix hack on NULL lastPtr.  Presently only works for single-lined fastas
	RecordPtr record(getRecord()), lastRecord(NULL), qualRecord(NULL), lastQualRecord(NULL);
	if (hasQuals()) {
	    qualRecord = getQualRecord();
		lastQualRecord = NULL;
	}
	SequenceRecordParser::parse(record, lastRecord, name, bases, quals, qualRecord, lastQualRecord);
}
Sequence::SequencePtr Sequence::readMmaped(bool usePreAllocation) const {
	std::string name, bases, quals;
	readMmaped(name, bases, quals);
	return SequencePtr(new Sequence(bases, usePreAllocation));
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

// Read constructors and operators

Read::Read(std::string name, std::string fasta, std::string qualBytes, bool usePreAllocation) {
	setRead(name, fasta, qualBytes, usePreAllocation);
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
Read Read::clone() const {
	return Read(getName(), getFasta(), getQuals());
}


Read::ReadPtr Read::readMmaped(bool usePreAllocation) const {
	BaseLocationVectorType markups = _getMarkups();
	std::string name, bases, quals;
	Sequence::readMmaped(name, bases, quals);
	TwoBitSequence::applyMarkup(bases, markups);
	return ReadPtr(new Read(name, bases, quals, usePreAllocation));
}

bool Read::recordHasQuals() const {
	assert(isMmaped());
	if (hasQuals())
		return true;
	else
		// TODO make this more general
	    return *getRecord() == '@'; // FASTQ
}

ProbabilityBases Read::getProbabilityBases() const {
	std::string fasta = getFasta();
	std::string quals = getQuals();
	ProbabilityBases probs(fasta.length());
	for(int i = 0; i < (int) fasta.length(); i++) {
		double prob = qualityToProbability[ (unsigned char) quals[i] ];
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

void Read::setRead(std::string name, std::string fasta, std::string qualBytes, bool usePreAllocation) {
	assert(!isMmaped());
	if (fasta.length() != qualBytes.length())
		throw std::invalid_argument(
				"fasta length != qual length for name = " + name + " " + fasta + " " + qualBytes);

	// do not store quals if it is a reference
	if (qualBytes.length() > 1 && qualBytes[0] == REF_QUAL) {
		qualBytes.clear();
	}

	Sequence::setSequence(fasta, qualBytes.length() + (name.length() + 1), usePreAllocation);

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
  if (isDiscarded())
	  return;

  string fasta = getFasta();
  SequenceLengthType len = getLength();
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
	if (hasQuals())
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
	if ( !isValid() ) {
		return string("");
	} else if (isMmaped()) {
		string name, bases, quals;
		Sequence::readMmaped(name, bases, quals);
		return name;
	} else {
		return string(_getName());
	}
}

string Read::getQuals(SequenceLengthType trimOffset, bool forPrinting, bool unmasked) const {
  if ((isDiscarded() || trimOffset <= 1) && !unmasked) {
	  // to support printing paired reads where 1 read is trimmed to 0
	  return string(1, FASTQ_START_CHAR+1);
  }
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
	} else if (_getData() == NULL) {
		// This is an invalid / empty sequence...
	    return string("");
    } else {
		// to support printing paired reads where 1 read is trimmed to 0
		return string(1, FASTQ_START_CHAR+1);
	}
  }
}

string Read::toFastq(SequenceLengthType trimOffset, std::string label, bool unmasked) const {
	std::string fasta;
	if (unmasked)
		fasta = getFastaNoMarkup(trimOffset);
	else
		fasta = getFasta(trimOffset);
	return string('@' + getName() + (label.length() > 0 ? " " + label : "")
			+ "\n" + fasta + "\n+\n"
			+ getQuals(trimOffset, true, unmasked) + "\n");
}
string Read::toFasta(SequenceLengthType trimOffset, std::string label, bool unmasked) const {
	std::string fasta;
	if (unmasked)
		fasta = getFastaNoMarkup(trimOffset);
	else
		fasta = getFasta(trimOffset);
	return string('>' + getName() + (label.length() > 0 ? " " + label : "")
			+ "\n" + fasta + "\n");
}
string Read::toQual(SequenceLengthType trimOffset, std::string label) const {
	return string('>' + getName() + (label.length() > 0 ? " " + label : "")
			+ "\n" + getFormattedQuals(trimOffset) + "\n");
}

// TODO format with linebreaks
string Read::getFormattedQuals(SequenceLengthType trimOffset) const {
	string quals = getQuals(trimOffset, true);
	stringstream ss;
	SequenceLengthType len = quals.length();
	if (trimOffset < len)
		len = trimOffset;
	for (SequenceLengthType i = 0; i < len; i++) {
		ss << (int) quals[i] - 64 << ' ';
	}
	return ss.str();
}

std::string Read::toString() const {
	return getFastaNoMarkup() + "\t" + getQuals(MAX_SEQUENCE_LENGTH, true, true) + "\t" + getName();
}


/*------------------------------------ BaseQual ----------------------------------------*/

char BaseQual::getQualChar(double prob, bool ignoreLow) {
	if (prob > 0.999) {
		return Sequence::FASTQ_START_CHAR + 30;
	} else if (prob > 0.99) {
		return Sequence::FASTQ_START_CHAR + 20;
	} else if (prob > 0.9) {
		return Sequence::FASTQ_START_CHAR + 10;
	} else if (!ignoreLow) {
		return Sequence::FASTQ_START_CHAR + 1;
	} else {
		if (prob>=.7)
			return '7';
		else if (prob>=.6)
			return '6';
		else if (prob>=.5)
			return '5';
		else if (prob >= .4)
			return '4';
		else if (prob >= .3)
			return '3';
		else if (prob >= .2)
			return '2';
		else if (prob >= .1)
			return '1';
		else
			return ' ';
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
	throw;
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
//
// $Log: Sequence.cpp,v $
// Revision 1.34  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.33.2.3  2010-05-10 22:40:20  regan
// minor refactor -- replaced invalid flag
//
// Revision 1.33.2.2  2010-05-10 20:44:35  regan
// minor refactor moved code into cpp
//
// Revision 1.33.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.33  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.32.4.1  2010-05-06 18:45:36  regan
// broke it...
//
// Revision 1.32  2010-05-05 06:28:35  regan
// merged changes from FixPairOutput-20100504
//
// Revision 1.31.4.1  2010-05-05 05:57:53  regan
// fixed pairing
// fixed name to exclude labels and comments after whitespace
// applied some performance optimizations from other branch
// created FixPair application
//
// Revision 1.31.2.1  2010-05-02 04:38:39  regan
// replaced mmap quals flag with paired
//
// Revision 1.31  2010-05-01 21:57:53  regan
// merged head with serial threaded build partitioning
//
// Revision 1.30  2010-04-22 23:41:32  regan
// fixed a few bugs
//
// Revision 1.29.4.8  2010-04-30 23:53:14  regan
// attempt to fix a bug.  clearing Sequence caches when it makes sense
//
// Revision 1.29.4.7  2010-04-30 23:33:37  regan
// replaced read cache with LRU 3rd party cache
//
// Revision 1.29.4.6  2010-04-30 22:29:58  regan
// added comments about dangling pointer
//
// Revision 1.29.4.5  2010-04-30 21:53:52  regan
// reuse memory efficiently for cache lookups
//
// Revision 1.29.4.4  2010-04-29 06:58:43  regan
// reworked output parameters to include option to print unmasked reads
//
// Revision 1.29.4.3  2010-04-28 22:28:11  regan
// refactored writing routines
//
// Revision 1.29.4.2  2010-04-26 22:53:45  regan
// bugfixes
//
// Revision 1.29.4.1  2010-04-23 17:46:20  regan
// merged bugfixes from head
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


