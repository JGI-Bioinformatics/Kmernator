// $Header: /repository/PI_annex/robsandbox/KoMer/src/Sequence.cpp,v 1.24 2010-03-03 17:10:49 regan Exp $
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
static const std::tr1::shared_ptr<TwoBitEncoding> nullSequence(
		new TwoBitEncoding[0]);

Sequence::Sequence() {
	reset();
}
;

Sequence::Sequence(std::string name) :
	_length(0) {
	setSequence(name);
}

Sequence::~Sequence() {
	reset();
}

void Sequence::setSequence(std::string fasta) {
	setSequence(fasta, 0);
}

void Sequence::setSequence(std::string fasta, unsigned int extraBytes) {
	reset();
	_length = fasta.length();

	unsigned long buffSize = getTwoBitEncodingSequenceLength();
	TwoBitEncoding *buffer = (TwoBitEncoding*) malloc(buffSize);
	if (buffer == NULL)
		throw std::runtime_error(
				"Could not allocate buffer memory in Sequence::setSequence");

	const char *f_str = fasta.c_str();
	BaseLocationVectorType markupBases = TwoBitSequence::compressSequence(
			f_str, buffer);
	SequenceLengthType markupBasesSize = markupBases.size();
	try {
		unsigned int size = getTwoBitEncodingSequenceLength()
				+ sizeof(SequenceLengthType) + markupBasesSize
				* sizeof(BaseLocationVectorType) + extraBytes;
		_data = std::tr1::shared_ptr<unsigned char>(new unsigned char[size]);
	} catch (...) {
		throw new std::runtime_error(
				"Cannot allocate memory in Sequence::setSequence()");
	}

	memcpy(getTwoBitSequence(), buffer, getTwoBitEncodingSequenceLength());
	free(buffer);

	memcpy(_getMarkupBasesCount(), &markupBasesSize, sizeof(markupBasesSize));
	BaseLocationType *ptr = _getMarkupBases();
	for (SequenceLengthType i = 0; i < markupBasesSize; i++) {
		memcpy(ptr++, &markupBases[i], sizeof(BaseLocationType));
	}

}

void Sequence::reset() {
	_length = 0;
	_data = nullSequence;
}

string Sequence::getFasta(SequenceLengthType trimOffset) const {
	if (_data == nullSequence)
		return string("");
	SequenceLengthType len = getLength();
	if (trimOffset < len)
		len = trimOffset;
	if (len == 0) {
		// to support printing paired reads where 1 read is trimmed to 0
		return string(1, 'N');
	}
	string fasta = TwoBitSequence::getFasta(getTwoBitSequence(), len);
	TwoBitSequence::applyMarkup(fasta, *_getMarkupBasesCount(),
			_getMarkupBases());

	return fasta;
}

const TwoBitEncoding *Sequence::getTwoBitSequence() const {
	return _data.get();
}
TwoBitEncoding *Sequence::getTwoBitSequence() {
	return const_cast<TwoBitEncoding*> (constThis().getTwoBitSequence());
}

const SequenceLengthType *Sequence::_getMarkupBasesCount() const {
	return (SequenceLengthType *) (getTwoBitSequence()
			+ getTwoBitEncodingSequenceLength());
}
SequenceLengthType *Sequence::_getMarkupBasesCount() {
	return const_cast<SequenceLengthType*> (constThis()._getMarkupBasesCount());
}

const BaseLocationType *Sequence::_getMarkupBases() const {
	return (BaseLocationType *) (_getMarkupBasesCount() + 1);
}
BaseLocationType *Sequence::_getMarkupBases() {
	return const_cast<BaseLocationType *> (constThis()._getMarkupBases());
}

SequenceLengthType Sequence::getLength() const {
	return _length;
}

SequenceLengthType Sequence::getTwoBitEncodingSequenceLength() const {
	return TwoBitSequence::fastaLengthToTwoBitLength(getLength());;
}

BaseLocationVectorType Sequence::getMarkups() const {
	SequenceLengthType size = *_getMarkupBasesCount();
	BaseLocationVectorType markups(size);
	if (size > 0) {
		const BaseLocationType *ptr = _getMarkupBases();
		for (unsigned int i = 0; i < size; i++)
			markups.push_back(*(ptr++));
	}
	return markups;
}

/*------------------------------------ READ ----------------------------------------*/

double Read::qualityToProbability[256];
int Read::initializeQualityToProbability(unsigned char minQualityScore) {
	for (int i = 0; i < 256; i++) {
		qualityToProbability[i] = 0;
	}
	int start = FASTQ_START_CHAR;
	for (int i = start + minQualityScore; i < FASTQ_START_CHAR + 100; i++)
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

const char * Read::_getQual() const {
	return (char *) (_getMarkupBases() + *_getMarkupBasesCount());
}
char * Read::_getQual() {
	return const_cast<char*> (constThis()._getQual());
}

const char * Read::_getName() const {
	return (char *) (_getQual() + _qualLength());
}
char * Read::_getName() {
	return const_cast<char*> (constThis()._getName());
}

void Read::setRead(std::string name, std::string fasta, std::string qualBytes) {
	if (fasta.length() != qualBytes.length())
		throw new std::invalid_argument(
				"fasta length != qual length for name = " + name);

	// set quality to one byte if this is a reference (or fasta without quality scores)
	if (qualBytes.length() > 1 && qualBytes[0] == REF_QUAL)
		qualBytes = string(1, REF_QUAL);

	Sequence::setSequence(fasta, qualBytes.length() + (name.length() + 1));

	memcpy(_getQual(), qualBytes.c_str(), qualBytes.length());
	strcpy(_getName(), name.c_str());
}

void Read::zeroQuals(SequenceLengthType offset, SequenceLengthType length) {
	char *qualPtr = _getQual();
	if (*qualPtr == REF_QUAL)
		return;
	for (unsigned int i = 0; i < offset; i++)
		qualPtr++;
	for (unsigned int i = 0; i < length; i++)
		*(qualPtr++) = FASTQ_START_CHAR;
}

string Read::getName() const {
	if (_data == nullSequence)
		return string("");
	else
		return string(_getName());
}

string Read::getQuals(SequenceLengthType trimOffset, bool forPrinting) const {
	const char * qualPtr = _getQual();
	SequenceLengthType len = (trimOffset <= _length ? trimOffset : _length);
	if (len > 0 && *qualPtr == REF_QUAL) {
		if (forPrinting)
			return string(len, PRINT_REF_QUAL);
		else
			return string(len, REF_QUAL);
	} else if (len > 0) {
		return string(qualPtr, len);
	} else {
		// to support printing paired reads where 1 read is trimmed to 0
		return string(1, FASTQ_START_CHAR);
	}
}

SequenceLengthType Read::_qualLength() const {
	if (_length == 0)
		return 0;
	else if (*(_getQual()) == REF_QUAL)
		return 1;
	else
		return _length;
}

string Read::toFastq(SequenceLengthType trimOffset, std::string label) const {
	return string('@' + getName() + (label.length() > 0 ? " " + label : "")
			+ "\n" + getFasta(trimOffset) + "\n+\n"
			+ getQuals(trimOffset, true) + "\n");
}
string Read::getFormattedQuals(SequenceLengthType trimOffset) const {
	string quals = getQuals();
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


