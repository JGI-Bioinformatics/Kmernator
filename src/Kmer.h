//
// Kmernator/src/Kmer.h
//
// Author: Rob Egan
//
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/


#ifndef _KMER_H
#define _KMER_H

#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sys/time.h>

#include <boost/functional/hash.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>

#include "config.h"
#include "Options.h"
#include "TwoBitSequence.h"
#include "MemoryUtils.h"
#include "KmerTrackingData.h"
#include "lookup3.h"
#include "MmapTempFile.h"
#include "Log.h"

using namespace TwoBitSequenceBase;

#ifndef MAX_KMER_SIZE
#define MAX_KMER_SIZE 1024
#endif

typedef boost::shared_ptr<TwoBitEncoding> KmerSharedPtr;

class KmerSizer {
public:
	typedef Kmernator::KmerIndexType  IndexType;

private:
	SequenceLengthType _sequenceLength;
	SequenceLengthType _twoBitLength;
	IndexType _totalSize;

	static KmerSizer singleton;

	KmerSizer() :
		_sequenceLength(0) {
		_set(_sequenceLength);
	}

	void _set(SequenceLengthType sequenceLength) {
		_sequenceLength = sequenceLength;
		_twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(_sequenceLength);
		_totalSize = _twoBitLength;
	}

public:

	static inline KmerSizer &getSingleton() {
		return singleton;
	}

	static void set(SequenceLengthType sequenceLength) {
		KmerSizer & singleton = getSingleton();
		singleton._set(sequenceLength);
		LOG_DEBUG_OPTIONAL(2, Logger::isMaster(), "Set Kmer Size: " << sequenceLength);
	}

	static inline SequenceLengthType getSequenceLength() {
		return getSingleton()._sequenceLength;
	}
	static inline SequenceLengthType getTwoBitLength() {
		return getSingleton()._twoBitLength;
	}
	static inline IndexType getByteSize() {
		return getSingleton()._totalSize;
	}

};

// this is the ONLY options class that is okay to extend, as there are no member variables
class _KmerBaseOptions : public OptionsBaseInterface {
public:
	_KmerBaseOptions(SequenceLengthType defaultKmerSize = KmerSizer::getSequenceLength()) : kmerSize(defaultKmerSize) {
	}
	~_KmerBaseOptions() {}
	void _resetOptions() {
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {

		po::options_description opts("Kmer Options");

		opts.add_options()

						("kmer-size", po::value<SequenceLengthType>()->default_value(kmerSize), "kmer size.  A size of 0 will skip k-mer calculations");

		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		if (vm.count("kmer-size") == 0) {
			LOG_WARN(1, "There was no kmer size specified!");
			ret = false;
		}
		setOpt("kmer-size", kmerSize);
		KmerSizer::set(kmerSize);

		return ret;
	}
	SequenceLengthType &getKmerSize()
	{
		return kmerSize;
	}
private:
	SequenceLengthType kmerSize;

};
typedef OptionsBaseTemplate< _KmerBaseOptions > KmerBaseOptions;


#define TEMP_KMER(name)  TwoBitEncoding _stack_##name[KmerSizer::getByteSize()]; Kmer &name = (Kmer &)(_stack_##name);

class Kmer {
public:
	typedef Kmernator::KmerNumberType NumberType;
	typedef Kmernator::KmerIndexType  IndexType;
	typedef Kmernator::KmerSizeType   SizeType;

private:
	inline static boost::hash<NumberType> &getHasher() {
		static boost::hash<NumberType> hasher;
		return hasher;
	}
	Kmer(); // never construct, just use as cast

#ifdef STRICT_MEM_CHECK
	TwoBitEncoding _someData[MAX_KMER_SIZE]; // need somedata to hold a pointer and a large amount to avoid memory warnings
public:
	const void *_data() const {return _someData;}
	void *_data() {return _someData;}
#else
	// No data for you!!!
public:
	const void *_data() const {
		return this;
	}
	void *_data() {
		return this;
	}
#endif
	// safely returns lowbits 64-bit numeric version of any sized kmer
	static NumberType toNumber(const Kmer &kmer) {
		NumberType val;
		switch (KmerSizer::getTwoBitLength()) {
		case 1:
			val = (NumberType) *((boost::uint8_t *) kmer.getTwoBitSequence());
			break;
		case 2:
			val = (NumberType) *((boost::uint16_t *) kmer.getTwoBitSequence());
			break;
		case 3:
			val = (((NumberType) *((boost::uint16_t *) kmer.getTwoBitSequence())))| (((NumberType) *((boost::uint8_t *) (kmer.getTwoBitSequence() + 2))) << 16);
			break;
		case 4:
			val = (NumberType) *((boost::uint32_t *) kmer.getTwoBitSequence());
			break;
		case 5:
			val = (((NumberType) *((boost::uint32_t *) kmer.getTwoBitSequence()))) | (((NumberType) *((boost::uint8_t *) (kmer.getTwoBitSequence() + 4))) << 32);
			break;
		case 6:
			val = (((NumberType) *((boost::uint32_t *) kmer.getTwoBitSequence()))) | (((NumberType) *((boost::uint16_t *) (kmer.getTwoBitSequence() + 4))) << 32);
			break;
		case 7:
			val = (((NumberType) *((boost::uint32_t *) kmer.getTwoBitSequence()))) | (((NumberType) *((boost::uint16_t *) (kmer.getTwoBitSequence() + 4))) << 32) | (((NumberType) *((boost::uint8_t *) (kmer.getTwoBitSequence() + 6))) << 48);
			break;
		default:
			val = (NumberType) *((boost::uint64_t *) kmer.getTwoBitSequence());
		}
		return val;
	}

public:

	Kmer(const Kmer &copy) {
		*this = copy;
	}

	inline int compare(const Kmer &other) const {
		return memcmp(_data(), other._data(), getTwoBitLength());
	}

	Kmer &operator=(const Kmer &other) {
		if (this == &other)
			return *this;

		memcpy(_data(), other._data(), getTwoBitLength());
		return *this;
	}

	inline Kmer *get() {
		return this;
	}
	inline const Kmer *get() const {
		return this;
	}
	inline bool operator ==(const Kmer &other) const {
		return compare(other) == 0;
	}
	inline bool operator !=(const Kmer &other) const {
		return compare(other) != 0;
	}
	inline bool operator <(const Kmer &other) const {
		return compare(other) < 0;
	}
	inline bool operator <=(const Kmer &other) const {
		return compare(other) <= 0;
	}
	inline bool operator >(const Kmer &other) const {
		return compare(other) > 0;
	}
	inline bool operator >=(const Kmer &other) const {
		return compare(other) >= 0;
	}

	void swap(Kmer &other) {
		TEMP_KMER(temp);
		temp = other;
		other = *this;
		*this = temp;
	}

	inline TwoBitEncoding *getTwoBitSequence() {
		return (TwoBitEncoding *) _data();
	}
	inline const TwoBitEncoding *getTwoBitSequence() const {
		return (const TwoBitEncoding *) _data();
	}
	inline SequenceLengthType getTwoBitLength() const {
		return KmerSizer::getTwoBitLength();
	}
	inline SequenceLengthType getByteSize() const {
		return KmerSizer::getByteSize();
	}
	inline SequenceLengthType getLength() const {
		return KmerSizer::getSequenceLength();
	}
	inline SequenceLengthType getGC() const {
		return TwoBitSequence::getGC(getTwoBitSequence(), getLength());
	}

	void buildReverseComplement(Kmer &output) const {
		TwoBitSequence::reverseComplement(getTwoBitSequence(),
				output.getTwoBitSequence(), getLength());
	}

	void set(std::string fasta, bool leastComplement = false) {
		assert(fasta.length() == KmerSizer::getSequenceLength());
		if (leastComplement) {
			TEMP_KMER(temp);
			TwoBitSequence::compressSequence(fasta, temp.getTwoBitSequence());
			temp.buildLeastComplement(*this);
		} else {
			TwoBitSequence::compressSequence(fasta, this->getTwoBitSequence());
		}
	}

	// returns true if this is the least complement, false otherwise (output is least)
	bool buildLeastComplement(Kmer &output) const {
		buildReverseComplement(output);
		if (*this <= output) {
			output = *this;
			return true;
		} else {
			return false;
		}
	}

	std::string toFasta() const {
		return TwoBitSequence::getFasta(getTwoBitSequence(), 0, getLength());
	}
	std::string toFastaFull() const {
		return TwoBitSequence::getFasta(getTwoBitSequence(), 0, getTwoBitLength()
				* 4);
	}
	inline NumberType toNumber() const {
		return Kmer::toNumber(*this);
	}


	inline NumberType hash() const {

		//		NumberType number = toNumber();
		//		return Lookup8::hash2(&number, 1, 0xDEADBEEF);
		uint64_t hash = 0xDEADBEEF;
		uint32_t *pc, *pb;
		pc = (uint32_t*) &hash;
		pb = pc+1;
		Lookup3::hashlittle2(this, getTwoBitLength(), pc, pb);
		return hash;

	}

};


template<typename Value>
class KmerArray {

public:
	typedef Kmer::NumberType    NumberType;
	typedef Kmer::IndexType     IndexType;
	typedef Kmer::SizeType      SizeType;

	typedef Value ValueType;
	static const IndexType MAX_INDEX = MAX_KMER_INDEX;
	static const float GROWTH_FACTOR = 1.10;
	static const IndexType MIN_GROWTH = 16;

	typedef std::vector< KmerArray > Vector;

	class ElementType {
	private:
		const KmerArray *_array;
		IndexType _idx;
		inline const ElementType &_constThis() const {return *this;}
	public:
		ElementType() :
			_array(NULL), _idx(MAX_INDEX) {
		}
		ElementType(IndexType idx,
				const KmerArray &array) :
					_array(&array), _idx(idx) {
		}
		ElementType(const ElementType &copy) :
			_array(NULL), _idx(MAX_INDEX)  {
			*this = copy;
		}
		~ElementType() {
			reset();
		}
		ElementType &operator=(const ElementType &copy) {
			reset();
			_idx = copy._idx;
			_array = copy._array;
			return *this;
		}
		void reset() {
			_idx = MAX_INDEX;
			_array = NULL;
		}
		inline bool isValid() const {
			return _array != NULL && _idx != MAX_INDEX;
		}

		const Kmer &key() const {
			assert(isValid());
			return _array->get(_idx);
		}
		Kmer &key() {
			return const_cast<Kmer&> (_constThis().key());
		}
		const ValueType &value() const {
			assert(isValid());
			return _array->valueAt(_idx);
		}
		ValueType &value() {
			return const_cast<ValueType&> (_constThis().value());
		}
		inline bool operator==(const ElementType &other) const {
			if (isValid() && other.isValid())
				return value() == other.value() && key() == other.key();
			else
				return false;
		}
		inline bool operator<(const ElementType &other) const {
			if (isValid() && other.isValid())
				if (value() == other.value())
					return key() < other.key();
				else
					return value() < other.value();
			else
				return false;
		}
		inline bool operator>(const ElementType &other) const {
			if (isValid() && other.isValid())
				if (value() == other.value())
					return key() > other.key();
				else
					return value() > other.value();
			else
				return false;
		}
	};

private:
	void *_begin; // safer than Kmer *: prevents incorrect ptr arithmetic: _begin+1 , use _add instead
	// TODO embed size & capacity into allocation -- save some memory because of padding
	IndexType _size;
	IndexType _capacity;
	IndexType _endSorted;

	// if capacity == MAX_INDEX, this is a constant mmaped instance
	inline bool isMmaped() const {
		return _capacity == MAX_INDEX;
	}

private:
	static inline const void *_add(const void *ptr, IndexType i) {
		return ((char *) ptr + i * KmerSizer::getByteSize());
	}
	static inline void *_add(void *ptr, IndexType i) {
		return ((char *) ptr + i * KmerSizer::getByteSize());
	}

	// these are never thread safe!
	inline const Kmer &get(IndexType index) const {
		return *((Kmer *) _add(_begin, index));
	}
	inline Kmer &get(IndexType index) {
		return *((Kmer *) _add(_begin, index));
	}

	KmerArray(void *begin, IndexType size, IndexType capacity) : _begin(begin), _size(size), _capacity(capacity), _endSorted(0) {
	}

public:

	KmerArray(IndexType size = 0) :
		_begin(NULL), _size(0), _capacity(0), _endSorted(0) {
		resize(size);
	}

	KmerArray(const TwoBitEncoding *twoBit, SequenceLengthType length, bool leastComplement = false, bool *bools = NULL) :
		_begin(NULL), _size(0), _capacity(0), _endSorted(0)  {
		SequenceLengthType kmerSize = KmerSizer::getSequenceLength();
		assert(kmerSize > 0);
		if (kmerSize <= length) {
			SequenceLengthType numKmers = length
					- kmerSize + 1;
			resize(numKmers, MAX_INDEX, false);
			build(twoBit, length, leastComplement, bools);
		} else {
			resize(0);
		}
	}

	KmerArray(const KmerArray &copy) :
		_begin(NULL), _size(0), _capacity(0), _endSorted(0)  {
		*this = copy;
	}

	KmerArray copyRange(SequenceLengthType offset, SequenceLengthType length) {
		assert(offset+length <= _size);
		if (offset == 0 && length == _size)
			return KmerArray(*this);
		KmerArray splice(length);
		splice._copyRange(_begin, getValueStart(), 0, offset, length, false);
		return splice;
	}

	~KmerArray() {
		reset();
	}

	KmerArray &operator=(const KmerArray &other) {
		
		if (this == &other)
			return *this;
		reset();
		if (other.isMmaped()) {
			_begin = other._begin;
			_size = other._size;
			_capacity = other._capacity;
			_endSorted = other._endSorted;
			assert(_endSorted <= size());
			return *this;
		}

		resize(other.size());
		if (size() == 0) {
			return *this;
		}
		if (_begin == NULL)
			LOG_THROW(
					"RuntimeError: KmerArray::operator=(): Could not allocate memory in KmerArray operator=()");

		_copyRange(other._begin, other.getValueStart(), 0, 0, _size, false);
		_endSorted = other._endSorted;
		assert(_endSorted <= size());

		return *this;
	}

	// restore a new array from a mmap, allocating new memory
	KmerArray(const void *src) : _begin(NULL), _size(0), _capacity(0), _endSorted(0) {
		IndexType *size = (IndexType *) src;
		resize(*size);
		void *ptr = ++size;
		assert(KmerSizer::getSequenceLength() > 0);
		if (_size > 0) {
			long kmerSize  = _size * KmerSizer::getByteSize();
			_copyRange(ptr, (ValueType *) (((char*)ptr)+kmerSize), 0, 0, _size, false);
		}
		setLastSorted();
	}

	// store an existing array to a mmap
	const void *store(void *dst) const {
		// store IndexType + _size * getElementByteSize()
		IndexType *sizePtr = (IndexType *) dst;
		*sizePtr = _size;
		char *ptr = (char*) ++sizePtr;
		long kmerSize  = _size * KmerSizer::getByteSize();
		long valueSize = _size * sizeof(Value);
		_copyRange((void*)ptr, (Value*)(ptr+kmerSize), _begin, getValueStart(), 0, 0, _size, false);
		return ptr + kmerSize + valueSize;
	}
	SizeType sizeToStore() const {
		return sizeToStore(_size);
	}
	static SizeType sizeToStore(IndexType size) {
		return sizeof(IndexType) + size * ( KmerSizer::getByteSize() + sizeof(Value) );
	}
	// create a new array using existing memory
	static const KmerArray restore(const void *src) {
		const IndexType *size = (const IndexType *) src;
		IndexType xsize = *(size++);
		const void *ptr = size;
		KmerArray array(const_cast<void*>(ptr), xsize, MAX_INDEX);
		array.setLastSorted();
		return array;
	}

protected:
	void setLastSorted() {
		_endSorted = size() > 0 ? 1: 0;
		for(IndexType i = 1; i < size(); i++) {
			if (get(i-1).compare(get(i)) <= 0) {
				_endSorted++;
			} else {
				break;
			}
		}
		LOG_DEBUG(5, "_endSorted = " << _endSorted << ", " << size() << " for " << this);
		assert(_endSorted <= size());
	}
	// never thread safe!
	const ValueType *getValueStart() const {
		if (capacity() > 0)
			return (ValueType*) _add(_begin, capacity());
		else
			return NULL;
	}

	// never thread safe!
	ValueType *getValueStart() {
		if (capacity() > 0)
			return (ValueType*) _add(_begin, capacity());
		else
			return NULL;
	}

public:
	// never thread safe!
	const Kmer &operator[](IndexType index) const {
		if (index >= size())
			LOG_THROW(
					"Invalid: Kmer::operator[](): attempt to access index greater than size in KmerArray operator[] const");
		return get(index);
	}

	// never thread safe!
	Kmer &operator[](IndexType index) {
		if (index >= size())
			LOG_THROW(
					"Invalid: Kmer::operator[](): attempt to access index greater than size in KmerArray operator[]");
		return get(index);
	}

	// never thread safe!
	const ValueType &valueAt(IndexType index) const {
		if (index >= _size) {
			LOG_THROW(
					"Invalid: Kmer::valueAt(): attempt to access index greater than size in KmerArray valueAt() const");
		}
		return *(getValueStart() + index);
	}

	// never thread safe!
	ValueType &valueAt(IndexType index) {
		if (index >= _size) {
			LOG_THROW(
					"Invalid: Kmer::valueAt(): attempt to access index greater than size in KmerArray valueAt()");
		}
		return *(getValueStart() + index);
	}

	const Value *beginValue() const {
		return getValueStart();
	}
	const Value *endValue() const {
		return getValueStart() + _size;
	}
	Value *beginValue() {
		return getValueStart();
	}
	Value *endValue() {
		return getValueStart() + _size;
	}
	ValueType sumAll() const {
		ValueType sum = ValueType();
		for(const Value *it = beginValue(); it != endValue(); it++)
			sum += *it;
		return sum;
	}

public:
	const ElementType getElement(IndexType idx) const {
		return ElementType(idx, *this);
	}
	ElementType getElement(IndexType idx) {
		return ElementType(idx, *this);
	}

	static inline IndexType getElementByteSize() {
		return (KmerSizer::getByteSize() + sizeof(ValueType));
	}

	void reset(bool releaseMemory = true) {
		if (isMmaped()) {
			_begin = NULL;
			_size = 0;
			_capacity = 0;
			_endSorted = 0;
			return;
		}


		if (_begin != NULL && releaseMemory) {
			// destruct old Values
			for (IndexType i = 0; i < _capacity; i++)
				(getValueStart() + i)->~Value();
			// free memory
			std::free(_begin);

			_begin = NULL;
			_capacity = 0;
		}
		_size = 0;
		_endSorted = 0;

	}

	inline IndexType size() const {
		return _size;
	}
	inline IndexType capacity() const {
		if (isMmaped())
			return _size;
		else
			return _capacity;
	}
	inline bool empty() const {
		return _size == 0;
	}
	void reserve(IndexType size, bool reserveExtra = false) {
		assert(!isMmaped()); // mmaped can not be modified!

		if (size < _capacity) {
			IndexType oldSize = _size;
			resize(size, MAX_INDEX, reserveExtra);
			_size = oldSize;
		}
	}

	void resize(IndexType size, bool reserveExtra = false) {
		resize(size, MAX_INDEX, reserveExtra);
	}
	void resize(IndexType size, IndexType idx, bool reserveExtra = true) {
		assert(!isMmaped()); // mmaped can not be modified!

		if (size == _size)
			return;
		IndexType oldSize = _size;

		// alloc / realloc memory
		_setMemory(size, idx, reserveExtra);

		if (_begin == NULL && size > 0) {
			LOG_THROW(
					"RuntimeError: KmerArray::resize(): Could not allocate memory in KmerArray resize()");
		}

		if (false && size > oldSize && idx == MAX_INDEX) {
			// zero fill remainder
			void *start = _add(_begin, oldSize);
			memset(start, 0, KmerSizer::getByteSize() * (size - oldSize));

			// Values should already have been constructed
		}
	}

	void _copyRange(const void * srcKmer, const ValueType *srcValue,
			IndexType idx, IndexType srcIdx, IndexType count,
			bool isOverlapped = false) {
		assert(!isMmaped()); // mmaped can not be modified!
		_copyRange(_begin, getValueStart(), srcKmer, srcValue, idx, srcIdx, count, isOverlapped);
	}
	static void _copyRange(void * dstKmer, ValueType *dstValue,
			const void * srcKmer, const ValueType *srcValue,
			IndexType idx, IndexType srcIdx, IndexType count,
			bool isOverlapped = false) {

		if (isOverlapped)
			memmove(_add(dstKmer, idx), _add(srcKmer, srcIdx), count
					* KmerSizer::getByteSize());
		else
			memcpy(_add(dstKmer, idx), _add(srcKmer, srcIdx), count
					* KmerSizer::getByteSize());

		//assignment copy Values
		if (isOverlapped) {// && dstValue + idx > srcValue + srcIdx) {
			bool needNew = count > 1024;
			Value _tmp[ needNew ? 0 : count];
			Value *tmp = _tmp;
			if (needNew) {
				tmp = new Value[count];
			}

			for (IndexType i = 0; i < count; i++)
				tmp[i] = *(srcValue + srcIdx + i);
			for (IndexType i = 0; i < count; i++)
				*(dstValue + idx + i) = tmp[i];

			if (needNew)
				delete [] tmp;

		} else {
			// copy forward, without a temp buffer, is safe
			for (IndexType i = 0; i < count; i++)
				*(dstValue + idx + i) = *(srcValue + srcIdx + i);
		}
	}

	void _setMemory(IndexType size, IndexType idx, bool reserveExtra = true) {
		assert(!isMmaped()); // mmaped can not be modified!
		LOG_DEBUG(4, "_setMemory(" << this << ", " << size << ", " << idx << ", " << reserveExtra << "): " << _begin << " size: " << _size << " capacity: " << _capacity);

		if (reserveExtra && idx == MAX_INDEX && size <= _capacity) {
			_size = size;
			return;
		}

		// preserve old pointers and metrics
		void *oldBegin = _begin;
		IndexType oldSize = _size;
		IndexType oldCapacity = _capacity;

		void *newBegin = oldBegin;
		IndexType newCapacity = oldCapacity;

		IndexType lesserSize = std::min(size, oldSize);
		ValueType *oldValueStart = NULL;
		if (oldCapacity > 0) {
			oldValueStart = getValueStart();
		}

		bool memChanged = false;
		if (size == 0) {
			// ignore reserveExtra for zero resize
			if (oldCapacity > 0) {
				newCapacity = 0;
				newBegin = NULL;
				memChanged = true;
			}
			_endSorted = 0;
		} else if ((size > oldCapacity) || (oldCapacity > size && !reserveExtra)) {
			// allocate new memory
			if (reserveExtra)
				newCapacity = std::max((IndexType) (size * GROWTH_FACTOR),
						(IndexType) (size + MIN_GROWTH));
			else
				newCapacity = size;

			newBegin = std::malloc(newCapacity * getElementByteSize());
			if (newBegin == NULL) {
				LOG_THROW("RuntimeError: Could not allocate memory in KmerArray _setMemory(): Attempt to malloc " << newCapacity
						* getElementByteSize() << " (" << newCapacity << " elements)");
			}

			memChanged = true;
		}

		ValueType *newValueStart = NULL;
		if (newCapacity > 0) {
			newValueStart = (ValueType*) _add(newBegin, newCapacity);
		}

		if (memChanged && newBegin != NULL) {
			// construct all values that are newly allocated
			for (IndexType i = 0; i < newCapacity; i++)
				::new ((void*) (newValueStart + i)) Value();
		}

		if (newBegin != NULL && oldBegin != NULL && oldValueStart != NULL
				&& lesserSize > 0) {
			// copy the old contents
			if (idx == MAX_INDEX || idx >= lesserSize) {
				if (memChanged) {
					// copy all records in order (default ; end is trimmed or expanded)
					_copyRange(newBegin, newValueStart, oldBegin, oldValueStart, 0, 0, lesserSize);
				} else {
					// noop
				}
			} else {
				// shrink or expand.

				if (memChanged && idx > 0) {
					// the first record(s) leading to idx will be copied
					_copyRange(newBegin, newValueStart, oldBegin, oldValueStart, 0, 0, idx);
				}

				if (lesserSize == size) {
					// shrink: skipping the old record at idx
					if (idx < size) {
						_copyRange(newBegin, newValueStart, oldBegin, oldValueStart, idx, idx + 1,
								lesserSize - idx, !memChanged);
					}
				} else {
					// expand: leaving new (uninitialized/unchanged) record at idx
					if (idx < oldSize) {
						_copyRange(newBegin, newValueStart, oldBegin, oldValueStart, idx + 1, idx,
								lesserSize - idx, !memChanged);
					}
				}
			}
		}

		if (memChanged) {
			if (newCapacity > oldCapacity) {
				// change begin first
				_begin = newBegin;
				_capacity = newCapacity;
			} else {
				// change size & capacity first
				_size = size;
				_capacity = newCapacity;
				_begin = newBegin;
			}
		}
		_size = size;
		if (memChanged && oldBegin != NULL) {
			// destruct old Values
			for (IndexType i = 0; i < oldCapacity; i++)
				(oldValueStart + i)->~Value();
			// free old memory
			std::free(oldBegin);
		}
		//LOG_DEBUG(5, "_setMemory(" << this << ", " << size << ", " << idx << ", " << reserveExtra << "): exited. memChanged:" << memChanged << " - " << _begin << " size: " << _size << " capacity: " << _capacity);
	}

	void build(const TwoBitEncoding *twoBit, SequenceLengthType length,
			bool leastComplement = false, bool *bools = NULL) {
		assert(!isMmaped()); // mmaped can not be modified!
		if (length < KmerSizer::getSequenceLength()) {
			resize(0, _capacity > 0);
			return;
		}
		SequenceLengthType numKmers = length - KmerSizer::getSequenceLength() + 1;
		if (numKmers != _size)
			resize(numKmers, MAX_INDEX, _capacity > 0);

		KmerArray &kmers = *this;
		long numBytes = (numKmers + 3) / 4;

#pragma omp parallel for if(numKmers >= 10000)
		for (long bytes = 0; bytes < numBytes; bytes++) {
			SequenceLengthType i = bytes * 4;
			const TwoBitEncoding *ref = twoBit + bytes;
			for (int bitShift = 0; bitShift < 4 && i + bitShift < numKmers; bitShift++) {
				TwoBitSequence::shiftLeft(ref, kmers[i + bitShift].get(),
						KmerSizer::getTwoBitLength(), bitShift, bytes < numBytes-1);
				TwoBitEncoding *lastByte =
						kmers[i + bitShift].getTwoBitSequence()
						+ KmerSizer::getTwoBitLength() - 1;
				switch (KmerSizer::getSequenceLength() & 0x03) {
				case 1:
					*lastByte &= 0xc0;
					break;
				case 2:
					*lastByte &= 0xf0;
					break;
				case 3:
					*lastByte &= 0xfc;
					break;
				}
			}
		}
		if (leastComplement) {
			bool isLeast;

			long longNumKmers = numKmers;
#pragma omp parallel for if(longNumKmers >= 10000)
			for (long i = 0; i < longNumKmers; i++) {
				TEMP_KMER(least);
				isLeast = kmers[i].buildLeastComplement(least);
				if (!isLeast)
					kmers[i] = least;
				if (bools != NULL)
					*(bools+i) = isLeast;
			}

		}
	}

	static SequenceLengthType _numPermutations(SequenceLengthType len, short editDistance) {
		SequenceLengthType s = 1;
		if ((SequenceLengthType) editDistance > len)
			editDistance = len;
		for(short e = 0; e < editDistance && (SequenceLengthType) e != len ; e++) {
			s *= 3*(len-e);
		}
		for (short e = editDistance ; e > 1 ; e--) {
			s /= e;
		}
		return s;
	}
	static void permuteBases(const Kmer &kmer, KmerArray &kmers, short editDistance, bool leastComplement = false) {
		if ((SequenceLengthType) editDistance > KmerSizer::getSequenceLength())
			editDistance = KmerSizer::getSequenceLength();

		SequenceLengthType size = 0;
		for(short e = 1 ; e <= editDistance; e++)
			size += _numPermutations(KmerSizer::getSequenceLength(), e);
		kmers.reset(false);
		kmers.resize(size);
		SequenceLengthType offset = 0;

		offset = __permuteBases(kmer, kmers, offset, 0, editDistance, leastComplement);
		if (offset != kmers.size()) {
			LOG_WARN(1, "Mismatching permute bases size: " << offset << " vs " << kmers.size());
			LOG_WARN(1, "permuteBases(" << kmer.toFasta() << ", " << editDistance << ", " << leastComplement << ")");
			//for(SequenceLengthType i = 0; i < offset; i++)
			//	LOG_WARN(1, kmers[i].toFasta());
			kmers.resize(offset);
		}
	}
	static SequenceLengthType __permuteBases(const Kmer &kmer, KmerArray &kmers, SequenceLengthType offset, SequenceLengthType startIdx, short editDistance, bool leastComplement = false) {
		if (editDistance == 0) {
			return offset;
		} else {
			for(SequenceLengthType baseIdx = startIdx; baseIdx < KmerSizer::getSequenceLength(); baseIdx++) {
				Kmer &v1 = kmers[offset++];
				Kmer &v2 = kmers[offset++];
				Kmer &v3 = kmers[offset++];
				TwoBitSequence::permuteBase(kmer.getTwoBitSequence(), v1.getTwoBitSequence(), v2.getTwoBitSequence(), v3.getTwoBitSequence(),
						KmerSizer::getSequenceLength(), baseIdx);
				if (editDistance > 1) {
					offset = __permuteBases(v1, kmers, offset, baseIdx+1, editDistance-1, leastComplement);
					offset = __permuteBases(v2, kmers, offset, baseIdx+1, editDistance-1, leastComplement);
					offset = __permuteBases(v3, kmers, offset, baseIdx+1, editDistance-1, leastComplement);
				}
			}
		}
		return offset;
	}
	// return a KmerArray that has one entry for each possible single-base substitution
	static KmerArray permuteBases(const Kmer &kmer, bool leastComplement = false) {
		KmerArray kmers(KmerSizer::getSequenceLength() * 3);
		_permuteBases(kmer, kmers, leastComplement);
		return kmers;
	}
	static KmerArray permuteBases(const Kmer &kmer, const ValueType defaultValue, bool leastComplement = false) {
		KmerArray kmers(KmerSizer::getSequenceLength() * 3);
		_permuteBases(kmer, kmers, leastComplement);
		for(SequenceLengthType idx = 0; idx < kmers.size(); idx++)
			kmers.valueAt(idx) = defaultValue;
		return kmers;
	}
	static void _permuteBases(const Kmer &kmer, KmerArray &kmers, bool leastComplement, SequenceLengthType offset = 0) {
		TEMP_KMER(tmp);
		for (SequenceLengthType baseIdx = 0; baseIdx < KmerSizer::getSequenceLength(); baseIdx++) {
			TwoBitSequence::permuteBase(kmer.getTwoBitSequence(), kmers[offset+baseIdx*3].getTwoBitSequence(), kmers[offset+baseIdx*3+1].getTwoBitSequence(), kmers[offset+baseIdx*3+2].getTwoBitSequence(),
					KmerSizer::getSequenceLength(), baseIdx);
			if (leastComplement) {
				for(int i = 0 ; i < 3; i++) {
					if (! kmers[offset+baseIdx*3+i].buildLeastComplement(tmp) ) {
						kmers[offset+baseIdx*3+i] = tmp;
					}
				}
			}
		}
	}

	static KmerArray extendKmer(const Kmer &kmer, bool toRight, bool leastComplement = false) {
		KmerArray kmers(4);
		std::string fasta = kmer.toFasta().substr(toRight ? 1: 0, KmerSizer::getSequenceLength() - 1);

		if (toRight) {
			TwoBitSequence::compressSequence(fasta + "A", kmers[0].getTwoBitSequence());
			TwoBitSequence::compressSequence(fasta + "C", kmers[1].getTwoBitSequence());
			TwoBitSequence::compressSequence(fasta + "G", kmers[2].getTwoBitSequence());
			TwoBitSequence::compressSequence(fasta + "T", kmers[3].getTwoBitSequence());
		} else {
			TwoBitSequence::compressSequence("A" + fasta, kmers[0].getTwoBitSequence());
			TwoBitSequence::compressSequence("C" + fasta, kmers[1].getTwoBitSequence());
			TwoBitSequence::compressSequence("G" + fasta, kmers[2].getTwoBitSequence());
			TwoBitSequence::compressSequence("T" + fasta, kmers[3].getTwoBitSequence());

		}
		if (leastComplement) {
			TEMP_KMER(tmp);
			for(int i = 0; i < 4; i++) {
				tmp = kmers[i];
				tmp.buildLeastComplement(kmers[i]);
			}
		}
		return kmers;
	}
protected:
	IndexType _find(const Kmer &target, IndexType start, IndexType end) const {
		for(IndexType i=start; i<end; i++) {
			if (target.compare(get(i)) == 0) {
				return i;
			}
		}
		return MAX_INDEX;
	}
public:
	IndexType find(const Kmer &target) const {
		bool targetIsFound;
		IndexType idx;
		if (_endSorted >= 16) {
			idx = _findSorted(target, targetIsFound, 0, _endSorted);
			if (targetIsFound)
				return idx;
		}
		idx = _find(target, _endSorted >= 16 ? _endSorted : 0, size());
		return idx;
	}
	IndexType find(const Kmer &target, bool &targetIsFound) const {
		IndexType idx = find(target);
		targetIsFound = idx != MAX_INDEX;
		return idx;
	}
protected:
	IndexType _findSorted(const Kmer &target, bool &targetIsFound, IndexType start, IndexType end) const {
		// binary search
		if (end <= start)
		{
			targetIsFound = false;
			return 0;
		}
		IndexType min = start;
		IndexType max = end;
		assert(min >= 0);
		assert(min < size());
		assert(max <= size());
		assert(max <= _endSorted);
		assert(max > 0);

		max--; // max must be a valid index... never cross the end boundary
		IndexType mid;
		int comp;
		do {
			mid = (min+max) / 2;
			assert(mid < size());
			comp = target.compare(get(mid));
			if (comp > 0)
				min = mid+1;
			else if (comp < 0)
				max = mid-1;
		} while (comp != 0 && max != MAX_INDEX && min <= max);
		if (comp == 0)
			targetIsFound = true;
		else
			targetIsFound = false;

		return mid + (comp>0 && end>mid?1:0);
	}
	IndexType _findSorted(const Kmer &target, bool &targetIsFound) const {
		return _findSorted(target, targetIsFound, 0, size());
	}
public:
	IndexType findSorted(const Kmer &target, bool &targetIsFound) const {
		IndexType idx = _findSorted(target, targetIsFound);
		return idx;
	}

protected:
	void _insertAt(IndexType idx, const Kmer &target) {
		assert(!isMmaped()); // mmaped can not be modified!
		if (idx > size())
			LOG_THROW("Invalid: attempt to access index greater than size in KmerArray insertAt");
		resize(size() + 1, idx);
		get(idx) = target;
	}
	void _insertAt(IndexType idx, const Kmer &target, const Value &value) {
		_insertAt(idx, target);
		valueAt(idx) = value;
	}

public:
	void insertAt(IndexType idx, const Kmer &target) {
		_insertAt(idx, target);
	}
	void insertAt(IndexType idx, const Kmer &target, const Value &value) {
		_insertAt(idx,target,value);
	}

	IndexType append(const Kmer &target) {
		IndexType idx = size();
		_insertAt(idx, target);
		return idx;
	}
	IndexType append(const Kmer &target, const Value &value) {
		IndexType idx = size();
		_insertAt(idx, target, value);
		return idx;
	}

protected:
	IndexType _insertSorted(const Kmer &target) {
		assert(isSorted());
		bool isFound;
		IndexType idx = findSorted(target, isFound);
		if (!isFound) {
			_insertAt(idx, target);
			_endSorted++;
			assert(_endSorted <= size());
		}
		return idx;
	}
	IndexType _insertSorted(const Kmer &target, const Value &value) {
		IndexType idx = _insertSorted(target);
		valueAt(idx) = value;
		return idx;
	}

public:
	IndexType insertSorted(const Kmer &target) {
		IndexType idx = _insertSorted(target);
		return idx;
	}
	IndexType insertSorted(const Kmer &target, const Value &value) {
		IndexType idx = _insertSorted(target,value);
		return idx;
	}

	void remove(const Kmer &target) {
		bool isFound;
		IndexType idx = find(target, isFound);
		if (isFound)
			remove(idx);
	}
	void remove(IndexType idx) {
		assert(!isMmaped()); // mmaped can not be modified!
		resize(size()-1,idx);
		if (_endSorted > idx)
			_endSorted--;
	}

	class CompareArrayIdx {
		const KmerArray &_kmerArray;
	public:
		CompareArrayIdx(const KmerArray &kmerArray): _kmerArray(kmerArray) {}
		bool operator()(IndexType i, IndexType j) {
			return _kmerArray.get(i).compare(_kmerArray.get(j)) <= 0;
		}
	};

	bool isSorted() const {
		return _endSorted == size();
	}

	IndexType numUnsorted() const {
		return _endSorted - size();
	}

	// sort the elements in the KmerArray
	// array will be resized
	void resort(bool reserveExtra = false) {
		if (isSorted())
			return;

		std::vector< IndexType > passingIndexes;
		passingIndexes.reserve(size());

		for(IndexType i = 0; i < size(); i++)
			passingIndexes.push_back(i);
		resort(passingIndexes, reserveExtra);
	}
private:
	void resort(std::vector< IndexType > &passingIndexes, bool reserveExtra = false) {
		std::sort(passingIndexes.begin(), passingIndexes.end(), CompareArrayIdx(*this));

		KmerArray s;
		IndexType passingSize = passingIndexes.size();
		s.reserve(passingSize + (reserveExtra ? std::max((IndexType) (passingSize * GROWTH_FACTOR + 0.5), (IndexType) (passingSize + MIN_GROWTH)) : 0));

		for(IndexType i = 0; i < passingSize; i++)
			s.append(get(passingIndexes[i]), valueAt(passingIndexes[i]));

		s._endSorted = passingSize;
		swap(s);
		LOG_DEBUG(5, "resort()" << _endSorted << ", " << size() << ", " << passingSize << " on " << this);
		assert(_endSorted <= size());
	}
public:
	void swap(IndexType idx1, IndexType idx2) {
		assert(!isMmaped()); // mmaped can not be modified!
		if (idx1 == idx2)
			return;
		if (idx1 >= size() || idx2 >= size())
			LOG_THROW("Invalid: attempt to access index greater than size in KmerArray swap()");

		get(idx1).swap(get(idx2));
		if (sizeof(ValueType) > 0) {
			ValueType tmp = valueAt(idx1);
			valueAt(idx1) = valueAt(idx2);
			valueAt(idx2) = tmp;
		}
	}
	void swap(KmerArray &other) {
		std::swap(_begin, other._begin);
		std::swap(_size, other._size);
		std::swap(_capacity, other._capacity);
		std::swap(_endSorted, other._endSorted);
		assert(_endSorted <= size());
	}

	// purge all element where the value is less than minimumCount
	// assumes that ValueType can be cast into long
	// resulting in a (potententially smaller) sorted KmerArray
	IndexType purgeMinCount(long minimumCount) {
		// scan values that pass, keep list and count
		IndexType maxSize = size();
		std::vector< IndexType > passingIndexes;
		passingIndexes.reserve(maxSize);
		IndexType affected;

		ValueType *valuePtr = getValueStart();
		for(IndexType i = 0; i < maxSize; i++) {
			if ( minimumCount <= (long) *(valuePtr++) ) {
				passingIndexes.push_back( i );
			}
		}
		if (passingIndexes.empty()) {
			reset(true);
			affected = maxSize;
		} else {
			affected = maxSize - passingIndexes.size();
			resort(passingIndexes);
		}

		return affected;
	}

	std::string toThis() const {
		std::stringstream ss;
		ss << "KmerArray " << this << "={" << _begin << ", " << _size << ", " << _capacity << "}";
		return ss.str();
	}

	std::string toString(bool inclMemPtr = false) const {
		std::stringstream ss;
		ss << "{";
		IndexType idx=0;
		for(idx=0; idx<size() && idx < 30; idx++) {
			ss << get(idx).toFasta() << ":" << valueAt(idx) << ", ";
		}
		if (idx < size())
			ss << " ... " << size() - idx << " more ";
		ss << "}";
		return ss.str();
	}

public:
	class Iterator : public std::iterator<std::forward_iterator_tag, KmerArray>
	{
	private:
		KmerArray *_tgt;
		IndexType _idx;
		ElementType thisElement;

		void setElement() {if ( !isEnd() ) thisElement = _tgt->getElement(_idx);}

	public:
		Iterator(KmerArray *target, IndexType idx = 0): _tgt(target), _idx(idx) {
			setElement();
		}
		Iterator(const Iterator &copy) {*this = copy;}
		Iterator() : _tgt(NULL), _idx(0), thisElement() {}

		~Iterator() {}
		Iterator& operator=(const Iterator& other) {
			_tgt = other._tgt;
			_idx = other._idx;
			thisElement = other.thisElement;
			return *this;
		}
		bool operator==(const Iterator& other) const {return _idx == other._idx && _tgt == other._tgt;}
		bool operator!=(const Iterator& other) const {return _idx != other._idx || _tgt != other._tgt;}
		Iterator& operator++() {if ( !isEnd() ) {++_idx; setElement();} return *this;}
		Iterator operator++(int unused) {Iterator tmp(*this); ++(*this); return tmp;}
		ElementType &operator*() {return thisElement;}
		Kmer &key() {return thisElement.key();}
		Value &value() {return thisElement.value();}
		ElementType *operator->() {return &thisElement;}
		bool isEnd() const {if (_tgt == NULL) { return true; } else { return _idx >= _tgt->size(); }}
	};
	class ConstIterator : public std::iterator<std::forward_iterator_tag, KmerArray>
	{
	private:
		Iterator _iterator;
	public:
		ConstIterator(KmerArray *target, IndexType idx = 0): _iterator(target,idx) {}
		ConstIterator(const ConstIterator &copy) {*this = copy;}
		ConstIterator() : _iterator() {}
		ConstIterator& operator=(const ConstIterator& other) {
			_iterator = other._iterator;
			return *this;
		}
		bool operator==(const ConstIterator& other) const {return _iterator == other._iterator;}
		bool operator!=(const ConstIterator& other) const {return _iterator != other._iterator;}
		ConstIterator& operator++() {++_iterator; return *this;}
		ConstIterator operator++(int unused) {return ConstIterator(_iterator++);}
		const ElementType &operator*() const {return _iterator.operator*();}
		const Kmer &key() const {return _iterator.key();}
		const Value &value() const {return _iterator.value();}
		const ElementType *operator->() {return _iterator.operator->();}
		bool isEnd() const {return _iterator.isEnd();}
	};

	Iterator begin() {return Iterator(this, 0);}
	Iterator end() {return Iterator(this, size());}

	ConstIterator begin() const {return ConstIterator(this,0);}
	ConstIterator end() const {return ConstIterator(this,size());}
};

template<typename Value>
class KmerMap {
public:
	static const bool defaultSort = false;
	typedef Kmer::NumberType    NumberType;
	typedef Kmer::IndexType     IndexType;
	typedef Kmer::SizeType      SizeType;

	typedef Kmer KeyType;
	typedef Value ValueType;
	typedef NumberType * NumberTypePtr;
	typedef KmerArray<Value> BucketType;
	typedef	typename BucketType::Iterator BucketTypeIterator;
	typedef typename BucketType::ElementType ElementType;

	typedef std::vector< BucketType > BucketsVector;
	typedef typename BucketsVector::iterator BucketsVectorIterator;

	static const IndexType MAX_UNSORTED = 24;

	static IndexType getMinPowerOf2(IndexType minBucketCount) {
		NumberType powerOf2 = minBucketCount;
		if (powerOf2 == 0) {
			powerOf2 = 1;
		} else if ( (powerOf2 & (powerOf2 -1)) == 0 ) {
			// argument is a power of 2
		} else {
			powerOf2--;
			for (size_t i = 1; i < sizeof(NumberType)*8; i<<=1)
				powerOf2 |= powerOf2 >> i;
			powerOf2++;
		}
		return powerOf2;
	}
private:
	BucketsVector _buckets;
	NumberType BUCKET_MASK;
	bool _isSorted;

	inline const KmerMap &_constThis() const {return *this;}

public:
	KmerMap() {
		BUCKET_MASK=0;
		_buckets.clear();
		_isSorted = defaultSort;
	}
	KmerMap(IndexType bucketCount) {

		// ensure buckets are a precise powers of two
		// with at least bucketCount buckets
		if (bucketCount > MAX_KMER_MAP_BUCKETS)
			bucketCount = MAX_KMER_MAP_BUCKETS;
		NumberType powerOf2 = getMinPowerOf2(bucketCount);
		BUCKET_MASK = powerOf2 - 1;
		_buckets.resize(powerOf2);
		_isSorted = defaultSort;

		LOG_DEBUG_OPTIONAL(3, Logger::isMaster(), "KmerMap(" << bucketCount << "): " << powerOf2);

	}
	~KmerMap()
	{
		clear();
	}
	KmerMap &operator=(const KmerMap &other) {
		_buckets = other._buckets;
		BUCKET_MASK = other.BUCKET_MASK;
		_isSorted = other._isSorted;
		return *this;
	}
	// restore new instance from mmap
	KmerMap(const void *src) {
		NumberType size(0), *offsetArray;
		const void *ptr;
		_getMmapSizes(src, size, BUCKET_MASK, offsetArray);
		_buckets.resize(size);
		_isSorted = true;
		for (NumberType idx = 0 ; idx < size ; idx++) {
			ptr = ((char*)src) + *(offsetArray + idx);
			_buckets[idx] = BucketType(ptr);
			_isSorted &= _buckets[idx].isSorted();
		}
	}
	inline bool isSorted() const {
		return _isSorted;
	}
	const Kmernator::MmapFile store(std::string permanentFile = "") const {
		long size = getSizeToStore();
		Kmernator::MmapFile mmap = MmapTempFile::buildNewMmap(size, permanentFile);
		store(mmap.data());
		return mmap;
	}
	// store to mmap
	const void *store(void *dst) const {
		NumberType size = (NumberType) _buckets.size();
		NumberType *numbers = (NumberType *) dst;
		*(numbers++) = size;
		*(numbers++) = BUCKET_MASK;
		NumberType *offsetArray = numbers;
		NumberType offset = sizeof(NumberType) * (2+size);

		for(NumberType idx = 0 ; idx < size; idx++) {
			*(offsetArray++) = offset;
			void *ptr = ((char*)dst) + offset;
			const char *newPtr = (const char *) _buckets[idx].store(ptr);
			NumberType newSize = (newPtr - (const char *) ptr);
			offset += newSize;
		}
		return ((char*)dst) + offset;
	}
	static const KmerMap restore(const void *src) {
		KmerMap map;
		NumberType size(0), *offsetArray;
		const void *ptr;
		_getMmapSizes(src, size, map.BUCKET_MASK, offsetArray);
		map._buckets.resize(size);
		map._isSorted = true;
		for (NumberType idx = 0 ; idx < size ; idx++) {
			ptr = ((char*)src) + *(offsetArray + idx);
			map._buckets[idx] = BucketType::restore(ptr);
			map._isSorted &= map._buckets[idx].isSorted();
		}
		return map;
	}

	void swap(KmerMap &other) {
		_buckets.swap(other._buckets);
		std::swap(BUCKET_MASK, other.BUCKET_MASK);
		std::swap(_isSorted, other._isSorted);
	}

	static const void _getMmapSizes(const void *src, NumberType &size, NumberType &mask, NumberTypePtr &offsetArray) {
		NumberType *numbers = (NumberType *) src;
		size = *(numbers++);
		mask = *(numbers++);
		offsetArray = numbers;
	}
	SizeType getSizeToStore() const {
		return getSizeToStoreCountsAndIndexes()
				+ getSizeToStoreElements();
	}
	SizeType getSizeToStoreCountsAndIndexes() const {
		return sizeof(NumberType)*(2+_buckets.size())
				+ _buckets.size()*sizeof(IndexType);
	}
	SizeType getSizeToStoreElements() const {
		return size()*BucketType::getElementByteSize();
	}
	NumberType getBucketMask() const {
		return BUCKET_MASK;
	}
	NumberType getBucketSize() const {
		return _buckets.size();
	}
	void reset(bool releaseMemory = true) {
		for(size_t i=0; i< _buckets.size(); i++) {
			_buckets[i].reset(releaseMemory);
		}
	}
	void clear(bool releaseMemory = true) {
		reset(releaseMemory);
		if (releaseMemory)
			_buckets.resize(0);
	}

	void resort() {
		if (!_isSorted) {
			LOG_DEBUG_OPTIONAL(1, Logger::isMaster(), "Sorting KmerMaps");
#pragma omp parallel for
			for(long i=0; i< (long) _buckets.size(); i++) {
				_buckets[i].resort();
			}
			_isSorted = true;
		}
	}

	void setReadOnlyOptimization() {
		for(size_t i = 0; i<_buckets.size(); i++) {
			_buckets[i].setReadOnlyOptimization( );
		}
	}
	void unsetReadOnlyOptimization() {
		for(size_t i = 0; i<_buckets.size(); i++) {
			_buckets[i].unsetReadOnlyOptimization();
		}
	}

	inline int getLocalThreadId(NumberType hash, int numThreads) const {
		// stripe across all buckets
		assert(numThreads > 0);
		if (_buckets.size() > 0)
			return (hash & BUCKET_MASK) % numThreads;
		else
			return 0;
	}
	inline int getLocalThreadId(const KeyType &key, int numThreads) const {
		return getLocalThreadId(key.hash(), numThreads);
	}
	inline int getDistributedThreadId(NumberType hash, NumberType numDistributedThreads) const {
		// partition in contiguous blocks of 'global' buckets
		if (_buckets.size() > 0)
			return numDistributedThreads * (hash & BUCKET_MASK) / _buckets.size();
		else
			return 0;
	}
	inline int getDistributedThreadId(const KeyType &key, NumberType numDistributedThreads) const {
		return getDistributedThreadId(key.hash(), numDistributedThreads);
	}

	// optimized to look at both possible thread partitions
	// if this is the correct distributed thread, return true and set the proper localThread
	// otherwise return false
	inline bool getLocalThreadId(NumberType hash, int &localThreadId, int numLocalThreads, int distributedThreadId, NumberType numDistributedThreads) const {
		if (numDistributedThreads == 1 || getDistributedThreadId(hash, numDistributedThreads) == distributedThreadId) {
			localThreadId = getLocalThreadId(hash, numLocalThreads);
			return true;
		} else {
			return false;
		}
	}
	inline bool getLocalThreadId(const KeyType &key, int &localThreadId, int numLocalThreads, int distributedThreadId, NumberType numDistributedThreads) const {
		NumberType hash = key.hash();
		return getLocalThreadId(hash, localThreadId, numLocalThreads, distributedThreadId, numDistributedThreads);
	}

	inline void getThreadIds(const KeyType &key, int &localThreadId, int numLocalThreads, int &distributedThreadId, NumberType numDistributedThreads) const {
		NumberType hash = key.hash();
		distributedThreadId = getDistributedThreadId(hash, numDistributedThreads);
		localThreadId = getLocalThreadId(hash, numLocalThreads);
	}

	// optimization to move the buckets with pre-allocated memory to the next DMP thread
	void rotateDMPBuffers(int numThreads) {
		IndexType block = _buckets.size() / numThreads;
		size_t i = 0;
		// skip to the first non-zero bucket
		while (i < _buckets.size() && _buckets[i].size() == 0)
			i++;
		for(size_t j = i; j < _buckets.size() - 1 && j < i+block; j++) {
			_buckets[j].swap(_buckets[ (j+block) % _buckets.size() ]);
		}
	}

	inline NumberType getBucketIdx(NumberType hash) const {
		NumberType bucketIdx = (hash & BUCKET_MASK);
		return bucketIdx;
	}
	inline NumberType getBucketIdx(const KeyType &key) const {
		return getBucketIdx(key.hash());
	}

	inline const BucketType &getBucket(NumberType hash) const {
		NumberType idx = getBucketIdx(hash);
		return _buckets[idx];
	}
	inline BucketType &getBucket(NumberType hash) {
		return const_cast<BucketType &>( _constThis().getBucket(hash) );
	}

	inline const BucketType &getBucket(const KeyType &key) const {
		NumberType idx = getBucketIdx(key);
		return getBucket(idx);
	}
	inline BucketType &getBucket(const KeyType &key) {
		return const_cast<BucketType &>( _constThis().getBucket(key) );
	}
	inline const BucketType &getBucketByIdx(IndexType idx) const {
		return _buckets[idx];
	}

	inline void setNumBuckets(IndexType numBuckets) {
		_buckets.clear();
		_buckets.resize(numBuckets);
	}
	inline IndexType getNumBuckets() const {
		return _buckets.size();
	}

	ElementType insert(const KeyType &key, const ValueType &value, BucketType &bucket) {
		IndexType idx;
		if (isSorted()) {
			idx = bucket.insertSorted(key,value);
		} else {
			idx = bucket.append(key,value);
			if (bucket.size() == bucket.capacity()) {
				bucket.resort(true);
				idx = bucket.find(key);
				assert(idx != BucketType::MAX_INDEX);
			}
		}

		ElementType element = bucket.getElement(idx);
		return element;
	}
	ElementType insert(const KeyType &key, const ValueType &value) {
		return insert(key,value, getBucket(key));
	}

	bool remove(const KeyType &key, BucketType &bucket) {
		bool isFound;
		IndexType idx = isSorted() ? bucket.findSorted(key, isFound) : bucket.find(key, isFound);
		if (isFound && idx != BucketType::MAX_INDEX)
			bucket.remove(idx);
		return isFound;
	}
	bool remove(const KeyType &key) {
		return remove(key, getBucket(key));
	}

	bool _exists(const KeyType &key, IndexType &existingIdx, const BucketType &bucket) const {
		bool isFound;
		existingIdx = isSorted() ? bucket.findSorted(key, isFound) : bucket.find(key, isFound);
		return isFound;
	}
	bool _exists(const KeyType &key, IndexType &existingIdx) const {
		return _exists(key,existingIdx, getBucket(key));
	}

	bool exists(const KeyType &key, const BucketType &bucket) const {
		IndexType dummy;
		return _exists(key, dummy, bucket);
	}
	bool exists(const KeyType &key) const {
		return exists(key, getBucket(key));
	}

	const bool getValueIfExists(const KeyType &key, ValueType &value) const {
		return getValueIfExists(key, value, getBucket(key));
	}
	const bool getValueIfExists(const KeyType &key, ValueType &value, const BucketType &bucket) const {
		bool isFound;
		IndexType existingIndex;
		existingIndex = isSorted() ? bucket.findSorted(key, isFound) : bucket.find(key, isFound);
		if (isFound)
			value = bucket.valueAt(existingIndex);
		return isFound;
	}

	const ElementType getElementIfExists(const KeyType &key, const BucketType &bucket) const {
		IndexType existingIdx;
		ElementType element;
		if (_exists(key, existingIdx, bucket)) {
			element = bucket.getElement(existingIdx);
		}
		return element;
	}
	ElementType getElementIfExists(const KeyType &key, BucketType &bucket) {
		IndexType existingIdx;
		ElementType element;
		if (_exists(key, existingIdx, bucket)) {
			element = bucket.getElement(existingIdx);
		}
		return element;
	}

	const ElementType getElementIfExists(const KeyType &key) const {
		return getElementIfExists(key, getBucket(key));
	}
	ElementType getElementIfExists(const KeyType &key) {
		return getElementIfExists(key, getBucket(key));
	}

	ElementType getOrSetElement(const KeyType &key, BucketType &bucket, ValueType value) {
		IndexType existingIdx;
		ElementType element;
		if (_exists(key, existingIdx, bucket)) {
			element = bucket.getElement(existingIdx);
		} else {
			element = insert(key, value, bucket);
		}
		return element;
	}
	ElementType getOrSetElement(const KeyType &key, ValueType value) {
		return getOrSetElement(key, getBucket(key), value);
	}
	ElementType getElement(const KeyType &key, BucketType &bucket) {
		ValueType value = ValueType();
		return getOrSetElement(key, bucket, value);
	}
	ElementType getElement(const KeyType &key) {
		ValueType value = ValueType();
		return getOrSetElement(key, value);
	}

	// not thread safe!
	ValueType &operator[](const KeyType &key) {
		BucketType &bucket = getBucket(key);
		IndexType existingIdx;
		if (_exists(key, existingIdx, bucket))
			return bucket.valueAt(existingIdx);
		else
			return insert(key, Value(), bucket).value();
	}

	SizeType size() const {
		SizeType size = 0;

		long bucketSize = _buckets.size();
#pragma omp parallel for reduction(+:size) if(_buckets.size()>1000000)
		for(long i = 0; i < bucketSize; i++) {
			size += _buckets[i].size();
		}
		return size;
	}

	std::string toString() const {
		std::stringstream ss;
		ss << this << "[";
		IndexType idx=0;
		for(; idx<_buckets.size() && idx < 30; idx++) {
			ss << "bucket:" << idx << ' ' << _buckets[idx].toString() << ", ";
		}
		if (idx < _buckets.size())
			ss << " ... " << _buckets.size() - idx << " more ";
		ss << "]";
		return ss.str();
	}

	IndexType maxBucket()
	{
		IndexType biggest = 0;
		IndexType imax;
		for(IndexType i = 0; i<_buckets.size(); i++)
			if (_buckets[i].size() > biggest)
			{
				imax = i;
				biggest = _buckets[i].size();
			}
		return biggest;
	}

	SizeType purgeMinCount(long minimumCount) {
		SizeType affected = 0;
		if (minimumCount <= 1)
			return affected;

		long bucketsSize = getNumBuckets();
#pragma omp parallel for reduction(+:affected)
		for(long idx = 0 ; idx < bucketsSize; idx++) {
			affected += _buckets[idx].purgeMinCount(minimumCount);
		}
		_isSorted = true; // purgeMinCount implicitly sorts...
		return affected;
	}

private:
	int _compare(const BucketType &a, const BucketType &b, IndexType &idxA, IndexType &idxB) {
		int cmp;
		if (idxA < a.size() && idxB < b.size()) {
			cmp = a[idxA].compare(b[idxB]);
		} else {
			if (idxA < a.size()) {
				cmp = -1;
			} else {
				cmp = 1;
			}
		}
		return cmp;
	}
public:

	static bool _mergeTrivial(BucketType &a, BucketType &b) {
		if (b.size() == 0) {
			// do nothing.  a is good as is
			return true;
		} else if (a.size() == 0) {
			// swap a and b
			a.swap(b);
			return true;
		}
		return false;
	}

	// optimized merge for DMP threaded (blocked) KmerMaps
	void merge(const KmerMap &src) const {
		if (getNumBuckets() != src.getNumBuckets()) {
			LOG_THROW("Invalid: Can not merge two KmerMaps of differing sizes!");
		}

		long bucketsSize = getNumBuckets();
#pragma omp parallel for
		for(long idx = 0 ; idx < bucketsSize; idx++) {
			BucketType &a = const_cast<BucketType&>(_buckets[idx]);
			BucketType &b = const_cast<BucketType&>(src._buckets[idx]);

			if (! _mergeTrivial(a,b) ) {
				LOG_THROW("Invalid: Expected one bucket to be zero in this optimized method: KmerMap::merge(const KmerMap &src) const");
			}
		}

	}

	void mergeAdd(KmerMap &src) {
		if (getNumBuckets() != src.getNumBuckets()) {
			LOG_THROW("Invalid: Can not merge two KmerMaps of differing sizes!");
		}
		BucketType merged;
		if (!isSorted())
			resort();
		if (!src.isSorted())
			src.resort();

		long bucketsSize = getNumBuckets();
#pragma omp parallel for private(merged)
		for(long idx = 0 ; idx < bucketsSize; idx++) {
			// buckets are sorted so perform a sorted merge by bucket
			BucketType &a = _buckets[idx];
			BucketType &b = src._buckets[idx];

			// short circuit if one of the buckets is empty
			if (_mergeTrivial(a,b))
				continue;

			IndexType capacity = std::max(a.size(), b.size());

			merged.reserve(capacity);
			IndexType idxA = 0;
			IndexType idxB = 0;
			while (idxA < a.size() || idxB < b.size()) {
				int cmp = _compare(a, b, idxA, idxB);
				if (cmp == 0) {
					merged.append(a[idxA], a.valueAt(idxA).add(b.valueAt(idxB)));
					idxA++; idxB++;
				} else {
					if (cmp < 0) {
						merged.append(a[idxA], a.valueAt(idxA));
						idxA++;
					} else {
						merged.append(b[idxB], b.valueAt(idxB));
						idxB++;
					}
				}
			}

			a.swap(merged);
			if (b.capacity() > merged.capacity()) {
				b.swap(merged);
			}
			merged.reset(false);
			b.reset(true);
		}
		src.clear();
	}

	template< typename OtherDataType >
	void mergePromote(KmerMap &src, KmerMap< OtherDataType > &mergeDest) {
		// TODO fixme
		LOG_THROW("Invalid: This method is broken for threaded execution somehow (even with pragmas disabled)");
		if (getNumBuckets() != src.getNumBuckets()) {
			LOG_THROW("Invalid: Can not merge two KmerMaps of differing sizes!");
		}
		// calculated chunk size to ensure thread safety of merge operation
		int chunkSize = mergeDest.getNumBuckets() / src.getNumBuckets();
		if (chunkSize <= 1)
			chunkSize = 2;
		BucketType merged;
		if (! src.isSorted() )
			src.resort();
		if (! mergeDest.isSorted() )
			mergeDest.resort();

		long bucketsSize = getNumBuckets();
		//#pragma omp parallel for private(merged) schedule(static, chunkSize)
		for(long idx = 0 ; idx < bucketsSize; idx++) {
			// buckets are sorted so perform a sorted merge by bucket
			BucketType &a = _buckets[idx];
			BucketType &b = src._buckets[idx];
			IndexType capacity = std::max(a.size(), b.size());
			merged.reserve(capacity);
			IndexType idxA = 0;
			IndexType idxB = 0;
			while (idxA < a.size() || idxB < b.size()) {
				int cmp = _compare(a, b, idxA, idxB);
				if (cmp == 0) {
					mergeDest[ a[idxA] ].add( a.valueAt(idxA) ).add( b.valueAt(idxB) );
					idxA++; idxB++;
				} else {
					if (cmp < 0) {
						merged.append(a[idxA], a.valueAt(idxA));
						idxA++;
					} else {
						merged.append(b[idxB], b.valueAt(idxB));
						idxB++;
					}
				}
			}

			a.swap(merged);
			if (b.capacity() > merged.capacity()) {
				b.swap(merged);
			}
			merged.reset(false);
			b.reset(true);
		}
		src.clear();
	}

public:
	class Iterator : public std::iterator<std::forward_iterator_tag, KmerMap>
	{
		friend class KmerMap;
	private:
		const KmerMap *_target;
		BucketsVectorIterator _iBucket;
		BucketTypeIterator _iElement;
		int _rank, _size;

	public:
		Iterator() : _target(NULL), _rank(0), _size(0) {}

		// iterator over rank/size will stripe across the buckets (modulus by size).
		Iterator(KmerMap *target, int rank = 0, int size = 1):
			_target(target),
			_iBucket(target->_buckets.begin()),
			_iElement(),
			_rank(rank), _size(size)
		{
			if (!isEnd())
				_iElement = _iBucket->begin();
			_moveToRank();
			_moveToNextValidElement();
		}

		Iterator(const Iterator &copy) {
			_target = copy._target;
			_iBucket = copy._iBucket;
			_iElement = copy._iElement;
			_rank = copy._rank;
			_size = copy._size;
		}

		Iterator(KmerMap *target, BucketsVectorIterator bucketPtr, int rank = 0, int size = 1):
			_target(target),
			_iBucket(bucketPtr),
			_iElement(),
			_rank(rank), _size(size)
		{
			if (!isEnd())
				_iElement = _iBucket->begin();
			_moveToRank();
			_moveToNextValidElement();
		}

		Iterator(KmerMap *target, BucketsVectorIterator bucketPtr, BucketTypeIterator elementPtr, int rank = 0, int size = 1):
			_target(target),
			_iBucket(bucketPtr),
			_iElement(elementPtr),
			_rank(rank), _size(size)
		{
			_moveToRank();
		}

	private:
		inline bool isEnd() const {
			return _iBucket == _target->_buckets.end();
		}

		void _moveToNextValidElement() {
			while (! isEnd() ) {
				if (_iElement == _iBucket->end()) {
					for(int i = 0 ; i < _size; i++) {
						++_iBucket;
						if (isEnd())
							break;
					}
					if (isEnd())
						break;
					_iElement = _iBucket->begin();
				} else
					break;
			}
		}
		void _moveToRank() {
			if (_size == 1)
				return;
			if (isEnd() || ((int) _target->_buckets.size() <= _rank)) {
				while (!isEnd())
					++_iBucket;
				return;
			}
			for(int i = 0 ; i < _rank ; i++) {
				++_iBucket;
				if (isEnd())
					return;
			}
			_iElement = _iBucket->begin();
		}

	public:

		bool operator==(const Iterator& other) const
		{	return _iBucket == other._iBucket && (isEnd() || _iElement == other._iElement);}

		bool operator!=(const Iterator& other) const
								{	return !(*this == other);}

		Iterator& operator++()
								{
			++_iElement;
			_moveToNextValidElement();
			return *this;
								}

		Iterator operator++(int unused)
								{	Iterator tmp(*this); ++(*this); return tmp;}

		ElementType &operator*() {return *_iElement;}
		const ElementType &operator*() const {return *_iElement;}
		ElementType *operator->() {return &(*_iElement);}
		const ElementType *operator->() const {return &(*_iElement);}

		Kmer &key() {return _iElement.key();}
		const Kmer &key() const {return _iElement.key();}
		Value &value() {return _iElement.value();}
		const Value &value() const {return _iElement.value();}

		BucketType &bucket() {return *_iBucket;}
		const BucketType &bucket() const {return *_iBucket;}

		IndexType bucketIndex() {return (_iBucket - _target->_buckets.begin());}
		const IndexType bucketIndex() const {return (_iBucket - _target->_buckets.begin());}

	};

	class ConstIterator : public std::iterator<std::forward_iterator_tag, KmerMap>
	{
	private:
		Iterator _iterator;
	public:
		ConstIterator(KmerMap *target, int rank = 0, int size = 1) : _iterator(target, rank, size) {}
		ConstIterator(const ConstIterator &copy) : _iterator(copy._iterator) {}
		ConstIterator(KmerMap *target, BucketsVectorIterator bucketPtr, int rank = 0, int size = 1) : _iterator(target,bucketPtr, rank, size) {}
		ConstIterator(KmerMap *target, BucketsVectorIterator bucketPtr,BucketTypeIterator elementPtr, int rank = 0, int size = 1): _iterator(target,bucketPtr,elementPtr, rank, size) {}
		bool operator==(const ConstIterator& other) const {return _iterator == other._iterator;}
		bool operator!=(const ConstIterator& other) const {return _iterator != other._iterator;}
		ConstIterator& operator++() {++_iterator; return *this;}
		ConstIterator operator++(int unused) {return ConstIterator(_iterator++);}
		const ElementType &operator*() const {return _iterator.operator*();}
		const ElementType *operator->() const {return _iterator.operator->();}
		const Kmer &key() const {return _iterator.key();}
		const Value &value() const {return _iterator.value();}
		const BucketType &bucket() const {return _iterator.bucket();}
		const IndexType bucketIndex() const {return _iterator.bucketIndex();}
	};

	Iterator begin(int rank = 0, int size = 1) {return Iterator(this, rank, size);}
	Iterator end() {return Iterator(this,_buckets.end());}

	ConstIterator begin(int rank = 0, int size = 1) const {return Iterator(this, rank, size);}
	ConstIterator end() const {return Iterator(this, _buckets.end());}

	Iterator beginThreaded(int rank = omp_get_thread_num(), int size = omp_get_num_threads()) {
		return begin(rank, size);
	}
	Iterator endThreaded() {
		return end();
	}

};

typedef KmerArray<char> Kmers;
typedef KmerArray<double> KmerWeights;
typedef KmerArray<Kmernator::UI32> KmerCounts;
typedef KmerArray< WeightedExtensionMessagePacket > KmerWeightedExtensions;

#endif

