//
// Kmernator/src/Kmer.h
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
#include "TwoBitSequence.h"
#include "MemoryUtils.h"
#include "KmerTrackingData.h"
#include "lookup8.h"
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
		_sequenceLength(21) {
		_set(_sequenceLength);
	}

	void _set(SequenceLengthType sequenceLength) {
		_verifyThreads();
		_sequenceLength = sequenceLength;
		_twoBitLength = TwoBitSequence::fastaLengthToTwoBitLength(_sequenceLength);
		_totalSize = _twoBitLength;
	}

	static void _verifyThreads() {
		int maxThreads = omp_get_max_threads();
		if ((maxThreads & (maxThreads-1)) != 0) {
			throw std::invalid_argument(
					(std::string("The maximum number of threads must be a power-of-two.\nPlease adjust the OMP_NUM_THREADS environmental variable. ")
					 + boost::lexical_cast<std::string>(maxThreads)).c_str());
		}
	}


public:

	static inline KmerSizer &getSingleton() {
		return singleton;
	}

	static void set(SequenceLengthType sequenceLength) {
		KmerSizer & singleton = getSingleton();
		singleton._set(sequenceLength);
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
			val
					= (((NumberType) *((boost::uint16_t *) kmer.getTwoBitSequence())))
							| (((NumberType) *((boost::uint8_t *) (kmer.getTwoBitSequence()
									+ 2))) << 16);
			break;
		case 4:
			val = (NumberType) *((boost::uint32_t *) kmer.getTwoBitSequence());
			break;
		case 5:
			val
					= (((NumberType) *((boost::uint32_t *) kmer.getTwoBitSequence())))
							| (((NumberType) *((boost::uint8_t *) (kmer.getTwoBitSequence()
									+ 4))) << 32);
			break;
		case 6:
			val
					= (((NumberType) *((boost::uint32_t *) kmer.getTwoBitSequence())))
							| (((NumberType) *((boost::uint16_t *) (kmer.getTwoBitSequence()
									+ 4))) << 32);
			break;
		case 7:
			val
					= (((NumberType) *((boost::uint32_t *) kmer.getTwoBitSequence())))
							| (((NumberType) *((boost::uint16_t *) (kmer.getTwoBitSequence()
									+ 4))) << 32)
							| (((NumberType) *((boost::uint8_t *) (kmer.getTwoBitSequence()
									+ 6))) << 48);
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
		return TwoBitSequence::getFasta(getTwoBitSequence(), getLength());
	}
	std::string toFastaFull() const {
		return TwoBitSequence::getFasta(getTwoBitSequence(), getTwoBitLength()
				* 4);
	}
	inline NumberType toNumber() const {
		return Kmer::toNumber(*this);
	}


	inline NumberType hash() const {

		NumberType number = toNumber();
		return Lookup8::hash2(&number, 1, 0xDEADBEEF);

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
			setLock();
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
			setLock();
			return *this;
		}
		void reset() {
			unsetLock();
			_idx = MAX_INDEX;
			_array = NULL;
		}
		inline bool isValid() const {
			return _array != NULL && _idx != MAX_INDEX;
		}
		void setLock() {
			if (isValid())
				_array->setSharedLock();
		}
		void unsetLock() {
			if (isValid())
				_array->unsetSharedLock();
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

	// if capacity == MAX_INDEX, this is a constant mmaped instance
	inline bool isMmaped() const {
		return _capacity == MAX_INDEX;
	}

#ifdef _USE_THREADSAFE_KMER

	typedef std::vector< omp_nest_lock_t > LockVector;
	mutable omp_nest_lock_t _lock;
	mutable LockVector _sharedLocks;

private:

	// thread locking routines
	void initLock() const {omp_init_nest_lock(&_lock); initSharedLocks();}
	void destroyLock() const {omp_destroy_nest_lock(&_lock); destroySharedLocks();}
	void initSharedLocks() const {
		_sharedLocks.clear();
		_sharedLocks.resize( omp_get_num_procs() );
		for(size_t i=0; i < _sharedLocks.size(); i++)
		omp_init_nest_lock( &(_sharedLocks[i]) );
	}
	void destroySharedLocks() const {
		for(size_t i=0; i < _sharedLocks.size(); i++)
		omp_destroy_nest_lock( &(_sharedLocks[i]) );
		_sharedLocks.clear();
	}

	inline void setLock() const {omp_set_nest_lock(&_lock);}
	inline bool testLock() const {return omp_test_nest_lock(&_lock);}
	inline void unsetLock() const {omp_unset_nest_lock(&_lock);}

public:
	// TODO make a true shared lock without a severe performance penalty

	void setSharedLock() const {
		if (_sharedLocks.empty())
		return;
		unsigned int myThread = omp_get_thread_num();
		bool gotShared = false;
		while (!gotShared) {

			bool gotExclusive = false;
			if (testLock()) {
				gotShared = omp_test_nest_lock( &( _sharedLocks[myThread] ) );
				unsetLock();
				gotExclusive = true;
			}
			if (! gotShared ) {

				usleep(1);

			}
		}

	}

	inline void unsetSharedLock() const {
		if (_sharedLocks.empty())
		return;
		unsigned int myThread = omp_get_thread_num();
		omp_unset_nest_lock( &( _sharedLocks[ myThread ] ) );
	}

	inline void setExclusiveLock() const {
		if (_sharedLocks.empty())
		return;
		bool gotExclusive = false;
		while(!gotExclusive) {

			gotExclusive = true;
			if (testLock()) {
				size_t i = 0;
				for(; i < _sharedLocks.size(); i++) {
					if (! omp_test_nest_lock( &( _sharedLocks[i] ) )) {
						gotExclusive = false;
						break;
					}
				}
				if (! gotExclusive) {
					// unset acquired locks
					unsetLock();
					for(size_t j = 0; j<i; j++)
					omp_unset_nest_lock( &( _sharedLocks[j] ) );
					usleep(1);
				}
			} else {
				gotExclusive = false;
			}

		}
	}
	inline void unsetExclusiveLock() const {
		if (_sharedLocks.empty())
		return;

		unsetLock();
		// unset acquired shared locks
		for(size_t j = 0; j < _sharedLocks.size(); j++)
		omp_unset_nest_lock( &( _sharedLocks[j] ) );

	}

#else // _USE_THREADSAFE_KMER
private:
	void initLock() const {
	}
	void destroyLock() const {
	}
	void initSharedLocks() const {
	}
	void destroySharedLocks() const {
	}
public:
	inline void setSharedLock() const {
	}
	inline void unsetSharedLock() const {
	}
	inline void setExclusiveLock() const {
	}
	inline void unsetExclusiveLock() const {
	}

#endif //  _USE_THREADSAFE_KMER
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

public:
	void setReadOnlyOptimization() const {
		destroyLock();
	}
	void unsetReadOnlyOptimization() {
		initLock();
	}

public:

	KmerArray(IndexType size = 0) :
		_begin(NULL), _size(0), _capacity(0) {
		initLock();
		resize(size);
	}

private:
	KmerArray(void *begin, IndexType size, IndexType capacity) : _begin(begin), _size(size), _capacity(capacity) {}

public:
	KmerArray(const TwoBitEncoding *twoBit, SequenceLengthType length, bool leastComplement = false, bool *bools = NULL) :
		_begin(NULL), _size(0), _capacity(0) {
		initLock();
		if (KmerSizer::getSequenceLength() <= length) {
			SequenceLengthType numKmers = length
					- KmerSizer::getSequenceLength() + 1;
			resize(numKmers, MAX_INDEX, false);
			build(twoBit, length, leastComplement, bools);
		} else {
			resize(0);
		}
	}

	KmerArray(const KmerArray &copy) :
		_begin(NULL), _size(0), _capacity(0) {
		initLock();
		*this = copy;
	}

	~KmerArray() {
		reset();
		destroyLock();
	}

	KmerArray &operator=(const KmerArray &other) {
		if (this == &other)
			return *this;
		reset();
		if (other.isMmaped()) {
			_begin = other._begin;
			_size = other._size;
			_capacity = other._capacity;
			return *this;
		}

		setExclusiveLock();
		resize(other.size());
		if (size() == 0) {
			unsetExclusiveLock();
			return *this;
		}
		if (_begin == NULL)
			throw std::runtime_error(
					"Could not allocate memory in KmerArray operator=()");

		_copyRange(other._begin, other.getValueStart(), 0, 0, _size, false);

		unsetExclusiveLock();
		return *this;
	}

	// restore a new array from a mmap
	KmerArray(const void *src) : _begin(NULL), _size(0), _capacity(0) {
		initLock();
		IndexType *size = (IndexType *) src;
	    resize(*size);
		void *ptr = ++size;
		if (_size > 0) {
			long kmerSize  = _size * KmerSizer::getByteSize();
			_copyRange(ptr, (ValueType *) (((char*)ptr)+kmerSize), 0, 0, _size, false);
		}
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
	static const KmerArray restore(const void *src) {
		const IndexType *size = (const IndexType *) src;
		IndexType xsize = *(size++);
		const void *ptr = size;
		KmerArray array(const_cast<void*>(ptr), xsize, MAX_INDEX);
		return array;
	}

protected:
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
			throw std::invalid_argument(
					"attempt to access index greater than size in KmerArray operator[] const");
		return get(index);
	}

	// never thread safe!
	Kmer &operator[](IndexType index) {
		if (index >= size())
			throw std::invalid_argument(
					"attempt to access index greater than size in KmerArray operator[]");
		return get(index);
	}

	// never thread safe!
	const ValueType &valueAt(IndexType index) const {
		if (index >= _size) {
			throw std::invalid_argument(
					"attempt to access index greater than size in KmerArray valueAt() const");
		}
		return *(getValueStart() + index);
	}

	// never thread safe!
	ValueType &valueAt(IndexType index) {
		if (index >= _size) {
			throw std::invalid_argument(
					"attempt to access index greater than size in KmerArray valueAt()");
		}
		return *(getValueStart() + index);
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
			return;
		}

		setExclusiveLock();

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

		unsetExclusiveLock();
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
	void reserve(IndexType size) {
		assert(!isMmaped()); // mmaped can not be modified!

		if (size < _capacity) {
		  IndexType oldSize = _size;
		  resize(size, MAX_INDEX, false);
		  _size = oldSize;
		}
	}

	void resize(IndexType size) {
		resize(size, MAX_INDEX, false);
	}
	void resize(IndexType size, IndexType idx, bool reserveExtra = true) {
		assert(!isMmaped()); // mmaped can not be modified!

		if (size == _size)
			return;
		IndexType oldSize = _size;

		setExclusiveLock();
		// alloc / realloc memory
		_setMemory(size, idx, reserveExtra);

		if (_begin == NULL && size > 0) {
			throw std::runtime_error(
					"Could not allocate memory in KmerArray resize()");
		}

		if (size > oldSize && idx == MAX_INDEX) {
			// zero fill remainder
			void *start = _add(_begin, oldSize);
			memset(start, 0, KmerSizer::getByteSize() * (size - oldSize));

			// Values should already have been constructed
		}
		unsetExclusiveLock();
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

		// preserve old pointers and metrics
		void *oldBegin = _begin;
		IndexType oldSize = _size;
		IndexType oldCapacity = _capacity;
		IndexType lesserSize = std::min(size, oldSize);
		ValueType *oldValueStart = NULL;
		if (oldCapacity > 0) {
			oldValueStart = getValueStart();
		}

		bool memChanged = false;
		if (size == 0) {
			// ignore reserveExtra for zero resize
			if (_capacity > 0) {
				_capacity = 0;
				_begin = NULL;
				memChanged = true;
			}
		} else if ((size > _capacity) || (_capacity > size && !reserveExtra)) {
			// allocate new memory
			if (reserveExtra)
				_capacity = std::max((IndexType) (size * 1.2),
						(IndexType) (size + 10));
			else
				_capacity = size;
			void *memory = std::malloc(_capacity * getElementByteSize());
			if (memory == NULL) {
				LOG_ERROR(0, "Attempt to malloc " << _capacity
						* getElementByteSize());
				throw std::runtime_error(
						"Could not allocate memory in KmerArray _setMemory()");
			}
			_begin = memory;
			memChanged = true;
		}

		_size = size;

		if (memChanged && _begin != NULL) {
			// construct all values that are newly allocated
			for (IndexType i = 0; i < _capacity; i++)
				::new ((void*) (getValueStart() + i)) Value();
		}

		if (_begin != NULL && oldBegin != NULL && oldValueStart != NULL
				&& lesserSize > 0) {
			// copy the old contents
			if (idx == MAX_INDEX || idx >= lesserSize) {
				if (memChanged) {
					// copy all records in order (default ; end is trimmed or expanded)
					_copyRange(oldBegin, oldValueStart, 0, 0, lesserSize);
				} else {
					// noop
				}
			} else {
				// shrink or expand.

				if (memChanged && idx > 0) {
					// the first record(s) leading to idx will be copied
					_copyRange(oldBegin, oldValueStart, 0, 0, idx);
				}

				if (lesserSize == size) {
					// shrink: skipping the old record at idx
					if (idx < size) {
						_copyRange(oldBegin, oldValueStart, idx, idx + 1,
								lesserSize - idx, !memChanged);
					}
				} else {
					// expand: leaving new (uninitialized/unchanged) record at idx
					if (idx < oldSize) {
						_copyRange(oldBegin, oldValueStart, idx + 1, idx,
								lesserSize - idx, !memChanged);
					}
				}
			}
		}

		if (memChanged && oldBegin != NULL) {
			// destruct old Values
			for (IndexType i = 0; i < oldCapacity; i++)
				(oldValueStart + i)->~Value();
			// free old memory
			std::free(oldBegin);
		}
	}

	void build(const TwoBitEncoding *twoBit, SequenceLengthType length,
			bool leastComplement = false, bool *bools = NULL) {
		assert(!isMmaped()); // mmaped can not be modified!
		setExclusiveLock();
		SequenceLengthType numKmers = length - KmerSizer::getSequenceLength()
				+ 1;
		if (_size != numKmers)
			throw std::invalid_argument(
					"attempt to build an incorrectly sized KmerArray in KmerArray build()");
		;

		KmerArray &kmers = *this;
		long numBytes = (numKmers + 3) / 4;

		#pragma omp parallel for if(numKmers >= 10000)
		for (long bytes = 0; bytes < numBytes; bytes++) {
			SequenceLengthType i = bytes * 4;
			const TwoBitEncoding *ref = twoBit + i / 4;
			for (int bitShift = 0; bitShift < 4 && i + bitShift < numKmers; bitShift++) {
				TwoBitSequence::shiftLeft(ref, kmers[i + bitShift].get(),
						KmerSizer::getTwoBitLength(), bitShift, bitShift != 0);
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
		unsetExclusiveLock();
	}

	// return a KmerArray that has one entry for each possible single-base substitution
	static KmerArray permuteBases(const Kmer &kmer, bool leastComplement = false) {
		KmerArray kmers(KmerSizer::getSequenceLength() * 3);
		TEMP_KMER(tmp);
		for (SequenceLengthType baseIdx = 0; baseIdx < KmerSizer::getSequenceLength(); baseIdx++) {
			TwoBitSequence::permuteBase(kmer.getTwoBitSequence(), kmers[baseIdx*3].getTwoBitSequence(), kmers[baseIdx*3+1].getTwoBitSequence(), kmers[baseIdx*3+2].getTwoBitSequence(),
					KmerSizer::getSequenceLength(), baseIdx);
			if (leastComplement) {
				for(int i = 0 ; i < 3; i++) {
				   if (! kmers[baseIdx*3+i].buildLeastComplement(tmp) ) {
				      kmers[baseIdx*3+i] = tmp;
				   }
				}
			}
		}
		return kmers;
	}
	static KmerArray permuteBases(const Kmer &kmer, const ValueType defaultValue, bool leastComplement = false) {
		KmerArray kmers = permuteBases(kmer, leastComplement);
		for(SequenceLengthType idx = 0; idx < kmers.size(); idx++)
			kmers.valueAt(idx) = defaultValue;
		return kmers;
	}

		protected:
			IndexType _find(const Kmer &target) const {
				for(IndexType i=0; i<_size; i++)
				if (target.compare(get(i)) == 0) {
					return i;
				}
				return MAX_INDEX;
			}
		public:
			IndexType find(const Kmer &target) const {
				setSharedLock();
				IndexType idx = _find(target);
				unsetSharedLock();
				return idx;
			}

		protected:
			IndexType _findSorted(const Kmer &target, bool &targetIsFound) const {
				// binary search
				IndexType min = 0;
				IndexType max = size();

				if (max == 0)
				{
					targetIsFound = false;
					return 0;
				}
				max--; // never let mid == size()
				IndexType mid;
				int comp;
				do {
					mid = (min+max) / 2;
					comp = target.compare(get(mid));
					if (comp > 0)
					min = mid+1;
					else if (comp < 0)
					max = mid-1;
				}while (comp != 0 && max != MAX_INDEX && min <= max);
				if (comp == 0)
				targetIsFound = true;
				else
				targetIsFound = false;

				return mid + (comp>0 && size()>mid?1:0);
			}

		public:
			IndexType findSorted(const Kmer &target, bool &targetIsFound) const {
				setSharedLock();
				IndexType idx = _findSorted(target, targetIsFound);
				unsetSharedLock();
				return idx;
			}

		protected:
			void _insertAt(IndexType idx, const Kmer &target) {
				assert(!isMmaped()); // mmaped can not be modified!
                if (idx > size())
				throw std::invalid_argument("attempt to access index greater than size in KmerArray insertAt");
				resize(size() + 1, idx);
				get(idx) = target;
			}
			void _insertAt(IndexType idx, const Kmer &target, const Value &value) {
				_insertAt(idx, target);
				valueAt(idx) = value;
			}

		public:
			void insertAt(IndexType idx, const Kmer &target) {
				setExclusiveLock();
				_insertAt(idx, target);
				unsetExclusiveLock();
			}
			void insertAt(IndexType idx, const Kmer &target, const Value &value) {
				setExclusiveLock();
				_insertAt(idx,target,value);
				unsetExclusiveLock();
			}

			IndexType append(const Kmer &target) {
				setExclusiveLock();
				IndexType idx = size();
				_insertAt(idx, target);
				unsetExclusiveLock();
				return idx;
			}
			IndexType append(const Kmer &target, const Value &value) {
				setExclusiveLock();
				IndexType idx = size();
				_insertAt(idx, target, value);
				unsetExclusiveLock();
				return idx;
			}

		protected:
			IndexType _insertSorted(const Kmer &target) {
				bool isFound;
				IndexType idx = findSorted(target, isFound);
				if (!isFound)
				_insertAt(idx, target);
				return idx;
			}
			IndexType _insertSorted(const Kmer &target, const Value &value) {
				IndexType idx = _insertSorted(target);
				valueAt(idx) = value;
				return idx;
			}

		public:
			IndexType insertSorted(const Kmer &target) {
				setExclusiveLock();
				IndexType idx = _insertSorted(target);
				unsetExclusiveLock();
				return idx;
			}
			IndexType insertSorted(const Kmer &target, const Value &value) {
				setExclusiveLock();
				IndexType idx = _insertSorted(target,value);
				unsetExclusiveLock();
				return idx;
			}

			void remove(const Kmer &target) {
				setExclusiveLock();
				bool isFound;
				IndexType idx = find(target, isFound);
				if (isFound)
				remove(idx);
				unsetExclusiveLock();
			}
			void remove(IndexType idx) {
				assert(!isMmaped()); // mmaped can not be modified!
				setExclusiveLock();
				resize(size()-1,idx);
				unsetExclusiveLock();
			}

			void swap(IndexType idx1, IndexType idx2) {
				assert(!isMmaped()); // mmaped can not be modified!
				if (idx1 == idx2)
				return;
				if (idx1 >= size() || idx2 >= size())
				throw std::invalid_argument("attempt to access index greater than size in KmerArray swap()");

				setExclusiveLock();
				get(idx1).swap(get(idx2));
				if (sizeof(ValueType) > 0) {
					ValueType tmp = valueAt(idx1);
					valueAt(idx1) = valueAt(idx2);
					valueAt(idx2) = tmp;
				}
				unsetExclusiveLock();
			}
			void swap(KmerArray &other) {
				std::swap(_begin, other._begin);
				std::swap(_size, other._size);
				std::swap(_capacity, other._capacity);
			}

			// purge all element where the value is less than minimumCount
			// assumes that ValueType can be cast into long
			IndexType purgeMinCount(long minimumCount) {
				setExclusiveLock();
				// scan values that pass, keep list and count
				IndexType maxSize = size();
				IndexType passing[maxSize];
				IndexType passed = 0;
				IndexType affected;

				ValueType *valuePtr = getValueStart();
				for(IndexType i = 0; i < maxSize; i++) {
					if ( minimumCount <= (long) *(valuePtr++) ) {
						passing[passed++] = i;
					}
				}
				if (passed == 0) {
					reset(true);
					affected = maxSize;
				} else if (passed == maxSize) {
					affected = 0;
				} else {
					affected = maxSize - passed;
					KmerArray tmp;
					tmp.reserve(passed);
					for(IndexType i = 0 ; i < passed ; i++) {
						ElementType elem = getElement(passing[i]);
						tmp.append( elem.key(), elem.value() );
					}
					swap(tmp);
				}

				unsetExclusiveLock();
				return affected;
			}

			std::string toString() {
				setSharedLock();
				std::stringstream ss;
				ss << "{";
				IndexType idx=0;
				for(idx=0; idx<size() && idx < 30; idx++) {
					ss << get(idx).toFasta() << ":" << valueAt(idx) << ", ";
				}
				if (idx < size())
				ss << " ... " << size() - idx << " more ";
				ss << "}";
				unsetSharedLock();
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
				Iterator& operator++() {if ( !isEnd() ) {++_idx; setElement();}return *this;}
				Iterator operator++(int unused) {Iterator tmp(*this); ++(*this); return tmp;}
				ElementType &operator*() {return thisElement;}
				Kmer &key() {return thisElement.key();}
				Value &value() {return thisElement.value();}
				ElementType *operator->() {return &thisElement;}
				bool isEnd() const {return _idx >= _tgt->size();}
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

private:
	BucketsVector _buckets;
	NumberType BUCKET_MASK;

	inline const KmerMap &_constThis() const {return *this;}

public:
	KmerMap(IndexType bucketCount = 1024) {

		// ensure buckets are a precicise power of two
		// with at least bucketCount buckets
		NumberType powerOf2 = bucketCount;
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

		BUCKET_MASK = powerOf2 - 1;
		_buckets.resize(powerOf2);
	}
	~KmerMap()
	{
		clear();
	}

	// restore new instance from mmap
	KmerMap(const void *src) {
		NumberType size(0), *offsetArray;
		const void *ptr;
		_getMmapSizes(src, size, BUCKET_MASK, offsetArray);
		_buckets.resize(size);
		for (NumberType idx = 0 ; idx < size ; idx++) {
			ptr = ((char*)src) + *(offsetArray + idx);
			_buckets[idx] = BucketType(ptr);
		}
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
		for (NumberType idx = 0 ; idx < size ; idx++) {
			ptr = ((char*)src) + *(offsetArray + idx);
			map._buckets[idx] = BucketType::restore(ptr);
		}
		return map;
	}

	void swap(KmerMap &other) {
		_buckets.swap(other._buckets);
		NumberType tmp = other.BUCKET_MASK;
		other.BUCKET_MASK = BUCKET_MASK;
		BUCKET_MASK = tmp;
	}

	static const void _getMmapSizes(const void *src, NumberType &size, NumberType &mask, NumberTypePtr &offsetArray) {
		NumberType *numbers = (NumberType *) src;
		size = *(numbers++);
		mask = *(numbers++);
		offsetArray = numbers;
	}
	SizeType getSizeToStore() const {
		return sizeof(NumberType)*(2+_buckets.size())
		    + _buckets.size()*sizeof(IndexType)
		    + size()*BucketType::getElementByteSize();
	}

	void reset(bool releaseMemory = true) {
		for(size_t i=0; i< _buckets.size(); i++) {
			_buckets[i].reset(releaseMemory);
		}
	}
	void clear(bool releaseMemory = true) {
		reset(releaseMemory);
		if (releaseMemory)
		    _buckets.resize(1); // iterators require at least 1
	}

	void setReadOnlyOptimization() {
		for(size_t i = 0; i<_buckets.size(); i++) {
			_buckets[i].setReadOnlyOptimization();
		}
	}
	void unsetReadOnlyOptimization() {
		for(size_t i = 0; i<_buckets.size(); i++) {
			_buckets[i].unsetReadOnlyOptimization();
		}
	}

	inline unsigned short getLocalThreadId(NumberType hash, unsigned short numThreads) const {
		// use the bottom bits of hash (which are used to sort by bucket)
		// partition by numThreads blocks
		// splits into numThread contiguous blocks
		return (hash & BUCKET_MASK) / (_buckets.size() / numThreads + 1);
	}
	inline unsigned short getLocalThreadId(const KeyType &key, unsigned short numThreads) const {
		return getLocalThreadId(key.hash(), numThreads);
	}
	inline unsigned short getDistributedThreadId(NumberType hash, NumberType threadBitMask) const {
		assert( threadBitMask != 0 && ((threadBitMask+1) & threadBitMask) == 0); // number of threads must be a power of 2

		// stripe within the contiguous blocks for the local thread
		return (hash & threadBitMask);

		// use top bits of hash (unused to sort by bucket)
		//return (hash >> 32 & BUCKET_MASK) & threadBitMask;
	}
	inline unsigned short getDistributedThreadId(const KeyType &key, NumberType threadBitMask) const {
		return getDistributedThreadId(key.hash(), threadBitMask);
	}

	// optimized to look at both possible thead partitions
	// if this is the correct distributed thread, return true and set the proper localThread
	// otherwise return false
	inline bool getLocalThreadId(const KeyType &key, unsigned short &localThreadId, unsigned short numLocalThreads, unsigned short distributedThreadId, NumberType distributedThreadBitMask) const {
		NumberType hash = key.hash();
		if (distributedThreadBitMask == 0 || getDistributedThreadId(hash, distributedThreadBitMask) == distributedThreadId) {
			localThreadId = getLocalThreadId(hash, numLocalThreads);
			return true;
		} else {
			return false;
		}
	}

	// optimization to move the buckets with pre-allocated memory to the next DMP thread
	void rotateDMPBuffers(unsigned short numThreads) {
		BucketType tmp;
		tmp.swap(_buckets[ 0 ]);
		for(size_t i = 0 ; i < _buckets.size() - 1; i++) {
			_buckets[i].swap(_buckets[i+1]);
		}
		tmp.swap(_buckets[ _buckets.size() - 1]);
	}


	inline NumberType getBucketIdx(NumberType hash) const {
		return hash & BUCKET_MASK;
	}
	inline NumberType getBucketIdx(const KeyType &key) const {
		return getBucketIdx(key.hash());
	}

	inline const BucketType &getBucket(NumberType hash) const {
		return _buckets[getBucketIdx(hash)];
	}
	inline BucketType &getBucket(NumberType hash) {
		return const_cast<BucketType &>( _constThis().getBucket(hash) );
	}

	inline const BucketType &getBucket(const KeyType &key) const {
		return getBucket(getBucketIdx(key));
	}
	inline BucketType &getBucket(const KeyType &key) {
		return const_cast<BucketType &>( _constThis().getBucket(key) );
	}

	inline void setNumBuckets(IndexType numBuckets) {
		_buckets.clear();
		_buckets.resize(numBuckets);
	}
	inline IndexType getNumBuckets() const {
		return _buckets.size();
	}

	ElementType insert(const KeyType &key, const ValueType &value, BucketType &bucket) {
		bucket.setExclusiveLock();
		IndexType idx = bucket.insertSorted(key,value);
		ElementType element = bucket.getElement(idx);
		bucket.unsetExclusiveLock();
		return element;
	}
	ElementType insert(const KeyType &key, const ValueType &value) {
		return insert(key,value, getBucket(key));
	}

	bool remove(const KeyType &key, BucketType &bucket) {
		bool isFound;
		bucket.setExclusiveLock();
		IndexType idx = bucket.findSorted(key, isFound);
		if (isFound && idx != BucketType::MAX_INDEX)
		bucket.remove(idx);
		bucket.unsetExclusiveLock();
		return isFound;
	}
	bool remove(const KeyType &key) {
		return remove(key, getBucket(key));
	}

	bool _exists(const KeyType &key, IndexType &existingIdx, const BucketType &bucket) const {
		bool isFound;
		existingIdx = bucket.findSorted(key, isFound);
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

	const ElementType getElementIfExists(const KeyType &key, const BucketType &bucket) const {
		IndexType existingIdx;
		ElementType element;
		bucket.setSharedLock();
		if (_exists(key, existingIdx, bucket)) {
		    element = bucket.getElement(existingIdx);
		}
		bucket.unsetSharedLock();
		return element;
	}
	ElementType getElementIfExists(const KeyType &key, BucketType &bucket) {
		IndexType existingIdx;
		ElementType element;
		bucket.setSharedLock();
		if (_exists(key, existingIdx, bucket)) {
		    element = bucket.getElement(existingIdx);
		}
		bucket.unsetSharedLock();
		return element;
	}

	const ElementType getElementIfExists(const KeyType &key) const {
		return getElementIfExists(key, getBucket(key));
	}
	ElementType getElementIfExists(const KeyType &key) {
		return getElementIfExists(key, getBucket(key));
	}

	ElementType getOrSetElement(const KeyType &key, BucketType &bucket, ValueType &value) {
		IndexType existingIdx;
		ElementType element;
		bucket.setExclusiveLock();
		if (_exists(key, existingIdx, bucket)) {
		    element = bucket.getElement(existingIdx);
		} else {
		    element = insert(key, value, bucket);
		}
		bucket.unsetExclusiveLock();
		return element;
	}
	ElementType getOrSetElement(const KeyType &key, ValueType &value) {
		return getOrSetElement(key, getBucket(key), value);
	}
	ElementType getElement(const KeyType &key, BucketType &bucket) {
		ValueType value = Value();
		return getOrSetElement(key, bucket, value);
	}
	ElementType getElement(const KeyType &key) {
		ValueType value = Value();
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

	std::string toString() {
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

		long bucketsSize = getNumBuckets();
 		#pragma omp parallel for reduction(+:affected)
		for(long idx = 0 ; idx < bucketsSize; idx++) {
			affected += _buckets[idx].purgeMinCount(minimumCount);
		}
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

	// optimized merge for DMP threaded (interlaced) KmerMaps
	void merge(const KmerMap &src) const {
	    if (getNumBuckets() != src.getNumBuckets()) {
	    	 throw std::invalid_argument("Can not merge two KmerMaps of differing sizes!");
	    }

	    long bucketsSize = getNumBuckets();
		#pragma omp parallel for
	  	for(long idx = 0 ; idx < bucketsSize; idx++) {
	    	   BucketType &a = const_cast<BucketType&>(_buckets[idx]);
	    	   BucketType &b = const_cast<BucketType&>(src._buckets[idx]);

	    	   if (! _mergeTrivial(a,b) ) {
	    		   throw std::invalid_argument("Expected one bucket to be zero in this optimized method: KmerMap::merge(const KmerMap &src) const");
	    	   }
	  	}

	}

	void mergeAdd(KmerMap &src) {
       if (getNumBuckets() != src.getNumBuckets()) {
    	   throw std::invalid_argument("Can not merge two KmerMaps of differing sizes!");
       }
       BucketType merged;

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
		throw std::invalid_argument("This method is broken for threaded execution somehow (even with pragmas disabled)");
	       if (getNumBuckets() != src.getNumBuckets()) {
	    	   throw std::invalid_argument("Can not merge two KmerMaps of differing sizes!");
	       }
	       // calculated chunk size to ensure thread safety of merge operation
	       int chunkSize = mergeDest.getNumBuckets() / src.getNumBuckets();
	       if (chunkSize <= 1)
	    	   chunkSize = 2;
	       BucketType merged;

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

	public:
		Iterator(KmerMap *target):
		_target(target),
		_iBucket(target->_buckets.begin()),
		_iElement(_iBucket->begin())
		{
			_moveToNextValidElement();
		}

		Iterator(const Iterator &copy) {
			_target = copy._target;
			_iBucket = copy._iBucket;
			_iElement = copy._iElement;
		}

		Iterator(KmerMap *target, BucketsVectorIterator bucketPtr):
		_target(target),
		_iBucket(bucketPtr),
		_iElement()
		{
			if (!isEnd())
			_iElement = _iBucket->begin();
			_moveToNextValidElement();
		}

		Iterator(KmerMap *target, BucketsVectorIterator bucketPtr,BucketTypeIterator elementPtr):
		_target(target),
		_iBucket(bucketPtr),
		_iElement(elementPtr)
		{
		}

	private:
		inline bool isEnd() const {
			return _iBucket == _target->_buckets.end();
		}

		void _moveToNextValidElement() {
			while( (!isEnd()) && _iElement == _iBucket->end() && ++_iBucket != _target->_buckets.end())
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
		ConstIterator(KmerMap *target) : _iterator(target) {}
		ConstIterator(const ConstIterator &copy) : _iterator(copy._iterator) {}
		ConstIterator(KmerMap *target, BucketsVectorIterator bucketPtr) : _iterator(target,bucketPtr) {}
		ConstIterator(KmerMap *target, BucketsVectorIterator bucketPtr,BucketTypeIterator elementPtr): _iterator(target,bucketPtr,elementPtr) {}
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

	Iterator begin() {return Iterator(this);}
	Iterator end() {return Iterator(this,_buckets.end());}

	ConstIterator begin() const {return Iterator(this);}
	ConstIterator end() const {return Iterator(this, _buckets.end());}

	typedef std::pair<BucketsVectorIterator, BucketsVectorIterator> ThreadedBuckets;
	ThreadedBuckets _getThreadedBuckets() {
		long threads = omp_get_num_threads();

		if (threads == 1) {
			return ThreadedBuckets(_buckets.begin(), _buckets.end());
		} else {
			long threadNum = omp_get_thread_num();
			long step = _buckets.size() / threads;
			if (step == 0)
				step = 1;
			BucketsVectorIterator begin = _buckets.begin() + step*threadNum;
			BucketsVectorIterator end = begin + step;
			if (threadNum + 1 == threads || end > _buckets.end())
				end = _buckets.end();
			if (begin > _buckets.end()) {
				begin = _buckets.end();
			}
			return ThreadedBuckets(begin,end);
		}
	}
	Iterator beginThreaded() {
		return Iterator(this, _getThreadedBuckets().first);
	}
	Iterator endThreaded() {
		return Iterator(this, _getThreadedBuckets().second);
	}

};

typedef KmerArray<char> Kmers;
typedef KmerArray<double> KmerWeights;
typedef KmerArray<Kmernator::UI32> KmerCounts;

#endif

//
// $Log: Kmer.h,v $
// Revision 1.88  2010-08-18 17:50:40  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.87.4.1  2010-07-19 18:43:51  regan
// fixed core dump in parallel build of large sequences
//
// Revision 1.87  2010-06-23 17:10:12  regan
// bugfix in boundary conditions of threaded iterator
//
// Revision 1.86  2010-06-22 23:06:31  regan
// merged changes in CorruptionBugfix-20100622 branch
//
// Revision 1.85.6.1  2010-06-22 22:59:55  regan
// named all critical sections
//
// Revision 1.85  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.84.2.5  2010-05-18 16:43:31  regan
// added count gc methods and lookup tables
//
// Revision 1.84.2.4  2010-05-13 20:29:30  regan
// new methods to support changes to CompareSpectrum
//
// Revision 1.84.2.3  2010-05-12 18:24:40  regan
// minor refactor
//
// Revision 1.84.2.2  2010-05-10 17:41:18  regan
// refactored kmer sizer.
// reordered member variables for better padding and alignment
// added test for # of threads divisable by 2
//
// Revision 1.84.2.1  2010-05-07 22:59:33  regan
// refactored base type declarations
//
// Revision 1.84  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.83  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
// Revision 1.82.2.1  2010-05-06 18:45:35  regan
// broke it...
//
// Revision 1.82  2010-05-06 16:43:56  regan
// merged changes from ConsensusTesting-20100505
//
// Revision 1.81.8.1  2010-05-05 23:46:22  regan
// checkpoint... seems to compile
//
// Revision 1.81.2.1  2010-05-04 19:49:51  regan
// minor rework on include headers
//
// Revision 1.81  2010-05-01 21:57:54  regan
// merged head with serial threaded build partitioning
//
// Revision 1.80.2.7  2010-04-29 20:33:52  regan
// *** empty log message ***
//
// Revision 1.80.2.6  2010-04-29 06:59:10  regan
// a few more methods
//
// Revision 1.80.2.5  2010-04-27 05:39:10  regan
// got distributed thread merge working
//
// Revision 1.80.2.4  2010-04-26 05:25:49  regan
// bugfix and optimization
//
// Revision 1.80.2.3  2010-04-24 04:56:20  regan
// bugfixes
//
// Revision 1.80.2.2  2010-04-23 23:39:41  regan
// a few changes in how lest complement is derived
//
// Revision 1.80.2.1  2010-04-22 22:55:23  regan
// checkpoint
//
// Revision 1.80  2010-04-21 23:39:37  regan
// got kmermap mmap store and restore working
//
// Revision 1.79  2010-04-21 00:33:20  regan
// merged with branch to detect duplicated fragment pairs with edit distance
//
// Revision 1.78.2.3  2010-04-20 23:55:37  regan
// added parallel/threaded iterator handles
//
// Revision 1.78.2.2  2010-04-19 18:20:51  regan
// refactored base permutation
//
// Revision 1.78.2.1  2010-04-16 23:46:06  regan
// checkpoint
//
// Revision 1.78  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.77.2.1  2010-04-04 11:59:41  regan
// minor formatting
//
// Revision 1.77  2010-03-15 14:59:50  regan
// minor refactor
// disabled mergePromote as there is some error/race condition
//
// Revision 1.76  2010-03-14 17:16:30  regan
// bugfix
//
// Revision 1.75  2010-03-14 16:57:28  regan
// checkpoint
//
// Revision 1.74  2010-03-10 13:17:27  regan
// removed singleton count
// fixed omp bug in discarded tracking
//
// Revision 1.73  2010-03-04 06:38:15  regan
// fixed compiler warnings
//
// Revision 1.72  2010-03-02 15:02:40  regan
// fixed bugs in tracking data and singleton counting
//
// Revision 1.71  2010-02-26 13:01:17  regan
// reformatted
//
// Revision 1.70  2010-01-14 03:25:34  cfurman
// fixed non-existent boost hashing
//
// Revision 1.69  2010-01-13 23:34:17  regan
// made const class modifications
//
// Revision 1.68  2010-01-13 00:23:31  regan
// added open mp to build kmers for long references
//
// Revision 1.67  2010-01-08 06:22:27  regan
// fixed toNumber to have the correct bit ordering (pointers go backwards in memory!)
// removed a thow when the kmer size is too short
//
// Revision 1.66  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.65  2009-12-24 00:55:57  regan
// made const iterators
// fixed some namespace issues
// added support to output trimmed reads
//
// Revision 1.64  2009-12-24 00:39:22  cfurman
// getAverageWeight() added
//
// Revision 1.63  2009-12-22 18:31:14  regan
// moved Open MP includes to a central config file
//
// Revision 1.62  2009-12-21 06:34:26  regan
// used openmp and clever partitioning to speed up building spectrum
//
// Revision 1.61  2009-12-18 19:05:09  regan
// added compile-time thread-safety option to kmer classes
//
// Revision 1.60  2009-12-14 05:31:35  regan
// optimized array resizing to malloc at logarithmic stepping
// fixed a bug in KmerArray<>::findSorted
//
// Revision 1.59  2009-11-29 19:04:45  regan
// optimized a bit, fixed a few bugs
//
// Revision 1.58  2009-11-28 01:00:07  regan
// fixed bugs and warnings
//
// Revision 1.57  2009-11-27 23:16:58  regan
// refactored and broke it
//
// Revision 1.56  2009-11-27 01:53:40  regan
// refactored and got first pass at error rate by position
//
// Revision 1.55  2009-11-26 09:03:29  regan
// refactored and stuff
//
// Revision 1.54  2009-11-25 18:39:08  regan
// fixed bug across versions of compiler
//
// Revision 1.53  2009-11-24 13:35:29  cfurman
// removed KmerPtr class.
//
// Revision 1.52  2009-11-22 08:16:41  regan
// some fixes some bugs... optimized vs debug vs deb4/5 give different results
//
// Revision 1.51  2009-11-21 18:46:53  regan
// added bugs
//
// Revision 1.50  2009-11-21 15:58:29  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.49  2009-11-12 17:01:51  regan
// checkpoint
//
// Revision 1.48  2009-11-11 17:23:24  regan
// fixed bugs in heap generation
// solid picking logic needs work
//
// Revision 1.47  2009-11-11 07:57:23  regan
// built framework for  (not working) - make_heap is broken
//
// Revision 1.46  2009-11-09 19:37:17  regan
// enhanced some debugging / analysis output
//
// Revision 1.45  2009-11-07 00:26:13  cfurman
// minor formatting
//
// Revision 1.44  2009-11-06 16:59:11  regan
// added base substitution/permutations table and build function
//
// Revision 1.43  2009-11-06 04:08:23  regan
// minor changes
//
// Revision 1.42  2009-11-04 18:24:25  regan
// reworked tracking data
//
// Revision 1.41  2009-11-03 17:15:40  regan
// minor refactor
//
// Revision 1.40  2009-11-02 21:19:25  regan
// fixed types and boundary tests
//
// Revision 1.39  2009-11-02 20:02:48  cfurman
// KmerMap::iterator refactor
//
// Revision 1.38  2009-11-02 18:48:18  regan
// minor refactor and performance tweaks
//
// Revision 1.37  2009-10-31 23:44:17  regan
// fixed bug in KmerArray::remove
// refactored memory pool out of KmerArray
//
// Revision 1.36  2009-10-31 00:16:35  regan
// minor changes and optimizations
//
// Revision 1.35  2009-10-30 20:56:27  regan
// fixed kmermap iterator
//
// Revision 1.34  2009-10-30 19:27:46  regan
// added iterator goodness, but KmerMap::Iterator still does not work
//
// Revision 1.33  2009-10-30 00:51:40  regan
// bug fix and working on executable
//
// Revision 1.32  2009-10-30 00:10:32  regan
// cleaned up a bit
//
// Revision 1.31  2009-10-30 00:07:59  regan
// bugfix on KmerArray.build trailing bits
//
// Revision 1.30  2009-10-29 23:30:01  regan
// checkpoint
//
// Revision 1.29  2009-10-29 23:04:49  regan
// works
//
// Revision 1.28  2009-10-29 20:59:23  cfurman
// fixed testing bugs
//
// Revision 1.27  2009-10-29 19:01:33  regan
// checkpoint
//
// Revision 1.26  2009-10-29 17:00:58  regan
// checkpoint (with bugs)
//
// Revision 1.25  2009-10-29 07:03:33  regan
// fixed some bugs , added others
// KmerArray is working, *Sorted methods are untested
//
// Revision 1.24  2009-10-28 18:50:57  regan
// made KmerArray behave properly and not like a KmerPtrArray
//
// Revision 1.23  2009-10-28 18:42:59  regan
// added debug flags, fixed tests, bugs
//
// Revision 1.22  2009-10-28 02:29:55  cfurman
// fixed KmerArray  bugs
//
// Revision 1.21  2009-10-28 00:00:41  regan
// added more bugs
//
// Revision 1.20  2009-10-27 07:16:09  regan
// checkpoint
// defined KmerMap and KmerArray lookup methods
//
// Revision 1.19  2009-10-26 23:04:33  regan
// checkpoint make Kmer private inner class
//
// Revision 1.18  2009-10-26 17:50:54  regan
// templated KmerArray; added boost pool allocation
//
// Revision 1.17  2009-10-24 00:32:46  regan
// added bugs
//
// Revision 1.16  2009-10-24 00:03:49  regan
// checkpoint
//
// Revision 1.15  2009-10-23 23:22:41  regan
// checkpoint
//
// Revision 1.14  2009-10-23 21:54:46  regan
// checkpoint
//
// Revision 1.13  2009-10-23 20:32:50  cfurman
// more kmer changes
//
// Revision 1.12  2009-10-23 17:22:39  regan
// added more tests
//
// Revision 1.11  2009-10-23 07:06:59  regan
// more unit testing
//   ReadSetTest
//   KmerTest
//
// Revision 1.10  2009-10-23 01:24:53  cfurman
// ReadSet test created
//
// Revision 1.9  2009-10-22 01:39:43  cfurman
// bug fix in kmer.h
//
// Revision 1.8  2009-10-22 00:07:43  cfurman
// more kmer related classes added
//
// Revision 1.7  2009-10-21 18:58:44  regan
// checkpoint
//
// Revision 1.6  2009-10-21 18:44:20  regan
// checkpoint
//
// Revision 1.5  2009-10-21 06:51:34  regan
// bug fixes
// build lookup tables for twobitsequence
//
//
