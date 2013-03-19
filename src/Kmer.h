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
#include <boost/thread.hpp>
#include <boost/lockfree/queue.hpp>

#include "config.h"
#include "Options.h"
#include "TwoBitSequence.h"
#include "MemoryUtils.h"
#include "KmerTrackingData.h"
#include "lookup3.h"
#include "MmapTempFile.h"
#include "Log.h"

using namespace TwoBitSequenceBase;

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
	_KmerBaseOptions(SequenceLengthType defaultKmerSize = KmerSizer::getSequenceLength())
	: kmerSize(defaultKmerSize), kmersPerBucket(32) {
	}
	~_KmerBaseOptions() {}
	void _resetOptions() {
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {

		po::options_description opts("Kmer Options");

		opts.add_options()

			("kmer-size", po::value<SequenceLengthType>()->default_value(kmerSize), "kmer size.  A size of 0 will skip k-mer calculations")

			("kmers-per-bucket", po::value<unsigned int>()->default_value(kmersPerBucket), "number of kmers to target per hash-bucket.  Lesser will use more memory, larger will be slower");


		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		if (vm.count("kmer-size") == 0) {
			LOG_WARN(1, "There was no kmer size specified!");
			ret = false;
		}
		setOpt("kmer-size", kmerSize);
		setOpt("kmers-per-bucket", getKmersPerBucket());

		KmerSizer::set(kmerSize);

		return ret;
	}
	SequenceLengthType &getKmerSize()
	{
		return kmerSize;
	}
	unsigned int &getKmersPerBucket() {
		return kmersPerBucket;
	}
private:
	SequenceLengthType kmerSize;
	unsigned int kmersPerBucket;

};
typedef OptionsBaseTemplate< _KmerBaseOptions > KmerBaseOptions;

class Kmer;
class KmerInstance;

class KmerHasher {
public:
	typedef Kmernator::KmerNumberType NumberType;

	// safely returns 64-bit numeric version of any sized kmer
	static NumberType toNumber(const void *ptr, int len) {
		NumberType val;
		if (len >= 8) {
			boost::uint64_t *x = (boost::uint64_t *) ptr;
			val = (NumberType) *x;
			if (len > 8)
				val += toNumber(++x, len-8);
		} else if (len >=4) {
			val = (NumberType) *((boost::uint32_t *) ptr);
		} else if (len >=2) {
			val = (NumberType) *((boost::uint16_t *) ptr);
		} else {
			val = (NumberType) *((boost::uint8_t *)  ptr);
		}
		return val;
	};
	static NumberType getHash(const void *ptr, int length) {
		// initialize it so something tasty
		uint64_t hash = 0xDEADBEEF;

		// Old hash algorithm, gave poor distribution, if I remember correctly
		//		NumberType number = toNumber();
		//		return Lookup8::hash2(&number, 1, 0xDEADBEEF);

		// Presently the act of caching the last hash takes longer than calculating it again
		//		NumberType val = toNumber(ptr, length);
		//		MicroCache &mc = getThreadCache();
		//		if (mc.isCached(ptr,length,val,hash))
		//			return hash;

		uint32_t *pc, *pb;
		pc = (uint32_t*) &hash;
		pb = pc+1;
		Lookup3::hashlittle2(ptr, length, pc, pb);

		// Do not cache it
		//		mc.set(ptr,length,val,hash);

		return hash; // The same thing as the recommended mixing: return *pc + (((uint64_t)*pb)<<32);
	}
	NumberType operator()(const Kmer& kmer) const;
	NumberType operator()(const KmerInstance& kmer) const;
protected:
	class MicroCache {
	public:
		MicroCache() {
			reset();
		}
		void reset() {
			_ptr = NULL;
			_len = _val = _hash = 0;
		}
		bool isCached(const void *ptr, int len, NumberType val, NumberType &hash) {
			if (ptr == _ptr && _len == len && _val == val) {
				hash = _hash;
				return true;
			}
			return false;
		}
		void set(const void *ptr, int len, NumberType val, const NumberType &hash) {
			_ptr = ptr;
			_len = len;
			_val = val;
			_hash = hash;
		}
	private:
		const void *_ptr;
		NumberType  _val, _hash;
		int _len;
	};
	static MicroCache &getThreadCache() {
		static boost::thread_specific_ptr<MicroCache> ptr;
		if (ptr.get() == NULL)
			ptr.reset(new MicroCache());
		return *ptr;
	}

};

#define TEMP_KMER(name)  TwoBitEncoding _stack_##name[KmerSizer::getByteSize()]; Kmer &name = (Kmer &)(_stack_##name);
class Kmer {
public:
	typedef Kmernator::KmerNumberType NumberType;
	typedef Kmernator::KmerIndexType  IndexType;
	typedef Kmernator::KmerSizeType   SizeType;

	inline static KmerHasher &getHasher() {
		static KmerHasher hasher;
		return hasher;
	}
	Kmer &operator=(const Kmer &other) {
		if (this == &other)
			return *this;

		set(other);
		return *this;
	}
	// return just the low bits of the kmer
	static NumberType toNumber(const Kmer &kmer) {
		return KmerHasher::toNumber(kmer.getTwoBitSequence(), std::min((SequenceLengthType)8, KmerSizer::getTwoBitLength()));
	}

	// safely returns lowbits 64-bit numeric version of any sized kmer
	inline int compare(const Kmer &other) const {
		return memcmp((const void*) getTwoBitSequence(), (const void*) other.getTwoBitSequence(), getTwoBitLength());
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


	void buildReverseComplement(Kmer &output) const {
		TwoBitSequence::reverseComplement(getTwoBitSequence(),
				output.getTwoBitSequence(), getLength());
	}

	// returns true if this is the least complement, false otherwise (output is least)
	bool buildLeastComplement(Kmer &output) const {
		buildReverseComplement(output);
		if (*this <= output) {
			output.set(*this);
			return true;
		} else {
			return false;
		}
	}
	void set(const Kmer &copy) {
		memcpy(getTwoBitSequence(), copy.getTwoBitSequence(), getTwoBitLength());
	}

	std::string toFasta() const {
		const TwoBitEncoding *tmp = this->getTwoBitSequence();
		return TwoBitSequence::getFasta(tmp, 0, getLength());
	}
	std::string toFastaFull() const {
		return TwoBitSequence::getFasta(getTwoBitSequence(), 0, getTwoBitLength() * 4);
	}

	inline NumberType hash() const { return getHasher()(*this); }

	TwoBitEncoding *getTwoBitSequence() {
		return (TwoBitEncoding *) _data();
	}
	const TwoBitEncoding *getTwoBitSequence() const {
		return (const TwoBitEncoding *) _data();
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
	void swap(Kmer &other) {
		TEMP_KMER(temp);
		temp = other;
		other = *this;
		*this = temp;
	}

protected:
	Kmer(); // never construct, just use as cast

	// No data for you!!!
	const void *_data() const {
		return this;
	}
	void *_data() {
		return this;
	}

};

#define MAX_KMER_INSTANCE_BYTES 24
class KmerInstance {
public:
	typedef Kmer Base;
	typedef Kmernator::KmerNumberType NumberType;
	typedef Kmernator::KmerIndexType  IndexType;
	typedef Kmernator::KmerSizeType   SizeType;

	inline static KmerHasher &getHasher() {
		return Kmer::getHasher();
	}
	KmerInstance() {
		init();
	}
	KmerInstance(const Kmer &copy) {
		init();
		*this = copy;
	}
	KmerInstance(const KmerInstance &copy) {
		init();
		*this = copy;
	}
	KmerInstance &operator=(const Kmer &other)  {
		memcpy(_data, other.getTwoBitSequence(), getByteSize());
		return *this;
	}
	KmerInstance &operator=(const KmerInstance &other) {
		if (this == &other)
			return *this;
		memcpy(_data, other._data, getByteSize());
		return *this;
	}
	~KmerInstance() {
		destroy();
	}

	// cast operator
	operator Kmer*() {
		return (Kmer*) getTwoBitSequence();
	}
	operator Kmer&() {
		return *((Kmer*) getTwoBitSequence());
	}
	operator const Kmer*() const {
		return (Kmer*) getTwoBitSequence();
	}
	operator const Kmer&() const {
		return *((Kmer*) getTwoBitSequence());
	}

	inline int compare(const Kmer &other) const {
		return memcmp((const void*) getTwoBitSequence(), (const void*) other.getTwoBitSequence(), getTwoBitLength());
	}
	inline int compare(const KmerInstance &other) const {
		return memcmp((const void*) getTwoBitSequence(), (const void*) other.getTwoBitSequence(), getTwoBitLength());
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
	inline bool operator ==(const Kmer &other) const {
		return compare(other) == 0;
	}
	inline bool operator ==(const KmerInstance &other) const {
		return compare(other) == 0;
	}
	inline bool operator !=(const Kmer &other) const {
		return compare(other) != 0;
	}
	inline bool operator !=(const KmerInstance &other) const {
		return compare(other) != 0;
	}
	inline bool operator <(const Kmer &other) const {
		return compare(other) < 0;
	}
	inline bool operator <(const KmerInstance &other) const {
		return compare(other) < 0;
	}
	inline bool operator <=(const Kmer &other) const {
		return compare(other) <= 0;
	}
	inline bool operator <=(const KmerInstance &other) const {
		return compare(other) <= 0;
	}
	inline bool operator >(const Kmer &other) const {
		return compare(other) > 0;
	}
	inline bool operator >(const KmerInstance &other) const {
		return compare(other) > 0;
	}
	inline bool operator >=(const Kmer &other) const {
		return compare(other) >= 0;
	}
	inline bool operator >=(const KmerInstance &other) const {
		return compare(other) >= 0;
	}


	void buildReverseComplement(Kmer &output) const {
		TwoBitSequence::reverseComplement(getTwoBitSequence(),
				output.getTwoBitSequence(), getLength());
	}
	void buildReverseComplement(KmerInstance &output) const {
		buildReverseComplement((Kmer&)output);
	}

	// returns true if this is the least complement, false otherwise (output is least)
	bool buildLeastComplement(Kmer &output) const {
		buildReverseComplement(output);
		if (*this <= output) {
			output.set(*this);
			return true;
		} else {
			return false;
		}
	}
	bool buildLeastComplement(KmerInstance &output) const {
		return buildLeastComplement((Kmer&) output);
	}
	void set(const Kmer &copy) {
		memcpy(getTwoBitSequence(), copy.getTwoBitSequence(), getTwoBitLength());
	}
	void set(const KmerInstance &copy) {
		memcpy(getTwoBitSequence(), copy.getTwoBitSequence(), getTwoBitLength());
	}

	std::string toFasta() const {
		const TwoBitEncoding *tmp = this->getTwoBitSequence();
		return TwoBitSequence::getFasta(tmp, 0, getLength());
	}
	std::string toFastaFull() const {
		return TwoBitSequence::getFasta(getTwoBitSequence(), 0, getTwoBitLength() * 4);
	}

	TwoBitEncoding *getTwoBitSequence() {
		return (TwoBitEncoding *) _data;
	}
	const TwoBitEncoding *getTwoBitSequence() const {
		return (const TwoBitEncoding *) _data;
	}

	void set(std::string fasta, bool leastComplement = false) {
		assert(fasta.length() == KmerSizer::getSequenceLength());
		if (leastComplement) {
			KmerInstance temp;
			TwoBitSequence::compressSequence(fasta, temp.getTwoBitSequence());
			temp.buildLeastComplement(*this);
		} else {
			TwoBitSequence::compressSequence(fasta, this->getTwoBitSequence());
		}
	}
	inline NumberType hash() const { return getHasher()(*this); }

private:
	inline void init() {
	}
	inline void destroy() {
	}
	char _data[MAX_KMER_INSTANCE_BYTES];
};

template<typename Value>
class KmerElementBaseRef : public std::pair<const Kmer &, Value&> {
public:
	typedef Value ValueType;
	typedef typename std::pair<const Kmer&, ValueType&> BaseRef;
	typedef typename std::pair<const Kmer&, const ValueType&> ConstBaseRef;

	KmerElementBaseRef(const Kmer &kmer, ValueType & value) : BaseRef(kmer, value) {}
	KmerElementBaseRef(BaseRef base): BaseRef(base) {}

	// cast to Base operators
	operator BaseRef&() { return (BaseRef&) *this; }
	operator BaseRef*() { return (BaseRef*) this; }

	virtual inline const Kmer &key() const { return this->first; }
	virtual const ValueType &value() const { return this->second; }
	virtual ValueType &value() { return this->second; };
};

template<typename Value>
class ConstKmerElementBaseRef : public std::pair<const Kmer &, const Value&> {
public:
	typedef Value ValueType;
	typedef typename std::pair<const Kmer&, ValueType&> BaseRef;
	typedef typename std::pair<const Kmer&, const ValueType&> ConstBaseRef;

	// cast to Base operators
	operator ConstBaseRef&() { return (ConstBaseRef&) *this; }
	operator ConstBaseRef*() { return (ConstBaseRef*) this; }

	ConstKmerElementBaseRef(const Kmer &kmer, const ValueType & value) : ConstBaseRef(kmer, value) {}
	ConstKmerElementBaseRef(BaseRef base) : BaseRef(base) {}
	virtual inline const Kmer &key() const { return this->first; }
	virtual const ValueType &value() const { return this->second; }
};

template<typename Value>
class KmerElementPair {
public:
	typedef Value ValueType;
	typedef typename std::pair<const Kmer*, ValueType*> Base;
	typedef typename Base::first_type _KmerType;
	typedef typename Base::second_type _ValueType;
	typedef KmerElementBaseRef<ValueType> BaseRef;
	typedef ConstKmerElementBaseRef<ValueType> ConstBaseRef;

	KmerElementPair() : _kmer(NULL), _value(NULL) {}
	KmerElementPair(const Base &pair) : _kmer(pair.first), _value(pair.second) {}
	KmerElementPair(const BaseRef &pair): _kmer(&pair.first), _value(&pair.second) {}
	KmerElementPair(const Kmer &kmer, ValueType &value) : _kmer(&kmer), _value(&value) {}
	KmerElementPair(const KmerElementPair &copy) : _kmer(copy._kmer), _value(copy._value) {}

	KmerElementPair &operator=(const KmerElementPair &copy) {
		_kmer = copy._kmer;
		_value = copy._value;
		return *this;
	}
	KmerElementPair &operator=(const Base &copy) {
		_kmer = copy.first;
		_value = copy.second;
		return *this;
	}

	inline bool isValid() const { return _kmer != NULL; }
	inline void reset() { _kmer = NULL; _value = NULL; }

	inline const Kmer &key() const { assert(isValid()); return *_kmer; }
	const ValueType &value() const { assert(isValid()); return *_value; }
	ValueType &value() { assert(isValid()); return *_value; };

	static const KmerElementPair getConst() {
		return KmerElementPair();
	}
	static const KmerElementPair getConst(const Kmer &kmer, const ValueType &value) {
		return KmerElementPair(kmer, const_cast<ValueType&>(value));
	}
protected:
	_KmerType _kmer;
	_ValueType _value;
};

template<typename Iterator>
class IteratorOfMapIterators {
public:
	typedef Iterator Base;
	typedef typename Base::value_type BaseValueType;
	typedef typename BaseValueType::iterator BaseIterator;
	typedef typename BaseValueType::key_type BaseIteratorKeyType;
	typedef typename BaseValueType::value_type BaseIteratorValueType;
	typedef typename BaseValueType::mapped_type BaseIteratorMappedType;
	typedef BaseIteratorKeyType key_type;
	typedef BaseIteratorValueType value_type;
	typedef BaseIteratorMappedType mapped_type;

	IteratorOfMapIterators() : _begin(), _end(), _it(), _rank(0), _size(1) {}
	IteratorOfMapIterators(Base begin, Base end, int rank = 0, int size = 1) :
		_begin(begin), _end(end), _it(), _rank(rank), _size(size) {
		assert(_rank < _size && _size > 0);
		incrBucket(_rank); // step to first bucket for this rank
		setNextValid();
	}
	IteratorOfMapIterators(const IteratorOfMapIterators &copy) {
		*this = copy;
	}
	IteratorOfMapIterators &operator=(const IteratorOfMapIterators &copy) {
		_begin = copy._begin;
		_end = copy._end;
		_it = copy._it;
		_rank = copy._rank;
		_size = copy._size;
		return *this;
	}

	bool operator==(const IteratorOfMapIterators &other) const {
		if (other._begin == other._end)
			return _begin == _end;
		else if (_begin == _end)
			return false;
		else
			return _it == other._it && _begin == other._begin && _end == other._end;
	}
	bool operator!=(const IteratorOfMapIterators &other) const {
		return !operator==(other);
	}
	IteratorOfMapIterators &operator++() {
		incr();
		return *this;
	}
	IteratorOfMapIterators operator++(int unused) {
		IteratorOfMapIterators tmp(*this); ++(*this); return tmp;
	}

	BaseIteratorValueType &operator*() { return *_it; }
	BaseIteratorValueType *operator->() { return &(*_it); }

	template<typename Container>
	static IteratorOfMapIterators getEnd(Container &cont, int rank = 0, int size = 1) {
		return IteratorOfMapIterators(cont.end(), cont.end(), rank, size);
	}

protected:
	void incr() {
		assert(_begin != _end && _it != _begin->end());
		++_it;
		setNextValid();
	}
	void incrBucket(int num) {
		for(int i = 0 ; i < num; i++) {
			if (_begin != _end)
				++_begin;
		}
		if (_begin != _end)
			_it = _begin->begin();
	}
	// jumps _it to next non-empty BaseIterator
	void setNextValid() {
		while (_begin != _end) {
			if (_it == _begin->end()) {
				incrBucket(_size); // step to next bucket for this rank
			} else
				break;
		}
	}

private:
	Base _begin, _end;
	BaseIterator _it;
	int _rank, _size;

};


template<typename Value>
class KmerArrayPair {

public:
	typedef Kmer::NumberType    NumberType;
	typedef Kmer::IndexType     IndexType;
	typedef Kmer::SizeType      SizeType;

	typedef Value ValueType;
	static const IndexType MAX_INDEX = MAX_KMER_INDEX;
	static const float GROWTH_FACTOR = 1.5;
	static const IndexType MIN_GROWTH = 48;
	static const IndexType MAX_UNSORTED = 16;

	typedef std::vector< KmerArrayPair > Vector;
	typedef KmerElementPair<ValueType> BaseElementType;

	class Iterator;
	class ConstIterator;

	class ElementType : public BaseElementType {
	private:
		KmerArrayPair *_array;
		IndexType _idx;
		inline const ElementType &_constThis() const {return *this;}

		static const Kmer &getKmer(KmerArrayPair &kmerArray, IndexType index) {
			return kmerArray.get(index);
		}
		static ValueType &getValue(KmerArrayPair &kmerArray, IndexType index) {
			return kmerArray.valueAt(index);
		}

	public:
		typedef BaseElementType Base;
		ElementType() : Base(), _array(NULL), _idx(MAX_INDEX) {	}
		ElementType(IndexType idx,	KmerArrayPair &array) :
			Base(getKmer(array,idx), getValue(array,idx)),	_array(&array), _idx(idx) {
		}
		ElementType(IndexType idx, const KmerArrayPair &array) :
			// HACK! const_cast ... need ConstElementType
			Base(getKmer(array,idx), getValue(array,idx)), _array(const_cast<KmerArrayPair*>(&array)), _idx(idx) {
		}

		ElementType(const ElementType &copy) :
			Base(),	_array(NULL), _idx(MAX_INDEX)  {
			*this = copy;
		}
		~ElementType() {
			reset();
		}
		ElementType &operator=(const ElementType &copy) {
			reset();
			*((Base*)this) = (const Base&) copy;
			_idx = copy._idx;
			_array = copy._array;
			return *this;
		}
		void reset() {
			Base::reset();
			_idx = MAX_INDEX;
			_array = NULL;
		}
		inline bool isValid() const {
			return Base::isValid() && _array != NULL && _idx != MAX_INDEX;
		}
		Base::key;
		Base::value;
		Base::getConst;

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

	KmerArrayPair(void *begin, IndexType size, IndexType capacity) : _begin(begin), _size(size), _capacity(capacity), _endSorted(0) {
	}

public:

	KmerArrayPair(IndexType size = 0) :
		_begin(NULL), _size(0), _capacity(0), _endSorted(0) {
		resize(size);
	}

	KmerArrayPair(const TwoBitEncoding *twoBit, SequenceLengthType length, bool leastComplement = false, bool *bools = NULL) :
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

	KmerArrayPair(const KmerArrayPair &copy) :
		_begin(NULL), _size(0), _capacity(0), _endSorted(0)  {
		*this = copy;
	}

	KmerArrayPair copyRange(SequenceLengthType offset, SequenceLengthType length) {
		assert(offset+length <= _size);
		if (offset == 0 && length == _size)
			return KmerArrayPair(*this);
		KmerArrayPair splice(length);
		splice._copyRange(_begin, getValueStart(), 0, offset, length, false);
		return splice;
	}

	~KmerArrayPair() {
		reset();
	}

	KmerArrayPair &operator=(const KmerArrayPair &other) {
		
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
					"RuntimeError: KmerArrayPair::operator=(): Could not allocate memory in KmerArrayPair operator=()");

		_copyRange(other._begin, other.getValueStart(), 0, 0, _size, false);
		_endSorted = other._endSorted;
		assert(_endSorted <= size());

		return *this;
	}

	// restore a new array from a mmap, allocating new memory
	KmerArrayPair(const void *src) : _begin(NULL), _size(0), _capacity(0), _endSorted(0) {
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
	static const KmerArrayPair restore(const void *src) {
		const IndexType *size = (const IndexType *) src;
		IndexType xsize = *(size++);
		const void *ptr = size;
		KmerArrayPair array(const_cast<void*>(ptr), xsize, MAX_INDEX);
		array.setLastSorted();
		return array;
	}
	// never thread safe!
	const ValueType *getValueStart() const {
		if (capacity() > 0)
			return (ValueType*) _add(_begin, capacity());
		else
			return NULL;
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
					"Invalid: Kmer::operator[](): attempt to access index greater than size in KmerArrayPair operator[] const");
		return get(index);
	}

	// never thread safe!
	Kmer &operator[](IndexType index) {
		if (index >= size())
			LOG_THROW(
					"Invalid: Kmer::operator[](): attempt to access index greater than size in KmerArrayPair operator[]");
		return get(index);
	}

	// never thread safe!
	const ValueType &valueAt(IndexType index) const {
		if (index >= _size) {
			LOG_THROW(
					"Invalid: Kmer::valueAt(): attempt to access index greater than size in KmerArrayPair valueAt() const");
		}
		return *(getValueStart() + index);
	}

	// never thread safe!
	ValueType &valueAt(IndexType index) {
		if (index >= _size) {
			LOG_THROW(
					"Invalid: Kmer::valueAt(): attempt to access index greater than size in KmerArrayPair valueAt()");
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
					"RuntimeError: KmerArrayPair::resize(): Could not allocate memory in KmerArrayPair resize()");
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
				LOG_THROW("RuntimeError: Could not allocate memory in KmerArrayPair _setMemory(): Attempt to malloc " << newCapacity
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

		KmerArrayPair &kmers = *this;
		long numBytes = (numKmers + 3) / 4;

#pragma omp parallel for if(numKmers >= 10000)
		for (long bytes = 0; bytes < numBytes; bytes++) {
			SequenceLengthType i = bytes * 4;
			const TwoBitEncoding *ref = twoBit + bytes;
			for (int bitShift = 0; bitShift < 4 && i + bitShift < numKmers; bitShift++) {
				TwoBitSequence::shiftLeft(ref, kmers[i + bitShift].getTwoBitSequence(),
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
	static void permuteBases(const Kmer &kmer, KmerArrayPair &kmers, short editDistance, bool leastComplement = false) {
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
	static SequenceLengthType __permuteBases(const Kmer &kmer, KmerArrayPair &kmers, SequenceLengthType offset, SequenceLengthType startIdx, short editDistance, bool leastComplement = false) {
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
	// return a KmerArrayPair that has one entry for each possible single-base substitution
	static KmerArrayPair permuteBases(const Kmer &kmer, bool leastComplement = false) {
		KmerArrayPair kmers(KmerSizer::getSequenceLength() * 3);
		_permuteBases(kmer, kmers, leastComplement);
		return kmers;
	}
	static KmerArrayPair permuteBases(const Kmer &kmer, const ValueType defaultValue, bool leastComplement = false) {
		KmerArrayPair kmers(KmerSizer::getSequenceLength() * 3);
		_permuteBases(kmer, kmers, leastComplement);
		for(SequenceLengthType idx = 0; idx < kmers.size(); idx++)
			kmers.valueAt(idx) = defaultValue;
		return kmers;
	}
	static void _permuteBases(const Kmer &kmer, KmerArrayPair &kmers, bool leastComplement, SequenceLengthType offset = 0) {
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

	static KmerArrayPair extendKmer(const Kmer &kmer, bool toRight, bool leastComplement = false) {
		KmerArrayPair kmers(4);
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
	IndexType _findIndex(const Kmer &target, IndexType start, IndexType end) const {
		for(IndexType i=start; i<end; i++) {
			if (target.compare(get(i)) == 0) {
				return i;
			}
		}
		return MAX_INDEX;
	}
public:
	IndexType findIndex(const Kmer &target) const {
		bool targetIsFound;
		IndexType idx;
		if (_endSorted >= 8) {
			idx = _findSortedIndex(target, targetIsFound, 0, _endSorted);
			if (targetIsFound)
				return idx;
		}
		idx = _findIndex(target, _endSorted >= 16 ? _endSorted : 0, size());
		return idx;
	}
	IndexType findIndex(const Kmer &target, bool &targetIsFound) const {
		IndexType idx = findIndex(target);
		targetIsFound = idx != MAX_INDEX;
		return idx;
	}
protected:
	IndexType _findSortedIndex(const Kmer &target, bool &targetIsFound, IndexType start, IndexType end) const {
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
	IndexType _findSortedIndex(const Kmer &target, bool &targetIsFound) const {
		return _findSortedIndex(target, targetIsFound, 0, size());
	}
public:
	IndexType findSortedIndex(const Kmer &target, bool &targetIsFound) const {
		IndexType idx = _findSortedIndex(target, targetIsFound);
		return idx;
	}

protected:
	void _insertAt(IndexType idx, const Kmer &target) {
		assert(!isMmaped()); // mmaped can not be modified!
		if (idx > size())
			LOG_THROW("Invalid: attempt to access index greater than size in KmerArrayPair insertAt");
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
		IndexType idx = findSortedIndex(target, isFound);
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
		IndexType idx = findIndex(target, isFound);
		if (isFound)
			remove(idx);
	}
	void remove(IndexType idx) {
		assert(!isMmaped()); // mmaped can not be modified!
		resize(size()-1,idx);
		if (_endSorted > idx)
			_endSorted--;
	}


	//
	// compatibility API with STL interfaces (should refactor eventually)
	//
	ConstIterator find(const Kmer &target) const {
		IndexType idx = findIndex(target);
		if (idx == MAX_INDEX)
			return end();
		else
			return ConstIterator(this, idx);
	}
	Iterator find(const Kmer &target) {
		IndexType idx = findIndex(target);
		if (idx == MAX_INDEX)
			return end();
		else
			return Iterator(this, idx);
	}
	Iterator erase(Iterator it) {
		if (it != end()) {
			Iterator it2 = it;
			it2++;
			remove(it.getIndex());
			return it2;
		} else
			return it;
	}
	int erase(const Kmer &target) const {
		IndexType idx = findIndex(target);
		if (idx == MAX_INDEX)
			return 0;
		else
			remove(idx);
		return 1;
	}
	Iterator insert(Iterator hint, BaseElementType element) {
		return insert(element.key(), element.value());
	}
	Iterator insert(const Kmer &target, const ValueType &value) {
		IndexType idx;
		if (isSorted())
			idx = insertSorted(target, value);
		else
			idx = append(target, value);
		return Iterator(this, idx);
	}

	void clear(bool releaseMemory = true) {
		reset(releaseMemory);
	}
	// for STL compatibility
	typedef Kmer& key_type;
	typedef BaseElementType value_type;
	typedef ValueType mapped_type;
	typedef Iterator iterator;
	typedef ConstIterator const_iterator;

	class CompareArrayIdx {
		const KmerArrayPair &_kmerArray;
	public:
		CompareArrayIdx(const KmerArrayPair &kmerArray): _kmerArray(kmerArray) {}
		bool operator()(IndexType i, IndexType j) {
			return _kmerArray.get(i).compare(_kmerArray.get(j)) <= 0;
		}
	};

	bool isSorted() const {
		return _endSorted == size();
	}

	IndexType getNumUnsorted() const {
		return _endSorted - size();
	}

	IndexType getNumSorted() const {
		return _endSorted;
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
	void resort(std::vector< IndexType > &passingIndexes, bool reserveExtra = false) {
		std::sort(passingIndexes.begin(), passingIndexes.end(), CompareArrayIdx(*this));

		KAPPtr s = popTmp();
		s->clear(false);
		IndexType passingSize = passingIndexes.size();
		IndexType newSize = passingSize + (reserveExtra ? std::max((IndexType) (passingSize * GROWTH_FACTOR + 0.5), (IndexType) (passingSize + MIN_GROWTH)) : 0);
		newSize = std::max(newSize, _capacity);
		s->reserve(newSize);

		for(IndexType i = 0; i < passingSize; i++)
			s->append(get(passingIndexes[i]), valueAt(passingIndexes[i]));

		s->_endSorted = passingSize;
		swap(*s);
		pushTmp(s);
		LOG_DEBUG(5, "resort()" << _endSorted << ", " << size() << ", " << passingSize << " on " << this);
		assert(_endSorted <= size());
	}
public:
	void swap(IndexType idx1, IndexType idx2) {
		assert(!isMmaped()); // mmaped can not be modified!
		if (idx1 == idx2)
			return;
		if (idx1 >= size() || idx2 >= size())
			LOG_THROW("Invalid: attempt to access index greater than size in KmerArrayPair swap()");

		get(idx1).swap(get(idx2));
		if (sizeof(ValueType) > 0) {
			ValueType tmp = valueAt(idx1);
			valueAt(idx1) = valueAt(idx2);
			valueAt(idx2) = tmp;
		}
	}
	void swap(KmerArrayPair &other) {
		std::swap(_begin, other._begin);
		std::swap(_size, other._size);
		std::swap(_capacity, other._capacity);
		std::swap(_endSorted, other._endSorted);
		assert(_endSorted <= size());
	}

	std::string toThis() const {
		std::stringstream ss;
		ss << "KmerArrayPair " << this << "={" << _begin << ", " << _size << ", " << _capacity << "}";
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
	class Iterator : public std::iterator<std::forward_iterator_tag, KmerArrayPair>
	{
	private:
		KmerArrayPair *_tgt;
		IndexType _idx;
		ElementType thisElement;

		void setElement() {if ( !isEnd() ) thisElement = _tgt->getElement(_idx);}
	protected:
		KmerArrayPair &getKmerArray() {
			assert(_tgt != NULL);
			return *_tgt;
		}
	public:
		typedef const Kmer key_type;
		typedef ElementType value_type;
		typedef typename ElementType::ValueType mapped_type;

		Iterator(KmerArrayPair *target, IndexType idx = 0): _tgt(target), _idx(idx) {
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
		const Kmer &key() {return thisElement.key();}
		const Kmer &key() const {return thisElement.key();}
		Value &value() {return thisElement.value();}
		const Value &value() const {return thisElement.value();}
		ElementType *operator->() {return &thisElement;}
		bool isEnd() const {if (_tgt == NULL) { return true; } else { return _idx >= _tgt->size(); }}

		IndexType getIndex() const {
			assert(_tgt != NULL);
			return _idx;
		}
		ElementType getElement() {
			return thisElement;
		}

	};
	class ConstIterator : public std::iterator<std::forward_iterator_tag, KmerArrayPair>
	{
	private:
		Iterator _iterator;
	public:
		typedef const Kmer key_type;
		typedef const ElementType value_type;
		typedef typename ElementType::ValueType mapped_type;

		ConstIterator(const KmerArrayPair *target, IndexType idx = 0): _iterator(const_cast<KmerArrayPair*>( target ),idx) {}
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

public:
	// temporary cache to accelerate reuse of memory
	typedef std::auto_ptr< KmerArrayPair > KAPPtr;
	static KAPPtr popTmp() {
		KmerArrayPair *tmp;
		if (!getTmpQueue().pop(tmp)) {
			tmp = new KmerArrayPair();
		}
		return KAPPtr(tmp);
	}
	static bool pushTmp(KAPPtr tmp) {
		if (tmp.get() == NULL)
			return false;
		KmerArrayPair *_tmp = tmp.release();
		if (getTmpQueue().bounded_push(_tmp)) {
			return true;
		} else {
			delete _tmp;
			return false;
		}
	}
	static void clearTmp() {
		KmerArrayPair *tmp;
		while (getTmpQueue().pop(tmp))
			if (tmp != NULL)
				delete tmp;
	}
	static void initTmp() {
		getTmpQueue();
	}
private:
	typedef boost::lockfree::queue< KmerArrayPair *, boost::lockfree::fixed_sized<true> > KAPPtrQueue;
	typedef std::auto_ptr< KAPPtrQueue > KAPPtrQueuePtr;
	static KAPPtrQueuePtr tmpQueue;
	static KAPPtrQueue &getTmpQueue() {
		if (tmpQueue.get() == NULL) {
			tmpQueue.reset( new KAPPtrQueue(64) );
			assert(tmpQueue.get() != NULL);
			assert(tmpQueue->empty());
			KmerArrayPair *x = NULL;
			assert(!tmpQueue->pop(x));
			assert(x == NULL);
		}
		return *tmpQueue;
	}
};

template<typename Value>
typename KmerArrayPair<Value>::KAPPtrQueuePtr KmerArrayPair<Value>::tmpQueue;

template<typename Iterator, typename Value>
class KmerPairIteratorWrapper {
public:
	typedef Iterator Base;
	typedef KmerElementPair<Value> KEP;
	typedef typename KEP::BaseRef BaseRef;
	typedef std::auto_ptr< BaseRef > BaseRefPtr;
	typedef Kmer& key_type;
	typedef KmerElementPair<Value> value_type;
	typedef Value mapped_type;
	typedef typename KmerArrayPair<Value>::Iterator KAPI;

	KmerPairIteratorWrapper() : _iter() {}
	KmerPairIteratorWrapper(const Base &iter): _iter(iter) {	}
	KmerPairIteratorWrapper(const KmerPairIteratorWrapper &copy): _iter(copy._iter) {	}
	~KmerPairIteratorWrapper() {
	}
	KmerPairIteratorWrapper &operator=(const Base &copy) {
		_iter = copy;
		brptr.reset();
		return *this;
	}
	KmerPairIteratorWrapper &operator=(const KmerPairIteratorWrapper &copy) {
		_iter = copy._iter;
		brptr.reset();
		return *this;
	}

	// cast operators
	operator Base&() {
		return _iter;
	}
	operator const Base&() const {
		return _iter;
	}
	bool operator==(const Base &other) const { return _iter == other; }
	bool operator!=(const Base &other) const { return _iter != other; }
	KmerPairIteratorWrapper &operator++() { _iter++; brptr.reset(); return *this; }
	KmerPairIteratorWrapper operator++(int unused) { KmerPairIteratorWrapper tmp(*this); ++(*this); return tmp; }
	BaseRef &operator*() { setRefs(); return *brptr; }
	BaseRef *operator->() { setRefs(); return brptr.get(); }

	inline const Kmer &key() {
		return _iter->first;
	}
	Value &value() {
		return _iter->second;
	}
private:
	Base _iter;
	// WARNING: only instantiate when required to dereference
	mutable BaseRefPtr brptr;
	void setRefs() {
		if (brptr.get() == NULL)
			brptr.reset(new BaseRef(key(), value()));
	}

};

// need this specialization for KmerArrayPair which does not have a std::pair
// as the return type of the dereferenced iterator..
template<typename Value>
class KmerPairIteratorWrapper<typename KmerArrayPair<Value>::Iterator, Value> {
public:
	typedef typename KmerArrayPair<Value>::Iterator Base;
	typedef KmerElementPair<Value> KEP;
	typedef typename KEP::BaseRef BaseRef;
	typedef std::auto_ptr< BaseRef > BaseRefPtr;
	typedef Kmer& key_type;
	typedef KmerElementPair<Value> value_type;
	typedef Value mapped_type;

	KmerPairIteratorWrapper() : _iter() {}
	KmerPairIteratorWrapper(const Base &iter): _iter(iter) {	}
	KmerPairIteratorWrapper(const KmerPairIteratorWrapper &copy): _iter(copy._iter) {	}
	~KmerPairIteratorWrapper() {
	}
	KmerPairIteratorWrapper &operator=(const Base &copy) {
		_iter = copy;
		brptr.reset();
		return *this;
	}
	KmerPairIteratorWrapper &operator=(const KmerPairIteratorWrapper &copy) {
		_iter = copy._iter;
		brptr.reset();
		return *this;
	}

	// cast operators
	operator Base&() {
		return _iter;
	}
	operator const Base&() const {
		return _iter;
	}
	bool operator==(const Base &other) const { return _iter == other; }
	bool operator!=(const Base &other) const { return _iter != other; }
	KmerPairIteratorWrapper &operator++() { _iter++; brptr.reset(); return *this; }
	KmerPairIteratorWrapper operator++(int unused) { KmerPairIteratorWrapper tmp(*this); ++(*this); return tmp; }
	BaseRef &operator*() { setRefs(); return *brptr; }
	BaseRef *operator->() { setRefs(); return brptr.get(); }

	inline const Kmer &key() {
		return _iter.key();
	}
	Value &value() {
		return _iter.value();
	}
private:
	Base _iter;
	// WARNING: only instantiate when required to dereference
	mutable BaseRefPtr brptr;
	void setRefs() {
		if (brptr.get() == NULL)
			brptr.reset(new BaseRef(key(), value()));
	}

};

template<typename Iterator, typename Value>
class ConstKmerPairIteratorWrapper {
public:
	typedef Iterator Base;
	typedef KmerElementPair<Value> KEP;
	typedef typename KEP::ConstBaseRef BaseRef;
	typedef std::auto_ptr< BaseRef > BaseRefPtr;
	typedef Kmer& key_type;
	typedef KmerElementPair<Value> value_type;
	typedef Value mapped_type;
	typedef typename KmerArrayPair<Value>::Iterator KAPI;

	ConstKmerPairIteratorWrapper() : _iter() {}
	ConstKmerPairIteratorWrapper(const Base &iter): _iter(iter) {	}
	ConstKmerPairIteratorWrapper(const ConstKmerPairIteratorWrapper &copy): _iter(copy._iter) {	}
	~ConstKmerPairIteratorWrapper() {
	}
	ConstKmerPairIteratorWrapper &operator=(const Base &copy) {
		_iter = copy;
		brptr.reset();
		return *this;
	}
	ConstKmerPairIteratorWrapper &operator=(const ConstKmerPairIteratorWrapper &copy) {
		_iter = copy._iter;
		brptr.reset();
		return *this;
	}

	// cast operators
	operator Base&() {
		return _iter;
	}
	operator const Base&() const {
		return _iter;
	}
	bool operator==(const Base &other) const { return _iter == other; }
	bool operator!=(const Base &other) const { return _iter != other; }
	ConstKmerPairIteratorWrapper &operator++() { _iter++; brptr.reset(); return *this; }
	ConstKmerPairIteratorWrapper operator++(int unused) { ConstKmerPairIteratorWrapper tmp(*this); ++(*this); return tmp; }
	BaseRef &operator*() { setRefs(); return *brptr; }
	BaseRef *operator->() { setRefs(); return brptr.get(); }

	inline const Kmer &key() {
		return _iter->first;
	}
	const Value &value() {
		return _iter->second;
	}
private:
	Base _iter;
	// WARNING: only instantiate when required to dereference
	mutable BaseRefPtr brptr;
	void setRefs() {
		if (brptr.get() == NULL)
			brptr.reset(new BaseRef(key(), value()));
	}

};

// need this specialization for KmerArrayPair which does not have a std::pair
// as the return type of the dereferenced iterator..
template<typename Value>
class ConstKmerPairIteratorWrapper<typename KmerArrayPair<Value>::ConstIterator, Value> {
public:
	typedef const typename KmerArrayPair<Value>::ConstIterator Base;
	typedef KmerElementPair<Value> KEP;
	typedef typename KEP::ConstBaseRef BaseRef;
	typedef std::auto_ptr< BaseRef > BaseRefPtr;
	typedef Kmer& key_type;
	typedef KmerElementPair<Value> value_type;
	typedef Value mapped_type;

	ConstKmerPairIteratorWrapper() : _iter() {}
	ConstKmerPairIteratorWrapper(const Base &iter): _iter(iter) {	}
	ConstKmerPairIteratorWrapper(const ConstKmerPairIteratorWrapper &copy): _iter(copy._iter) {	}
	~ConstKmerPairIteratorWrapper() {
	}
	ConstKmerPairIteratorWrapper &operator=(const Base &copy) {
		_iter = copy;
		brptr.reset();
		return *this;
	}
	ConstKmerPairIteratorWrapper &operator=(const ConstKmerPairIteratorWrapper &copy) {
		_iter = copy._iter;
		brptr.reset();
		return *this;
	}

	// cast operators
	operator Base&() {
		return _iter;
	}
	operator const Base&() const {
		return _iter;
	}
	bool operator==(const Base &other) const { return _iter == other; }
	bool operator!=(const Base &other) const { return _iter != other; }
	ConstKmerPairIteratorWrapper &operator++() { _iter++; brptr.reset(); return *this; }
	ConstKmerPairIteratorWrapper operator++(int unused) { ConstKmerPairIteratorWrapper tmp(*this); ++(*this); return tmp; }
	BaseRef &operator*() { setRefs(); return *brptr; }
	BaseRef *operator->() { setRefs(); return brptr.get(); }

	inline const Kmer &key() {
		return _iter.key();
	}
	const Value &value() {
		return _iter.value();
	}
private:
	Base _iter;
	// WARNING: only instantiate when required to dereference
	mutable BaseRefPtr brptr;
	void setRefs() {
		if (brptr.get() == NULL)
			brptr.reset(new BaseRef(key(), value()));
	}

};

template<typename _KeyType, typename _ValueType, typename _BucketType, typename _Hasher>
class BucketExposedMapLogic {
public:
	typedef _KeyType KeyType;
	typedef _ValueType ValueType;
	typedef _BucketType BucketType;
	typedef _Hasher HasherType;

	typedef Kmer::NumberType    NumberType;
	typedef Kmer::IndexType     IndexType;
	typedef Kmer::SizeType      SizeType;

	typedef typename BucketType::iterator BaseBucketTypeIterator;
	typedef typename BucketType::const_iterator ConstBaseBucketTypeIterator;
	typedef typename BucketType::value_type BaseBucketTypeElementType;

	typedef std::vector< BucketType > BucketsVector;
	typedef BucketsVector BaseBucketsVector;
	typedef typename BucketsVector::iterator BaseBucketsVectorIterator;
	typedef typename BucketsVector::const_iterator ConstBaseBucketsVectorIterator;

	typedef IteratorOfMapIterators< BaseBucketsVectorIterator > IofMI;
	typedef KmerPairIteratorWrapper< IofMI, ValueType> Iterator;
	typedef KmerPairIteratorWrapper< BaseBucketsVectorIterator, ValueType > BucketsVectorIterator;
	typedef KmerPairIteratorWrapper< BaseBucketTypeIterator, ValueType > BucketTypeIterator;

	typedef IteratorOfMapIterators< ConstBaseBucketsVectorIterator > ConstIofMI;
	typedef ConstKmerPairIteratorWrapper< ConstIofMI, ValueType> ConstIterator;
	typedef ConstKmerPairIteratorWrapper< ConstBaseBucketsVectorIterator, ValueType > ConstBucketsVectorIterator;
	typedef ConstKmerPairIteratorWrapper< ConstBaseBucketTypeIterator, ValueType > ConstBucketTypeIterator;

	// for STL compatibility
	typedef Iterator iterator;
	typedef ConstIterator const_iterator;
	typedef typename Iterator::key_type key_type;
	typedef typename Iterator::value_type value_type;
	typedef typename Iterator::mapped_type mapped_type;

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

	BucketExposedMapLogic(int numBuckets = 0, int elementsPerBucket = 0) {
		// ensure buckets are a precise powers of two
		// with at least bucketCount buckets
		resizeBuckets(numBuckets, elementsPerBucket);
	}
	BucketExposedMapLogic &operator=(const BucketExposedMapLogic &copy) {
		_buckets = copy._buckets;
		BUCKET_MASK = copy.BUCKET_MASK;
		return *this;
	}
	void resizeBuckets(int bucketCount, int elementsPerBucket) {
		if (bucketCount > MAX_KMER_MAP_BUCKETS)
			bucketCount = MAX_KMER_MAP_BUCKETS;
		NumberType powerOf2 = getMinPowerOf2(bucketCount);
        BUCKET_MASK = powerOf2 - 1;
        _buckets.resize(powerOf2);
        if (elementsPerBucket > 16) {
        	for(int i = 0; i < (int) powerOf2; i++)
        		_buckets[i].reserve(elementsPerBucket);
        } // else lazy allocate the buckets.
		LOG_VERBOSE_OPTIONAL(1, bucketCount > 0 && Logger::isMaster(), "BucketExposedMap::resizeBuckets(" << bucketCount << " ): " << powerOf2 << " of " << elementsPerBucket);
	}
	void clear(bool releaseMemory = true) { // release memory ignored by default
		if (releaseMemory) {
			_buckets.clear();
			resizeBuckets(0,0);
		} else {
			for(BaseBucketsVectorIterator it = _buckets.begin(); it != _buckets.end(); it++)
				it->clear();
		}
	}
	void reset(bool releaseMemory = true) { // release memory ignored by default
		for(size_t i=0; i< _buckets.size(); i++) {
			_buckets[i].clear();
		}
	}
	void swap(BucketExposedMapLogic &other) {
		_buckets.swap(other._buckets);
		std::swap(BUCKET_MASK, other.BUCKET_MASK);
	}

	NumberType &getBucketMask() {
		return BUCKET_MASK;
	}
	const NumberType getBucketMask() const {
		return BUCKET_MASK;
	}
	NumberType getNumBuckets() const {
		return _buckets.size();
	}

	// Optionally organize buckets by Local thread and Distributed Rank
	inline int getLocalThreadId(NumberType hash, int numThreads) const {
		// stripe across all buckets
		assert(numThreads > 0);
		if (getBuckets().size() > 0)
			return (hash & getBucketMask()) % numThreads;
		else
			return 0;
	}
	inline int getLocalThreadId(const KeyType &key, int numThreads) const {
		return getLocalThreadId(HasherType()(key), numThreads);
	}
	inline int getDistributedThreadId(NumberType hash, NumberType numDistributedThreads) const {
		// partition in contiguous blocks of 'global' buckets
		if (getBuckets().size() > 0)
			return numDistributedThreads * (hash & getBucketMask()) / getBuckets().size();
		else
			return 0;
	}
	inline int getDistributedThreadId(const KeyType &key, NumberType numDistributedThreads) const {
		return getDistributedThreadId(HasherType()(key), numDistributedThreads);
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
		NumberType hash = HasherType()(key);
		return getLocalThreadId(hash, localThreadId, numLocalThreads, distributedThreadId, numDistributedThreads);
	}

	inline void getThreadIds(const KeyType &key, int &localThreadId, int numLocalThreads, int &distributedThreadId, NumberType numDistributedThreads) const {
		NumberType hash = HasherType()(key);
		distributedThreadId = getDistributedThreadId(hash, numDistributedThreads);
		localThreadId = getLocalThreadId(hash, numLocalThreads);
	}

	// Exposed BucketType interface

	inline NumberType getBucketIdx(NumberType hash) const {
		NumberType bucketIdx = (hash & getBucketMask());
		return bucketIdx;
	}
	inline NumberType getBucketIdx(const KeyType &key) const {
		return getBucketIdx(HasherType()(key));
	}

	SizeType size() const {
		SizeType size = 0;

		long bucketSize = _buckets.size();
#pragma omp parallel for reduction(+:size) if(getBuckets().size()>1000000)
		for(long i = 0; i < bucketSize; i++) {
			size += _buckets[i].size();
		}
		return size;
	}

	IndexType maxBucketSize() const	{
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
	bool empty() const {
		for(ConstBaseBucketsVectorIterator it = _buckets.begin(); it != _buckets.end(); it++)
			if (! it->empty() )
				return false;
		return true;
	}

	inline BucketsVector &getBuckets() { return _buckets; }
	inline const BucketsVector &getBuckets() const { return _buckets; }
	inline BucketType &getBucketByIdx(IndexType index) {
		assert(index < _buckets.size());
		return _buckets[index];
	}
	inline const BucketType &getBucketByIdx(IndexType index) const {
		assert(index < _buckets.size());
	        return _buckets[index];
	}

	inline BucketType &getBucket(NumberType hash) {
		NumberType idx = getBucketIdx(hash);
		return getBucketByIdx(idx);
	}
	inline const BucketType &getBucket(NumberType hash) const {
		NumberType idx = getBucketIdx(hash);
		return getBucketByIdx(idx);
	}
	inline BucketType &getBucket(const KeyType &key) {
		NumberType idx = getBucketIdx(key);
		return getBucketByIdx(idx);
	}
	inline const BucketType &getBucket(const KeyType &key) const {
		NumberType idx = getBucketIdx(key);
		return getBucketByIdx(idx);
	}

	std::string toString() {
		return "TODO implement toString()";
	}


	Iterator begin(int rank = 0, int size = 1) { return Iterator( IofMI( _buckets.begin(), _buckets.end(), rank, size) ); }
	Iterator end() { return Iterator( IofMI::getEnd(_buckets) ); }

	ConstIterator begin(int rank = 0, int size = 1) const { return ConstIterator( IofMI( _buckets.begin(), _buckets.end(), rank, size) ); }
	ConstIterator end() const { return ConstIterator( IofMI::getEnd(_buckets) ); }

	Iterator beginThreaded(int rank = omp_get_thread_num(), int size = omp_get_num_threads()) {
		return begin(rank, size);
	}
	Iterator endThreaded() {
		return end();
	}

private:
	BucketsVector _buckets;
	NumberType BUCKET_MASK;
};
template<typename _KeyType, typename _ValueType, typename _BucketType, typename _Hasher>
class BucketExposedMap : public BucketExposedMapLogic<_KeyType, _ValueType, _BucketType, _Hasher> {
public:
	typedef BucketExposedMapLogic<_KeyType, _ValueType,_BucketType, _Hasher> BEML;
	typedef typename BEML::KeyType KeyType;
	typedef typename BEML::ValueType ValueType;
	typedef typename BEML::BucketType BucketType;
	typedef typename BEML::HasherType HasherType;

	typedef Kmer::NumberType    NumberType;
	typedef Kmer::IndexType     IndexType;
	typedef Kmer::SizeType      SizeType;

	typedef typename BEML::Iterator Iterator;
	typedef typename BEML::ConstIterator ConstIterator;
	typedef typename BEML::key_type key_type;
	typedef typename BEML::value_type value_type;
	typedef typename BEML::mapped_type mapped_type;

	typedef typename BEML::BaseBucketsVector BaseBucketsVector;
	typedef typename BEML::BaseBucketsVectorIterator BaseBucketsVectorIterator;

	typedef typename BEML::BaseBucketTypeIterator BaseBucketTypeIterator;
	typedef typename BEML::BaseBucketTypeElementType BaseBucketTypeElementType;

	typedef typename BEML::BucketTypeIterator BucketTypeIterator;
	typedef typename BEML::ConstBucketTypeIterator ConstBucketTypeIterator;

	typedef typename BEML::BucketsVector BucketsVector;
	typedef typename BEML::BucketsVectorIterator BucketsVectorIterator;

	typedef KmerElementPair<ValueType> BaseElementType;
	typedef BaseElementType ElementType;

	BEML::getMinPowerOf2;
	BEML::clear;
	BEML::reset;
	BEML::getNumBuckets;
	BEML::getBucketMask;
	BEML::getBucketIdx;
	BEML::getThreadIds;
	BEML::getLocalThreadId;
	BEML::getDistributedThreadId;

	BEML::getBuckets;
	BEML::getBucket;
	BEML::getBucketByIdx;

	BEML::begin;
	BEML::end;
	BEML::beginThreaded;
	BEML::endThreaded;



public:
	BucketExposedMap() : BEML() {
		clear();
	}
	BucketExposedMap(int numBuckets, int elementsPerBucket = 0) : BEML(numBuckets, elementsPerBucket) {
	}
	BucketExposedMap &operator=(const BucketExposedMap &other) {
		(BEML&)*this = (BEML&) other;
		return *this;
	}
	virtual ~BucketExposedMap() {
		clear();
	}

	void swap(BucketExposedMap &other) {
		BEML::swap(other);
	}
	// optimization to move the buckets with pre-allocated memory to the next DMP thread
	void rotateDMPBuffers(int numThreads) {
		IndexType block = getBuckets().size() / numThreads;
		size_t i = 0;
		// skip to the first non-zero bucket
		while (i < getBuckets().size() && getBuckets()[i].size() == 0)
			i++;
		for(size_t j = i; j < getBuckets().size() - 1 && j < i+block; j++) {
			getBuckets()[j].swap(getBuckets()[ (j+block) % getBuckets().size() ]);
		}
	}

	// insert/add, remove, exists methods
	ElementType insert(const KeyType &key, const ValueType &value) {
		return insert(key, value, getBucket(key));
	}
	ElementType insert(const KeyType &key, const ValueType &value, NumberType hash) {
		assert(HasherType()(key) == hash);
		return insert(key,value, getBucket(hash));
	}

	bool remove(const KeyType &key) {
		return remove(key, getBucket(key));
	}
	bool remove(const KeyType &key, NumberType hash) {
		assert(HasherType()(key) == hash);
		return remove(key, getBucket(hash));
	}

	Iterator find(const KeyType &key) {
		return find(key, getBucket(key));
	}
	Iterator find(const KeyType &key, NumberType hash) {
		assert(HasherType()(key) == hash);
		return find(key, getBuckets(hash));
	}

	ValueType &operator[](const KeyType &key) {
		NumberType hash = HasherType()(key);
		BucketType &bucket = getBucket(hash);
		BucketTypeIterator it = bucket.find(key);
		if (it != bucket.end()) {
			return it->value();
		} else {
			ValueType value;
			return insert(key, value, bucket).value();
		}
	}

	bool exists(const KeyType &key) const {
		return exists(key, getBucket(key));
	}
	bool exists(const KeyType &key, NumberType hash) const {
		assert(HasherType()(key) == hash);
		return exists(key, getBuckets(hash));
	}

	const bool getValueIfExists(const KeyType &key, ValueType &value) const {
		return getValueIfExists(key, value, getBucket(key));
	}

	const bool getValueIfExists(const KeyType &key, ValueType &value, NumberType hash) const {
		assert(HasherType()(key) == hash);
		return getValueIfExists(key, value, getBucket(hash));
	}

	const ElementType getElementIfExists(const KeyType &key) const {
		return getElementIfExists(key, getBucket(key));
	}
	ElementType getElementIfExists(const KeyType &key) {
		return getElementIfExists(key, getBucket(key));
	}

	const ElementType getElementIfExists(const KeyType &key, NumberType hash) const {
		assert(HasherType()(key) == hash);
		return getElementIfExists(key, getBucket(hash));
	}
	ElementType getElementIfExists(const KeyType &key, NumberType hash) {
		assert(HasherType()(key) == hash);
		return getElementIfExists(key, getBucket(hash));
	}

	ElementType getOrSetElement(const KeyType &key, ValueType value) {
		return getOrSetElement(key, value, getBucket(key));
	}
	ElementType getElement(const KeyType &key) {
		ValueType value = ValueType();
		return getOrSetElement(key, value);
	}
	ElementType getOrSetElement(const KeyType &key, ValueType value, NumberType hash) {
		assert(HasherType()(key) == hash);
		return getOrSetElement(key, value, getBucket(hash));
	}
	ElementType getElement(const KeyType &key, NumberType hash) {
		assert(HasherType()(key) == hash);
		ValueType value = ValueType();
		return getOrSetElement(key, value);
	}


	std::string toString() const {
		std::stringstream ss;
		ss << this << "[";
		IndexType idx=0;
		for(; idx<getBuckets().size() && idx < 30; idx++) {
			ss << "bucket:" << idx << " (" << getBuckets()[idx].size() << "), ";
		}
		if (idx < getBuckets().size())
			ss << " ... " << getBuckets().size() - idx << " more ";
		ss << "]";
		return ss.str();
	}

	virtual void optimize() {}

	bool exists(const KeyType &key, const BucketType &bucket) const {
		return bucket.find(key) != bucket.end();
	}
	const bool getValueIfExists(const KeyType &key, ValueType &value, const BucketType &bucket) const {
		ConstBucketTypeIterator it = bucket.find(key);
		if (it == bucket.end())
			return false;
		else {
			value = it->value();
			return true;
		}
	}
	const ElementType getElementIfExists(const KeyType &key, const BucketType &bucket) const {
		ConstBucketTypeIterator it = bucket.find(key);
		if (it == bucket.end())
			return ElementType::getConst();
		else {
			return ElementType::getConst(it->key(), it->value());
		}
	}

	// optimized merge for DMP threaded (i.e. blocked where only one bucket per map is populated)
	void mergeStripedBuckets(BucketExposedMap &src) {
		if (getNumBuckets() != src.getNumBuckets()) {
			LOG_THROW("Invalid: Can not merge two BucketExposedMaps of differing sizes!");
		}

		long bucketsSize = getNumBuckets();
		#pragma omp parallel for
		for(long idx = 0 ; idx < bucketsSize; idx++) {
			BucketType &a = getBucketByIdx(idx);
			BucketType &b = src.getBucketByIdx(idx);

			if (! mergeTriviallyInterleavedBuckets(a,b) ) {
				LOG_THROW("Invalid: Expected one bucket to be zero in this optimized method: KmerMapByKmerArrayPair::merge(const KmerMapByKmerArrayPair &src) const");
			}
		}
	}

	void mergeAdd(BucketExposedMap &src) {
		if (getNumBuckets() != src.getNumBuckets()) {
			LOG_THROW("Invalid: Can not merge two BucketExposedMaps of differing sizes!");
		}
		BucketType merged;

		long bucketsSize = getNumBuckets();

		#pragma omp parallel for private(merged)
		for(long idx = 0 ; idx < bucketsSize; idx++) {
			// buckets are sorted so perform a sorted merge by bucket
			BucketType &a = getBucketByIdx(idx);
			BucketType &b = src.getBucketByIdx(idx);

			// short circuit if one of the buckets is empty
			if (mergeTriviallyInterleavedBuckets(a,b))
				continue;

			for(Iterator ib = b.begin(); ib != b.end(); ib++) {
				Iterator ia = a.find(ib->key());
				if (ia == a.end())
					a.insert(ib->key(), ib->value());
				else
					ia->value() = ia->value + ib->value();
			}
			b.reset(true);
		}
		src.clear();
	}

	template< typename OtherBucketExposedMap >
	void mergePromote(BucketExposedMap &src, OtherBucketExposedMap &mergeDest) {
		// TODO fixme
		LOG_THROW("Invalid: This method is broken for threaded execution somehow (even with pragmas disabled)");
		if (getNumBuckets() != src.getNumBuckets()) {
			LOG_THROW("Invalid: Can not merge two KmerMapByKmerArrayPair of differing sizes!");
		}
		/*
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
			BucketType &a = getBuckets()[idx];
			BucketType &b = src.getBuckets()[idx];
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
		*/
	}


protected:
	static bool mergeTriviallyInterleavedBuckets(BucketType &a, BucketType &b) {
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




protected:

	ElementType insert(const KeyType &key, ValueType value, BucketType &bucket) {
		BucketTypeIterator it = bucket.find(key);
		if (it == bucket.end()) {
			BaseBucketTypeElementType elem(key,value);
			it = bucket.insert((BaseBucketTypeIterator)it, elem);
		}
		return ElementType(it->key(), it->value());
	}

	bool remove(const KeyType &key, BucketType &bucket) {
		BucketTypeIterator it = bucket.find(key);
		if (it != bucket.end()) {
			bucket.erase((BaseBucketTypeIterator)it);
			return true;
		} else {
			return false;
		}
	}
	bool exists(const KeyType &key, BucketType &bucket) {
		return bucket.find(key) != bucket.end();
	}
	Iterator find(const KeyType &key, BucketType &bucket) {
		return bucket.find(key);
	}

	ElementType getElementIfExists(const KeyType &key, BucketType &bucket) {
		BucketTypeIterator it = bucket.find(key);
		if (it == bucket.end())
			return ElementType();
		else {
			return ElementType(it->key(), it->value());
		}
	}

	ElementType getOrSetElement(const KeyType &key, ValueType value, BucketType &bucket) {
		BucketTypeIterator it = bucket.find(key);
		if (it == bucket.end())
			return insert(key, value, bucket);
		else
			return ElementType(it->key(), it->value());
	}
	ElementType getElement(const KeyType &key, BucketType &bucket) {
		ValueType value = ValueType();
		return getOrSetElement(key, bucket, value);
	}

protected:
	BucketsVector _buckets;
	NumberType BUCKET_MASK;

private:
	inline const BucketExposedMap &_constThis() const {return *this;}

};

template<typename Value>
class KmerMapByKmerArrayPair : public BucketExposedMap<Kmer, Value, KmerArrayPair<Value>, KmerHasher > {
public:

	// typedef propagation

	typedef KmerArrayPair<Value> KAP;
	typedef BucketExposedMap<Kmer, Value, KAP, KmerHasher > Base;
	typedef typename Kmer::NumberType    NumberType;
	typedef typename Kmer::IndexType     IndexType;
	typedef typename Kmer::SizeType      SizeType;

	typedef typename Base::KeyType KeyType;
	typedef typename Base::ValueType ValueType;
	typedef typename Base::BucketType BucketType;
	typedef typename Base::HasherType HasherType;
	typedef	typename Base::BucketTypeIterator BucketTypeIterator;
	typedef typename Base::BaseElementType BaseElementType;
	typedef typename Base::ElementType ElementType;

	typedef typename Base::BaseBucketsVector BaseBucketsVector;
	typedef typename Base::BaseBucketsVectorIterator BaseBucketsVectorIterator;

	typedef typename Base::BucketsVector BucketsVector;
	typedef typename Base::BucketsVectorIterator BucketsVectorIterator;

	typedef typename Base::BaseBucketTypeIterator BaseBucketTypeIterator;

	typedef NumberType * NumberTypePtr;

private:
	inline const KmerMapByKmerArrayPair &_constThis() const {return *this;}

public:
	KmerMapByKmerArrayPair() : Base() {
		_isSorted = defaultSort;
		KAP::initTmp();
	}
	KmerMapByKmerArrayPair(long estimatedRawKmers) : Base(estimatedRawKmers / KmerBaseOptions::getOptions().getKmersPerBucket() + 1, KmerBaseOptions::getOptions().getKmersPerBucket()) {
		_isSorted = defaultSort;
		KAP::initTmp();
	}
	virtual ~KmerMapByKmerArrayPair()
	{
		clear();
		KAP::clearTmp();
	}
	KmerMapByKmerArrayPair &operator=(const KmerMapByKmerArrayPair &other) {
		*((Base*)this) = (const Base&) other;
		_isSorted = other._isSorted;
		return *this;
	}
	Base::getBuckets;
	Base::getBucketMask;
	Base::getBucket;
	Base::getBucketByIdx;
	Base::insert;
	Base::remove;

	// Base LocalThread / DistributedRank methods
	Base::getLocalThreadId;
	Base::getDistributedThreadId;
	Base::getThreadIds;
	Base::getBucketIdx;
	Base::getNumBuckets;


	// specific overloading, taking advantage of KmerArrayPair memory release directives
	void reset(bool releaseMemory = true) {
		for(size_t i=0; i< getBuckets().size(); i++) {
			getBuckets()[i].reset(releaseMemory);
		}
	}
	void clear(bool releaseMemory = true) {
		reset(releaseMemory);
		if (releaseMemory) {
			this->resizeBuckets(0,0);
		}
	}

	void swap(KmerMapByKmerArrayPair &other) {
		Base::swap((Base&) other);
		std::swap(_isSorted, other._isSorted);
	}

public:
	class Iterator : public std::iterator<std::forward_iterator_tag, KmerMapByKmerArrayPair>
	{
		// Provides an iterator that can optionaly partition the buckets by
		// into rank of size slices, where each rank will interleave responsibility
		// for a given bucket by modulus of the size.
		friend class KmerMapByKmerArrayPair;
	private:
		const KmerMapByKmerArrayPair *_target;
		BaseBucketsVectorIterator _iBucket;
		BaseBucketTypeIterator _iElement;
		int _rank, _size;

	public:
		Iterator() : _target(NULL), _rank(0), _size(1) {}

		// iterator over rank/size will stripe across the buckets (modulus by size).
		Iterator(KmerMapByKmerArrayPair *target, int rank = 0, int size = 1):
			_target(target),
			_iBucket(target->getBuckets().begin()),
			_iElement(),
			_rank(rank), _size(size)
		{
			if (!isEnd())
				_iElement = _iBucket->begin();
			_moveToRank();
			_moveToNextValidElement();
		}

		Iterator(const Iterator &copy) {
			*this = copy;
		}

		Iterator(KmerMapByKmerArrayPair *target, BucketsVectorIterator bucketPtr, int rank = 0, int size = 1):
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

		Iterator(KmerMapByKmerArrayPair *target, BucketsVectorIterator bucketPtr, BucketTypeIterator elementPtr, int rank = 0, int size = 1):
			_target(target),
			_iBucket(bucketPtr),
			_iElement(elementPtr),
			_rank(rank), _size(size)
		{
			_moveToRank();
		}

		Iterator &operator=(const Iterator &copy) {
			_target = copy._target;
			_iBucket = copy._iBucket;
			_iElement = copy._iElement;
			_rank = copy._rank;
			_size = copy._size;
			return *this;
		}

	private:
		inline bool isEnd() const {
			return _iBucket == _target->getBuckets().end();
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
			if (isEnd() || ((int) _target->getBuckets().size() <= _rank)) {
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

		bool operator==(const Iterator& other) const {
			return _iBucket == other._iBucket
					&& (isEnd() || _iElement == other._iElement);
		}

		bool operator!=(const Iterator& other) const {
			return !(*this == other);
		}

		Iterator& operator++() {
			++_iElement;
			_moveToNextValidElement();
			return *this;
		}

		Iterator operator++(int unused) {
			Iterator tmp(*this); ++(*this); return tmp;
		}

		ElementType &operator*() {return *_iElement;}
		const ElementType &operator*() const {return *_iElement;}
		ElementType *operator->() {return &(*_iElement);}
		const ElementType *operator->() const {return &(*_iElement);}

		const Kmer &key() const {return _iElement.key();}
		const Kmer &key() {return _iElement.key();}
		ValueType &value() { return _iElement.value();}
		const ValueType &value() const {return _iElement.value();}

		BucketType &bucket() {return *_iBucket;}
		const BucketType &bucket() const {return *_iBucket;}

		IndexType bucketIndex() {return (_iBucket - _target->getBuckets().begin());}
		const IndexType bucketIndex() const {return (_iBucket - _target->getBuckets().begin());}

	};

	class ConstIterator : public std::iterator<std::forward_iterator_tag, KmerMapByKmerArrayPair>
	{
	private:
		Iterator _iterator;
	public:
		ConstIterator(KmerMapByKmerArrayPair *target, int rank = 0, int size = 1) : _iterator(target, rank, size) {}
		ConstIterator(const ConstIterator &copy) : _iterator(copy._iterator) {}
		ConstIterator(KmerMapByKmerArrayPair *target, BucketsVectorIterator bucketPtr, int rank = 0, int size = 1) : _iterator(target,bucketPtr, rank, size) {}
		ConstIterator(KmerMapByKmerArrayPair *target, BucketsVectorIterator bucketPtr,BucketTypeIterator elementPtr, int rank = 0, int size = 1): _iterator(target,bucketPtr,elementPtr, rank, size) {}
		ConstIterator &operator=(const ConstIterator &copy) { _iterator = copy._iterator; return *this; }
		bool operator==(const ConstIterator& other) const {return _iterator == other._iterator;}
		bool operator!=(const ConstIterator& other) const {return _iterator != other._iterator;}
		ConstIterator& operator++() {++_iterator; return *this;}
		ConstIterator operator++(int unused) {return ConstIterator(_iterator++);}
		const ElementType &operator*() const {return _iterator.operator*();}
		const ElementType *operator->() const {return _iterator.operator->();}
		const Kmer &key() const {return _iterator.key();}
		const ValueType &value() const {return _iterator.value();}
		const BucketType &bucket() const {return _iterator.bucket();}
		const IndexType bucketIndex() const {return _iterator.bucketIndex();}
	};

	typedef Iterator iterator;
	typedef ConstIterator const_iterator;

	Iterator begin(int rank = 0, int size = 1) {return Iterator(this, rank, size);}
	Iterator end() {return Iterator(this,getBuckets().end());}

	ConstIterator begin(int rank = 0, int size = 1) const {return Iterator(this, rank, size);}
	ConstIterator end() const {return Iterator(this, getBuckets().end());}

	Iterator beginThreaded(int rank = omp_get_thread_num(), int size = omp_get_num_threads()) {
		return begin(rank, size);
	}
	Iterator endThreaded() {
		return end();
	}

	// other methods
	Base::size;

	// Sorted KmerArrayPair specializations
public:
	static const bool defaultSort = false;

private:
	bool _isSorted;

public:
	// sorting buckets specificity
	inline bool isSorted() const {
		return _isSorted;
	}
	void optimize() {
		resort();
	}
	void resort() {
		if (!_isSorted) {
			LOG_DEBUG_OPTIONAL(1, Logger::isMaster(), "Sorting KmerMapByKmerArrayPair");
#pragma omp parallel for
			for(long i=0; i< (long) getBuckets().size(); i++) {
				getBuckets()[i].resort();
			}
			_isSorted = true;
		}
	}
	void _setIsSorted(bool sorted) {
		_isSorted = sorted;
	}

protected:
	// insert specialiazation, using sorted KmerArrayPair
	ElementType insert(const KeyType &key, const ValueType &value, BucketType &bucket) {
		IndexType idx;
		if (isSorted()) {
			idx = bucket.insertSorted(key,value);
		} else {
			idx = bucket.append(key,value);
			if (bucket.getNumUnsorted() >= bucket.MAX_UNSORTED || bucket.size() == bucket.capacity()) {
				bucket.resort( bucket.capacity() - bucket.size() < bucket.MIN_GROWTH / 2 );
				idx = bucket.findIndex(key);
				assert(idx != BucketType::MAX_INDEX);
			}
		}

		ElementType element = bucket.getElement(idx);
		return element;
	}

	// remove specialization, using sorted KmerArrayPair
	bool remove(const KeyType &key, BucketType &bucket) {
		bool isFound;
		IndexType idx = isSorted() ? bucket.findSortedIndex(key, isFound) : bucket.findIndex(key, isFound);
		if (isFound && idx != BucketType::MAX_INDEX)
			bucket.remove(idx);
		return isFound;
	}

	// store/restore specialiazation
public:
	// restore new instance from mmap
	KmerMapByKmerArrayPair(const void *src) {
		NumberType size(0), *offsetArray;
		const void *ptr;
		_getMmapSizes(src, size, getBucketMask(), offsetArray);
		getBuckets().resize(size);
		_isSorted = true;
		for (NumberType idx = 0 ; idx < size ; idx++) {
			ptr = ((char*)src) + *(offsetArray + idx);
			getBuckets()[idx] = BucketType(ptr);
			_isSorted &= getBuckets()[idx].isSorted();
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
		NumberType size = (NumberType) getBuckets().size();
		NumberType *numbers = (NumberType *) dst;
		*(numbers++) = size;
		*(numbers++) = getBucketMask();
		NumberType *offsetArray = numbers;
		NumberType offset = sizeof(NumberType) * (2+size);

		for(NumberType idx = 0 ; idx < size; idx++) {
			*(offsetArray++) = offset;
			void *ptr = ((char*)dst) + offset;
			const char *newPtr = (const char *) getBuckets()[idx].store(ptr);
			NumberType newSize = (newPtr - (const char *) ptr);
			offset += newSize;
		}
		return ((char*)dst) + offset;
	}
	static const KmerMapByKmerArrayPair restore(const void *src) {
		KmerMapByKmerArrayPair map;
		NumberType size(0), *offsetArray;
		const void *ptr;
		_getMmapSizes(src, size, map.getBucketMask(), offsetArray);
		map.getBuckets().resize(size);
		map._isSorted = true;
		for (NumberType idx = 0 ; idx < size ; idx++) {
			ptr = ((char*)src) + *(offsetArray + idx);
			map.getBuckets()[idx] = BucketType::restore(ptr);
			map._isSorted &= map.getBuckets()[idx].isSorted();
		}
		return map;
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
		return sizeof(NumberType)*(2+getBuckets().size())
				+ getBuckets().size()*sizeof(IndexType);
	}
	SizeType getSizeToStoreElements() const {
		return size()*(KmerSizer::getByteSize() + sizeof(ValueType));
	}

	std::string toStringDetailed() const {
		std::stringstream ss;
		ss << this << "[";
		IndexType idx=0;
		for(; idx<getBuckets().size() && idx < 30; idx++) {
			ss << "bucket:" << idx << ' ' << getBuckets()[idx].toString() << ", ";
		}
		if (idx < getBuckets().size())
			ss << " ... " << getBuckets().size() - idx << " more ";
		ss << "]";
		return ss.str();
	}

	Base::mergeTriviallyInterleavedBuckets;

	// specialized merge for KmerArrayPair
	void mergeAdd(KmerMapByKmerArrayPair &src) {
		if (getNumBuckets() != src.getNumBuckets()) {
			LOG_THROW("Invalid: Can not merge two BucketExposedMaps of differing sizes!");
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
			BucketType &a = getBucketByIdx(idx);
			BucketType &b = src.getBucketByIdx(idx);

			// short circuit if one of the buckets is empty
			if (mergeTriviallyInterleavedBuckets(a,b))
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


};



template<typename Value, typename BucketType>
class KmerMapBySTLMap : public BucketExposedMap< KmerInstance, Value, BucketType, KmerHasher >{
public:
	typedef Kmer::NumberType    NumberType;
	typedef Kmer::IndexType     IndexType;
	typedef Kmer::SizeType      SizeType;

	typedef BucketExposedMap< KmerInstance, Value, BucketType, KmerHasher > Base;
	typedef Value ValueType;
	//typedef	typename BucketType::iterator BucketTypeIterator;
	//typedef typename BucketType::value_type BucketElementType;
	//typedef std::vector< BucketType > BucketsVector;
	//typedef typename BucketsVector::iterator BucketsVectorIterator;
	//typedef typename BucketsVector::const_iterator ConstBucketsVectorIterator;
	typedef typename Base::Iterator Iterator;
	//typedef typename Iterator::value_type ElementType;

	Base::getBuckets;
	Base::getBucketByIdx;
	KmerMapBySTLMap() : Base() {}
	KmerMapBySTLMap(long estimatedRawKmers) : Base(omp_get_max_threads(), estimatedRawKmers / omp_get_max_threads()) {
	}
	KmerMapBySTLMap(const void *src) {
		LOG_THROW("Unimplemented restore");
	}
	KmerMapBySTLMap(const KmerMapBySTLMap &copy) {
		*this = copy;
	}
	KmerMapBySTLMap &operator=(const KmerMapBySTLMap &copy) {
		getBuckets().resize(copy.getBuckets().size());
		for(int i = 0; i < (int) getBuckets().size(); i++) {
			getBucketByIdx(i).clear();
			getBucketByIdx(i).insert(copy.getBucketByIdx(i).begin(), copy.getBucketByIdx(i).end());
		}
		return *this;
	}

};

// For Example:
// typedef KmerMapBySTLMap<DataType, boost::unordered_map<KmerInstance, DataType, KmerHasher> > KmerMapBoostType;
//
#include <boost/unordered_map.hpp>
template<typename Value>
class KmerMapBoost : public KmerMapBySTLMap<Value, boost::unordered_map<KmerInstance, Value, KmerHasher> > {
public:
	typedef KmerMapBySTLMap<Value, boost::unordered_map<KmerInstance, Value, KmerHasher> > Base;
	KmerMapBoost() : Base() {}
	KmerMapBoost(long estimatedRawKmers) : Base(4*estimatedRawKmers) {} // allocate 4x the buckets we think we need...
	KmerMapBoost(const void *src) : Base(src) {}
	KmerMapBoost(const Base &copy) : Base( (const Base&) copy ) {}
	KmerMapBoost &operator=(const Base &other) {
		*((Base*)this) = other;
		return *this;
	}
	typedef typename Base::BucketType BucketType;
};

#include <sparsehash/sparse_hash_map>
template<typename Value>
class GSHWrapper : public google::sparse_hash_map<KmerInstance, Value, KmerHasher> {
public:
	typedef google::sparse_hash_map<KmerInstance, Value, KmerHasher> Base;
	typedef typename Base::mapped_type mapped_type;
	typedef typename Base::iterator iterator;
	typedef typename Base::size_type size_type;
	GSHWrapper(size_type n = 0): Base(n) {}
	// GSH has no reserve method, but resize acts the same
	void reserve(size_type n) { this->resize(n); }
	// erase is not supported by GSH without a deleted_key...
	void erase(iterator pos) { pos->second = mapped_type(); }

};
template<typename Value>
class KmerMapGoogleSparse : public KmerMapBySTLMap<Value, GSHWrapper<Value> > {
public:
	typedef KmerMapBySTLMap<Value, GSHWrapper<Value> > Base;
	KmerMapGoogleSparse(): Base() {}
	KmerMapGoogleSparse(long estimatedRawKmers): Base(8*estimatedRawKmers) {} // allocate 8x the buckets we think we need
	KmerMapGoogleSparse(const void *src) : Base(src) {}
	KmerMapGoogleSparse(const Base &copy) : Base( (const Base&) copy ) {}
	KmerMapGoogleSparse &operator=(const Base &other) {
		*((Base*)this) = other;
		return *this;
	}
	typedef typename Base::BucketType BucketType;
};


template<typename Value>
class KmerMap : public KmerMapByKmerArrayPair<Value> {
//class KmerMap : public KmerMapBySTLMap<Value, GSHWrapper<Value> > {
public:
	typedef KmerMapByKmerArrayPair<Value> Base;
//	typedef KmerMapBySTLMap<Value, GSHWrapper<Value> > Base;
	KmerMap(): Base() {}
	KmerMap(long estimatedRawKmers) : Base(estimatedRawKmers) {}
	KmerMap(const void *src) : Base(src) {}
	KmerMap(const Base &copy) : Base( (const Base&) copy ) {}
	KmerMap &operator=(const Base &other) {
		*((Base*)this) = other;
		return *this;
	}
	typedef typename Base::BucketType BucketType;
};

typedef KmerArrayPair<char> Kmers;
typedef KmerArrayPair<double> KmerWeights;
typedef KmerArrayPair<Kmernator::UI32> KmerCounts;
typedef KmerArrayPair< WeightedExtensionMessagePacket > KmerWeightedExtensions;

#endif

