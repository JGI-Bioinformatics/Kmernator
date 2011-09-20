//
// Kmernator/src/KmerTrackingData.h
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

#ifndef _KMER_TRACKING_DATA_H
#define _KMER_TRACKING_DATA_H

#include "config.h"
#include "Log.h"

class Extension {
public:
	enum ExtensionType {A = 0, C, G, T, N, X, MAX_EXTENSIONS};

	static char getBase(ExtensionType e) {
		static const char ExtensionToBase[MAX_EXTENSIONS] = {'A', 'C', 'G', 'T', 'N', 'X'};
		return ExtensionToBase[e];
	}
	Extension() : _base(MAX_EXTENSIONS), _quality(0) {}
	Extension(char c, int quality) : _base(MAX_EXTENSIONS), _quality((unsigned char) quality) {
		switch(c) {
		case 'A' : case 'a' : _base = A; break;
		case 'C' : case 'c' : _base = C; break;
		case 'G' : case 'g' : _base = G; break;
		case 'T' : case 't' : _base = T; break;
		case 'N' : case 'n' : _base = N; break;
		case 'X' : case 'x' : _base = X; break;
		default : LOG_THROW("Invalid Extension type: " << c << " q: " << quality);
		}
		if (quality > 255 || quality < 0)
			LOG_THROW("Invalid Extension quality: " << quality << " type: " << c);
	}
	bool isValid() const {
		return _base < MAX_EXTENSIONS;
	}
	bool isBase() const {
		return _base <= T;
	}
	ExtensionType getExtension() const {
		return _base;
	}
	char getBase() const {
		return getBase(_base);
	}
	unsigned char getQuality() const {
		return _quality;
	}
	Extension getReverseComplement() const {
		Extension revComp = *this;
		if (_base ==  A)
			revComp._base = T;
		else if (_base == C)
			revComp._base = G;
		else if (_base == G)
			revComp._base = C;
		else if (_base == T)
			revComp._base = A;

		return revComp;
	}

private:
	ExtensionType _base;
	unsigned char _quality;
};


class ExtensionTracking {
public:
	enum DirectionType {Left = 0, Right, MAX_DIRECTIONS};

	static unsigned char &getMinQuality() {
		static unsigned char _minQuality = 20;
		return _minQuality;
	}
	static void setMinQuality(unsigned char minQuality) {
		getMinQuality() = minQuality;
	}

	ExtensionTracking() {
		for(int j = 0; j < MAX_DIRECTIONS; j++)
			for(int i = 0; i < Extension::MAX_EXTENSIONS; i++)
				_extensionCounts[j][i] = 0;
	}
	ExtensionTracking(const ExtensionTracking &copy) {
		*this = copy;
	}
	ExtensionTracking &operator=(const ExtensionTracking &copy) {
		for(int j = 0; j < MAX_DIRECTIONS; j++)
			for(int i = 0; i < Extension::MAX_EXTENSIONS; i++)
				_extensionCounts[j][i] = copy._extensionCounts[j][i];
		return *this;
	}
	~ExtensionTracking() {}

	ExtensionTracking &operator+(const ExtensionTracking &other) {
		for(int j = 0; j < MAX_DIRECTIONS; j++)
			for(int i = 0; i < Extension::MAX_EXTENSIONS; i++)
				_extensionCounts[j][i] += other._extensionCounts[j][i];
		return *this;
	}

	void trackExtensions(Extension left, Extension right) {
		trackExtension(left, Left);
		trackExtension(right, Right);
	}
	void trackExtension(Extension ext, DirectionType dir) {
		if (ext.isValid() && ( ext.getQuality() >= getMinQuality() || !ext.isBase())) {
			getExtensionCount(ext.getExtension(), dir)++;
		}
	}
	unsigned int &getExtensionCount(Extension::ExtensionType ext, DirectionType dir) {
		assert(dir < MAX_DIRECTIONS);
		return _extensionCounts[dir][ext];
	}

private:
	unsigned int _extensionCounts[MAX_DIRECTIONS][Extension::MAX_EXTENSIONS];
};

class ExtensionMessagePacket {
public:
	ExtensionMessagePacket() : _leftB('\0'), _rightB('\0'), _leftQ(0), _rightQ(0) {}
	ExtensionMessagePacket(const ExtensionMessagePacket &copy) {
		*this = copy;
	}
	~ExtensionMessagePacket() {}
	ExtensionMessagePacket &operator=(const ExtensionMessagePacket &copy) {
		_leftB = copy._leftB;
		_rightB = copy._rightB;
		_leftQ = copy._leftQ;
		_rightQ = copy._rightQ;
		return *this;
	}

	Extension getLeft() const {
		return Extension(_leftB, _leftQ);
	}
	Extension getRight() const {
		return Extension(_rightB, _rightQ);
	}
	void setExtensions(Extension left, Extension right) {
		_leftB = left.getBase();
		_leftQ = left.getQuality();
		_rightB = right.getBase();
		_rightQ = right.getQuality();
	}
	void setExtensions(ExtensionTracking &extTrack) {
		for(int dir = 0; dir < ExtensionTracking::MAX_DIRECTIONS; dir++) {
			for(int ext = 0; ext < Extension::MAX_EXTENSIONS; ext++) {
				Extension::ExtensionType e = (Extension::ExtensionType) ext;
				ExtensionTracking::DirectionType d = (ExtensionTracking::DirectionType) dir;
				if (extTrack.getExtensionCount(e, d) > 0) {
					if (d == ExtensionTracking::Left) {
						_leftB = Extension::getBase(e);
					    _leftQ = ExtensionTracking::getMinQuality();
					} else {
						_rightB = Extension::getBase(e);
					    _rightQ = ExtensionTracking::getMinQuality();
					}
				}
			}
		}
	}
protected:
	char _leftB, _rightB;
	unsigned char _leftQ, _rightQ;

};

class ExtensionTrackingInterface {
public:
	void trackExtensions(Extension left, Extension right) {}
	ExtensionTracking getExtensionTracking() const { return ExtensionTracking(); }
};

class TrackingDataSingleton;
class TrackingDataWithAllReads;

class TrackingData : public ExtensionTrackingInterface {
public:
	typedef Kmernator::UI16 CountType;
	typedef Kmernator::ReadSetSizeType ReadIdType;
	typedef Kmernator::SequenceLengthType PositionType;
	typedef float WeightType;

	static const CountType    MAX_COUNT    = MAX_UI16;
	static const ReadIdType   MAX_READ_ID  = MAX_READ_SET_SIZE;
	static const PositionType MAX_POSITION = MAX_SEQUENCE_LENGTH;

	class ReadPosition {
	public:
		ReadIdType readId;
		PositionType position;

		ReadPosition(ReadIdType _read = 0, PositionType _pos = 0) :
			readId(_read), position(_pos) {
		}
	};

	class ReadPositionWeight: public ReadPosition {
	public:
		WeightType weight;

		ReadPositionWeight(ReadIdType _read = 0, PositionType _pos = 0,
				WeightType _weight = 0.0) :
			ReadPosition(_read, _pos), weight(_weight) {
		}
	};
	typedef std::vector<ReadPositionWeight> ReadPositionWeightVector;


	static void resetGlobalCounters() {
		discarded = 0;
		maxCount = 0;
		maxWeightedCount = 0.0;
	}
	static void resetForGlobals(CountType count) {
	}
	static void setGlobals(CountType count, WeightType weightedCount) {
		if (count > maxCount) {
			maxCount = count;
		}
		if (weightedCount > maxWeightedCount) {
			maxWeightedCount = weightedCount;
		}
	}
	static inline bool isDiscard(WeightType weight) {
		if (weight < minimumWeight) {

			#pragma omp atomic
			discarded++;

			return true;
		} else
			return false;
	}
	static inline bool useWeighted() {
		return useWeightedByDefault;
	}
	static void setUseWeighted(bool _useWeighted = true) {
		useWeightedByDefault = _useWeighted;
	}
	static void setMinimumWeight(WeightType _minimumWeight) {
		minimumWeight = _minimumWeight;
	}
	static inline WeightType getMinimumWeight() {
		return minimumWeight;
	}
	static void setMinimumDepth(CountType _minimumDepth) {
		minimumDepth = _minimumDepth;
	}
	static inline CountType getMinimumDepth() {
		return minimumDepth;
	}
	static void discard() {
		discarded++;
	}
	static unsigned long getDiscarded() {
		return discarded;
	}

private:
	static WeightType minimumWeight;
	static CountType minimumDepth;
	static unsigned long discarded;

	static CountType maxCount;
	static WeightType maxWeightedCount;
	static bool useWeightedByDefault;


protected:
	CountType count;
	WeightType weightedCount;

public:
	TrackingData() :
		count(0),  weightedCount(0.0) {
	}

	void reset() {
		resetForGlobals(getCount());
		count = 0;
		weightedCount = 0.0;
	}
	inline bool operator==(const TrackingData &other) const {
		return getCount() == other.getCount();// && weightedCount == other.weightedCount;
	}
	inline bool operator<(const TrackingData &other) const {
		return getCount() < other.getCount();// || (count == other.count && weightedCount < other.weightedCount);
	}

	bool track(double weight, bool forward, ReadIdType readIdx = 0,
			PositionType readPos = 0) {
		assert(weight <= 1.0);

		if (isDiscard(weight))
			return false;

		if (count < MAX_COUNT) {
#ifdef _USE_THREADSAFE_KMER
#pragma omp critical (TrackingData)
#endif
			{
				count++;
				weightedCount += weight;
			}

			setGlobals(getCount(), getWeightedCount());

			return true;
		} else
			return false;

	}

	inline unsigned long getCount() const {
		return count;
	}
	inline unsigned long getDirectionBias() const {
		return count / 2;
	}
	inline double getWeightedCount() const {
		return weightedCount;
	}
	inline double getAverageWeight() const {
		return weightedCount / count;
	}
	inline double getNormalizedDirectionBias() const {
		return 0.5;
	}

	ReadPositionWeightVector getEachInstance() const {
		//returns count entries from read -1, position 0
		ReadPositionWeight dummy1((ReadIdType) -1, 0, getWeightedCount()
				/ getCount());
		ReadPositionWeightVector dummy(getCount(), dummy1);
		return dummy;
	}

	std::string toString() const {
		std::stringstream ss;
		ss << count << ":" << std::fixed << std::setprecision(2)
				<< getNormalizedDirectionBias();
		ss << ':' << std::fixed << std::setprecision(2)
				<< ((double) weightedCount / (double) count);
		return ss.str();
	}

	// cast operator
	operator unsigned long() const {
		return getCount();
	}

	template<typename U>
	TrackingData &add(const U &other) {
			count += other.getCount();
			weightedCount += other.getWeightedCount();
			return *this;
	}
	template<typename U>
	TrackingData &operator=(const U &other) {
		count = other.getCount();
		weightedCount = other.getWeightedCount();
		return *this;
	};
};

class TrackingDataWithDirection: public TrackingData {
public:

protected:
	CountType directionBias;

public:
	TrackingDataWithDirection(): TrackingData(), directionBias(0) {}
	void reset() {
		TrackingData::reset();
		directionBias = 0;
	}

	bool track(double weight, bool forward, ReadIdType readIdx, PositionType readPos) {
		assert(weight <= 1.0);

		bool ret = TrackingData::track(weight,forward);

		if (ret) {
			if (forward)
				directionBias++;
		}
		return ret;
	}
	inline unsigned long getDirectionBias() const {
		return directionBias;
	}
	inline double getNormalizedDirectionBias() const {
		return (double) directionBias / (double) count;
	}

	template<typename U>
	TrackingDataWithDirection &add(const U &other) {
		TrackingData::add(other);
		directionBias += other.getDirectionBias();
		return *this;
	};
	template<typename U>
	TrackingDataWithDirection &operator=(const U &other) {
		count = other.getCount();
		weightedCount = other.getWeightedCount();
		directionBias = other.getDirectionBias();
		return *this;
	}

};

class TrackingDataWithLastRead: public TrackingDataWithDirection {
public:

protected:
	ReadPosition readPosition;

public:
	TrackingDataWithLastRead() :
		TrackingDataWithDirection(), readPosition() {
	}

	void reset() {
		TrackingDataWithDirection::reset();
		readPosition = ReadPosition();
	}
	bool track(double weight, bool forward, ReadIdType readIdx,
			PositionType readPos) {
		assert(weight <= 1.0);

		bool ret = TrackingDataWithDirection::track(weight, forward, readIdx, readPos);
		if (ret)
			readPosition = ReadPosition(readIdx, readPos);
		return ret;
	}
	ReadPositionWeightVector getEachInstance() const {
		//returns count entries from read -1, position 0
		ReadPositionWeight dummy1(readPosition.readId, readPosition.position,
				getWeightedCount() / getCount());
		ReadPositionWeightVector dummy(getCount(), dummy1);
		return dummy;
	}

	template<typename U>
	TrackingDataWithLastRead &add(const U &other) {
		TrackingDataWithDirection::add(other);
		ReadPositionWeightVector v = other.getEachInstance();
		if (v.size() > 0) {
			ReadPositionWeight &rpw = v[v.size()-1];
			readPosition.position = rpw.position;
			readPosition.readId = rpw.readId;
		}
		return *this;
	}
	template<typename U>
	TrackingDataWithDirection &operator=(const U &other) {
		count = other.getCount();
		weightedCount = other.getWeightedCount();
		directionBias = other.getDirectionBias();
		ReadPositionWeight v = other.getEachInstance()[0];
		readPosition.position = v.position;
		readPosition.readId = v.readId;
		return *this;
	}


};

class TrackingDataSingleton : public ExtensionTrackingInterface {
public:
	typedef TrackingData::PositionType PositionType;
	typedef TrackingData::ReadPositionWeight ReadPositionWeight;
	typedef TrackingData::ReadPositionWeightVector ReadPositionWeightVector;
	typedef TrackingData::CountType CountType;
	typedef TrackingData::WeightType WeightType;
	typedef TrackingData::ReadIdType ReadIdType;

protected:
	unsigned char _weight;

public:
	TrackingDataSingleton() : _weight(0) {}
	void reset() {
		TrackingData::resetForGlobals(getCount());
		_weight = 0;
	}
	inline bool operator==(const TrackingDataSingleton &other) const {
		return getCount() == other.getCount();
	}
	inline bool operator<(const TrackingDataSingleton &other) const {
		return getCount() < other.getCount();
	}
	bool track(double weight, bool forward, ReadIdType readIdx = 0, PositionType readPos = 0) {
		assert(weight <= 1.0);

		if (TrackingData::isDiscard(weight))
			return false;
		_weight = (unsigned char) ((weight * 254.0)) + 1;
		TrackingData::setGlobals(1, weight);
		return true;
	}

	inline unsigned long getCount() const {
		return _weight == 0 ? 0 : 1;
	}
	inline unsigned long getDirectionBias() const {
		return 0;
	}
	inline double getWeightedCount() const {
		return _weight == 0 ? 0.0 : (_weight - 1 ) / 254.0;
	}
	inline double getNormalizedDirectionBias() const {
		return 0.5;
	}

	ReadPositionWeightVector getEachInstance() const {
		//returns count entries from read -1, position 0
		ReadPositionWeight dummy1((ReadIdType) -1, 0, getWeightedCount()
				/ getCount());
		ReadPositionWeightVector dummy(getCount(), dummy1);
		return dummy;
	}

	// cast operator
	operator unsigned long() const {
		return getCount();
	}

	std::string toString() const {
		std::stringstream ss;
		ss << getCount() << ":" << std::fixed << std::setprecision(2)
				<< getNormalizedDirectionBias();
		ss << ':' << std::fixed << std::setprecision(2)
				<< ((double) getWeightedCount() / (double) getCount());
		return ss.str();
	}

};

class TrackingDataSingletonWithReadPosition : public ExtensionTrackingInterface {
public:
	typedef TrackingData::ReadPositionWeight ReadPositionWeight;
	typedef TrackingData::ReadPositionWeightVector ReadPositionWeightVector;
	typedef TrackingData::CountType CountType;
	typedef TrackingData::WeightType WeightType;
	typedef TrackingData::ReadIdType ReadIdType;
	typedef TrackingData::PositionType PositionType;

protected:
	ReadPositionWeight instance; // save space and store direction within sign of weight.

public:
	TrackingDataSingletonWithReadPosition() :
		instance(0, 0, 0.0) {
	}

	void reset() {
		// destructor should *not* call this
		TrackingData::resetForGlobals(getCount());
		instance = ReadPositionWeight(0, 0, 0.0);
	}
	inline bool operator==(const TrackingDataSingletonWithReadPosition &other) const {
		return getCount() == other.getCount();// && weightedCount == other.weightedCount;
	}
	inline bool operator<(const TrackingDataSingletonWithReadPosition &other) const {
		return getCount() < other.getCount();// || (count == other.count && weightedCount < other.weightedCount);
	}

	bool track(double weight, bool forward, ReadIdType readIdx = 0,
			PositionType readPos = 0) {
		assert(weight <= 1.0);

		if (TrackingData::isDiscard(weight))
			return false;
#ifdef _USE_THREADSAFE_KMER
#pragma omp critical (TrackingDataSingletonWithReadPosition)
#endif
		instance = ReadPositionWeight(readIdx, readPos, forward ? weight : -1.0
				* weight);

		TrackingData::setGlobals(getCount(), getWeightedCount());

		return true;
	}

	inline unsigned long getCount() const {
		return instance.weight == 0.0 ? 0 : 1;
	}
	inline unsigned long getDirectionBias() const {
		return instance.weight > 0 ? 1 : 0;
	}
	inline double getWeightedCount() const {
		return instance.weight > 0 ? instance.weight : -1.0 * instance.weight;
	}
	inline double getNormalizedDirectionBias() const {
		return (double) getDirectionBias() / (double) getCount();
	}

	inline ReadPositionWeightVector getEachInstance() const {
		return ReadPositionWeightVector(1, instance);
	}

	inline ReadIdType getReadId() const {
		return instance.readId;
	}
	inline PositionType getPosition() const {
		return instance.position;
	}

	// cast operator
	operator unsigned long() const {
		return getCount();
	}

	std::string toString() const {
		std::stringstream ss;
		ss << getCount() << ":" << std::fixed << std::setprecision(2)
				<< getNormalizedDirectionBias();
		ss << ':' << std::fixed << std::setprecision(2)
				<< ((double) getWeightedCount() / (double) getCount());
		return ss.str();
	}


};

class TrackingDataWithAllReads : public ExtensionTrackingInterface {
public:
	typedef TrackingData::ReadPositionWeight ReadPositionWeight;
	typedef TrackingData::ReadPositionWeightVector ReadPositionWeightVector;
	typedef TrackingData::CountType CountType;
	typedef TrackingData::WeightType WeightType;
	typedef TrackingData::ReadIdType ReadIdType;
	typedef TrackingData::PositionType PositionType;

protected:
	ReadPositionWeightVector instances;
	unsigned int directionBias;
	WeightType weightedCount; // for performance reasons

public:
	TrackingDataWithAllReads() :
		instances(), directionBias(0), weightedCount(0.0) {
	}

	void reset() {
		TrackingData::resetForGlobals(getCount());
		instances.clear();
		directionBias = 0;
		weightedCount = 0.0;
	}
	inline bool operator==(const TrackingDataWithAllReads &other) const {
		return getCount() == other.getCount();// && weightedCount == other.weightedCount;
	}
	inline bool operator<(const TrackingDataWithAllReads &other) const {
		return getCount() < other.getCount();// || (count == other.count && weightedCount < other.weightedCount);
	}

	bool track(double weight, bool forward, ReadIdType readIdx = 0,
			PositionType readPos = 0) {
		assert(weight <= 1.0);

		if (TrackingData::isDiscard(weight))
			return false;
#ifdef _USE_THREADSAFE_KMER
#pragma omp critical (TrackingDataWithAllReads)
#endif
		{
			if (forward)
				directionBias++;

			ReadPositionWeight rpw(readIdx, readPos, weight);

			instances.push_back(rpw);
			weightedCount += weight;
		}
		TrackingData::setGlobals(getCount(), getWeightedCount());

		return true;
	}

	inline unsigned long getCount() const {
		return instances.size();
	}
	inline unsigned long getDirectionBias() const {
		return directionBias;
	}
	inline double getWeightedCount() const {
		return weightedCount;
	}
	double _getWeightedCount() const {
		double weightedCount = 0.0;
		for (ReadPositionWeightVector::const_iterator it = instances.begin(); it
				!= instances.end(); it++)
			weightedCount += it->weight;
		return weightedCount;
	}
	inline double getNormalizedDirectionBias() const {
		return (double) getDirectionBias() / (double) getCount();
	}

	inline ReadPositionWeightVector getEachInstance() const {
		return instances;
	}

	// cast operator
	operator unsigned long() const {
		return getCount();
	}

	std::string toString() const {
		std::stringstream ss;
		ss << getCount() << ":" << std::fixed << std::setprecision(2)
				<< getNormalizedDirectionBias();
		ss << ':' << std::fixed << std::setprecision(2)
				<< ((double) getWeightedCount() / (double) getCount());
		return ss.str();
	}
	TrackingDataWithAllReads &operator=(const TrackingDataSingleton &other) {
		this->reset();
		if (other.getCount() > 0) {
			directionBias = other.getDirectionBias();
			instances = other.getEachInstance();
			weightedCount = _getWeightedCount();
		}
		return *this;
	}
	TrackingDataWithAllReads &add(const TrackingDataWithAllReads &other) {
		instances.insert(instances.end(), other.instances.begin(), other.instances.end());
		directionBias += other.getDirectionBias();
		weightedCount += other.getWeightedCount();
		return *this;
	}

	TrackingDataWithAllReads &add(const TrackingDataSingleton &other) {
		if (other.getCount() > 0) {
	       instances.push_back(ReadPositionWeight((ReadIdType) -1, 0, other.getWeightedCount()));
	       weightedCount += other.getWeightedCount();
		}
		return *this;
	}

};

std::ostream &operator<<(std::ostream &stream, TrackingData &ob);
std::ostream &operator<<(std::ostream &stream, TrackingDataSingleton &ob);
std::ostream &operator<<(std::ostream &stream, TrackingDataWithAllReads &ob);

class WeightedExtensionMessagePacket : public ExtensionMessagePacket {
public:
	typedef TrackingData::WeightType WeightType;

	WeightedExtensionMessagePacket() : ExtensionMessagePacket(), _weight(0) {}
	WeightedExtensionMessagePacket(const WeightedExtensionMessagePacket &copy) {
		*this = copy;
	}
	~WeightedExtensionMessagePacket() {}
	WeightedExtensionMessagePacket &operator=(const WeightedExtensionMessagePacket &copy) {
		*((ExtensionMessagePacket*)this) = (ExtensionMessagePacket) copy;
		_weight = copy._weight;
		return *this;
	}
	WeightType getWeight() const {
		return _weight;
	}
	void setWeight(WeightType weight) {
		_weight = weight;
	}

protected:
	WeightType _weight;
};

class ExtensionTrackingData : public TrackingDataWithDirection {
public:
protected:
	ExtensionTracking _extensionTracking;

public:
	ExtensionTrackingData(): TrackingDataWithDirection(), _extensionTracking() {}
	ExtensionTrackingData(const ExtensionTrackingData &copy) {
		*this = copy;
	}
	~ExtensionTrackingData() {}

	void reset() {
		TrackingDataWithDirection::reset();
	}

	bool track(double weight, bool forward, ReadIdType readIdx, PositionType readPos) {
		assert(weight <= 1.0);

		bool ret = TrackingDataWithDirection::track(weight,forward,readIdx,readPos);

		return ret;
	}
	void trackExtensions(Extension left, Extension right) {
		_extensionTracking.trackExtensions(left, right);
	}
	ExtensionTracking getExtensionTracking() const {
		return _extensionTracking;
	}

	template<typename U>
	ExtensionTrackingData &add(const U &other) {
		TrackingDataWithDirection::add(other);
		_extensionTracking = _extensionTracking + other.getExtensionTracking();

		return *this;
	};

	template<typename U>
	ExtensionTrackingData &operator=(const U &other) {
		*((TrackingDataWithDirection*) this) = other;
		_extensionTracking = other.getExtensionTracking();

		return *this;
	};

};
std::ostream &operator<<(std::ostream &stream, ExtensionTrackingData &ob);


class ExtensionTrackingDataSingleton : public TrackingDataSingleton {
public:
	ExtensionTrackingDataSingleton() : TrackingDataSingleton(), extensionMsgPacket() {}
	ExtensionTrackingDataSingleton(ExtensionTrackingDataSingleton &copy) {
		*this = copy;
	}
	~ExtensionTrackingDataSingleton() {}

	void reset() {
		TrackingDataSingleton::reset();
	}

	bool track(double weight, bool forward, ReadIdType readIdx, PositionType readPos) {
		assert(weight <= 1.0);

		bool ret = TrackingDataSingleton::track(weight,forward,readIdx,readPos);

		return ret;
	}
	void trackExtensions(Extension left, Extension right) {
		extensionMsgPacket.setExtensions(left, right);
	}
	ExtensionTracking getExtentionTracking() const {
		ExtensionTracking extTrack;
		extTrack.trackExtensions(extensionMsgPacket.getLeft(), extensionMsgPacket.getRight());
		return extTrack;
	}

	template<typename U>
	ExtensionTrackingDataSingleton &add(const U &other) {
		TrackingDataWithDirection::add(other);
		ExtensionTracking extensionTracking = other.getExtensionTracking();
		extensionMsgPacket.setExtensions(extensionTracking);
		return *this;
	};

	template<typename U>
	ExtensionTrackingDataSingleton &operator=(const U &other) {
		*((TrackingDataWithDirection*) this) = other;
		ExtensionTracking extensionTracking = other.getExtensionTracking();
		extensionMsgPacket.setExtensions(extensionTracking);
		return *this;
	};

protected:
	ExtensionMessagePacket extensionMsgPacket;
};

std::ostream &operator<<(std::ostream &stream, ExtensionTrackingDataSingleton &ob);


template<typename T>
class TrackingDataMinimal : public ExtensionTrackingInterface {
public:
	typedef T DataType;
    typedef	typename TrackingData::CountType CountType;
	typedef typename TrackingData::WeightType WeightType;
	typedef typename TrackingData::ReadIdType ReadIdType;
	typedef typename TrackingData::PositionType PositionType;
	typedef typename TrackingData::ReadPositionWeightVector ReadPositionWeightVector;

private:
	DataType count;

public:
	TrackingDataMinimal() : count(0) {}
	void reset() {
		TrackingData::resetForGlobals(getCount());
	    count = 0;
	}

	bool track(double weight, bool forward, ReadIdType readIdx, PositionType readPos)
	{
		assert(weight <= 1.0);
		if (TrackingData::isDiscard(weight)) {
			return false;
		}
		DataType test = (DataType) weight;
		// handle both integers and floats "properly"
		if (test < weight) {
			count += 1;
		} else {
			count += (DataType) weight;
		}
		TrackingData::setGlobals(getCount(), getWeightedCount());
		return true;
	}
	inline unsigned long getCount() const { unsigned long tmp = count; return (tmp<count) ? tmp+1 : tmp; }
	inline unsigned long getDirectionBias() const {return count / 2;}
	inline double getWeightedCount() const {return count;}
	inline double getNormalizedDirectionBias() const {return 0.5;}
	ReadPositionWeightVector getEachInstance() const
	{
		return ReadPositionWeightVector(0);
	}

	// cast operator
	operator unsigned long() const {
		return getCount();
	}

	std::string toString() const {
		std::stringstream ss;
		ss << getCount();
		return ss.str();
	}
	template<typename U>
	TrackingDataMinimal &operator=(const U &other) {
		double weight = other.getWeightedCount();
		DataType weightD = weight;
		if (weightD < weight)
			count = other.getCount();
		else
		    count = weightD;
		return *this;
	};

	template<typename U>
	TrackingDataMinimal &add(const U &other) {
		double weight = other.getWeightedCount();
		DataType weightD = weight;
		if (weightD < weight)
			count += other.getCount();
		else
		    count += weightD;
		return *this;
	};

};

typedef TrackingDataMinimal<unsigned char> TrackingDataMinimal1;
typedef TrackingDataMinimal<unsigned short> TrackingDataMinimal2;
typedef TrackingDataMinimal<unsigned int> TrackingDataMinimal4;
typedef TrackingDataMinimal<float> TrackingDataMinimal4f;
typedef TrackingDataMinimal<double> TrackingDataMinimal8;

template<typename T>
std::ostream &operator<<(std::ostream &stream, TrackingDataMinimal<T> &ob) {
	stream << ob.toString();
	return stream;
};


#endif

// $Log: KmerTrackingData.h,v $
// Revision 1.6  2010-08-18 17:50:40  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.5.4.1  2010-07-13 19:45:51  regan
// bugfix in casting TrackingData*
//
// Revision 1.5  2010-06-22 23:06:31  regan
// merged changes in CorruptionBugfix-20100622 branch
//
// Revision 1.4.6.1  2010-06-22 23:01:55  regan
// named all critical sections
//
// Revision 1.4  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.3.2.1  2010-05-07 22:59:33  regan
// refactored base type declarations
//
// Revision 1.3  2010-05-06 22:55:05  regan
// merged changes from CodeCleanup-20100506
//
// Revision 1.2.10.1  2010-05-06 18:45:35  regan
// broke it...
//
// Revision 1.2  2010-05-01 21:57:54  regan
// merged head with serial threaded build partitioning
//
// Revision 1.1.2.2  2010-04-27 05:38:32  regan
// added count cast operator
//
// Revision 1.1.2.1  2010-04-26 04:59:46  regan
// bugfix and templated some tracking methods
//
// Revision 1.1  2010-04-21 23:39:37  regan
// got kmermap mmap store and restore working
//
//
