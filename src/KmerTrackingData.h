// $Header: /repository/PI_annex/robsandbox/KoMer/src/KmerTrackingData.h,v 1.1 2010-04-21 23:39:37 regan Exp $
//

#ifndef _KMER_TRACKING_DATA_H
#define _KMER_TRACKING_DATA_H

#include "config.h"

class TrackingDataSingleton;
class TrackingDataWithAllReads;

class TrackingData {
public:
	typedef unsigned short CountType;
	typedef float WeightType;
	typedef unsigned int ReadIdType;
	typedef unsigned short PositionType;

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

	static WeightType minimumWeight;
	static CountType minimumDepth;
	static const CountType MAX_COUNT = (CountType) -1;
	static unsigned long discarded;

	static CountType maxCount;
	static WeightType maxWeightedCount;

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
#ifdef _USE_OPENMP
#pragma omp atomic
#endif
			discarded++;
			return true;
		} else
			return false;
	}

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
		if (isDiscard(weight))
			return false;

		if (count < MAX_COUNT) {
#ifdef _USE_THREADSAFE_KMER
#pragma omp critical
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

	TrackingData &add(const TrackingData &other) {
		count += other.count;
		weightedCount += other.weightedCount;
		return *this;
	}
	TrackingData &operator=(const TrackingDataSingleton &other);
	TrackingData &operator=(const TrackingDataWithAllReads &other);
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

	TrackingDataWithDirection &add(const TrackingDataWithDirection &other) {
		TrackingData::add((TrackingData) other);
		directionBias += other.directionBias;
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

	TrackingDataWithLastRead &add(const TrackingDataWithLastRead &other) {
		TrackingDataWithDirection::add((TrackingDataWithDirection) other);
		readPosition = other.readPosition;
		return *this;
	}


};

class TrackingDataSingleton {
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

	std::string toString() const {
		std::stringstream ss;
		ss << getCount() << ":" << std::fixed << std::setprecision(2)
				<< getNormalizedDirectionBias();
		ss << ':' << std::fixed << std::setprecision(2)
				<< ((double) getWeightedCount() / (double) getCount());
		return ss.str();
	}
};

class TrackingDataSingletonWithReadPosition {
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
		if (TrackingData::isDiscard(weight))
			return false;
#ifdef _USE_THREADSAFE_KMER
#pragma omp critical
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

	std::string toString() const {
		std::stringstream ss;
		ss << getCount() << ":" << std::fixed << std::setprecision(2)
				<< getNormalizedDirectionBias();
		ss << ':' << std::fixed << std::setprecision(2)
				<< ((double) getWeightedCount() / (double) getCount());
		return ss.str();
	}

};

class TrackingDataWithAllReads {
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
		if (TrackingData::isDiscard(weight))
			return false;
#ifdef _USE_THREADSAFE_KMER
#pragma omp critical
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
		directionBias += other.directionBias;
		weightedCount += other.weightedCount;
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

template<typename T>
class TrackingDataMinimal {
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
		if (weight > 1.0) {
			std::stringstream ss;
			ss << "How can the weight be greater than one? " << weight << "," << forward << "," << readIdx << "," << readPos;
			throw std::invalid_argument(ss.str());
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

	std::string toString() const {
		std::stringstream ss;
		ss << getCount();
		return ss.str();
	}
	TrackingDataMinimal &operator=(const TrackingData &other) {
		count = (DataType) other.getWeightedCount();
		return *this;
	}TrackingDataMinimal &operator=(const TrackingDataWithAllReads &other) {
		count = (DataType) other.getWeightedCount();
		return *this;
	}
	TrackingDataMinimal &operator=(const TrackingDataSingleton &other) {
		count = (DataType) other.getWeightedCount();
		return *this;
	}
	TrackingDataMinimal &add(const TrackingDataMinimal &other) {
		count += other.count;
		return *this;
	}

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
// Revision 1.1  2010-04-21 23:39:37  regan
// got kmermap mmap store and restore working
//
//
