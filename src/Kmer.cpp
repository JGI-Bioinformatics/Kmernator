#include "Kmer.h"
#include "MemoryUtils.h"

KmerSizer KmerSizer::singleton;

TrackingData::WeightType TrackingData::minimumWeight = 0.01;
TrackingData::CountType TrackingData::minimumDepth = 10;

TrackingData::CountType TrackingData::maxCount = 0;
TrackingData::WeightType TrackingData::maxWeightedCount = 0;
unsigned long TrackingData::discarded = 0;

std::ostream &operator<<(std::ostream &stream, TrackingData &ob) {
	stream << ob.toString();
	return stream;
};

std::ostream &operator<<(std::ostream &stream, TrackingDataSingleton &ob) {
	stream << ob.toString();
	return stream;
};

std::ostream &operator<<(std::ostream &stream, TrackingDataWithAllReads &ob) {
	stream << ob.toString();
	return stream;
};


ClassicMemory ClassicMemory::singleton;
BoostPoolManager BoostPoolManager::singleton;
