#include "Kmer.h"
#include "MemoryUtils.h"

KmerSizer KmerSizer::singleton;

TrackingData::WeightType TrackingData::minimumWeight = 0.01;
TrackingData::CountType  TrackingData::minimumDepth = 10;

TrackingData::CountType  TrackingData::maxCount = 0;
TrackingData::WeightType TrackingData::maxWeightedCount = 0;
unsigned long            TrackingData::discarded = 0;
unsigned long            TrackingData::singletonCount = 0;

std::ostream &operator<<(std::ostream &stream, TrackingData &ob)
{
  stream << ob.toString();
  return stream;
};

std::ostream &operator<<(std::ostream &stream, TrackingDataWithAllReads &ob)
{
  stream << ob.toString();
  return stream;
};

TrackingData &TrackingData::operator=(const TrackingDataWithAllReads &copy) {
	count = copy.getCount();
	weightedCount = copy.getWeightedCount();
	directionBias = copy.getDirectionBias();
	return *this;
}

ClassicMemory ClassicMemory::singleton;
BoostPoolManager BoostPoolManager::singleton;
