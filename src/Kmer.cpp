#include "Kmer.h"
#include "MemoryUtils.h"

KmerSizer KmerSizer::singleton;

double SolidTrackingData::minimumWeight = 0.1;
unsigned short SolidTrackingData::maxCount = 0;
unsigned short SolidTrackingData::minimumDepth = 10;

std::ostream &operator<<(std::ostream &stream, SolidTrackingData &ob)
{
  stream << ob.toString();
  return stream;
};


ClassicMemory ClassicMemory::singleton;
BoostPoolManager BoostPoolManager::singleton;
