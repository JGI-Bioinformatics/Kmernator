// $Header: /repository/PI_annex/robsandbox/KoMer/src/Kmer.cpp,v 1.18 2010-05-18 20:50:24 regan Exp $

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


// $Log: Kmer.cpp,v $
// Revision 1.18  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.17.2.1  2010-05-10 17:57:34  regan
// added header
//
// Revision 1.17  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
//

