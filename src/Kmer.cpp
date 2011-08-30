//
// Kmernator/src/Kmer.cpp
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

#include "Kmer.h"
#include "MemoryUtils.h"

KmerSizer KmerSizer::singleton;

TrackingData::WeightType TrackingData::minimumWeight = 0.01;
TrackingData::CountType TrackingData::minimumDepth = 10;
unsigned long TrackingData::discarded = 0;

TrackingData::CountType TrackingData::maxCount = 0;
TrackingData::WeightType TrackingData::maxWeightedCount = 0;
bool TrackingData::useWeightedByDefault = true;

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

