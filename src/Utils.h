// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.30 2010-03-03 17:10:05 regan Exp $
//

#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <tr1/memory>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "config.h"

class OfstreamMap {
public:
	typedef std::tr1::shared_ptr< std::ofstream > OStreamPtr;
	typedef boost::unordered_map< std::string, OStreamPtr > Map;
	typedef Map::iterator Iterator;
private:
    Map _map;
    std::string _outputFilePathPrefix;
    std::string _suffix;

public:
	OfstreamMap(std::string outputFilePathPrefix, std::string suffix = "")
	 : _outputFilePathPrefix(outputFilePathPrefix), _suffix(suffix) {}
	~OfstreamMap() {
        clear();
	}
	void clear() {
		for(Iterator it = _map.begin() ; it != _map.end(); it++) {
			if (Options::getDebug()) {
				std::cerr << "Closing " << it->first << std::endl;
			}
			it->second->close();
		}
		_map.clear();
	}

	std::ofstream &getOfstream(std::string key) {
		std::string filename = _outputFilePathPrefix + key + _suffix;
		Iterator it = _map.find(filename);
		if (it == _map.end()) {
			if (Options::getDebug()) {
			  std::cerr << "Opening " << filename << std::endl;
		    }

			OStreamPtr osp(new std::ofstream(filename.c_str()));
			it = _map.insert( it, Map::value_type(filename, osp) );
		}
		return *(it->second);
	}
};

template<typename S>
class PartitioningData {
public:
	typedef S DataType;
	typedef std::vector< DataType > Partitions;

private:
	Partitions _partitions;

public:
	PartitioningData() : _partitions() {}

    inline bool hasPartitions() const {
    	return ! _partitions.empty();
    }
    inline int getPartitionNum(DataType score) const {
		// TODO binary search?
		for(int i = 0; i < _partitions.size(); i++)
			if (score < _partitions[i])
				return i+1;
		return _partitions.size() + 1;
	}
    inline Partitions getPartitions() const {
    	return _partitions;
    }

    inline int addPartition(DataType partition) {
    	_partitions.push_back(partition);
    	std::sort(_partitions.begin(), _partitions.end());
    }


};

template<typename Raw, typename Store>
class BucketedData {
public:
	typedef Store StoreType;
	typedef Raw RawType;

private:
	RawType _minValue, _maxValue;
	StoreType steps;

public:
	// TODO inclusive permutations: (), [), (], []
	//      presently it truncates fraction, so [)
	// TODO implement log scale
	BucketedData(RawType minValue, RawType maxValue) :
		_minValue(minValue), _maxValue(maxValue) {
		// assumes unsigned Store...
		steps = ((StoreType) -1);
	}
	~BucketedData() {
	}

	StoreType getStore(RawType value) {
		if (value < _minValue || value > _maxValue)
			throw std::invalid_argument("getStore() out of range");
		return ((RawType) steps) * (value - _minValue)
				/ (_maxValue - _minValue);
	}

	RawType getValue(StoreType store) {
		return (_maxValue - _minValue) * ((RawType) store) / ((RawType) steps);
	}
};

typedef unsigned char OneByte;
typedef unsigned short TwoByte;
typedef unsigned int FourByte;
typedef unsigned long EightByte;

typedef BucketedData<double, OneByte> DoubleToOneByte;
typedef BucketedData<double, TwoByte> DoubleToTwoByte;
typedef BucketedData<double, FourByte> DoubleToFourByte;
typedef BucketedData<double, EightByte> DoubleToEightByte;

#endif

//
// $Log: Utils.h,v $
// Revision 1.30  2010-03-03 17:10:05  regan
// added two helper classes
// partitioning data and ofstream mapper
//
// Revision 1.29  2010-02-26 13:01:17  regan
// reformatted
//
// Revision 1.28  2010-01-13 23:48:51  regan
// refactored
//
// Revision 1.27  2010-01-13 07:20:08  regan
// refactored filter
// checkpoint on read picker
//
// Revision 1.26  2010-01-13 00:25:43  regan
// use less memory for reference sequences and those without quality
//
// Revision 1.25  2010-01-11 19:14:10  regan
// minor performance enhancements
//
// Revision 1.24  2010-01-11 05:40:09  regan
// believe that FilterKnownOddities is working
//
// Revision 1.23  2010-01-08 06:25:18  regan
// refactored some code
// still working on FilterKnownOddities
//
// Revision 1.22  2010-01-06 15:20:24  regan
// code to screen out primers
//
// Revision 1.21  2010-01-05 06:43:47  regan
// bugfix
//
// Revision 1.20  2009-12-24 00:56:11  regan
// started class to output picked reads
//
// Revision 1.19  2009-12-21 06:34:26  regan
// used openmp and clever partitioning to speed up building spectrum
//
// Revision 1.18  2009-11-28 01:00:07  regan
// fixed bugs and warnings
//
// Revision 1.17  2009-11-26 09:03:29  regan
// refactored and stuff
//
// Revision 1.16  2009-11-24 13:35:29  cfurman
// removed KmerPtr class.
//
// Revision 1.15  2009-11-22 08:16:41  regan
// some fixes some bugs... optimized vs debug vs deb4/5 give different results
//
// Revision 1.14  2009-11-21 18:46:53  regan
// added bugs
//
// Revision 1.13  2009-11-21 15:58:29  regan
// changed some types
// bugfix in reading and using qual files
//
// Revision 1.12  2009-11-12 17:01:51  regan
// checkpoint
//
// Revision 1.11  2009-11-12 01:29:22  cfurman
// Solid picking logic bug fixed, params tweaked.
//
// Revision 1.10  2009-11-11 17:23:23  regan
// fixed bugs in heap generation
// solid picking logic needs work
//
// Revision 1.9  2009-11-11 07:57:23  regan
// built framework for autoPromote (not working) - make_heap is broken
//
// Revision 1.8  2009-11-10 07:05:37  regan
// changes for debugging
//
// Revision 1.7  2009-11-09 19:37:17  regan
// enhanced some debugging / analysis output
//
// Revision 1.6  2009-11-06 16:59:11  regan
// added base substitution/permutations table and build function
//
// Revision 1.5  2009-11-06 04:10:21  regan
// refactor of cmd line option handling
// added methods to evaluate spectrums
//
// Revision 1.4  2009-11-04 18:26:17  regan
// refactored
// added statistics calculations and histograms
//
// Revision 1.3  2009-11-03 17:15:40  regan
// minor refactor
//
// Revision 1.2  2009-11-02 21:19:25  regan
// fixed types and boundary tests
//
// Revision 1.1  2009-11-02 18:47:34  regan
// added some code without a permanent home
//
//
