//
// Kmernator/src/Utils.h
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

#ifndef _UTILS_H
#define _UTILS_H

#include <ostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <set>
#include <cstring>
#include <cmath>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#define foreach BOOST_FOREACH

#include "config.h"
#include "Options.h"
#include "Log.h"


class FormatOutput
{
private:
	int _type;

public:
	FormatOutput(int type) : _type(type) {
		if (type < 0 || type > 3)
			throw std::invalid_argument("Format can only be 0-3");
	}
	static const int FASTQ = 0;
	static const int FASTA = 1;
	static const int FASTQ_UNMASKED = 2;
	static const int FASTA_UNMASKED = 3;

	static inline FormatOutput getDefault() {
		return FormatOutput(Options::getFormatOutput());
	}
	static std::string getDefaultSuffix() {
		return getDefault().getSuffix();
	}

	inline int getType() const {
		return _type;
	}
	bool inline operator==(const FormatOutput &other) const {
		return _type == other._type;
	}
	inline std::string getSuffix() const {
		return getSuffix(*this);
	}
	static inline std::string getSuffix(const FormatOutput &format) {
		return getSuffix(format._type);
	}
	static inline std::string getSuffix(int type) {
		std::string ret;
		switch(type) {
		case 0:
		case 2: ret = std::string(".fastq"); break;
		case 1:
		case 3: ret = std::string(".fasta"); break;
		default : ret = std::string(".txt");
		}
		return ret;
	}
};


class OfstreamMap {
public:
	typedef boost::shared_ptr< std::ofstream > OStreamPtr;
	typedef boost::unordered_map< std::string, OStreamPtr > Map;
	typedef Map::iterator Iterator;
	typedef boost::shared_ptr< Map > MapPtr;
protected:
    MapPtr _map;
    std::string _outputFilePathPrefix;
    std::string _suffix;
    bool _append;

public:
	static bool &getDefaultAppend() {
		static bool _defaultAppend = false;
		return _defaultAppend;
	}

	OfstreamMap(std::string outputFilePathPrefix = Options::getOutputFile(), std::string suffix = FormatOutput::getDefaultSuffix())
	 : _map(new Map()), _outputFilePathPrefix(outputFilePathPrefix), _suffix(suffix) {
		_append = getDefaultAppend();
	}
	virtual ~OfstreamMap() {
        clear();
	}
	std::set<std::string> getFiles() {
		std::set<std::string> files;
		std::string rank = getRank();
		for(Iterator it = _map->begin() ; it != _map->end(); it++) {
			std::string file = it->first;
			if (!rank.empty()) {
				file = file.substr(0, file.find(rank));
			}
			files.insert(file);
		}
		return files;
	}
	bool &getAppend() {
		return _append;
	}
	const bool &getAppend() const {
		return _append;
	}

	void clear() {
		close();
		_map->clear();
	}
	virtual void close() {
		for(Iterator it = _map->begin() ; it != _map->end(); it++) {
			LOG_VERBOSE_OPTIONAL(1, true, "Closing " << it->first);
			it->second->close();
		}
	}
	virtual std::string getRank() const {
		return std::string();
	}
	std::string getFilename(std::string key) const {
		return _outputFilePathPrefix + key + _suffix + getRank();
	}
	std::ofstream &getOfstream(std::string key) {
		std::string filename = getFilename(key);
		// lockless lookup
		MapPtr thisMap = _map;
		Iterator it = thisMap->find(filename);
		if (it == thisMap->end()) {

			// lock if not found and map needs to be updated
			#pragma omp critical (ofStreamMap)
			{
				// recheck map
				thisMap = _map;
				it = thisMap->find(filename);
				if (it == thisMap->end()) {
					LOG_VERBOSE_OPTIONAL(1, true, "Writing to " << filename);
					std::ios_base::openmode mode = std::ios_base::out;
					if (getAppend())
						mode |= std::ios_base::app;
					else
						mode |= std::ios_base::trunc;
					OStreamPtr osp(new std::ofstream(filename.c_str(), mode));
					if( osp->fail() )
						throw std::runtime_error((std::string("Could not open file for writing: ") + filename).c_str());

					MapPtr copy = MapPtr(new Map(*thisMap));
					it = copy->insert( it, Map::value_type(filename, osp) );
					_map = copy;
				}
			}
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
    inline int getPartitionIdx(DataType score) const {
		// TODO binary search?
		for(unsigned int i = 0; i < _partitions.size(); i++)
			if (score < _partitions[i])
				return i;
		return _partitions.size();
	}
    inline Partitions getPartitions() const {
    	return _partitions;
    }

    inline int addPartition(DataType partition) {
    	_partitions.push_back(partition);
    	std::sort(_partitions.begin(), _partitions.end());
    	return _partitions.size();
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


typedef BucketedData<double, Kmernator::UI8> DoubleToOneByte;
typedef BucketedData<double, Kmernator::UI16> DoubleToTwoByte;

class SequenceRecordParser
{
public:
	static inline std::string &nextLine(std::string &buffer, Kmernator::RecordPtr &recordPtr) {
		Kmernator::RecordPtr nextPtr = strchr(recordPtr, '\n');
		long len = nextPtr - recordPtr;
		if (len > 0) {
		  buffer.assign(recordPtr, len);
		} else {
		  buffer.clear();
		}
	    recordPtr = ++nextPtr;
	    return buffer;
	}
	static std::string &trimName(std::string &nameLine) {
		if (nameLine.length() == 0)
			return nameLine;

		char &marker = nameLine[0];
		if (marker != '>' && marker != '@') {
			throw std::invalid_argument( (std::string("Can not parse name without a marker: ") + nameLine).c_str() );
		} else {
		    // remove marker
		    nameLine.erase(0,1);
		}

		// trim at first whitespace
		size_t pos = nameLine.find_first_of(" \t\r\n");
		if (pos != std::string::npos) {
			nameLine.erase(pos);
		}
		return nameLine;
	}
	static void parse(Kmernator::RecordPtr record, Kmernator::RecordPtr lastRecord,
			          std::string &name, std::string &bases, std::string &quals,
			          Kmernator::RecordPtr qualRecord = NULL, Kmernator::RecordPtr lastQualRecord = NULL,
			          char fastqStartChar = Kmernator::FASTQ_START_CHAR_ILLUMINA) {
		std::string buf;
		if (*record == '@') {
			// FASTQ
			nextLine(name,  record);
			trimName(name);
			nextLine(bases, record);
			nextLine(buf,   record);
			nextLine(quals, record);
			if (buf[0] != '+') {
				throw "Invalid FASTQ record!";
			}
		} else if (*record == '>') {
			// FASTA
			nextLine(name, record);
			trimName(name);

			// TODO FIXME HACK!
			if (lastRecord == NULL) // only read one line if no last record is given
				lastRecord = record+1;

			while (record < lastRecord && *record != '>') {
				bases += nextLine(buf, record);
			}
			if (qualRecord != NULL) {
				nextLine(buf, qualRecord);
				trimName(buf);
				if (buf != name) {
					throw "fasta and qual do not match names!";
				}

                // TODO FIXME HACK!
				if (lastQualRecord == NULL) // only read one line if no last record is given
					lastQualRecord = qualRecord+1;
				while (qualRecord < lastQualRecord && *qualRecord != '>') {
					quals += nextLine(buf, qualRecord);
				}
				quals = convertQualIntsToChars(quals, fastqStartChar);
			} else {
				quals.assign(bases.length(), Kmernator::REF_QUAL);
			}
		} else {
			throw "Do not know how to parse this file!";
		}
	}
	static std::string convertQualIntsToChars(const std::string &qualInts, char fastqStartChar) {
		std::istringstream ss(qualInts);
		std::ostringstream oss;
		while (!ss.eof()) {
			int qVal;
			ss >> qVal;
			if (ss.fail())
				break;
			qVal += fastqStartChar;
			oss << (char) qVal;
		}
        return oss.str();
	}

	static std::string commonName(const std::string &readName) {
		return readName.substr(0, readName.length() - 1);
	}
	static int readNum(const std::string &readName) {
		int retVal = 0;
		int len = readName.length();
		char c = readName[len - 1];
		switch (c) {
		case '1':
			if (readName[len - 2] != '/')
				break;
		case 'A':
		case 'F':
			retVal = 1;
			break;
		case '2':
			if (readName[len - 2] != '/')
				break;
		case 'B':
		case 'R':
			retVal = 2;
			break;
		}
		return retVal;
	}

	static bool isPair(const std::string &readNameA, const std::string readNameB) {
		std::string commonA = commonName(readNameA);
		std::string commonB = commonName(readNameB);
		if (commonA == commonB) {
			int readNumA = readNum(readNameA);
			int readNumB = readNum(readNameB);
			if (readNumA != 0 && readNumB != 0 && readNumA != readNumB) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}
};

class FileUtils
{
public:
	static unsigned long getFileSize(std::string &filePath) {
		std::ifstream ifs(filePath.c_str());
		return getFileSize(ifs);
	}
	static unsigned long getFileSize(std::ifstream &ifs) {
		assert( !ifs.eof() );
		assert( ifs.is_open() && ifs.good() );
		assert( !ifs.fail() );
		std::ifstream::streampos current = ifs.tellg();
		ifs.seekg(0, std::ios_base::end);
		unsigned long size = ifs.tellg();
		ifs.seekg(current);
		return size;
	}
	static bool fileExists(std::string &filePath) {
		std::ifstream ifs(filePath.c_str());
		return fileExists(ifs);
	}
	static bool fileExists(std::ifstream &ifs) {
		if (ifs.fail())
			return false;
		else
			return true;
	}
};

class Statistics
{
public:
	class MeanStdCount
	{
	public:
		double mean;
		double stdDev;
		long count;
		MeanStdCount() : mean(0.0), stdDev(0.0), count(0) {}

		template<typename forwardIterator>
		MeanStdCount(forwardIterator begin, forwardIterator end, bool isPoisson = true) : mean(0.0), stdDev(0.0), count(0) {
			for(forwardIterator it = begin; it != end ; it++) {
				count++;
				mean += *it;
			}
			if (count > 1) {
				mean /= (double) count;
				for(forwardIterator it = begin; it != end; it++) {
					double diff = ((double)*it) - mean;
					stdDev += diff*diff;
				}
				stdDev /= (double) (count-1);
			}
			if (isPoisson && mean > 0.0) {
				double poissonStdDev = sqrt(mean);
				if (stdDev < poissonStdDev) {
					stdDev = poissonStdDev;
				}
			}
		}
	};
	
	template<typename forwardIterator>
	static forwardIterator findBimodalPartition(double numSigmas, MeanStdCount &firstMSC, MeanStdCount &secondMSC, forwardIterator begin, forwardIterator end, bool isPoisson = true) {
		if (begin == end)
			return end;
		forwardIterator bestPart = end;
		double bestDiff = 0.0;
		forwardIterator part = begin;
		while (++part != end) {
			MeanStdCount first(begin, part, isPoisson);
			MeanStdCount second(part, end, isPoisson);
			if (first.count == 1 && second.count == 1)
				continue; // can not evaluate unless one sample is >1
			double diff = abs(first.mean - second.mean);
			// use the larger of the two standard deviations
			double stdDev = std::max(first.stdDev, second.stdDev);
			if (diff > numSigmas * stdDev) {	
				if (diff > bestDiff) {
					bestDiff = diff;
					bestPart = part;
					firstMSC = first;
					secondMSC = second;
				}
			}
		}
		return bestPart;
	};
	template<typename forwardIterator>
	static forwardIterator findBimodalPartition(double numSigmas, forwardIterator begin, forwardIterator end, bool isPoisson = true) {
		MeanStdCount f, s;
		return findBimodalPartition(numSigmas, f, s, begin, end, isPoisson);
	};
};

#endif

//
// $Log: Utils.h,v $
// Revision 1.40  2010-08-18 17:50:40  regan
// merged changes from branch FeaturesAndFixes-20100712
//
// Revision 1.39.4.1  2010-07-20 20:02:56  regan
// autodetect fastq quality range
//
// Revision 1.39  2010-06-22 23:06:30  regan
// merged changes in CorruptionBugfix-20100622 branch
//
// Revision 1.38.4.1  2010-06-22 23:02:20  regan
// named all critical sections
// added a critical section when modifying ostreamap data
//
// Revision 1.38  2010-05-24 21:48:46  regan
// merged changes from RNADedupMods-20100518
//
// Revision 1.37.2.1  2010-05-19 00:20:46  regan
// refactored fomat output options
// added options to fastq2fasta
//
// Revision 1.37  2010-05-18 20:50:24  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.36.2.2  2010-05-12 18:25:20  regan
// refactored
//
// Revision 1.36.2.1  2010-05-07 22:59:32  regan
// refactored base type declarations
//
// Revision 1.36  2010-05-06 21:46:54  regan
// merged changes from PerformanceTuning-20100501
//
// Revision 1.34.2.1  2010-05-04 19:49:51  regan
// minor rework on include headers
//
// Revision 1.35  2010-05-05 06:28:35  regan
// merged changes from FixPairOutput-20100504
//
// Revision 1.34.4.1  2010-05-05 05:57:53  regan
// fixed pairing
// fixed name to exclude labels and comments after whitespace
// applied some performance optimizations from other branch
// created FixPair application
//
// Revision 1.34  2010-05-01 21:57:54  regan
// merged head with serial threaded build partitioning
//
// Revision 1.33.2.5  2010-05-01 21:28:44  regan
// fixed constructor
//
// Revision 1.33.2.4  2010-04-28 22:28:10  regan
// refactored writing routines
//
// Revision 1.33.2.3  2010-04-27 18:25:20  regan
// bugfix in temp directory usage
//
// Revision 1.33.2.2  2010-04-26 22:56:36  regan
// honor temp-dir option
//
// Revision 1.33.2.1  2010-04-26 04:59:19  regan
// more output
//
// Revision 1.33  2010-04-21 23:39:17  regan
// added tempfile util
//
// Revision 1.32  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.31.2.3  2010-04-14 20:53:49  regan
// checkpoint and passes unit tests!
//
// Revision 1.31.2.2  2010-04-14 17:51:43  regan
// checkpoint
//
// Revision 1.31.2.1  2010-04-14 03:51:19  regan
// checkpoint. compiles but segfaults
//
// Revision 1.31  2010-03-04 06:37:42  regan
// fixed compiler warnings
//
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
