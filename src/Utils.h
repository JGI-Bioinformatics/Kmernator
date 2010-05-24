// $Header: /repository/PI_annex/robsandbox/KoMer/src/Utils.h,v 1.38 2010-05-24 21:48:46 regan Exp $
//

#ifndef _UTILS_H
#define _UTILS_H

#include <ostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cstring>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#define foreach BOOST_FOREACH

#include "config.h"
#include "Options.h"


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
private:
    Map _map;
    std::string _outputFilePathPrefix;
    std::string _suffix;

public:
	OfstreamMap(std::string outputFilePathPrefix = Options::getOutputFile(), std::string suffix = FormatOutput::getDefaultSuffix())
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
			if (Options::getVerbosity()) {
			  std::cerr << "Writing to " << filename << std::endl;
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


typedef BucketedData<double, KoMer::UI8> DoubleToOneByte;
typedef BucketedData<double, KoMer::UI16> DoubleToTwoByte;

class SequenceRecordParser
{
public:
	static inline std::string &nextLine(std::string &buffer, KoMer::RecordPtr &recordPtr) {
		KoMer::RecordPtr nextPtr = strchr(recordPtr, '\n');
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
	static void parse(KoMer::RecordPtr record, KoMer::RecordPtr lastRecord,
			          std::string &name, std::string &bases, std::string &quals,
			          KoMer::RecordPtr qualRecord = NULL, KoMer::RecordPtr lastQualRecord = NULL) {
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
				quals = convertQualIntsToChars(quals);
			} else {
				quals.assign(bases.length(), KoMer::REF_QUAL);
			}
		} else {
			throw "Do not know how to parse this file!";
		}
	}
	static std::string convertQualIntsToChars(const std::string &qualInts) {
		std::istringstream ss(qualInts);
		std::ostringstream oss;
		while (!ss.eof()) {
			int qVal;
			ss >> qVal;
			if (ss.fail())
				break;
			qVal += KoMer::FASTQ_START_CHAR;
			oss << (char) qVal;
		}
        return oss.str();
	}
};

#endif

//
// $Log: Utils.h,v $
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
