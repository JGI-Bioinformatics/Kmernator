// $Header: /repository/PI_annex/robsandbox/KoMer/apps/CompareSpectrums.cpp,v 1.10 2010-05-18 20:50:18 regan Exp $
//

#include <iostream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "config.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "KmerReadUtils.h"

typedef TrackingDataMinimal4 DataType;
typedef KmerMap<DataType> KmerSolidMap;
typedef KmerSpectrum<DataType, DataType> KS;

using namespace std;


class CS_Options : public Options {

	// cache of variables (for inline lookup and defaults)

public:

	static bool getCircularReference() { return  getVarMap()["circular-reference"].as<unsigned int>() != 0; }
	static bool getPerRead() { return getVarMap()["per-read"].as<unsigned int>() != 0; }

	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("kmer-size", 1);
		getPosDesc().add("reference-file", 1);
		getPosDesc().add("input-file", -1);

		getDesc().add_options()
		 ("circular-reference", po::value<unsigned int>()->default_value(0),
				 "reference file should be treated as circular")
		 ("per-read", po::value<unsigned int>()->default_value(0),
				 "if set, each read in readset1 will be compared to the entire readset2 separately")
	    ;

		bool ret = Options::parseOpts(argc, argv);

		if (ret) {
		}
		return ret;
	}

};

typedef std::vector<unsigned long> NumbersVector;

// returns:
//  0: the common count of unique kmers
//  1: the cumulative count of the common kmers from map 1
//  2: the cumulative count of the common kmers from map 2
//  3: the total cumulative count of all kmers from map 1
//  4: the total cumulative count of all kmers from map 2
static NumbersVector countCommonKmers(KmerSolidMap &m1, KmerSolidMap &m2) {
	NumbersVector ret(5);

	for (KmerSolidMap::Iterator it(m1.begin()), itEnd(m1.end()); it != itEnd; it++) {
		ret[3] += it->value().getCount();
		KmerSolidMap::ElementType element = m2.getElementIfExists(it->key());
		if (element.isValid()) {
			ret[0]++;
			ret[1] += it->value().getCount();
			ret[2] += element.value().getCount();
		}
	}
	for (KmerSolidMap::Iterator it(m2.begin()), itEnd(m2.end()); it != itEnd; it++) {
		ret[4] += it->value().getCount();
	}
	return ret;
}

void evaluate(KS &ks1, KS &ks2, std::string label = "");
void evaluatePerRead(KS &ks1, KS &ks2, ReadSet &readSet1);

int main(int argc, char *argv[]) {

	if (!CS_Options::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	MemoryUtils::getMemoryUsage();
	ReadSet readSet1, readSet2;

	KmerSizer::set(CS_Options::getKmerSize());
	CS_Options::FileListType fileList1 = CS_Options::getReferenceFiles();
	CS_Options::FileListType fileList2 = CS_Options::getInputFiles();

	cerr << "Reading 1st file set:";
	readSet1.appendAllFiles(fileList1);

	cerr << " loaded " << readSet1.getSize() << " Reads, "
			<< readSet1.getBaseCount() << " Bases " << endl;

	if (CS_Options::getCircularReference())
		readSet1.circularize(KmerSizer::getSequenceLength());

	cerr << "Reading 2nd file set:";
	readSet2.appendAllFiles(fileList2);
	cerr << " loaded " << readSet2.getSize() << " Reads, "
			<< readSet2.getBaseCount() << " Bases " << endl;

	long buckets = std::max(KS::estimateWeakKmerBucketSize(readSet1),
			KS::estimateWeakKmerBucketSize(readSet2));

	cerr << "Estimated bucket size: " << buckets << std::endl;

	KS ks1;
	KS ks2;
	ks1.setSolidOnly(buckets);
	ks2.setSolidOnly(buckets);

	cerr << "Building map 2\n";
	ks2.buildKmerSpectrum(readSet2, true);

	cout << endl;
	cout << "Set 1\tSet 2\tCommon\t%Uniq1\t%Tot1\t%Uniq2\t%Tot2\n";


	if (CS_Options::getPerRead()) {
		evaluatePerRead(ks1, ks2, readSet1);
	} else {
		cerr << "Building map 1\n";
		ks1.buildKmerSpectrum(readSet1, true);
		evaluate(ks1, ks2);
	}
}

void evaluate(KS &ks1, KS &ks2, std::string label) {
	cerr << "Counting common Kmers\n";
	KmerSolidMap &m1 = ks1.solid;
	KmerSolidMap &m2 = ks2.solid;

	NumbersVector common = countCommonKmers(m1, m2);

	// TODO fix check common to iterater through same-bucketed KmerMaps in sorted order (fast)
	// TODO if perRead1, iterate through common matches and output non-zero % matches

	cout << m1.size() << '\t' << m2.size() << '\t' << common[0] << '\t'
			<< setprecision(4) << (common[0] * 100.0) / m1.size() << '\t'
			<< (common[1] * 100.0 / common[3]) << "\t" << (common[0]
			* 100.0) / m2.size() << "\t" << (common[2] * 100.0 / common[4])
			<< '\t' << label
			<< endl;

}
void evaluatePerRead(KS &ks1, KS &ks2, ReadSet &readSet1) {
	for(ReadSet::ReadSetSizeType readIdx = 0; readIdx < readSet1.getSize(); readIdx++ ) {

		TrackingData::resetGlobalCounters();
		ks1.reset(false);

		ReadSet a;
		Read &read = readSet1.getRead(readIdx);
		a.append(read);
		ks1.buildKmerSpectrum(a, true);

		evaluate(ks1,ks2,read.getName());
	}
}
//
// $Log: CompareSpectrums.cpp,v $
// Revision 1.10  2010-05-18 20:50:18  regan
// merged changes from PerformanceTuning-20100506
//
// Revision 1.9.20.4  2010-05-13 20:29:13  regan
// minor refactor
//
// Revision 1.9.20.3  2010-05-12 22:44:36  regan
// reworked option handling
//
// Revision 1.9.20.2  2010-05-10 17:57:11  regan
// fixing types
//
// Revision 1.9.20.1  2010-05-07 22:59:29  regan
// refactored base type declarations
//
// Revision 1.9  2010-03-04 06:36:36  regan
// fixed compiler warnings for non-openmp compilers
//
// Revision 1.8  2010-02-26 13:01:15  regan
// reformatted
//
// Revision 1.7  2010-02-22 14:41:31  regan
// checkpoint
//
// Revision 1.6  2010-01-14 19:27:43  regan
// bugfixes
//
// Revision 1.5  2010-01-14 00:46:51  regan
// refactor and other changes
//
// Revision 1.4  2010-01-08 18:34:51  regan
// refactored a bit
// enabled openmp parallelizations
//
// Revision 1.3  2010-01-05 06:44:37  regan
// fixed warnings
//
// Revision 1.2  2009-12-24 00:54:08  regan
// fixed reading of fasta files
// parallelized reading of multiple files
//
// Revision 1.1  2009-11-27 23:28:07  cfurman
// CompareSpectrum application added
//

