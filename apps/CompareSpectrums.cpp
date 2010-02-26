// $Header: /repository/PI_annex/robsandbox/KoMer/apps/CompareSpectrums.cpp,v 1.8 2010-02-26 13:01:15 regan Exp $
//

#include <iostream>

#include "config.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "KmerReadUtils.h"

//#include "Options.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

typedef TrackingDataMinimal4 DataType;
typedef KmerMap<DataType> KmerSolidMap;
typedef KmerSpectrum<DataType, DataType> KS;

class CS_Options {
public:
	typedef std::vector<std::string> FileListType;

private:
	static inline CS_Options &getOptions() {
		static CS_Options singleton;
		return singleton;
	}

	po::options_description desc;
	po::positional_options_description p;
	po::variables_map vm;

	// cache of variables (for inline lookup and defaults)
	FileListType referenceFiles;
	FileListType inputFiles;
	unsigned int kmerSize;
	bool perRead1;

public:

	static inline FileListType &getFileSet1() {
		return getOptions().referenceFiles;
	}
	static inline FileListType &getFileSet2() {
		return getOptions().inputFiles;
	}
	static inline unsigned int &getKmerSize() {
		return getOptions().kmerSize;
	}
	static inline bool &getPerRead1() {
		return getOptions().perRead1;
	}

	static bool parseOpts(int argc, char *argv[]) {
		try {
			po::options_description & desc = getOptions().desc;

			desc.add_options()("help", "produce help message")

			("file-set-1", po::value<FileListType>(), "1st file(s)")(
					"kmer-size", po::value<unsigned int>(), "kmer size")(
					"file-set-2", po::value<FileListType>(), "2nd file(s)")(
					"per-read1", po::value<bool>(),
					"compare for each read in set 1");

			po::positional_options_description & p = getOptions().p;
			p.add("kmer-size", 1);
			p.add("file-set-1", 1);
			p.add("file-set-2", -1);

			po::variables_map & vm = getOptions().vm;
			po::store(
					po::command_line_parser(argc, argv). options(desc).positional(
							p).run(), vm);
			po::notify(vm);

			if (vm.count("help")) {
				std::cerr << desc << std::endl;
				return false;
			}

			if (vm.count("file-set-1")) {
				std::cerr << "1st file set: ";
				FileListType & referenceFiles = getFileSet1()
						= vm["file-set-1"].as<FileListType> ();
				for (FileListType::iterator it = referenceFiles.begin(); it
						!= referenceFiles.end(); it++)
					std::cerr << *it << ", ";
				std::cerr << std::endl;
			} else {
				std::cerr << "file set 1 not specified" << std::endl;
				return false;
			}

			if (vm.count("file-set-2")) {
				std::cerr << "2nd file set: ";
				FileListType inputs = getFileSet2() = vm["file-set-2"].as<
						FileListType> ();
				for (FileListType::iterator it = inputs.begin(); it
						!= inputs.end(); it++)
					std::cerr << *it << ", ";
				std::cerr << std::endl;
			} else {
				std::cerr << desc << "file set 2 not specified!" << std::endl;
				return false;
			}

			if (vm.count("kmer-size")) {
				getKmerSize() = vm["kmer-size"].as<unsigned int> ();
				std::cerr << "Kmer size is: " << getKmerSize() << std::endl;
			} else {
				std::cerr << desc << "There was no kmer size specified!"
						<< std::endl;
				return false;
			}

			if (vm.count("per-read1")) {
				getPerRead1() = vm["per-read1"].as<bool> ();
				std::cerr << "per-read1: " << getPerRead1() << std::endl;
			} else {
				getPerRead1() = false;
			}
		} catch (std::exception& e) {
			std::cerr << "error: " << e.what() << std::endl;
			return false;
		} catch (...) {
			std::cerr << "Exception of unknown type!" << std::endl;
			return false;
		}

		return true;
	}
};

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;

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

int main(int argc, char *argv[]) {

	if (!CS_Options::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	ReadSet readSet1, readSet2;

	KmerSizer::set(CS_Options::getKmerSize());
	CS_Options::FileListType fileList1 = CS_Options::getFileSet1();
	CS_Options::FileListType fileList2 = CS_Options::getFileSet2();

	cerr << "Reading 1st file set:";
	readSet1.appendAllFiles(fileList1);

	cerr << " loaded " << readSet1.getSize() << " Reads, "
			<< readSet1.getBaseCount() << " Bases " << endl;

	cerr << "Reading 2nd file set:";
	readSet2.appendAllFiles(fileList2);
	cerr << " loaded " << readSet2.getSize() << " Reads, "
			<< readSet2.getBaseCount() << " Bases " << endl;

	long buckets = std::max(KS::estimateWeakKmerBucketSize(readSet1),
			KS::estimateWeakKmerBucketSize(readSet2));

	KS ks1(buckets);
	KS ks2(buckets);

	cerr << "Building map 2\n";
	ks2.buildKmerSpectrum(readSet2, true);

	cout << endl;
	cout << "Set 1\tSet 2\tCommon\t%Uniq1\t%Tot1\t%Uniq2\t%Tot2\n";

	if (CS_Options::getPerRead1()) {
#pragma omp parallel for schedule(dynamic)
		for (long i = 0; i < readSet1.getSize(); i++) {
			ReadSet tmpSet;
			Read &read = readSet1.getRead(i);
			tmpSet.append(read);
			KS tmpKs1(buckets);
			tmpKs1.buildKmerSpectrum(tmpSet, true);
			KmerSolidMap &m1 = tmpKs1.solid;
			KmerSolidMap &m2 = ks2.solid;
			NumbersVector common = countCommonKmers(m1, m2);
			if (common[0] > 0) {
#pragma omp critical
				{
					cout << m1.size() << '\t' << m2.size() << '\t' << common[0]
							<< '\t' << setprecision(4) << (common[0] * 100.0)
							/ m1.size() << '\t' << (common[1] * 100.0
							/ common[3]) << "\t" << (common[0] * 100.0)
							/ m2.size() << "\t" << (common[2] * 100.0
							/ common[4]) << "\t" << read.getName() << endl;
				}
			}
		}
	} else {
		cerr << "Building map 1\n";
		ks1.buildKmerSpectrum(readSet1, true);
		TrackingData::resetGlobalCounters();

		KmerSolidMap &m1 = ks1.solid;
		KmerSolidMap &m2 = ks2.solid;

		cerr << "Counting common Kmers\n";
		NumbersVector common = countCommonKmers(m1, m2);

		// TODO fix check common to iterater through same-bucketed KmerMaps in sorted order (fast)
		// TODO if perRead1, iterate through common matches and output non-zero % matches

		cout << m1.size() << '\t' << m2.size() << '\t' << common[0] << '\t'
				<< setprecision(4) << (common[0] * 100.0) / m1.size() << '\t'
				<< (common[1] * 100.0 / common[3]) << "\t" << (common[0]
				* 100.0) / m2.size() << "\t" << (common[2] * 100.0 / common[4])
				<< endl;
	}
}

//
// $Log: CompareSpectrums.cpp,v $
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

