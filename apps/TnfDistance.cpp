// $Header: /repository/PI_annex/robsandbox/KoMer/apps/TnfDistance.cpp,v 1.2 2010-08-18 17:47:49 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>

#include "config.h"
#include "Options.h"
#include "ReadSet.h"
#include "KmerSpectrum.h"
#include "Utils.h"


using namespace std;
typedef TrackingDataMinimal4f DataType;
typedef KmerMap<DataType> KM;
typedef KmerSpectrum<DataType, DataType> KS;

class TnfDistanceOptions : public Options {
public:
	static string getInterFile() {
		return getVarMap()["inter-distance-file"].as<string> ();
	}
	static string getIntraFile() {
		return getVarMap()["intra-distance-file"].as<string> ();
	}
	static int getIntraWindow() {
		return getVarMap()["intra-window-size"].as<int> ();
	}
	static string getClusterFile() {
		return getVarMap()["cluster-file"].as<string>();
	}
	static double getClusterThreshold() {
		return getVarMap()["cluster-threshold-distance"].as<double>();
	}
	static bool parseOpts(int argc, char *argv[]) {
		// set options specific to this program
		getPosDesc().add("input-file", -1);
		getDesc().add_options()

		("inter-distance-file", po::value<string>()->default_value(""),
				"output inter-distance LT matrix to this filename")

		("intra-distance-file", po::value<string>()->default_value(""),
				"output intra-distance matrix (one line per read) to this filename")

		("intra-window-size", po::value<int>()->default_value(25000),
				"size of adjacent intra-distance windows")

		("cluster-file", po::value<string>()->default_value(""),
				"cluster output filename")

		("cluster-threshold-distance", po::value<double>()->default_value(0.175),
				"Euclidean distance threshold for clusters")
		;

		bool ret = Options::parseOpts(argc, argv);

		return ret;
	}

};

class TNF {
public:
	typedef vector<double> Vector;
	static KM stdMap;
	static long stdSize;
	static void initStdMap() {
	  for(int i=0; i < 256 ; i++) {
		char c = (char) i;
		Kmer &kmer = (Kmer&) c;
		TEMP_KMER(lc);
		kmer.buildLeastComplement(lc);
		if ( !stdMap.exists(lc) ) {
			if (Options::getDebug() > 0) cerr << lc.toFasta() << "\t" << (int) i << std::endl;
		    stdMap.insert(lc, DataType());
		}
	  }
	  stdSize = stdMap.size();
	}
	static std::string toHeader() {
		std::stringstream ss;
		bool printTab = false;
		for(KM::Iterator it = stdMap.begin(); it != stdMap.end(); it++) {
			if (printTab) ss << "\t";
			ss << it->key().toFasta();
			printTab = true;
		}
		return ss.str();
	}
private:
	double length;
	Vector tnfValues;

public:
	TNF() : length(1.0) {
		tnfValues.resize(stdSize);
	}
	TNF(const TNF &copy) {
		*this = copy;
	}
	TNF(const KM &map, bool normalize = false) : length(1.0) {
		buildTnf(map);
		if (normalize)
			buildLength();
	}
	TNF &operator=(const TNF &copy) {
		if (this == &copy)
			return *this;
		length = copy.length;
		tnfValues = copy.tnfValues;
		return *this;
	}
	TNF &operator*(double factor) {
		bool wasNormalized = isNormalized();
		for(unsigned int i = 0; i < tnfValues.size(); i++) {
			tnfValues[i] *= factor;
		}
		if (wasNormalized) {
			buildLength();
		}
		return *this;
	}
	TNF &operator+(const TNF &rh) {
		bool wasNormalized = isNormalized();

		for(unsigned int i = 0; i < tnfValues.size(); i++) {
			tnfValues[i] = tnfValues[i] + rh.tnfValues[i];
		}
		if (wasNormalized) {
			buildLength();
		}
		return *this;
	}
	std::string toString() const {
		std::stringstream ss;
		for(unsigned int i = 0 ; i < tnfValues.size(); i++) {
			if (i != 0) ss << "\t";
			ss << tnfValues[i] / length;
		}
		return ss.str();
	}
	bool isNormalized() const {
		return length != 1.0;
	}
	void buildLength() {
		length = 0.0;
	    bool debug = Options::getDebug() > 1;
	    if (debug) cerr << "buildLength()\t";
	    for(unsigned int i = 0 ; i < tnfValues.size(); i++) {
	    	length += tnfValues[i] * tnfValues[i];
	    	if (debug) cerr << tnfValues[i] << "\t";
	    }
	    length = sqrt(length);
	    if (length == 0.0)
	    	length = 1.0;
	    if (debug) cerr << length << endl;
	}
	void buildTnf(const KM & map) {
		tnfValues.empty();
		tnfValues.reserve(stdSize);
	    for(KM::Iterator it = stdMap.begin(); it != stdMap.end(); it++) {
	    	KM::ElementType elem = map.getElementIfExists(it->key());
	    	double t = 0.0;
	    	if (elem.isValid()) {
	    		t = elem.value() / length;
	    	}
	    	tnfValues.push_back(t);
	    }
	}
	double getDistance(const TNF &other) const {
		bool debug = Options::getDebug() > 1;
		double dist = 0.0;

		if (debug) cerr << "getDistance()\t";
		for(unsigned int i = 0 ; i < tnfValues.size(); i++) {
			double d = tnfValues[i] / length - other.tnfValues[i] / other.length;
			dist += d*d;
			if (debug) cerr << tnfValues[i] / length << "\t" << other.tnfValues[i] / other.length << "\t" << d*d << endl;
		}
		dist = sqrt(dist);
		if (debug) cerr << dist << endl;
		return dist;
	}
	double getDistance(const KM &query) const {
		TNF other(query, isNormalized());
		return getDistance(other);
	}
	static double getDistance(const KM &target, const KM &query, bool normalize = true) {
		TNF tgt(target, normalize), qry(query, normalize);
		return tgt.getDistance(qry);
	}
	static double getDistance(KM &kmers, KM &target, KM &query) {
		double distance = 0.0;

	    double targetSum = 0.0;
	    double querySum = 0.0;
	    bool debug = Options::getDebug() > 1;

	    if (debug) cerr << "getDistance(...)" << endl;
	    for(KM::Iterator it = kmers.begin(); it != kmers.end(); it++) {
	    	double t = target[it->key()];
	    	targetSum += t*t;
	    	double q = query[it->key()];
	    	querySum += q*q;
			if (debug) cerr << t << "\t" << q << endl;
	    }
	    if (debug) cerr << endl;
	    targetSum = sqrt(targetSum);
	    querySum = sqrt(querySum);
	    if (debug) cerr << targetSum << "\t" << querySum << endl;

	    for(KM::Iterator it = kmers.begin(); it != kmers.end(); it++) {
	    	double t = target[it->key()] / targetSum;
	    	double q = query[it->key()] / querySum;
	    	double d = t-q;
	    	distance += d*d;
	    	if (debug) cerr << t << "\t" << q << "\t" << d*d << endl;
	    }
	    distance = sqrt(distance);
	    if (debug) cerr << distance << endl;

		return distance;
	}

};
KM TNF::stdMap;
long TNF::stdSize;

typedef std::vector<TNF> TNFS;
TNFS buildTnfs(const ReadSet &reads, bool normalize = false) {
	TNFS tnfs;
	long size = reads.getSize();
	tnfs.resize(size);

	#pragma omp parallel for
	for(long readIdx = 0; readIdx < size; readIdx++) {
		KS ksReads;
		ksReads.setSolidOnly();
		ReadSet singleRead;
		Read read = reads.getRead(readIdx);
		if (Options::getDebug()) {
			cerr << "buildTnfs() for: " << read.getName() << endl;
		}
		singleRead.append(read);

		ksReads.buildKmerSpectrum(singleRead, true);

		TNF tnf(ksReads.solid, normalize);
		tnfs[readIdx] = tnf;
	}
	return tnfs;
}

TNFS buildIntraTNFs(const Read &read, long window, long step, bool normalize = false) {
	string fullFasta = read.getFasta();
	long length = fullFasta.length();
	ReadSet reads;
	for(long i = 0 ; i < length - window; i+= step) {
		string fasta = fullFasta.substr(i, window);
		Read newRead(read.getName() + ":" + boost::lexical_cast<string>(i) + "-" + boost::lexical_cast<string>(i+fasta.length()), fasta, string(fasta.length(), KoMer::PRINT_REF_QUAL));
		reads.append(newRead);
	}
	return buildTnfs(reads, normalize);
}

class Result {
public:
  double distance;
  string label;
  TNF tnf;
  Result() : distance(0.0) {}
  Result(double _dist, string _label) : distance(_dist), label(_label) {}
  Result(double _dist, string _label, TNF _tnf) :  distance(_dist), label(_label), tnf(_tnf) {}
  Result(const Result &copy) {
	  *this = copy;
  }
  Result &operator=(const Result &copy) {
	  if (this == &copy)
		  return *this;
	  distance = copy.distance;
	  label = copy.label;
	  tnf = copy.tnf;
	  return *this;
  }
};
typedef std::vector< Result > Results;
class ResultCompare {
public:
	bool operator() (const Result &a, const Result &b) { return (a.distance < b.distance); }
};


Results calculateDistances(const TNF &refTnf, const ReadSet &reads, const TNFS &tnfs) {
	Results results;
	long size = reads.getSize();
	results.resize(size);

	#pragma omp parallel for
	for(long readIdx = 0; readIdx < size; readIdx++) {

		ReadSet singleRead;
		Read read = reads.getRead(readIdx);
		string name = read.getName();

		double dist = refTnf.getDistance(tnfs[readIdx]);
		if (dist <= 0.20) {
			string fullFasta = read.getFasta();
			long len = fullFasta.length();
			int divs = 5;
			TNFS intraTnfs = buildIntraTNFs(read,len/divs, len/divs, true);

			stringstream ss;
			ss.precision(3);
			ss << name;
			for(long i = 0 ; i < (long) intraTnfs.size() ; i++) {
				double distPart = refTnf.getDistance(intraTnfs[i]);
				ss << "\t" << fixed << distPart;
			}
			name = ss.str();
		}
		results[readIdx] = Result( dist, name );
	}
	return results;
}
Results calculateDistances(const TNF &refTNF, const ReadSet &reads) {
	return calculateDistances(refTNF, reads, buildTnfs(reads));
}

typedef vector<double> DVector;
typedef vector< DVector > DMatrix;
class MinVecOper {
public:
	bool operator()( const vector<double>::iterator &a, const vector<double>::iterator &b) const {
		return (*a < *b);
	}
} mvo;

int main(int argc, char *argv[]) {
	if (!TnfDistanceOptions::parseOpts(argc, argv))
		throw std::invalid_argument("Please fix the command line arguments");

	MemoryUtils::getMemoryUsage();

	ReadSet refs;
	ReadSet reads;
	KmerSizer::set(4);
	TNF::initStdMap();

	Options::FileListType inputs = Options::getInputFiles();
	reads.appendAllFiles(inputs);

	KS ksRef;
	ksRef.setSolidOnly();

	TNFS readTnfs = buildTnfs(reads, true);

	Options::FileListType referenceInputs = Options::getReferenceFiles();
	if (!referenceInputs.empty()) {
		// compare distances from reference to each read in the input
		refs.appendAllFiles(referenceInputs);
		ksRef.buildKmerSpectrum(refs, true);
		TNF refTnf(ksRef.solid, true);

		Results results = calculateDistances(refTnf, reads, readTnfs);
		std::sort(results.begin(), results.end(), ResultCompare());
		for(Results::iterator it = results.begin(); it != results.end(); it++) {
			cout << it->distance << "\t" << it->label << endl;
		}
	} else {
		// output TNF matrix for each read in the input

		cout << "Label\t" << TNF::toHeader() << endl;
		for(ReadSet::ReadSetSizeType readIdx = 0; readIdx < reads.getSize(); readIdx++) {

			Read read = reads.getRead(readIdx);

			cout << read.getName() << "\t" << readTnfs[readIdx].toString() << endl;
		}

	}

	string interFile = TnfDistanceOptions::getInterFile();
	if (!interFile.empty()) {
		OfstreamMap om(interFile, "");
		ostream &os = om.getOfstream("");
		for(ReadSet::ReadSetSizeType readIdxi = 0; readIdxi < reads.getSize(); readIdxi++) {
			os << reads.getRead(readIdxi).getName();
			for(ReadSet::ReadSetSizeType readIdxj = 0; readIdxj < readIdxi; readIdxj++) {
				double dist = readTnfs[readIdxi].getDistance(readTnfs[readIdxj]);
				os << "\t" << dist;
			}
			os << endl;
		}
	}

	string intraFile = TnfDistanceOptions::getIntraFile();
	if (!intraFile.empty()) {
		OfstreamMap om(intraFile, "");
		ostream &os = om.getOfstream("");
		long window = TnfDistanceOptions::getIntraWindow();
		long step = window / 10;
		for(ReadSet::ReadSetSizeType readIdx = 0; readIdx < reads.getSize(); readIdx++) {
			Read read = reads.getRead(readIdx);
			TNFS intraTnfs = buildIntraTNFs(read, window, step, true);
			os << read.getName();
			for(unsigned long i = 0 ; i < intraTnfs.size(); i++) {
				for(unsigned long j = 0 ; j < i ; j++) {
					double dist = intraTnfs[i].getDistance(intraTnfs[j]);
					os << "\t" << dist;
				}
			}
			os << endl;
		}
	}

	string clusterFile = TnfDistanceOptions::getClusterFile();
	if (!clusterFile.empty()) {
		bool debug = Options::getDebug() > 0;
		long size = readTnfs.size();

		cerr << "Building clusters" << endl;
		OfstreamMap om(clusterFile, "");
		ostream &os = om.getOfstream("");
		DMatrix distMatrix;
		distMatrix.resize( size );
		vector< vector< string > > clusterNames;
		clusterNames.resize( size );

		//TODO optimize this and the the while loop (use LT, and directed updates to minVec)
		vector< vector<double>::iterator > minVec;
		minVec.resize( size );
		float clusterThreshold = TnfDistanceOptions::getClusterThreshold();


        #pragma omp parallel for schedule(dynamic)
		for(long i = 0; i < size; i++) {
			clusterNames[i].push_back( reads.getRead(i).getName() );
			DVector &vect = distMatrix[i];
			vect.resize(i+1);
			for(long j = 0 ; j <= i; j++) {
				if (j==i)
					vect[j] = clusterThreshold+1.0;
				else
					vect[j] = readTnfs[i].getDistance(readTnfs[j]);
			}
			minVec[i] = std::min_element( vect.begin(), vect.end() );
		}

		cerr << "Finished dist matrix" << endl;
		if (debug) {
			for(long i = 0 ; i < size ; i++) {
				cerr << "distMatrix: " << i << "\t" << clusterNames[i][0];
				for(long j = 0 ; j <= i ; j++)
					cerr << "\t" << distMatrix[i][j];
				cerr << "\n";
			}
		}


		TNFS readTnfs2 = TNFS(readTnfs);

		while( true ) {
			// set the thresholds for each remaining read

			vector< vector<double>::iterator >::iterator minElem = std::min_element(minVec.begin(), minVec.end(), mvo);
			if (**minElem > clusterThreshold)
				break;
			if (debug) cerr << **minElem << endl;

			unsigned long tgtj = minElem - minVec.begin();
			unsigned long tgti = *minElem - distMatrix[tgtj].begin();

			unsigned long lastIdx = minVec.size() - 1;
			if (debug) {
				for(unsigned long i2 = 0 ; i2 < minVec.size(); i2++) {
					unsigned long j = i2 <= tgti ? i2 : tgti;
					unsigned long i = i2 <= tgti ? tgti : i2;

					cerr << distMatrix[i][j] << "\t";
				}
				cerr << endl;
			}
			clusterNames[tgti].insert( clusterNames[tgti].end(), clusterNames[tgtj].begin(), clusterNames[tgtj].end() );
			if (Options::getVerbosity()) {
				cerr << "Cluster merged " << tgti << " with " << tgtj << " " << distMatrix[tgtj][tgti] << " " << **minElem << endl;
			}
			readTnfs2[tgti] = readTnfs2[tgti] + readTnfs2[tgtj];

			// remove tgtj from vectors and matrix (1 row and 1 column);
			std::swap(minVec[tgtj], minVec[lastIdx]); minVec.pop_back();
			std::swap(clusterNames[tgtj], clusterNames[lastIdx]); clusterNames.pop_back();
			std::swap(distMatrix[tgtj], distMatrix[lastIdx]); distMatrix.pop_back();
			std::swap(readTnfs2[tgtj], readTnfs2[lastIdx]); readTnfs2.pop_back();

			for(unsigned long i = tgtj + 1; i < minVec.size(); i++) {
				float dist = distMatrix[tgtj][i];
				distMatrix[i][tgtj] = dist;

				if (*minVec[i] > dist) {
					minVec[i] = distMatrix[i].begin() + tgtj;
				}
			}
			distMatrix[tgtj].resize(tgtj+1);
			distMatrix[tgtj][tgtj] = clusterThreshold + 1.0;
			minVec[tgtj] = std::min_element( distMatrix[tgtj].begin(), distMatrix[tgtj].end() );

			// update row & column tgti and affected minVects
			for(unsigned long i = 0 ; i < minVec.size(); i++) {
				if (i == tgti)
					distMatrix[tgti][tgti] = clusterThreshold+1.0;
				else {
					float dist = readTnfs2[tgti].getDistance(readTnfs2[i]);
					if (i < tgti) {
						distMatrix[tgti][i] = dist;
					} else {
						distMatrix[i][tgti] = dist;
						if (*minVec[i] > dist) {
							minVec[i] = distMatrix[i].begin() + tgti;
						}
					}
				}
			}
			minVec[tgti] = std::min_element( distMatrix[tgti].begin(), distMatrix[tgti].end() );

		}

		for(unsigned long i = 0 ; i < minVec.size(); i++) {
			for(unsigned long j = 0 ; j < clusterNames[i].size(); j++)
				os << clusterNames[i][j] << "\t";
			if (!clusterNames[i].empty())
				os << endl;
		}

	}

	return 0;
}

// $Log: TnfDistance.cpp,v $
// Revision 1.2  2010-08-18 17:47:49  regan
// added TNF distance application
//
// Revision 1.1.2.10  2010-08-02 22:10:50  regan
// refactored and optimized/parallelized cluster building
//
// Revision 1.1.2.9  2010-07-26 20:44:00  regan
// parallelized some parts
//
// Revision 1.1.2.8  2010-07-16 17:05:56  regan
// fixed debug output levels
//
// Revision 1.1.2.7  2010-07-16 16:38:33  regan
// added clustering option and output
//
// Revision 1.1.2.6  2010-07-16 00:25:42  regan
// refactored and refeatured.. added inter, intra options
//
// Revision 1.1.2.5  2010-07-14 21:13:50  regan
// added mode to print TNF
//
// Revision 1.1.2.4  2010-07-14 21:02:13  regan
// refactor
//
// Revision 1.1.2.3  2010-07-13 21:39:34  regan
// changed options
//
// Revision 1.1.2.2  2010-07-13 19:45:37  regan
// bugfix
//
// Revision 1.1.2.1  2010-07-13 17:39:11  regan
// created TnfDistance application
//
//
