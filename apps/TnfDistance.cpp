//
// Kmernator/apps/TnfDistance.cpp
//
// Author: Rob Egan
//
/*****************

Kmernator Copyright (c) 2012, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of any
required approvals from the U.S. Dept. of Energy).  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley
National Laboratory, U.S. Dept. of Energy nor the names of its contributors may
be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to Lawrence Berkeley National
Laboratory, without imposing a separate written license agreement for such
Enhancements, then you hereby grant the following license: a  non-exclusive,
royalty-free perpetual license to install, use, modify, prepare derivative
works, incorporate into other computer software, distribute, and sublicense
such enhancements or derivative works thereof, in binary and source code form.

*****************/

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
#include "Log.h"


using namespace std;
typedef TrackingDataMinimal4f DataType;
typedef KmerMap<DataType> KM;
typedef KmerSpectrum<KM, KM> KS;

struct DistanceFormula {
	enum Enum{EUCLIDEAN = 0, SPEARMAN, MAX};
};

class _TnfDistanceBaseOptions : public OptionsBaseInterface {
public:
	_TnfDistanceBaseOptions() :
		kmerSize(4),
		interDistanceFile(),
		intraDistanceFile(),
		clusterFile(),
		intraWindowSize(25000),
		clusterThresholdDistance(.175),
		referenceFiles(),
		distanceFormula(DistanceFormula::EUCLIDEAN) {}
	virtual ~_TnfDistanceBaseOptions() {}

	int &getKmerSize() {
		return kmerSize;
	}

	string &getInterFile() {
		return interDistanceFile;
	}
	string &getIntraFile() {
		return intraDistanceFile;
	}
	string &getClusterFile() {
		return clusterFile;
	}
	int &getIntraWindow() {
		return intraWindowSize;
	}
	float &getClusterThreshold() {
		return clusterThresholdDistance;
	}

	enum DistanceFormula::Enum &getDistanceFormula() {
		return distanceFormula;
	}

	FileListType &getReferenceFiles()
	{
		return referenceFiles;
	}


	void _resetDefaults() {
		GeneralOptions::_resetDefaults();
		GeneralOptions::getOptions().getVerbose() = 0;
		GeneralOptions::getOptions().getMmapInput() = false;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		// set options specific to this program
		po::options_description opts("Tetra Nucleotide Frequency (TNF) Distance Options");
		opts.add_options()

		        ("kmer-size", po::value<int>()->default_value(kmerSize), "set to 4 for tetra-mer, 5 penta-mer, 1 for AT/GC ...")

				("reference-file", po::value<FileListType>(), "set reference file(s).  If set calculate distance of each input to this reference TNF vector")

				("inter-distance-file", po::value<string>()->default_value(interDistanceFile), "output inter-distance LT matrix of inputs to this filename")

				("intra-distance-file", po::value<string>()->default_value(intraDistanceFile), "output intra-distance matrix (one line per input) to this filename")

				("intra-window-size", po::value<int>()->default_value(intraWindowSize), "size of adjacent intra-distance windows")

				("cluster-file", po::value<string>()->default_value(clusterFile), "cluster output filename")

				("cluster-threshold-distance", po::value<float>()->default_value(clusterThresholdDistance), "Euclidean distance threshold for clusters")

				("distance-formula", po::value<int>()->default_value(distanceFormula), "0 - Euclidean, 1 - Spearman")

				;


		desc.add(opts);

		GeneralOptions::_setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		setOpt("kmer-size", kmerSize);
		setOpt("inter-distance-file", interDistanceFile);
		setOpt("intra-distance-file", intraDistanceFile);
		setOpt("intra-window-size", intraWindowSize);
		setOpt("cluster-file", clusterFile);
		setOpt("cluster-threshold-distance", clusterThresholdDistance);
		setOpt2("reference-file", referenceFiles);
		int df;
		setOpt("distance-formula", df);
		if (df < 0 || df >= DistanceFormula::MAX) {
			setOptionsErrorMsg("Invalid --distance-formula.  Please choose either 0 or 1");
			ret = false;
		} else {
			distanceFormula = (DistanceFormula::Enum) df;
		}

		ret |= GeneralOptions::_parseOptions(vm);
		return ret;
	}
private:
	int kmerSize;
	string interDistanceFile, intraDistanceFile, clusterFile;
	int intraWindowSize;
	float clusterThresholdDistance;
	FileListType referenceFiles;
	DistanceFormula::Enum distanceFormula;

};
typedef OptionsBaseTemplate< _TnfDistanceBaseOptions > TnfDistanceBaseOptions;

class _TnfDistanceOptions : public OptionsBaseInterface {
public:
	_TnfDistanceOptions() {}
	virtual ~_TnfDistanceOptions() {}
	void _resetDefaults() {
		TnfDistanceBaseOptions::_resetDefaults();
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		p.add("input-file", -1);
		po::options_description opts("TnfDistance <options> input.fa");
		desc.add(opts);

		TnfDistanceBaseOptions::_setOptions(desc,p);

	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		ret |= TnfDistanceBaseOptions::_parseOptions(vm);
		if (GeneralOptions::getOptions().getInputFiles().empty()) {
			ret = false;
			setOptionsErrorMsg("You must specify at least one input!");
		}
		return ret;
	}
};

typedef OptionsBaseTemplate< _TnfDistanceOptions > TnfDistanceOptions;
typedef std::vector<std::string> FastaVector;

class TNF {
public:
	typedef vector<float> Vector;
	typedef RankVector<float> RVector;

	static KM stdMap;
	static long stdSize;
	static DistanceFormula::Enum distanceFormula;

	static void init(int kmerSize, DistanceFormula::Enum df = DistanceFormula::EUCLIDEAN) {
		distanceFormula = df;
		KmerSizer::set(kmerSize);
		TEMP_KMER(kmer);
		TEMP_KMER(lc);
		kmer.clear();
		char fasta[kmerSize];
		strcpy(fasta, kmer.toFasta().c_str());

		do {
			kmer.set(fasta);
			kmer.buildLeastComplement(lc);

			bool exists = stdMap.exists(lc);
			LOG_DEBUG_OPTIONAL(2, exists, lc.toFasta() << "\t" << "\t" << kmer.toFasta());
			if ( !exists ) {

				stdMap.insert(lc, DataType());
			}
		} while(TwoBitSequence::permuteFasta(fasta, kmerSize));

		stdSize = stdMap.size();
		LOG_DEBUG(1, "Found " << stdSize << " unique " << kmerSize << "-mers");
	}
	static std::string toHeader() {
		assert(stdSize > 0);
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
	float length;
	Vector tnfValues;

public:
	TNF() : length(1.0) {
		assert(stdSize > 0);
		tnfValues.resize(stdSize);
	}
	TNF(const TNF &copy) {
		assert(stdSize > 0);
		*this = copy;
	}
	TNF(const KM &map, bool normalize = false) : length(1.0) {
		assert(stdSize > 0);
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
	TNF &operator*(float factor) {
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
		bool debug = Log::isDebug(1);
		ostream *debugPtr = NULL;
		if (debug) {
			debugPtr = Log::Debug("buildLength()").getOstreamPtr();
		}
		for(unsigned int i = 0 ; i < tnfValues.size(); i++) {
			length += tnfValues[i] * tnfValues[i];
			if (debug) *debugPtr << tnfValues[i] << "\t";
		}
		length = sqrt(length);
		if (length == 0.0)
			length = 1.0;
		if (debug) *debugPtr << length << endl;
	}
	void buildTnf(const KM & map) {
		tnfValues.empty();
		tnfValues.reserve(stdSize);
		for(KM::Iterator it = stdMap.begin(); it != stdMap.end(); it++) {
			KM::ElementType elem = map.getElementIfExists(it->key());
			float t = 0.0;
			if (elem.isValid()) {
				t = elem.value() / length;
			}
			tnfValues.push_back(t);
		}
	}
	float getDistance(const TNF &other) const {
		float dist = 0.0;
		if (distanceFormula == DistanceFormula::EUCLIDEAN) {
			for(unsigned int i = 0 ; i < tnfValues.size(); i++) {
				float d = tnfValues[i] / length - other.tnfValues[i] / other.length;
				dist += d*d;
			}
			dist = sqrt(dist);
		} else if (distanceFormula == DistanceFormula::SPEARMAN) {
			dist = RVector::getSpearmanDistance(RVector(tnfValues), RVector(other.tnfValues));
		}
		return dist;
	}
	float getDistance(const KM &query) const {
		TNF other(query, isNormalized());
		return getDistance(other);
	}
	static float getDistance(const KM &target, const KM &query, bool normalize = true) {
		TNF tgt(target, normalize), qry(query, normalize);
		return tgt.getDistance(qry);
	}

};
KM TNF::stdMap(2112); // initialize room for hexamers
long TNF::stdSize = 0;
DistanceFormula::Enum TNF::distanceFormula;

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
		LOG_DEBUG(1, "buildTnfs() for: " << read.getName());
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
		Read newRead(read.getName() + ":" + boost::lexical_cast<string>(i) + "-" + boost::lexical_cast<string>(i+fasta.length()), fasta, string(fasta.length(), Kmernator::PRINT_REF_QUAL), "");
		reads.append(newRead);
	}
	return buildTnfs(reads, normalize);
}

class Result {
public:
	float distance;
	string label;
	TNF tnf;
	Result() : distance(0.0) {}
	Result(float _dist, string _label) : distance(_dist), label(_label) {}
	Result(float _dist, string _label, TNF _tnf) :  distance(_dist), label(_label), tnf(_tnf) {}
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

		float dist = refTnf.getDistance(tnfs[readIdx]);
		if (dist <= 0.20) {
			string fullFasta = read.getFasta();
			long len = fullFasta.length();
			int divs = 5;
			TNFS intraTnfs = buildIntraTNFs(read,len/divs, len/divs, true);

			stringstream ss;
			ss.precision(3);
			ss << name;
			for(long i = 0 ; i < (long) intraTnfs.size() ; i++) {
				float distPart = refTnf.getDistance(intraTnfs[i]);
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

typedef vector<float> DVector;
typedef vector< DVector > DMatrix;
class MinVecOper {
public:
	bool operator()( const DVector::iterator &a, const DVector::iterator &b) const {
		return (*a < *b);
	}
} mvo;

int main(int argc, char *argv[]) {

	if (!TnfDistanceOptions::parseOpts(argc, argv)) exit(1);

	Cleanup::prepare();

	ReadSet refs;
	ReadSet reads;
	TNF::init(TnfDistanceBaseOptions::getOptions().getKmerSize(),
			TnfDistanceBaseOptions::getOptions().getDistanceFormula());

	OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
	reads.appendAllFiles(inputs);

	KS ksRef;
	ksRef.setSolidOnly();

	TNFS readTnfs = buildTnfs(reads, true);

	ostream *out = &cout;
	OfstreamMap om(Options::getOptions().getOutputFile(), "");
	if (!Options::getOptions().getOutputFile().empty()) {
		out = &om.getOfstream("");
	}
	OptionsBaseInterface::FileListType referenceInputs = TnfDistanceBaseOptions::getOptions().getReferenceFiles();
	if (!referenceInputs.empty()) {
		// compare distances from reference to each read in the input
		refs.appendAllFiles(referenceInputs);
		ksRef.buildKmerSpectrum(refs, true);
		TNF refTnf(ksRef.solid, true);

		Results results = calculateDistances(refTnf, reads, readTnfs);
		std::sort(results.begin(), results.end(), ResultCompare());
		for(Results::iterator it = results.begin(); it != results.end(); it++) {
			*out << it->distance << "\t" << it->label << endl;
		}
	} else {
		// output TNF matrix for each read in the input

		*out << "Label\t" << TNF::toHeader() << endl;
		for(ReadSet::ReadSetSizeType readIdx = 0; readIdx < reads.getSize(); readIdx++) {

			Read read = reads.getRead(readIdx);

			*out << read.getName() << "\t" << readTnfs[readIdx].toString() << endl;
		}

	}

	string interFile = TnfDistanceBaseOptions::getOptions().getInterFile();
	if (!interFile.empty()) {
		OfstreamMap om(interFile, "");
		ostream &os = om.getOfstream("");
		for(ReadSet::ReadSetSizeType readIdxi = 0; readIdxi < reads.getSize(); readIdxi++) {
			os << reads.getRead(readIdxi).getName();
			for(ReadSet::ReadSetSizeType readIdxj = 0; readIdxj < readIdxi; readIdxj++) {
				float dist = readTnfs[readIdxi].getDistance(readTnfs[readIdxj]);
				os << "\t" << dist;
			}
			os << endl;
		}
	}

	string intraFile = TnfDistanceBaseOptions::getOptions().getIntraFile();
	if (!intraFile.empty()) {
		OfstreamMap om(intraFile, "");
		ostream &os = om.getOfstream("");
		long window = TnfDistanceBaseOptions::getOptions().getIntraWindow();
		long step = window / 10;
		for(ReadSet::ReadSetSizeType readIdx = 0; readIdx < reads.getSize(); readIdx++) {
			Read read = reads.getRead(readIdx);
			TNFS intraTnfs = buildIntraTNFs(read, window, step, true);
			os << read.getName();
			for(unsigned long i = 0 ; i < intraTnfs.size(); i++) {
				for(unsigned long j = 0 ; j < i ; j++) {
					float dist = intraTnfs[i].getDistance(intraTnfs[j]);
					os << "\t" << dist;
				}
			}
			os << endl;
		}
	}

	string clusterFile = TnfDistanceBaseOptions::getOptions().getClusterFile();
	if (!clusterFile.empty()) {
		bool debug = Log::isDebug(1);
		long size = readTnfs.size();

		LOG_VERBOSE(1, "Building clusters");
		OfstreamMap om(clusterFile, "");
		ostream &os = om.getOfstream("");
		DMatrix distMatrix;
		distMatrix.resize( size );
		vector< vector< string > > clusterNames;
		clusterNames.resize( size );

		//TODO optimize this and the the while loop (use LT, and directed updates to minVec)
		vector< DVector::iterator > minVec;
		minVec.resize( size );
		float clusterThreshold = TnfDistanceBaseOptions::getOptions().getClusterThreshold();

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

		LOG_VERBOSE(1, "Finished dist matrix");
		if (debug) {
			ostream &debugOut = Log::Debug("DistMatrix");
			for(long i = 0 ; i < size ; i++) {
				debugOut << "distMatrix: " << i << "\t" << clusterNames[i][0];
				for(long j = 0 ; j <= i ; j++)
					debugOut << "\t" << distMatrix[i][j];
				debugOut << std::endl;
			}
		}


		TNFS readTnfs2 = TNFS(readTnfs);

		while( true ) {
			// set the thresholds for each remaining read

			vector< DVector::iterator >::iterator minElem = std::min_element(minVec.begin(), minVec.end(), mvo);
			if (**minElem > clusterThreshold)
				break;
			LOG_DEBUG(1, **minElem );

			unsigned long tgtj = minElem - minVec.begin();
			unsigned long tgti = *minElem - distMatrix[tgtj].begin();

			unsigned long lastIdx = minVec.size() - 1;
			if (debug) {
				ostream &debugOut = Log::Debug("MinValue");
				for(unsigned long i2 = 0 ; i2 < minVec.size(); i2++) {
					unsigned long j = i2 <= tgti ? i2 : tgti;
					unsigned long i = i2 <= tgti ? tgti : i2;

					debugOut << distMatrix[i][j] << "\t";
				}
				debugOut << endl;
			}
			clusterNames[tgti].insert( clusterNames[tgti].end(), clusterNames[tgtj].begin(), clusterNames[tgtj].end() );
			LOG_VERBOSE(1, "Cluster merged " << tgti << " with " << tgtj << " " << **minElem);
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
		clusterNames.clear();

	}

	return 0;
}

