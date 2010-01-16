// $Header: /repository/PI_annex/robsandbox/KoMer/apps/FilterReads.cpp,v 1.4 2010-01-16 01:12:45 regan Exp $
//

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "config.h"
#include "Options.h"
#include "ReadSet.h"
#include "FilterKnownOddities.h"
#include "KmerSpectrum.h"
#include "ReadSelector.h"

using namespace std;

typedef TrackingDataMinimal8 DataType;
typedef KmerSpectrum<DataType,DataType> KS;
typedef ReadSelector<DataType> RS;

class FilterReadsOptions : Options
{
public:
	static bool parseOpts(int argc, char *argv[])
	{
		getPosDesc().add("kmer-size", 1);
        getPosDesc().add("input-file", -1);
		return Options::parseOpts(argc, argv);
	}
};

int main(int argc, char *argv[]) {
   if (!FilterReadsOptions::parseOpts(argc,argv))
      throw std::invalid_argument("Please fix the command line arguments");

    cerr << MemoryUtils::getMemoryUsage() << endl;
	
	ReadSet reads;
    FilterKnownOddities filter;
    KmerSizer::set(Options::getKmerSize());
    
    Options::FileListType inputs     = Options::getInputFiles();
    cerr << "Reading Input Files" << endl;
    reads.appendAllFiles(inputs);
    cerr << "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount() << " Bases " << endl;
    cerr << MemoryUtils::getMemoryUsage() << endl;    
    
    cerr << "Identifying Pairs: ";
    long numPairs =  reads.identifyPairs();
    cerr << numPairs << endl;
    cerr << MemoryUtils::getMemoryUsage() << endl;
    
    cerr << "Applying filter to Input Files" << endl;
    unsigned long filtered = filter.applyFilter(reads);
    cerr << "filter affected " << filtered << " Reads " << endl;
    cerr << MemoryUtils::getMemoryUsage() << endl;

    long numBuckets = KS::estimateWeakKmerBucketSize( reads, 64 );
    cerr << "targetting " << numBuckets << " buckets for reads " << endl;

    KS spectrum(numBuckets);
    cerr << MemoryUtils::getMemoryUsage() << endl;

    TrackingData::minimumWeight = 0.25;

    spectrum.buildKmerSpectrum( reads );
    cerr << MemoryUtils::getMemoryUsage() << endl;

    cerr << "Picking reads: " << endl;
    RS selector(reads, spectrum.weak, Options::getMinDepth());
    long oldPicked = 0;
    long picked = 0;
    for(unsigned int depth = 4096 ; depth >= 2; depth /= 2) {
      if (reads.hasPairs())
        picked = selector.pickAllPassingPairs(depth);
      else
        picked = selector.pickAllPassingReads(depth);
      cerr << "At or above coverage: " << depth << " Picked " << picked << " / " << reads.getSize() << " reads" << endl;
    
      std::string outputFile = Options::getOutputFile();
      if (picked > 0 && !outputFile.empty()) {
        outputFile += "-" + boost::lexical_cast<std::string>( depth );
    	cerr << "Writing reads to output file: " << outputFile << endl;
    	ofstream out;
    	out.open(outputFile.c_str());
    	selector.writePicks(out, oldPicked);
    	out.close();
      }
      oldPicked += picked;
    }
}
