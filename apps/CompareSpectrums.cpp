// $Header: /repository/PI_annex/robsandbox/KoMer/apps/CompareSpectrums.cpp,v 1.1 2009-11-27 23:28:07 cfurman Exp $
//

#include <iostream>
 

#include "ReadSet.h"
#include "Kmer.h"
#include "Utils.h"
 

//#include "Options.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

typedef KmerMap<TrackingData> KmerSolidMap;

class CS_Options
{
public:
  typedef std::vector<std::string> FileListType;
  
private:
  static inline CS_Options &getOptions() { static CS_Options singleton; return singleton; }
  
  po::options_description desc;
  po::positional_options_description p;
  po::variables_map vm;
  
  // cache of variables (for inline lookup and defaults)
  FileListType referenceFiles;
  FileListType inputFiles;
  unsigned int kmerSize;
 
  
public:

  static inline FileListType        &getFileSet1()   { return getOptions().referenceFiles; }
  static inline FileListType        &getFileSet2()       { return getOptions().inputFiles; }
  static inline unsigned int        &getKmerSize()         { return getOptions().kmerSize; }
 
    
  static bool parseOpts(int argc, char *argv[]) {
    try {
        po::options_description &desc = getOptions().desc;
        
        desc.add_options()
            ("help", "produce help message")
 
            ("file-set-1", po::value< FileListType >(), "1st file(s)")
 
            ("kmer-size", po::value< unsigned int >(), "kmer size")
            ("file-set-2", po::value< FileListType >(), "2nd file(s)")
         ;

        po::positional_options_description &p = getOptions().p;
        p.add("kmer-size", 1);
        p.add("file-set-1", 1);
        p.add("file-set-2", -1);
       
        po::variables_map &vm = getOptions().vm;        
        po::store(po::command_line_parser(argc, argv).
                 options(desc).positional(p).run(), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            std::cerr << desc << std::endl;
            return false; 
        }
 
        if (vm.count("file-set-1")) {
            std::cerr << "1st file set: ";
            FileListType &referenceFiles = getFileSet1() = vm["file-set-1"].as< FileListType >();
            for(FileListType::iterator it = referenceFiles.begin(); it != referenceFiles.end() ; it++ )
                std::cerr << *it << ", ";
            std::cerr << std::endl;
        } else {
            std::cerr << "file set 1 not specified" << std::endl;
            return false;
        }

        if (vm.count("file-set-2")) {
            std::cerr << "2nd file set: ";
            FileListType inputs = getFileSet2() = vm["file-set-2"].as< FileListType >();
            for(FileListType::iterator it = inputs.begin(); it != inputs.end() ; it++ )
                std::cerr << *it << ", ";
            std::cerr << std::endl;
        } else {
            std::cerr << desc << "file set 2 not specified!" << std::endl;
            return false;
        }
        
        if (vm.count("kmer-size")) {
            getKmerSize() = vm["kmer-size"].as< unsigned int >();
            std::cerr << "Kmer size is: " << getKmerSize() << std::endl;
        } else {
            std::cerr << desc << "There was no kmer size specified!" << std::endl;
            return false;
        }        
    }
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << std::endl;
        return false;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!" << std::endl;
        return false;
    }
    
    return true;
  }
};



#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using namespace std;


 
static unsigned long  countCommonKmers(KmerSolidMap &m1, KmerSolidMap &m2)
{
    unsigned long common = 0;
    
    for(KmerSolidMap::Iterator it(m1.begin()), itEnd(m1.end()); it != itEnd; it++)
    {
        if (m2.exists(it->key()))
          common++;
    }

    return common;
}

void buildKmerMap(KmerSolidMap &m, ReadSet &store)
{
   for (unsigned int i=0 ; i < store.getSize(); i++) {
    KmerWeights kmers = buildWeightedKmers(store.getRead(i));

    for (int j=0; j < kmers.size(); j++)
    {
        TEMP_KMER (least);
        bool keepDirection = kmers[j].buildLeastComplement( least);
        m[ least ].track( kmers.valueAt(j), keepDirection );
    }
  }
}
 
int main(int argc, char *argv[])
{
 
    if (!CS_Options::parseOpts(argc,argv))
      throw std::invalid_argument("Please fix the command line arguments");
      
    ReadSet readSet1, readSet2;

    KmerSizer::set(CS_Options::getKmerSize());
    CS_Options::FileListType fileList1  = CS_Options::getFileSet1();
    CS_Options::FileListType fileList2  = CS_Options::getFileSet2();
   
    cerr << "Reading 1st file set:";
    foreach( string file, fileList1 ) {
        readSet1.appendFastq( file  );
    }
    cerr << " loaded " << readSet1.getSize() << " Reads, " << readSet1.getBaseCount() << " Bases " << endl;

    
    cerr << "Reading 2nd file set:" ;
    foreach( string file, fileList2 ) {
       readSet2.appendFastq(file);
     }
    cerr << " loaded " << readSet2.getSize() << " Reads, " << readSet2.getBaseCount() << " Bases " << endl;
 
   KmerSolidMap m1,m2;

   cerr << "Building map 1\n";
   buildKmerMap(m1,readSet1);
    
   TrackingData::resetGlobalCounters();

   cerr << "Building map 2\n";
   buildKmerMap(m2,readSet2);

   cerr << "Counting common Kmers\n";
   unsigned long common = countCommonKmers(m1,m2);

   cout << endl;
   cout << "Set 1\tSet 2\tCommon\t% Set1\t% Set2\n";
   cout << m1.size() << '\t' << m2.size()<< '\t'<< common << '\t'  <<  setprecision(4) << (common*100.0)/m1.size() << '\t' <<  (common*100.0)/m2.size() << endl;
   
}







//
// $Log: CompareSpectrums.cpp,v $
// Revision 1.1  2009-11-27 23:28:07  cfurman
// CompareSpectrum application added
//
 
