// $Header: /repository/PI_annex/robsandbox/KoMer/src/Options.h,v 1.3 2010-01-14 00:47:44 regan Exp $
//

#ifndef _OPTIONS_H
#define _OPTIONS_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

class Options
{
public:
  typedef std::vector<std::string> FileListType;
  
private:
  static inline Options &getOptions() { static Options singleton; return singleton; }
  
  po::options_description desc;
  po::positional_options_description p;
  po::variables_map vm;
  
  // cache of variables (for inline lookup and defaults)
  FileListType  referenceFiles;
  FileListType  inputFiles;
  std::string   outputFile;
  unsigned int  kmerSize;
  double        solidQuantile;
  double        minKmerQuality;
  unsigned int  verbosity;
  double        firstOrderWeight;
  double        secondOrderWeight;
  unsigned int  minQuality;
  
public:

  static inline FileListType        &getReferenceFiles()   { return getOptions().referenceFiles; }
  static inline FileListType        &getInputFiles()       { return getOptions().inputFiles; }
  static inline std::string         &getOutputFile()       { return getOptions().outputFile; }
  static inline unsigned int        &getKmerSize()         { return getOptions().kmerSize; }
  static inline double              &getSolidQuantile()    { return getOptions().solidQuantile; }
  static inline double              &getMinKmerQuality()   { return getOptions().minKmerQuality; }
  static inline unsigned int        &getVerbosity()        { return getOptions().verbosity; }
  static inline double              &getFirstOrderWeight() { return getOptions().firstOrderWeight; }
  static inline double              &getSecondOrderWeight(){ return getOptions().secondOrderWeight; }
  static inline unsigned int        &getMinQuality()       { return getOptions().minQuality; }
    
  static bool parseOpts(int argc, char *argv[]) {
    try {
    	po::options_description &desc = getOptions().desc;
    	
        desc.add_options()
            ("help", "produce help message")
            ("verbose", po::value< unsigned int >()->default_value(0), "level of verbosity (0+)")
            ("reference-file", po::value< FileListType >(), "set reference file(s)")
            ("solid-quantile", po::value< double >()->default_value(0.05), "quantile threshold for solid kmers (0-1)")
            ("kmer-size", po::value< unsigned int >(), "kmer size")
            ("input-file", po::value< FileListType >(), "input file(s)")
            ("output-file", po::value< std::string >(), "output file or dir")
            ("min-kmer-quality", po::value< double >()->default_value(0.25), "minimum quality-adjusted kmer probability (0-1)")
            ("min-quality-score", po::value< unsigned int >()->default_value(10), "minimum quality score over entire kmer")
            ("first-order-weight", po::value< double >()->default_value(0.10), "first order permuted bases weight")
            ("second-order-weight", po::value< double >()->default_value(0.00), "second order permuted bases weight")
        ;

        po::positional_options_description &p = getOptions().p;
        p.add("kmer-size", 1);
        p.add("input-file", -1);
        
        po::variables_map &vm = getOptions().vm;        
        po::store(po::command_line_parser(argc, argv).
                 options(desc).positional(p).run(), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            std::cerr << desc << std::endl;
            return false;
        }
        getVerbosity() = vm["verbose"].as< unsigned int >();

        if (vm.count("reference-file")) {
        	std::cerr << "Reference files are: ";
            FileListType &referenceFiles = getReferenceFiles() = vm["reference-file"].as< FileListType >();
        	for(FileListType::iterator it = referenceFiles.begin(); it != referenceFiles.end() ; it++ )
                std::cerr << *it << ", ";
            std::cerr << std::endl;
        } else {
            std::cerr << "reference was not set." << std::endl;
        }
        
        if (vm.count("kmer-size")) {
        	getKmerSize() = vm["kmer-size"].as< unsigned int >();
        	std::cerr << "Kmer size is: " << getKmerSize() << std::endl;
        } else {
        	std::cerr << desc << "There was no kmer size specified!" << std::endl;
        	return false;
        }
        if (vm.count("input-file")) {	
        	std::cerr << "Input files are: ";
        	FileListType inputs = getInputFiles() = vm["input-file"].as< FileListType >();
        	for(FileListType::iterator it = inputs.begin(); it != inputs.end() ; it++ )
                std::cerr << *it << ", ";
            std::cerr << std::endl;
        } else {
        	std::cerr << desc << "There were no input files specified!" << std::endl;
        	return false;
        }
        
        
     
        // set solid-quantile
      	getSolidQuantile() = vm["solid-quantile"].as< double >();
        std::cerr << "solid-quantile is: " << getSolidQuantile() << std::endl;
        
        // set kmer quality
        getMinKmerQuality() = vm["min-kmer-quality"].as< double >();
        std::cerr << "min-kmer-quality is: " << getMinKmerQuality() << std::endl;
        
        // set minimum quality score
        getMinQuality() = vm["min-quality-score"].as< unsigned int >();
        std::cerr << "min-quality-score is: " << getMinQuality() << std::endl;
        
        // set permuted weights
        getFirstOrderWeight() = vm["first-order-weight"].as< double >();
        getSecondOrderWeight() = vm["second-order-weight"].as< double>();
        
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

#endif

//
// $Log: Options.h,v $
// Revision 1.3  2010-01-14 00:47:44  regan
// added two common options
//
// Revision 1.2  2009-11-09 19:37:17  regan
// enhanced some debugging / analysis output
//
// Revision 1.1  2009-11-06 04:10:21  regan
// refactor of cmd line option handling
// added methods to evaluate spectrums
//
//
