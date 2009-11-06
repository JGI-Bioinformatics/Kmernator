// $Header: /repository/PI_annex/robsandbox/KoMer/src/Options.h,v 1.1 2009-11-06 04:10:21 regan Exp $
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
private:
  static Options &getOptions() { static Options singleton; return singleton; }
  static po::variables_map &getVM() { return getOptions().vm; }
  
  po::options_description desc;
  po::positional_options_description p;
  po::variables_map vm;
  
public:

  typedef std::vector<std::string> FileListType;

  static FileListType        getReferenceFiles() { if (getVM().count("reference-file"))
  	                                                  return getVM()["reference-file"].as< FileListType >();
  	                                               FileListType empty(0);
  	                                               return empty;
  	                                             }
  static FileListType        getInputFiles()     { return getVM()["input-file"]    .as< FileListType >(); }
  static unsigned int        getKmerSize()       { return getVM()["kmer-size"]     .as< unsigned int >(); }
  static double              getSolidQuantile()  { return getVM()["solid-quantile"].as< double       >(); }
  
  static bool parseOpts(int argc, char *argv[]) {
    try {
    	po::options_description &desc = getOptions().desc;
    	
        desc.add_options()
            ("help", "produce help message")
            ("reference-file", po::value< FileListType >(), "set reference file(s)")
            ("solid-quantile", po::value< double >()->default_value(0.05), "quantile threshold for solid kmers (0-1)")
            ("kmer-size", po::value< unsigned int >(), "kmer size")
            ("input-file", po::value< FileListType >(), "input file(s)")
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

        if (vm.count("reference-file")) {
        	std::cerr << "Reference files are: ";
            FileListType referenceFiles = getReferenceFiles();
        	for(FileListType::iterator it = referenceFiles.begin(); it != referenceFiles.end() ; it++ )
                std::cerr << *it << ", ";
            std::cerr << std::endl;
        } else {
            std::cerr << "reference was not set." << std::endl;
        }
        
        if (vm.count("kmer-size")) {
        	std::cerr << "Kmer size is: " << getKmerSize() << std::endl;
        } else {
        	std::cerr << desc << "There was no kmer size specified!" << std::endl;
        	return false;
        }
        if (vm.count("input-file")) {	
        	std::cerr << "Input files are: ";
        	FileListType inputs = getInputFiles();
        	for(FileListType::iterator it = inputs.begin(); it != inputs.end() ; it++ )
                std::cerr << *it << ", ";
            std::cerr << std::endl;
        } else {
        	std::cerr << desc << "There were no input files specified!" << std::endl;
        	return false;
        }
        std::cerr << "solid-quantile is: " << getSolidQuantile() << std::endl;
        
        
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
// Revision 1.1  2009-11-06 04:10:21  regan
// refactor of cmd line option handling
// added methods to evaluate spectrums
//
//
