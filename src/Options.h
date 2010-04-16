// $Header: /repository/PI_annex/robsandbox/KoMer/src/Options.h,v 1.12 2010-04-16 22:44:18 regan Exp $
//

#ifndef _OPTIONS_H
#define _OPTIONS_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// put common, universal options in this class
// extend the class for specific options for each application
class Options {
public:
	typedef std::vector<std::string> FileListType;

protected:
	static inline Options &getOptions() {
		static Options singleton;
		return singleton;
	}

private:
	po::options_description desc;
	po::positional_options_description p;
	po::variables_map vm;

	// cache of variables (for inline lookup and defaults)
	FileListType referenceFiles;
	FileListType inputFiles;
	std::string outputFile;
	unsigned int kmerSize;
	double solidQuantile;
	double minKmerQuality;
	unsigned int verbosity;
	unsigned int debug;
	double firstOrderWeight;
	double secondOrderWeight;
	unsigned int minQuality;
	unsigned int minDepth;
	unsigned int minReadLength;
	unsigned int ignoreQual;
	unsigned int periodicSingletonPurge;
	unsigned int skipArtifactFilter;
	unsigned int mmapInput;

	Options() {
		setOptions();
	}

public:

	static inline FileListType &getReferenceFiles() {
		return getOptions().referenceFiles;
	}
	static inline FileListType &getInputFiles() {
		return getOptions().inputFiles;
	}
	static inline std::string &getOutputFile() {
		return getOptions().outputFile;
	}
	static inline unsigned int &getKmerSize() {
		return getOptions().kmerSize;
	}
	static inline double &getSolidQuantile() {
		return getOptions().solidQuantile;
	}
	static inline double &getMinKmerQuality() {
		return getOptions().minKmerQuality;
	}
	static inline unsigned int &getVerbosity() {
		return getOptions().verbosity;
	}
	static inline unsigned int &getDebug() {
		return getOptions().debug;
	}
	static inline double &getFirstOrderWeight() {
		return getOptions().firstOrderWeight;
	}
	static inline double &getSecondOrderWeight() {
		return getOptions().secondOrderWeight;
	}
	static inline unsigned int &getMinQuality() {
		return getOptions().minQuality;
	}
	static inline unsigned int &getMinDepth() {
		return getOptions().minDepth;
	}
	static inline unsigned int &getMinReadLength() {
		return getOptions().minReadLength;
	}
	static inline unsigned int &getIgnoreQual() {
		return getOptions().ignoreQual;
	}
	static inline unsigned int &getPeriodicSingletonPurge() {
		return getOptions().periodicSingletonPurge;
	}
	static inline unsigned int &getSkipArtifactFilter() {
		return getOptions().skipArtifactFilter;
	}
	static inline unsigned int &getMmapInput() {
		return getOptions().mmapInput;
	}
	const static unsigned int MAX_INT = (unsigned int) -1;

private:
	void setOptions() {

		desc.add_options()("help", "produce help message")

		("verbose", po::value<unsigned int>()->default_value(0),
				"level of verbosity (0+)")

		("debug", po::value<unsigned int>()->default_value(0),
				"level of debug verbosity (0+)")

		("reference-file", po::value<FileListType>(), "set reference file(s)")

		("kmer-size", po::value<unsigned int>(), "kmer size.  A size of 0 will skip k-mer calculations")

		("input-file", po::value<FileListType>(), "input file(s)")

		("output-file", po::value<std::string>(), "output file or dir")

		("min-read-length",
				po::value<int>()->default_value(0),
				"minimum (trimmed) read length of selected reads.  0: (default) no minimum, -1: full read length")

		("min-kmer-quality", po::value<double>()->default_value(0.10),
				"minimum quality-adjusted kmer probability (0-1)")

		("min-quality-score", po::value<unsigned int>()->default_value(10),
				"minimum quality score over entire kmer")

		("min-depth", po::value<unsigned int>()->default_value(10),
				"minimum depth for a solid kmer")

		("solid-quantile", po::value<double>()->default_value(0.05),
				"quantile threshold for solid kmers (0-1)")

		("first-order-weight", po::value<double>()->default_value(0.10),
				"first order permuted bases weight")

		("second-order-weight", po::value<double>()->default_value(0.00),
				"second order permuted bases weight")

		("ignore-quality", po::value<unsigned int>()->default_value(0),
				"ignore the quality score, to save memory or if they are untrusted")

		("periodic-singleton-purge", po::value<unsigned int>()->default_value(0),
				"Purge singleton memory structure every # of reads")

		("skip-artifact-filter", po::value<unsigned int>()->default_value(0),
				"Skip homo-polymer, primer-dimer and duplicated fragment pair filtering")

		("mmap-input", po::value<unsigned int>()->default_value(1),
				"If set to 0, prevents input files from being mmaped, instead import reads into memory (somewhat faster if memory is abundant)");

	}

public:

	static po::options_description &getDesc() {
		return getOptions().desc;
	}
	static po::positional_options_description &getPosDesc() {
		return getOptions().p;
	}
	static po::variables_map &getVarMap() {
		return getOptions().vm;
	}

	static bool parseOpts(int argc, char *argv[]) {
		try {

			po::options_description & desc = getDesc();
			po::positional_options_description & p = getPosDesc();

			po::variables_map & vm = getVarMap();
			po::store(
					po::command_line_parser(argc, argv). options(desc).positional(
							p).run(), vm);
			po::notify(vm);

			if (vm.count("help")) {
				std::cerr << desc << std::endl;
				return false;
			}
			getVerbosity() = vm["verbose"].as<unsigned int> ();
			getDebug() = vm["debug"].as<unsigned int> ();

			if (vm.count("reference-file")) {
				std::cerr << "Reference files are: ";
				FileListType & referenceFiles = getReferenceFiles()
						= vm["reference-file"].as<FileListType> ();
				for (FileListType::iterator it = referenceFiles.begin(); it
						!= referenceFiles.end(); it++)
					std::cerr << *it << ", ";
				std::cerr << std::endl;
			} else {
				std::cerr << "reference was not set." << std::endl;
			}

			if (vm.count("kmer-size")) {
				getKmerSize() = vm["kmer-size"].as<unsigned int> ();
				std::cerr << "Kmer size is: " << getKmerSize() << std::endl;
			} else {
				std::cerr << desc << "There was no kmer size specified!"
						<< std::endl;
				return false;
			}
			if (vm.count("input-file")) {
				std::cerr << "Input files are: ";
				FileListType inputs = getInputFiles() = vm["input-file"].as<
						FileListType> ();
				for (FileListType::iterator it = inputs.begin(); it
						!= inputs.end(); it++)
					std::cerr << *it << ", ";
				std::cerr << std::endl;
			} else {
				std::cerr << desc << "There were no input files specified!"
						<< std::endl;
				return false;
			}
			if (vm.count("output-file")) {
				getOutputFile() = vm["output-file"].as<std::string> ();
				std::cerr << "Output file (or dir) is: " << getOutputFile()
						<< std::endl;
			}

			// set kmer quality
			getMinKmerQuality() = vm["min-kmer-quality"].as<double> ();

			std::cerr << "min-kmer-quality is: " << getMinKmerQuality()
					<< std::endl;

			// set minimum quality score
			getMinQuality() = vm["min-quality-score"].as<unsigned int> ();
			std::cerr << "min-quality-score is: " << getMinQuality()
					<< std::endl;

			// set minimum depth
			getMinDepth() = vm["min-depth"].as<unsigned int> ();
			std::cerr << "min-depth is: " << getMinDepth() << std::endl;

			getMinReadLength() = vm["min-read-length"].as<int> ();
			std::cerr << "min-read-length is: " << getMinReadLength()
					<< std::endl;

			// set solid-quantile
			getSolidQuantile() = vm["solid-quantile"].as<double> ();
			std::cerr << "solid-quantile is: " << getSolidQuantile()
					<< std::endl;

			// set permuted weights
			getFirstOrderWeight() = vm["first-order-weight"].as<double> ();
			getSecondOrderWeight() = vm["second-order-weight"].as<double> ();

			// set the ignore quality value
			getIgnoreQual() = vm["ignore-quality"].as<unsigned int> ();

			// set periodic singleton purge value
			getPeriodicSingletonPurge() = vm["periodic-singleton-purge"].as<unsigned int> ();

			// set skipArtifactFilterin
			getSkipArtifactFilter() = vm["skip-artifact-filter"].as<unsigned int> ();

			// set mmapInput
			getMmapInput() = vm["mmap-input"].as<unsigned int>();

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

#endif

//
// $Log: Options.h,v $
// Revision 1.12  2010-04-16 22:44:18  regan
// merged HEAD with changes for mmap and intrusive pointer
//
// Revision 1.11.2.2  2010-04-15 17:59:52  regan
// made mmap optional
//
// Revision 1.11.2.1  2010-04-15 17:47:19  regan
// added option to skip artifact filter
//
// Revision 1.11  2010-03-15 07:42:39  regan
// added kmer of 0 to skip kmer calculations
//
// Revision 1.10  2010-03-10 13:18:19  regan
// added quality ignore and singleton purge to save memory
//
// Revision 1.9  2010-03-03 17:08:52  regan
// fixed options to support -1 in min read length
//
// Revision 1.8  2010-03-02 15:03:33  regan
// added debug option and reformatted
//
// Revision 1.7  2010-02-26 13:01:16  regan
// reformatted
//
// Revision 1.6  2010-02-22 14:40:05  regan
// added optoin
//
// Revision 1.5  2010-01-16 01:06:28  regan
// refactored
//
// Revision 1.4  2010-01-14 01:17:48  regan
// working out options
//
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
