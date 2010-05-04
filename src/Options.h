// $Header: /repository/PI_annex/robsandbox/KoMer/src/Options.h,v 1.14 2010-05-04 23:49:21 regan Exp $
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
	FileListType inputFilePrefixes;
	std::string outputFile;
	std::string tmpDir;
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
	unsigned int maskSimpleRepeats;
	unsigned int phiXOutput;
	unsigned int filterOutput;
	unsigned int deDupMode;
	unsigned int deDupSingle;
	unsigned int deDupEditDistance;
	unsigned int mmapInput;
	unsigned int buildPartitions;

	Options() : tmpDir("/tmp"), kmerSize(31), minKmerQuality(0.10), verbosity(1), debug(0),
	minQuality(10), minDepth(10), minReadLength(20), ignoreQual(0),
	periodicSingletonPurge(0), skipArtifactFilter(0), maskSimpleRepeats(1), phiXOutput(0), filterOutput(0),
	deDupMode(1), deDupSingle(0), deDupEditDistance(0), mmapInput(1), buildPartitions(0) {
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
	static inline std::string &getTmpDir() {
		return getOptions().tmpDir;
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
	static inline unsigned int &getMaskSimpleRepeats() {
		return getOptions().maskSimpleRepeats;
	}
	static inline unsigned int &getPhiXOutput(){
		return getOptions().phiXOutput;
	}
	static inline unsigned int &getFilterOutput(){
		return getOptions().filterOutput;
	}
	static inline unsigned int &getDeDupMode() {
		return getOptions().deDupMode;
	}
	static inline unsigned int &getDeDupSingle() {
		return getOptions().deDupSingle;
	}
	static inline unsigned int &getDeDupEditDistance() {
		return getOptions().deDupEditDistance;
	}
	static inline unsigned int &getMmapInput() {
		return getOptions().mmapInput;
	}
	static inline unsigned int &getBuildPartitions() {
		return getOptions().buildPartitions;
	}
	const static unsigned int MAX_INT = (unsigned int) -1;

	static std::string &getInputFileSubstring(unsigned int fileIdx) {
		if (getOptions().inputFilePrefixes.empty()) {
			// populate all the file indexes
			for(FileListType::iterator it = getOptions().inputFiles.begin(); it != getOptions().inputFiles.end(); it++) {
				size_t start = it->find_last_of('/');
				if (start == std::string::npos)
					start = 0;
				else
					start++;
				size_t end = it->find_last_of('.');
				if (end == std::string::npos)
					end = it->length() - 1;

				std::string fileprefix = it->substr( start, end-start );
				if (getDebug() > 0) {
					std::cerr << "InputFilePrefix: " << fileprefix << std::endl;
				}
				getOptions().inputFilePrefixes.push_back( fileprefix );
			}
		}
		return getOptions().inputFilePrefixes[fileIdx];
	}
private:
	void setOptions() {

		desc.add_options()("help", "produce help message")

		("verbose", po::value<unsigned int>()->default_value(verbosity),
				"level of verbosity (0+)")

		("debug", po::value<unsigned int>()->default_value(debug),
				"level of debug verbosity (0+)")

		("reference-file", po::value<FileListType>(), "set reference file(s)")

		("kmer-size", po::value<unsigned int>()->default_value(kmerSize), "kmer size.  A size of 0 will skip k-mer calculations")

		("input-file", po::value<FileListType>(), "input file(s)")

		("output-file", po::value<std::string>(), "output file pattern")

        ("phix-output", po::value<unsigned int>()->default_value(phiXOutput),
		        "if set, artifact filter also screens for PhiX174, and any matching reads will be output into a separate file (requires --output-file set)")

		("filter-output", po::value<unsigned int>()->default_value(filterOutput),
				"if set, artifact filter reads will be output into a separate file. If not set, then affected reads will be trimmed and then output normally.  (requires --output-file set)")

		("temp-dir", po::value<std::string>()->default_value(tmpDir), "temporary directory to deposit mmap file")

		("min-read-length",
				po::value<int>()->default_value(minReadLength),
				"minimum (trimmed) read length of selected reads.  0: (default) no minimum, -1: full read length")

		("min-kmer-quality", po::value<double>()->default_value(minKmerQuality),
				"minimum quality-adjusted kmer probability (0-1)")

		("min-quality-score", po::value<unsigned int>()->default_value(minQuality),
				"minimum quality score over entire kmer")

		("min-depth", po::value<unsigned int>()->default_value(minDepth),
				"minimum depth for a solid kmer")

		("solid-quantile", po::value<double>()->default_value(0.05),
				"quantile threshold for solid kmers (0-1)")

		("first-order-weight", po::value<double>()->default_value(0.10),
				"first order permuted bases weight")

		("second-order-weight", po::value<double>()->default_value(0.00),
				"second order permuted bases weight")

		("ignore-quality", po::value<unsigned int>()->default_value(ignoreQual),
				"ignore the quality score, to save memory or if they are untrusted")

		("periodic-singleton-purge", po::value<unsigned int>()->default_value(periodicSingletonPurge),
				"Purge singleton memory structure every # of reads")

		("skip-artifact-filter", po::value<unsigned int>()->default_value(skipArtifactFilter),
				"Skip homo-polymer, primer-dimer and duplicated fragment pair filtering")

		("mask-simple-repeats", po::value<unsigned int>()->default_value(maskSimpleRepeats),
				"if filtering artifacts, also mask simple repeats")

		("dedup-mode", po::value<unsigned int>()->default_value(deDupMode),
				"if 0, no fragment de-duplication will occur.  if 1, single orientation (AB and BA are separated) will collapse to consensus. if 2, both orientations (AB and BA are the same) will collapse")

        ("dedup-single", po::value<unsigned int>()->default_value(deDupSingle),
		    "if 0, no single read de-duplication will occur.  if 1, then single read deduplication will occur")

		("dedup-edit-distance", po::value<unsigned int>()->default_value(deDupEditDistance),
				"if -1, no fragment de-duplication will occur, if 0, only exact match, ...")

		("mmap-input", po::value<unsigned int>()->default_value(mmapInput),
				"If set to 0, prevents input files from being mmaped, instead import reads into memory (somewhat faster if memory is abundant)")

		("build-partitions", po::value<unsigned int>()->default_value(buildPartitions),
				"If set, kmer spectrum will be computed in stages and then combined in mmaped files on disk.  Must be a power of 2");

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
				std::cerr << "Output file pattern: " << getOutputFile()
						<< std::endl;
			}

			if (vm.count("temp-dir")) {
				getTmpDir() = vm["temp-dir"].as<std::string> ();
				std::cerr << "Tmp dir is: " << getTmpDir() << std::endl;
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

			// set skipArtifactFiltering
			getSkipArtifactFilter() = vm["skip-artifact-filter"].as<unsigned int> ();

			// set simple repeat masking
			getMaskSimpleRepeats() = vm["mask-simple-repeats"].as<unsigned int> ();
			// set phix masking
			getPhiXOutput() = vm["phix-output"].as<unsigned int> ();
			// set simple repeat masking
			getFilterOutput() = vm["filter-output"].as<unsigned int> ();

			// set dedup mode
			getDeDupMode() = vm["dedup-mode"].as<unsigned int>();
			// set dedup single
			getDeDupSingle() = vm["dedup-single"].as<unsigned int>();
			if (getDeDupSingle() > 0) {
			    std::cerr << "Unsupported options dedup-single (unimplemented)" << std::endl;
			    return false;
			}
			// set dedup edit distance
			getDeDupEditDistance() = vm["dedup-edit-distance"].as<unsigned int>();
			if (getDeDupEditDistance() > 1) {
				std::cerr <<"Unsupported option dedup-edit-distance > 1" << std::endl;
			    return false;
			}

			// set mmapInput
			getMmapInput() = vm["mmap-input"].as<unsigned int>();

			// set buildPartitions
			getBuildPartitions() = vm["build-partitions"].as<unsigned int>();

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
// Revision 1.14  2010-05-04 23:49:21  regan
// merged changes for FixConsensusOutput-201000504
//
// Revision 1.13.6.2  2010-05-04 23:41:30  regan
// checkpoint
//
// Revision 1.13.6.1  2010-05-04 21:51:47  regan
// checkpoint
//
// Revision 1.13  2010-05-01 21:57:53  regan
// merged head with serial threaded build partitioning
//
// Revision 1.12.4.14  2010-05-01 05:57:40  regan
// made edit distance for dedup fragment pairs optional
//
// Revision 1.12.4.13  2010-04-30 16:33:21  regan
// better option descriptions
//
// Revision 1.12.4.12  2010-04-29 20:33:29  regan
// bugfix
//
// Revision 1.12.4.11  2010-04-29 16:54:07  regan
// bugfixes
//
// Revision 1.12.4.10  2010-04-29 06:58:55  regan
// modified default optoins
//
// Revision 1.12.4.9  2010-04-29 04:26:44  regan
// phix changes
//
// Revision 1.12.4.8  2010-04-28 22:28:10  regan
// refactored writing routines
//
// Revision 1.12.4.7  2010-04-28 16:57:00  regan
// bugfix in output filenames
//
// Revision 1.12.4.6  2010-04-27 23:17:50  regan
// fixed naming of output files
//
// Revision 1.12.4.5  2010-04-27 22:53:20  regan
// reworked defaults
//
// Revision 1.12.4.4  2010-04-27 18:25:20  regan
// bugfix in temp directory usage
//
// Revision 1.12.4.3  2010-04-27 05:48:46  regan
// added bulild in parts to main code and options
//
// Revision 1.12.4.2  2010-04-26 22:56:26  regan
// bugfix
//
// Revision 1.12.4.1  2010-04-26 22:53:08  regan
// added more global options
//
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
