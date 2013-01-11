//
// DistributedNucleatingAssembler.cpp
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

#include "config.h"
#include "ReadSet.h"
#include "Options.h"
#include "Kmer.h"
#include "KmerSpectrum.h"
#include "FilterKnownOddities.h"
#include "DuplicateFragmentFilter.h"
#include "ContigExtender.h"
#include "DistributedFunctions.h"
#include "MemoryUtils.h"
#include "Utils.h"
#include "Log.h"
#include "Vmatch.h"
#include "Cap3.h"
#include "KmerMatch.h"
#include "MatcherInterface.h"

using namespace std;
typedef TrackingDataMinimal4f DataType;
typedef KmerSpectrum<DataType, DataType> KS;

typedef KmerSpectrum< TrackingDataWithAllReads, TrackingDataWithAllReads > KS2;


class _DistributedNucleatingAssemblerOptions: public OptionsBaseInterface {
public:
	_DistributedNucleatingAssemblerOptions(): maxIterations(1000), maxContigLength(3000) {}
	virtual ~_DistributedNucleatingAssemblerOptions() {}

	int &getMaxIterations() {
		return maxIterations; // getVarMap()["max-iterations"].as<int> ();
	}
	int &getMaxContigLength() {
		return maxContigLength; // getVarMap()["max-contig-length"].as<int>();
	}
	void _resetDefaults() {
		Cap3Options::_resetDefaults();
		ContigExtenderBaseOptions::_resetDefaults();
		MatcherInterfaceOptions::_resetDefaults();
		VmatchOptions::_resetDefaults();
		KmerMatchOptions::_resetDefaults();
		KmerBaseOptions::_resetDefaults();
		KmerSpectrumOptions::_resetDefaults();
		MPIOptions::_resetDefaults();
		FilterKnownOdditiesOptions::_resetDefaults();

		GeneralOptions::_resetDefaults();
		FilterKnownOdditiesOptions::getOptions().getSkipArtifactFilter() = 0;
		KmerSpectrumOptions::getOptions().getMinDepth() = 1;

		// override the default output format!
		GeneralOptions::getOptions().getFormatOutput() = 3;
		GeneralOptions::getOptions().getMmapInput() = false;
		GeneralOptions::getOptions().getVerbose() = 1;
		GeneralOptions::getOptions().getMaxThreads() = 1;
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {

		po::options_description opts("DistributedNucleatingAssembly <options> [input-file ...]\n\tNote: --input-file can be specified either as a positional argument or within the options\n\nDistributedNucleatingAssembler Options");

		p.add("input-file", -1);

		opts.add_options()
						("max-iterations", po::value<int>()->default_value(maxIterations), "the maximum number of rounds to extend the set of contigs")

						("max-contig-length", po::value<int>()->default_value(maxContigLength), "the maximum size of a contig to continue extending");

		desc.add(opts);

		MatcherInterfaceOptions::_setOptions(desc,p);
		KmerMatchOptions::_setOptions(desc,p);
		KmerBaseOptions::_setOptions(desc,p);
		KmerSpectrumOptions::_setOptions(desc,p);
		VmatchOptions::_setOptions(desc,p);
		ContigExtenderBaseOptions::_setOptions(desc,p);
		Cap3Options::_setOptions(desc,p);
		MPIOptions::_setOptions(desc,p);
		FilterKnownOdditiesOptions::_setOptions(desc,p);
		GeneralOptions::_setOptions(desc,p);

	};
	bool _parseOptions(po::variables_map &vm) {

		bool ret = true;
		ret &= GeneralOptions::_parseOptions(vm);
		ret &= MatcherInterfaceOptions::_parseOptions(vm);
		ret &= KmerMatchOptions::_parseOptions(vm);
		ret &= KmerBaseOptions::_parseOptions(vm);
		ret &= KmerSpectrumOptions::_parseOptions(vm);
		ret &= VmatchOptions::_parseOptions(vm);
		ret &= ContigExtenderBaseOptions::_parseOptions(vm);
		ret &= Cap3Options::_parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);
		ret &= FilterKnownOdditiesOptions::_parseOptions(vm);

		setOpt("max-iterations", maxIterations);
		setOpt("max-contig-length", maxContigLength);

		if (Options::getOptions().getOutputFile().empty()) {
			LOG_ERROR(1, "You must specify an --output");
			ret = false;
		}

		return ret;
	}
protected:
	int maxIterations;
	int maxContigLength;

};
typedef OptionsBaseTemplate<_DistributedNucleatingAssemblerOptions>
DistributedNucleatingAssemblerOptions;


std::string extendContigsWithCap3(ReadSet & contigs,
		ReadSet::ReadSetVector &contigReadSet, ReadSet & changedContigs,
		ReadSet & finalContigs, ReadSet::ReadSetSizeType minimumCoverage) {
	std::stringstream extendLog;
	std::string cap3Path = Cap3Options::getOptions().getCap3Path();

	int poolsWithoutMinimumCoverage = 0;

	#pragma omp parallel for
	for (long i = 0; i < (long) contigs.getSize(); i++) {
		const Read &oldRead = contigs.getRead(i);
		Read newRead = oldRead;
		SequenceLengthType oldLen = oldRead.getLength(), newLen = 0;

		ReadSet::ReadSetSizeType poolSize = contigReadSet[i].getSize();

		double extTime = MPI_Wtime();
		if (poolSize > minimumCoverage) {
			LOG_VERBOSE_OPTIONAL(2, true, "Extending " << oldRead.getName() << " with " << poolSize << " pool of reads");
			newRead = Cap3::extendContig(oldRead, contigReadSet[i]);
			newLen = newRead.getLength();
		} else {
			poolsWithoutMinimumCoverage++;
		}
		extTime = MPI_Wtime() - extTime;
		long deltaLen = (long)newLen - (long)oldLen;
		if (deltaLen > 0) {
			extendLog << std::endl << "Cap3 Extended " << oldRead.getName() << " "
					<< deltaLen << " bases to " << newRead.getLength() << ": "
					<< newRead.getName() << " with " << poolSize
					<< " reads in the pool, in " << extTime << " sec";
			//#pragma omp critical
			changedContigs.append(newRead);
		} else {
			extendLog << std::endl << "Did not extend " << oldRead.getName() << " with " << poolSize << " reads in the pool, in " << extTime << " sec";
			//#pragma omp critical
			finalContigs.append(oldRead);
		}
	}

	LOG_VERBOSE_OPTIONAL(2, true, "Extended " << contigs.getSize() - poolsWithoutMinimumCoverage << " contigs out of " << contigs.getSize());

	return extendLog.str();
}
std::string extendContigsWithContigExtender(ReadSet & contigs,
		ReadSet::ReadSetVector &contigReadSet, ReadSet & changedContigs,
		ReadSet & finalContigs, SequenceLengthType minKmerSize,
		double minimumCoverage, SequenceLengthType maxKmerSize,
		SequenceLengthType maxExtend, SequenceLengthType kmerStep) {

	std::stringstream extendLog;
	//#pragma omp parallel for
	for (ReadSet::ReadSetSizeType i = 0; i < contigs.getSize(); i++) {
		const Read &oldRead = contigs.getRead(i);
		Read newRead;
		SequenceLengthType oldLen = oldRead.getLength(), newLen = 0;
		ReadSet::ReadSetSizeType poolSize = contigReadSet[i].getSize();
		SequenceLengthType myKmerSize = minKmerSize;
		if (poolSize > minimumCoverage) {
			LOG_VERBOSE_OPTIONAL(2, true, "kmer-Extending " << oldRead.getName() << " with " << poolSize << " pool of reads");
			ReadSet myContig;
			myContig.append(oldRead);
			ReadSet newContig;

			while (newLen <= oldLen && myKmerSize <= maxKmerSize) {
				newContig = ContigExtender<KS>::extendContigs(myContig,
						contigReadSet[i], maxExtend, myKmerSize, myKmerSize);
				newLen = newContig.getRead(0).getLength();
				myKmerSize += kmerStep;
			}
			newRead = newContig.getRead(0);
		} else {
			newRead = oldRead;
		}
		long deltaLen = (long) newLen - (long) oldLen;
		if (deltaLen > 0) {
			extendLog << std::endl << "Kmer Extended " << oldRead.getName() << " "
					<< deltaLen << " bases to " << newRead.getLength() << ": "
					<< newRead.getName() << " with " << poolSize
					<< " reads in the pool K " << (myKmerSize - kmerStep);
			//#pragma omp critical
			changedContigs.append(newRead);
		} else {
			extendLog << std::endl << "Did not extend " << oldRead.getName() << " with " << poolSize << " reads in the pool";
			//#pragma omp critical
			finalContigs.append(oldRead);
		}
	}
	return extendLog.str();
}

void finishLongContigs(long maxContigLength, ReadSet &changedContigs, ReadSet &finalContigs) {
	ReadSet keepContigs;
	for(long i = 0; i < (long) changedContigs.getSize(); i++) {
		const Read &read = changedContigs.getRead(i);
		if ((long) read.getLength() >= maxContigLength) {
			LOG_VERBOSE_OPTIONAL(1, true, read.getName() << " (" << read.getLength() << ") has exceeded maxContiglength, terminating extension");
			finalContigs.append(read);
		} else
			keepContigs.append(read);
	}
	changedContigs.swap(keepContigs);
}

int main(int argc, char *argv[]) {

	ScopedMPIComm< DistributedNucleatingAssemblerOptions > world(argc, argv);

	Cleanup::prepare();

	try {

		double timing1, timing2;

		timing1 = MPI_Wtime();

		OptionsBaseInterface::FileListType &inputFiles =
				Options::getOptions().getInputFiles();
		std::string contigFile =
				ContigExtenderBaseOptions::getOptions().getContigFile();
		std::string finalContigFile;
		double minimumCoverage =
				ContigExtenderBaseOptions::getOptions().getMinimumCoverage();
		long maxIterations =
				DistributedNucleatingAssemblerOptions::getOptions().getMaxIterations();

		ReadSet reads;
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files" );
		reads.appendAllFiles(inputFiles, world.rank(), world.size());
		reads.identifyPairs();
		setGlobalReadSetOffsets(world, reads);

		timing2 = MPI_Wtime();

		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "loaded " << reads.getGlobalSize() << " Reads, (local:" << reads.getSize() << " pair:" << reads.getPairSize() << ") in " << (timing2-timing1) << " seconds" );

		if (FilterKnownOdditiesOptions::getOptions().getSkipArtifactFilter() == 0) {

			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Preparing artifact filter: ");

			FilterKnownOddities filter;
			LOG_VERBOSE_OPTIONAL(2, world.rank() == 0, "Applying sequence artifact filter to Input Files");

			unsigned long filtered = filter.applyFilter(reads);

			LOG_VERBOSE(2, "local filter affected (trimmed/removed) " << filtered << " Reads ");

			unsigned long allFiltered;
			mpi::reduce(world, filtered, allFiltered, std::plus<unsigned long>(), 0);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "distributed filter (trimmed/removed) " << allFiltered << " Reads ");

		}

		boost::shared_ptr< MatcherInterface > matcher;
		if (KmerBaseOptions::getOptions().getKmerSize() == 0) {
			matcher.reset( new Vmatch(world, UniqueName::generateHashName(inputFiles), reads) );
		} else {
			matcher.reset( new KmerMatch(world, reads) );
		}

		SequenceLengthType minKmerSize, maxKmerSize, kmerStep, maxExtend;
		ContigExtender<KS>::getMinMaxKmerSize(reads, minKmerSize, maxKmerSize,
				kmerStep);
		maxKmerSize = boost::mpi::all_reduce(world, maxKmerSize, mpi::minimum<
				SequenceLengthType>());
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Kmer size ranges: " << minKmerSize << "\t" << maxKmerSize << "\t" << kmerStep);
		maxExtend = maxKmerSize;

		timing1 = timing2;
		timing2 = MPI_Wtime();
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Prepared Matcher indexes in " << (timing2-timing1) << " seconds");

		ReadSet finalContigs;
		ReadSet contigs;
		contigs.appendFastaFile(contigFile, world.rank(), world.size());

		short iteration = 0;
		while (++iteration <= maxIterations) {

			matcher->resetTimes("Start Iteration", MPI_Wtime());

			setGlobalReadSetOffsets(world, contigs);

			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Iteration: " << iteration << ". Contig File: " << contigFile << ". contains " << contigs.getGlobalSize() << " Reads");
			if (contigs.getGlobalSize() == 0) {
				LOG_VERBOSE_OPTIONAL(1, true, "There are no contigs to extend in " << contigFile);
				break;
			}

			MatcherInterface::MatchReadResults contigReadSet = matcher->match(contigs, contigFile);
			assert(contigs.getSize() == contigReadSet.size());

			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Iteration: " << iteration << ". Matches made");

			ReadSet changedContigs;
			std::string extendLog;

			if (!Cap3Options::getOptions().getCap3Path().empty()) {
				extendLog = extendContigsWithCap3(contigs, contigReadSet, changedContigs, finalContigs, minimumCoverage);
			} else {
				extendLog = extendContigsWithContigExtender(contigs, contigReadSet,
						changedContigs, finalContigs,
						minKmerSize, minimumCoverage, maxKmerSize, maxExtend, kmerStep);
			}

			matcher->recordTime("extendContigs", MPI_Wtime());
			LOG_DEBUG(1, (extendLog));

			finishLongContigs(DistributedNucleatingAssemblerOptions::getOptions().getMaxContigLength(), changedContigs, finalContigs);

			LOG_DEBUG(1, "Changed contigs: " << changedContigs.getSize() << " finalContigs: " << finalContigs.getSize());
			setGlobalReadSetOffsets(world, changedContigs);
			setGlobalReadSetOffsets(world, finalContigs);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Changed contigs: " << changedContigs.getGlobalSize() << " finalContigs: " << finalContigs.getGlobalSize());

			std::string oldFinalContigFile = finalContigFile;
			std::string oldContigFile = contigFile;
			{
				// write out the state of the contig files (so far) so we do not loose them
				DistributedOfstreamMap om(world,
						Options::getOptions().getOutputFile(), "");
				om.setBuildInMemory();
				if (finalContigs.getGlobalSize() > 0) {
					std::string fileKey = "final-" + boost::lexical_cast<
							std::string>(iteration);
					finalContigs.writeAll(om.getOfstream(fileKey),
							FormatOutput::Fasta());
					finalContigFile = om.getRealFilePath(fileKey);
				}
				if (changedContigs.getGlobalSize() > 0) {
					std::string filekey = "-inputcontigs-" + boost::lexical_cast<
							std::string>(iteration) + ".fasta";
					changedContigs.writeAll(om.getOfstream(filekey),
							FormatOutput::Fasta());
					contigFile = om.getRealFilePath(filekey);
				}
				contigs = changedContigs;
			}

			if (world.rank() == 0) {
				// preserve the final contigs in case of crash
				unlink(Options::getOptions().getOutputFile().c_str());
				link(finalContigFile.c_str(), Options::getOptions().getOutputFile().c_str());
			}

			matcher->recordTime("writeFinalTime", MPI_Wtime());

			if (!Log::isDebug(1) && world.rank() == 0) {
				// remove most recent contig files (if not debugging)
				if (!oldFinalContigFile.empty()) {
					LOG_VERBOSE_OPTIONAL(1, true, "Removing " << oldFinalContigFile);
					unlink(oldFinalContigFile.c_str());
				}

				if (ContigExtenderBaseOptions::getOptions().getContigFile().compare(
						oldContigFile) != 0) {
					LOG_VERBOSE_OPTIONAL(1, true, "Removing " << oldContigFile);
					unlink(oldContigFile.c_str());
				}
			}

			if (changedContigs.getGlobalSize() == 0) {
				LOG_VERBOSE_OPTIONAL(1, world.rank() == 1, "No more contigs to extend " << changedContigs.getSize());
				break;
			}

			matcher->recordTime("finishIteration", MPI_Wtime());
			LOG_DEBUG(1, matcher->getTimes("") + " " + MemoryUtils::getMemoryUsage());

		}

		matcher.reset(); // release the matcher interface

		if (world.rank() == 0 && !Log::isDebug(1)) {
			if (ContigExtenderBaseOptions::getOptions().getContigFile().compare(
					contigFile) != 0) {
				LOG_DEBUG_OPTIONAL(1, true, "Removing " << contigFile);
				unlink(contigFile.c_str());
			}
		}

		// write final contigs (and any unfinished contigs still remaining)
		finalContigs.append(contigs);
		std::string tmpFinalFile = DistributedOfstreamMap::writeGlobalReadSet(world, finalContigs, Options::getOptions().getOutputFile(), ".tmp", FormatOutput::Fasta());
		if (world.rank() == 0 && !finalContigFile.empty()) {
			LOG_DEBUG_OPTIONAL(1, true, "Removing " << finalContigFile);
			unlink(finalContigFile.c_str());
		}
		finalContigFile = tmpFinalFile;
		if (world.rank() == 0) {
			unlink(Options::getOptions().getOutputFile().c_str());
			rename(finalContigFile.c_str(), Options::getOptions().getOutputFile().c_str());
		}
		finalContigFile = Options::getOptions().getOutputFile();
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Final contigs are in: " << finalContigFile);

		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");

	} catch (std::exception &e) {
		LOG_ERROR(1, "DistributedNucleatingAssembler threw an exception! Aborting..." << e.what());
		world.abort(1);
	} catch (...) {
		LOG_ERROR(1, "DistributedNucleatingAssembler threw an error!" );
		world.abort(1);
	}

	return 0;
}

