//
// Kmernator/apps/FilterReads-P.cpp
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
//

#include <mpi.h>
#include <boost/mpi.hpp>

#include "FilterReads.h"
#include "DistributedFunctions.h"

typedef TrackingDataWithDirection DataType;
typedef TrackingDataSingleton SDataType;
typedef KmerMapByKmerArrayPair< DataType > MapType;
typedef KmerMapByKmerArrayPair< SDataType > SMapType;
typedef DistributedKmerSpectrum<MapType, MapType, SMapType> KS;
typedef DistributedReadSelector< MapType > RS;

class _MPIFilterReadsOptions : public OptionsBaseInterface {
public:
	void _resetDefaults() {
		FilterReadsBaseOptions::_resetDefaults();
		MPIOptions::_resetDefaults();

		// assign defaults
		GeneralOptions::getOptions().getMmapInput() = false;
		GeneralOptions::getOptions().getVerbose() = 2;
		KmerSpectrumOptions::getOptions().getSaveKmerMmap() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		po::options_description opts("Usage: mpirun FilterReads-P <options> [[--kmer-size] 51] [[--input-file] ...]\n\tNote: --kmer-size and --input-file can either be specified as positional arguments at the end or within <options>");

		p.add("kmer-size", 1);
		p.add("input-file", -1);
		desc.add(opts);

		FilterReadsBaseOptions::_setOptions(desc, p);
		MPIOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		ret &= FilterReadsBaseOptions::_parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);

		if (KmerSpectrumOptions::getOptions().getSaveKmerMmap() || !KmerSpectrumOptions::getOptions().getLoadKmerMmap().empty()) {
			if (Logger::isMaster())
				LOG_WARN(1, "FilterReads-P can no longer load and save the kmer-map to disk.");
			KmerSpectrumOptions::getOptions().getSaveKmerMmap() = false;
			KmerSpectrumOptions::getOptions().getLoadKmerMmap().clear();
		}
		return ret;
	}
};
typedef OptionsBaseTemplate< _MPIFilterReadsOptions > MPIFilterReadsOptions;

template<typename KS, typename RS>
void BuildSpectrumAndFilter(ScopedMPIComm< MPIFilterReadsOptions > &world, ReadSet &reads, std::string &outputFilename, boost::shared_ptr< KS > subtractingSpectrum = NULL)
{
	long rawKmers = 0;
	if(KmerBaseOptions::getOptions().getKmerSize() > 0){
		rawKmers = KS::estimateRawKmers(world, reads);
	}
	unsigned int minDepth = KmerSpectrumOptions::getOptions().getMinDepth();
	KS spectrum(world, rawKmers);
	Kmernator::MmapFileVector spectrumMmaps;
	if (KmerBaseOptions::getOptions().getKmerSize() > 0 && !KmerSpectrumOptions::getOptions().getLoadKmerMmap().empty()) {
		spectrum.restoreMmap(KmerSpectrumOptions::getOptions().getLoadKmerMmap());
	} else if (KmerBaseOptions::getOptions().getKmerSize() > 0) {
		LOG_DEBUG_GATHER(1, MemoryUtils::getMemoryUsage());

		spectrum.subtractReference(subtractingSpectrum);
		subtractingSpectrum.reset(); // spectrum now has its own copy

		spectrum.buildKmerSpectrum(reads);
		spectrum.optimize();
		spectrum.trackSpectrum(true);
		typename KS::SizeTracker reducedSizeTracker = spectrum.reduceSizeTracker(world);

		std::string sizeHistoryFile = FilterReadsBaseOptions::getOptions().getSizeHistoryFile();
		if (!sizeHistoryFile.empty()) {
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Writing size history file to: " << sizeHistoryFile);
			if (world.rank() == 0) {
				OfstreamMap ofm(sizeHistoryFile, "");
				ofm.getOfstream("") << reducedSizeTracker.toString();
			}
		} else {
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Kmer Size History:" << std::endl << reducedSizeTracker.toString())
		}

		if (Log::isVerbose(1)) {
			std::string hist = spectrum.getHistogram(false);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Kmer Histogram\n" << hist);
		}

		if (!FilterReadsBaseOptions::getOptions().getHistogramFile().empty()) {
			std::string hist = spectrum.getHistogram(false);
			if (world.rank() == 0) {
				ofstream of(FilterReadsBaseOptions::getOptions().getHistogramFile().c_str());
				of << hist;
			}
		}

	}
	if (KmerBaseOptions::getOptions().getKmerSize() > 0) {

		if (KmerSpectrumOptions::getOptions().getVariantSigmas() > 0.0) {
			long purgedVariants = spectrum.purgeVariants();
			long totalPurgedVariants = mpi::all_reduce(world, purgedVariants, std::plus<long>());
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Distributed Purged " << totalPurgedVariants << " kmer variants");

			std::string hist = spectrum.getHistogram(false);

			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Collective Variant Purged Kmer Histogram\n" << hist);
			world.barrier();

		}
		/*
		 * TODO Reimplement the ability to store and restore GoogleSparse maps...
			if (!outputFilename.empty() && KmerSpectrumOptions::getOptions().getSaveKmerMmap() > 0) {
				spectrumMmaps = spectrum.writeKmerMaps(outputFilename + "-mmap");
				LOG_DEBUG_GATHER(1, MemoryUtils::getMemoryUsage());
			}
		 */

		if (minDepth > 1) {
			LOG_DEBUG_GATHER(1, "Clearing singletons from memory");
			spectrum.purgeMinDepth(minDepth, true);
			LOG_DEBUG_GATHER(1, MemoryUtils::getMemoryUsage());
		} else {
			spectrum.optimize(true);
		}
	}
	if (!outputFilename.empty()) {
		if (KmerBaseOptions::getOptions().getKmerSize() > 0) {
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Trimming reads with minDepth: " << minDepth);
		} else {
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Trimming reads that pass Artifact Filter with length: " << ReadSelectorOptions::getOptions().getMinReadLength());
		}
		RS selector(world, reads, spectrum.weak);
		selector.scoreAndTrimReads(minDepth);

		// TODO implement a more efficient algorithm to output data in order

		// rank 0 will overwrite, all others will append
		if (world.rank() != 0)
			OfstreamMap::getDefaultAppend() = true;

		// let only one rank at a time write to the files
		LOG_VERBOSE_GATHER(1, "Writing Files");

		selectReads(minDepth, reads, selector, outputFilename);
	}
	spectrum.reset();
	LOG_DEBUG_GATHER(1, "Finished, waiting for rest of collective");
}

void importAndProcessReadset(ScopedMPIComm<MPIFilterReadsOptions>& world,
		ReadSet& reads, OptionsBaseInterface::FileListType& inputs) {
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files");
	reads.appendAllFiles(inputs, world.rank(), world.size());
	LOG_DEBUG_GATHER(1, MemoryUtils::getMemoryUsage());
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Identifying Pairs: ");
	setGlobalReadSetConstants(world, reads);
	LOG_DEBUG_GATHER(1, MemoryUtils::getMemoryUsage());
	if (FilterKnownOdditiesOptions::getOptions().getSkipArtifactFilter() == 0) {

		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0,
				"Preparing artifact filter: ");

		FilterKnownOddities filter;
		LOG_DEBUG_GATHER(1, MemoryUtils::getMemoryUsage());

		LOG_VERBOSE_OPTIONAL(2, world.rank() == 0,
				"Applying sequence artifact filter to Input Files");

		unsigned long filtered = filter.applyFilter(reads);

		LOG_VERBOSE_GATHER(2,
				"local filter affected (trimmed/removed) " << filtered
				<< " Reads ");
		LOG_DEBUG_GATHER(1, MemoryUtils::getMemoryUsage());

		unsigned long allFiltered;
		mpi::reduce(world, filtered, allFiltered, std::plus<unsigned long>(),
				0);
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0,
				"distributed filter (trimmed/removed) " << allFiltered << " Reads.");

	}
	if (DuplicateFragmentFilterOptions::getOptions().getDeDupMode() > 0
			&& DuplicateFragmentFilterOptions::getOptions().getDeDupEditDistance()
			>= 0) {
		if (world.size() == 1) {
			LOG_VERBOSE_GATHER(2,
					"Applying DuplicateFragmentPair Filter to Input Files");
			unsigned long duplicateFragments =
					DuplicateFragmentFilter::filterDuplicateFragments(reads);

			LOG_VERBOSE(2,
					"filter removed duplicate fragment pair reads: " << duplicateFragments);
			LOG_DEBUG_GATHER(1, MemoryUtils::getMemoryUsage());

			unsigned long allDuplicateFragments;
			mpi::reduce(world, duplicateFragments, allDuplicateFragments,
					std::plus<unsigned long>(), 0);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0,
					"distributed removed duplicate fragment pair reads: " << allDuplicateFragments);
		} else {
			if (world.rank() == 0)
				LOG_WARN(1,
						"Distributed DuplicateFragmentPair Filter is not supported (yet)." << std::endl << "If you want this feature please run the non-MPI FilterReads");
		}

	}
}

int main(int argc, char *argv[]) {

	ScopedMPIComm< MPIFilterReadsOptions > world (argc, argv);

	Cleanup::prepare();

	if (ReadSelectorOptions::getOptions().getMaxKmerDepth() > 0 && ReadSelectorOptions::getOptions().getNormalizationMethod() == "OPTIMAL" && world.size() > 1) {
		if (Logger::isMaster())
			LOG_WARN(1, "Setting --normalization-method to RANDOM, as Distributed version does not support 'OPTIMAL' option");
		ReadSelectorOptions::getOptions().getNormalizationMethod() = "RANDOM";
	}

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	ReadSet reads;

	try {
		boost::shared_ptr< KS > subtractingSpectrum;
		OptionsBaseInterface::FileListType &subtractFiles = FilterReadsBaseOptions::getOptions().getSubtractFiles();
		OptionsBaseInterface::FileListType &referenceFiles = FilterReadsBaseOptions::getOptions().getReferenceFiles();
		OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
		importAndProcessReadset(world, reads, inputs);
		long rawKmers = KS::estimateRawKmers(world, reads);
		if (!referenceFiles.empty()) {
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Subtracting reference-file set: " << OptionsBaseInterface::toString(referenceFiles));
			if (subtractingSpectrum.get() == NULL)
				subtractingSpectrum.reset( new KS(world, rawKmers) );
			ReadSetStream rss(referenceFiles);
			subtractingSpectrum->buildKmerSpectrum(rss, true, 0);
		}
		if (!subtractFiles.empty()) {

			ReadSet subtractReads;
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Subtracting abundant kmers in subtract-file set: " << OptionsBaseInterface::toString(subtractFiles));
			importAndProcessReadset(world, subtractReads, subtractFiles);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Building subtraction abundant kmers");
			if (subtractingSpectrum.get() == NULL)
				subtractingSpectrum.reset( new KS(world, rawKmers) );
			subtractingSpectrum->buildKmerSpectrum(subtractReads);
			subtractReads.clear();
			unsigned int minDepth = KmerSpectrumOptions::getOptions().getMinDepth();
			if (minDepth > 1)
				subtractingSpectrum->purgeMinDepth(minDepth, true);
			subtractingSpectrum->optimize();
		}

		BuildSpectrumAndFilter<KS, RS>(world, reads, outputFilename, subtractingSpectrum);

	} catch (std::exception &e) {
		LOG_ERROR(1, "FilterReads-P caught an exception!\n\t" << e.what());
		world.abort(1);
	} catch (...) {
		LOG_ERROR(1, "FilterReads-P caught an error!\n\t");
		world.abort(1);
	}
	LOG_DEBUG(2, "Clearing spectrum");

	world.barrier();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");

	return 0;
}
