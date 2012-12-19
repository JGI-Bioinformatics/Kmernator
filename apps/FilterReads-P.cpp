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


#include "FilterReads.h"
#include "DistributedFunctions.h"

typedef TrackingDataWithDirection DataType;
typedef DistributedKmerSpectrum<DataType, DataType> KS;
typedef DistributedReadSelector<DataType> RS;

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
		FilterReadsBaseOptions::_setOptions(desc, p);
		MPIOptions::_setOptions(desc,p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		ret &= FilterReadsBaseOptions::_parseOptions(vm);
		ret &= MPIOptions::_parseOptions(vm);

		return ret;
	}
};
typedef OptionsBaseTemplate< _MPIFilterReadsOptions > MPIFilterReadsOptions;

int main(int argc, char *argv[]) {

	ScopedMPIComm< MPIFilterReadsOptions > world (argc, argv);

	if (ReadSelectorOptions::getOptions().getMaxKmerDepth() > 0 && ReadSelectorOptions::getOptions().getNormalizationMethod() == "OPTIMAL" && world.size() > 1) {
		if (Logger::isMaster())
			LOG_WARN(1, "Setting --normalization-method to RANDOM, as Distributed version does not support 'OPTIMAL' option");
		ReadSelectorOptions::getOptions().getNormalizationMethod() = "RANDOM";
	}

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	ReadSet reads;

	try {
		OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Reading Input Files");

		reads.appendAllFiles(inputs, world.rank(), world.size());

		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Identifying Pairs: ");

		unsigned long counts[3], totalCounts[3];
		unsigned long &readCount = counts[0] = reads.getSize();
		unsigned long &numPairs  = counts[1] = reads.identifyPairs();
		unsigned long &baseCount = counts[2] = reads.getBaseCount();
		LOG_VERBOSE(2, "loaded " << readCount << " Reads, " << baseCount << " Bases ");
		LOG_VERBOSE(2, "Pairs + single = " << numPairs);

		mpi::all_reduce(world, (unsigned long*) counts, 3, (unsigned long*) totalCounts, std::plus<unsigned long>());
		LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Loaded " << totalCounts[0] << " distributed reads, " << totalCounts[1] << " distributed pairs, " << totalCounts[2] << " distributed bases");

		setGlobalReadSetOffsets(world, reads);

		if (FilterKnownOdditiesOptions::getOptions().getSkipArtifactFilter() == 0) {

			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Preparing artifact filter: ");

			FilterKnownOddities filter;
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			LOG_VERBOSE_OPTIONAL(2, world.rank() == 0, "Applying sequence artifact filter to Input Files");

			unsigned long filtered = filter.applyFilter(reads);

			LOG_VERBOSE(2, "local filter affected (trimmed/removed) " << filtered << " Reads ");
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			unsigned long allFiltered;
			mpi::reduce(world, filtered, allFiltered, std::plus<unsigned long>(), 0);
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "distributed filter (trimmed/removed) " << allFiltered << " Reads ");

		}

		if ( DuplicateFragmentFilterOptions::getOptions().getDeDupMode() > 0 && DuplicateFragmentFilterOptions::getOptions().getDeDupEditDistance() >= 0) {
			if (world.size() == 1) {
				LOG_VERBOSE(2, "Applying DuplicateFragmentPair Filter to Input Files");
				unsigned long duplicateFragments = DuplicateFragmentFilter::filterDuplicateFragments(reads);

				LOG_VERBOSE(2, "filter removed duplicate fragment pair reads: " << duplicateFragments);
				LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

				unsigned long allDuplicateFragments;
				mpi::reduce(world, duplicateFragments, allDuplicateFragments, std::plus<unsigned long>(), 0);
				LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "distributed removed duplicate fragment pair reads: " << allDuplicateFragments);
			} else {
				if (world.rank() == 0)
					LOG_WARN(1, "Distributed DuplicateFragmentPair Filter is not supported (yet)." << std::endl
							<< "If you want this feature please run the non-MPI FilterReads");
			}

		}

		long numBuckets = 0;
		if (KmerBaseOptions::getOptions().getKmerSize() > 0) {

			numBuckets = KS::estimateWeakKmerBucketSize(reads);

			numBuckets = all_reduce(world, numBuckets, mpi::maximum<int>());
			LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "targeting " << numBuckets << " buckets for reads");
		}
		KS spectrum(world, numBuckets);
		Kmernator::MmapFileVector spectrumMmaps;
		if (KmerBaseOptions::getOptions().getKmerSize() > 0 && !KmerSpectrumOptions::getOptions().getLoadKmerMmap().empty()) {
			spectrum.restoreMmap(KmerSpectrumOptions::getOptions().getLoadKmerMmap());
		} else if (KmerBaseOptions::getOptions().getKmerSize() > 0) {
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			spectrum.buildKmerSpectrum(reads);
			spectrum.optimize();
			spectrum.trackSpectrum(true);
			KS::SizeTracker reducedSizeTracker = spectrum.reduceSizeTracker(world);

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
				ofstream of(FilterReadsBaseOptions::getOptions().getHistogramFile().c_str());
				spectrum.printHistograms(of);
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

			if (!outputFilename.empty() && KmerSpectrumOptions::getOptions().getSaveKmerMmap() > 0) {
				spectrumMmaps = spectrum.writeKmerMaps(outputFilename + "-mmap");
				LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
			}

			if (KmerSpectrumOptions::getOptions().getMinDepth() > 1) {
				LOG_DEBUG(1, "Clearing singletons from memory");
				spectrum.singleton.clear();
				LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
			} else {
				spectrum.optimize(true);
			}
		}


		unsigned int minDepth = KmerSpectrumOptions::getOptions().getMinDepth();

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
			LOG_VERBOSE(1, "Writing Files");

			selectReads(minDepth, reads, selector, outputFilename);
		}
		spectrum.reset();
		LOG_DEBUG(1, "Finished, waiting for rest of collective");

	} catch (...) {
		LOG_ERROR(1, "caught an error!" << StackTrace::getStackTrace());
	}
	LOG_DEBUG(2, "Clearing spectrum");

	world.barrier();
	LOG_VERBOSE_OPTIONAL(1, world.rank() == 0, "Finished");

	return 0;
}
