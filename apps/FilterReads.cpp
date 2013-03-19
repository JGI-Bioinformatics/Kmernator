//
// Kmernator/apps/FilterReads.cpp
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

#include "FilterReads.h"

typedef TrackingDataWithDirection DataType;
typedef KmerMap< DataType > MapType;
typedef KmerSpectrum<MapType, MapType> KS;
typedef ReadSelector< MapType > RS;
class _FilterReadsOptions : public OptionsBaseInterface {
public:
	void _resetDefaults() {
		FilterReadsBaseOptions::_resetDefaults();

		KmerSpectrumOptions::getOptions().getSaveKmerMmap() = 0;
	}
	void _setOptions(po::options_description &desc, po::positional_options_description &p) {
		FilterReadsBaseOptions::_setOptions(desc, p);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;

		ret &= FilterReadsBaseOptions::_parseOptions(vm);
		return ret;
	}
};
typedef OptionsBaseTemplate< _FilterReadsOptions > FilterReadsOptions;

int main(int argc, char *argv[]) {

	FilterReadsOptions::parseOpts(argc, argv);

	Cleanup::prepare();

	MemoryUtils::getMemoryUsage();
	std::string outputFilename = Options::getOptions().getOutputFile();

	ReadSet reads;

	try {
		OptionsBaseInterface::FileListType &inputs = Options::getOptions().getInputFiles();
		LOG_VERBOSE(1, "Reading Input Files");
		reads.appendAllFiles(inputs);
		LOG_VERBOSE(1, "loaded " << reads.getSize() << " Reads, " << reads.getBaseCount()
				<< " Bases ");
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		LOG_VERBOSE(1, "Identifying Pairs: ");
		long numPairs = reads.identifyPairs();
		LOG_VERBOSE(1, "Pairs + single = " << numPairs);
		LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		if (FilterKnownOdditiesOptions::getOptions().getSkipArtifactFilter() == 0) {

			LOG_VERBOSE(1, "Preparing artifact filter: ");
			FilterKnownOddities filter;
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			LOG_VERBOSE(2, "Applying sequence artifact filter to Input Files");
			unsigned long filtered = filter.applyFilter(reads);
			LOG_VERBOSE(1, "filter affected (trimmed/removed) " << filtered << " Reads ");;
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

		}
		if (DuplicateFragmentFilterOptions::getOptions().getDeDupMode() > 0 && DuplicateFragmentFilterOptions::getOptions().getDeDupEditDistance() >= 0) {
			LOG_VERBOSE(2, "Applying DuplicateFragmentPair Filter to Input Files");
			unsigned long duplicateFragments = DuplicateFragmentFilter::filterDuplicateFragments(reads);
			LOG_VERBOSE(1, "filter removed duplicate fragment pair reads: " << duplicateFragments);
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());
		}

		KS spectrum(0);

		Kmernator::MmapFileVector spectrumMmaps;
		if (KmerBaseOptions::getOptions().getKmerSize() > 0 && !KmerSpectrumOptions::getOptions().getLoadKmerMmap().empty()) {
			spectrum.restoreMmap(KmerSpectrumOptions::getOptions().getLoadKmerMmap());
		} else if (KmerBaseOptions::getOptions().getKmerSize() > 0) {

			long rawKmers = KS::estimateRawKmers(reads);
			LOG_DEBUG(1, "targeting " << rawKmers << " raw kmers for reads ");

			spectrum = KS(rawKmers);
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			spectrumMmaps = spectrum.buildKmerSpectrumInParts(reads, KmerSpectrumOptions::getOptions().getBuildPartitions(), outputFilename.empty() ? "" : outputFilename + "-mmap");
			spectrum.optimize();
			spectrum.trackSpectrum(true);
			std::string sizeHistoryFile = FilterReadsBaseOptions::getOptions().getSizeHistoryFile();
			if (!sizeHistoryFile.empty()) {
				LOG_VERBOSE(1, "Writing size history file to: " << sizeHistoryFile);
				OfstreamMap ofm(sizeHistoryFile, "");
				ofm.getOfstream("") << spectrum.getSizeTracker().toString();
			} else {
				LOG_VERBOSE(1, "Kmer Size History:" << std::endl << spectrum.getSizeTracker().toString());
			}

			if (Log::isVerbose(1))
				spectrum.printHistograms(Log::Verbose("Kmer Histogram"));

			if (!FilterReadsBaseOptions::getOptions().getHistogramFile().empty()) {
				ofstream of(FilterReadsBaseOptions::getOptions().getHistogramFile().c_str());
				spectrum.printHistograms(of);
			}

			if (KmerSpectrumOptions::getOptions().getVariantSigmas() > 0.0) {
				spectrum.purgeVariants();
				if (Log::isVerbose(1)) {
					spectrum.printHistograms(Log::Verbose("Variant-Removed Kmer Histogram"));
				}
			}
		}

		if (KmerBaseOptions::getOptions().getKmerSize() > 0) {
			LOG_DEBUG(1, MemoryUtils::getMemoryUsage());

			if (KmerSpectrumOptions::getOptions().getGCHeatMap() && ! outputFilename.empty()) {
				LOG_VERBOSE(1, "Creating GC Heat Map ");
				LOG_DEBUG(1,  MemoryUtils::getMemoryUsage());
				OfstreamMap ofmap(outputFilename + "-GC", ".txt");
				spectrum.printGC(ofmap.getOfstream(""));
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
				LOG_VERBOSE(1, "Trimming reads with minDepth: " << minDepth);
			} else {
				LOG_VERBOSE(1, "Trimming reads that pass Artifact Filter with length: " << ReadSelectorOptions::getOptions().getMinReadLength());
			}

			RS selector(reads, spectrum.weak);
			selector.scoreAndTrimReads(minDepth);

			selectReads(minDepth, reads, selector, outputFilename);
		}
		LOG_DEBUG(1, "Clearing spectrum");
		spectrum.reset();

	} catch (std::exception &e) {
		LOG_ERROR(1, "FilterReads threw an exception!\n\t" << e.what());
		return 1;
	} catch (...) {
		LOG_ERROR(1, "FilterReads threw an error!");
		return 1;
	}

	LOG_VERBOSE(1, "Finished");

	return 0;
}


