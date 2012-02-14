/*
 * Cap3.h
 *
 *  Created on: Sep 7, 2011
 *      Author: regan
 */

#ifndef CAP3_H_
#define CAP3_H_

#include "Options.h"
#include "Utils.h"
#include "ReadSet.h"
#include "Kmer.h"
#include "KmerSpectrum.h"

class _Cap3Options : public OptionsBaseInterface {
public:
	_Cap3Options() : cap3Path() {}
	virtual ~_Cap3Options() {}
	std::string &getCap3Path() {
		return cap3Path;
	}
	void _resetDefaults() {
		Options::getOptions().getMmapInput() = 0;
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {
		po::options_description opts("Cap3 Options");

		opts.add_options()
		("cap3-path", po::value<std::string>()->default_value(cap3Path),
				"if set, cap3 will be used to extend contigs")

		;

		desc.add(opts);
	}
	bool _parseOptions(po::variables_map &vm) {
		bool ret = true;
		setOpt<std::string>("cap3-path", cap3Path);
		return ret;
	}
protected:
	std::string cap3Path;
};
typedef OptionsBaseTemplate< _Cap3Options > Cap3Options;

class Cap3 {
public:
	typedef ReadSet::ReadSetSizeType ReadSetSizeType;
	static const int repeatContig = 1;
	static const double minOverlapFraction = 0.9;
	static std::string getNewName(std::string oldName, int deltaLen) {
		std::string newName;
		size_t pos = oldName.find_last_of("+");
		if (pos != oldName.npos) {
			newName = oldName.substr(0, pos );
			deltaLen += atoi(oldName.substr(pos).c_str());
		} else {
			newName = oldName;
		}

		return newName + "+" + boost::lexical_cast<std::string>(deltaLen);
	}
	static Read selectBestContig(const ReadSet &candidateContigs, const Read &targetContig) {
		Read bestRead = targetContig;
		SequenceLengthType oldSize = KmerSizer::getSequenceLength();
		KmerSizer::set(16);
		KmerSpectrum<> tgtSpectrum(16);
		tgtSpectrum.setSolidOnly();
		tgtSpectrum.buildKmerSpectrum(targetContig);
		long tgtKmers = tgtSpectrum.solid.size();
		LOG_DEBUG(4, "Cap3::selectBestContig(): tgtContig: " << targetContig.getLength() << ", " << tgtKmers);

		for(ReadSetSizeType i = 0; i < candidateContigs.getSize(); i++) {
			const Read &tgtRead = candidateContigs.getRead(i);
			KmerWeights readKmers(tgtRead.getTwoBitSequence(), tgtRead.getLength(), true);
			tgtSpectrum.getCounts(readKmers, false);
			long readOverlap = readKmers.sumAll();
			if (readOverlap >= minOverlapFraction * tgtKmers && tgtRead.getLength() >= targetContig.getLength() && bestRead.getLength() < tgtRead.getLength()) {
				bestRead = tgtRead;
				bestRead.setName(getNewName(targetContig.getName(), tgtRead.getLength() - targetContig.getLength()));
			}
			LOG_DEBUG(4, "Cap3::selectBestContig(): tgtRead: " << tgtRead.getLength() << ", " << readOverlap << ": " << bestRead.getLength());
		}

		if (bestRead.getLength() > targetContig.getLength()) {
			LOG_DEBUG_OPTIONAL(2, true, "Cap3 new contig: " << bestRead.getName() << " from " << targetContig.getLength() << " to " << bestRead.getLength());
		} else {
			LOG_DEBUG_OPTIONAL(2, true, "Cap3 failed to extend: " << bestRead.getName());
		}
		KmerSizer::set(oldSize);
		return bestRead;
	}
	static Read extendContig(const Read &oldContig, const ReadSet &_inputReads) {
		ReadSet inputReads = _inputReads;

		for(int i = 0; i < repeatContig; i++)
			inputReads.append(oldContig);

		FormatOutput format = FormatOutput::FastaUnmasked();
		std::string outputDir = Options::getOptions().getTmpDir() + UniqueName::generateUniqueName("/.cap3-assembly");
		int status = mkdir(outputDir.c_str(), 0777);

		if (status != 0) {
			LOG_WARN(1, "Could not mkdir: " << outputDir << " bailing...");
			return Read();
		}
		std::string outputName = outputDir + "/input" + format.getSuffix();
		{
			OfstreamMap ofm(outputName, "");
			inputReads.writeAll(ofm.getOfstream(""), format);
		}
		std::string log = outputName + ".log";
		std::string cmd = Cap3Options::getOptions().getCap3Path() + " " + outputName + " > " + log + " 2>&1";
		LOG_DEBUG_OPTIONAL(1, true, "Executing: " << cmd);
		status = system(cmd.c_str());
		if (status == 0) {
			std::string newContigFile = outputName + ".cap.contigs";
			long fileSize = FileUtils::getFileSize(newContigFile);
			if (fileSize > 0) {
				ReadSet newContig;
				newContig.appendFastaFile(newContigFile);
				Read bestRead = selectBestContig(newContig, oldContig);
				clean(outputDir);
				return bestRead;
			}
		}
		LOG_WARN(1, "Could not assemble " << oldContig.getName() << " with pool of " << _inputReads.getSize() << " reads: " << FileUtils::dumpFile(log));

		clean(outputDir);
		return Read();
	}
	static void clean(std::string fileDir) {
		std::string cmd;
		cmd = "/bin/rm -rf " + fileDir;
		int status = system(cmd.c_str());
		if (status != 0)
			LOG_WARN(1, "Could not clean up directory: " + fileDir);
	}
};

#endif /* CAP3_H_ */
