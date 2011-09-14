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

class _Cap3Options : public OptionsBaseInterface {
public:
	static std::string getCap3Path() {
		return getVarMap()["cap3-path"].as<std::string>();
	}
	static int getCap3MaxReads() {
		return getVarMap()["cap3-max-reads"].as<int>();
	}
	void _resetDefaults() {
		Options::getOptions().getMmapInput() = 0;
	}
	void _setOptions(po::options_description &desc,
			po::positional_options_description &p) {
		po::options_description opts("Cap3 Options");

		opts.add_options()
		("cap3-path", po::value<std::string>()->default_value(""),
				"if set, cap3 will be used to extend contigs")

		("cap3-max-reads", po::value<int>()->default_value(500),
				"maximum number of reads to send to cap 3 per contig extension (randomly sampled, of course)");

		desc.add(opts);
	}
};
typedef OptionsBaseTemplate< _Cap3Options > Cap3Options;

class Cap3 {
public:
	static const int repeatContig = 10;
	static ReadSet extendContig(const Read &oldContig, const ReadSet &_inputReads) {
		ReadSet inputReads = _inputReads.randomlySample(Cap3Options::getOptions().getCap3MaxReads());

		for(int i = 0; i < repeatContig; i++)
			inputReads.append(oldContig);

		FormatOutput format = FormatOutput::FastaUnmasked();
		std::string outputDir = Options::getOptions().getTmpDir() + UniqueName::generateUniqueName("/.cap3-assembly");
		int status = mkdir(outputDir.c_str(), 0x770);
		if (status != 0) {
			LOG_WARN(1, "Could not mkdir: " << outputDir << " bailing...");
			return ReadSet();
		}
		std::string outputName = outputDir + "input" + format.getSuffix();
		{
			OfstreamMap ofm(outputName, "");
			inputReads.writeAll(ofm.getOfstream(""), format);
		}
		std::string cmd = Cap3Options::getOptions().getCap3Path() + " " + outputName + " > " + outputName + ".log" + " 2>&1";
		LOG_DEBUG_OPTIONAL(1, true, "Executing: " << cmd);
		status = system(cmd.c_str());
		if (status == 0) {
			std::string newContigFile = outputName + ".cap3.contigs.fa";
			long fileSize = FileUtils::getFileSize(newContigFile);
			if (fileSize > 0) {
				ReadSet newContig;
				newContig.appendFastaFile(newContigFile);
				clean(outputDir);
				return newContig;
			}
		}
		LOG_WARN(1, "Could not assemble " << oldContig.getName());
		clean(outputDir);
		return ReadSet();
	}
	static void clean(std::string fileDir) {
		std::string cmd;
		cmd = "rm -rf " + fileDir;
		int status = system(cmd.c_str());
		if (status != 0)
			LOG_WARN(1, "Could not clean up directory: " + fileDir);
	}
};

#endif /* CAP3_H_ */
