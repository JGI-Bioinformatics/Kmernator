/*
 * FilterReadsTemplate.h
 *
 *  Created on: Oct 19, 2010
 *      Author: regan
 */

#ifndef FILTER_READS_TEMPLATE_H_
#define FILTER_READS_TEMPLATE_H_

long selectReads(unsigned int minDepth, ReadSet &reads, KS &spectrum, RS &selector, std::string outputFilename)
{

	LOG_VERBOSE_OPTIONAL(1, true, "selectReads with minDepth " << minDepth << ": " << reads.getSize());
	LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());

	long oldPicked = 0;
	long picked = 0;

	int maximumKmerDepth = FilterReadsOptions::getMaxKmerDepth();

	OfstreamMap ofmap(outputFilename, ".fastq");

	if (maximumKmerDepth > 0) {
		for (int depth = 1; depth < maximumKmerDepth; depth++) {
			LOG_VERBOSE_OPTIONAL(2, true, "Picking depth " << depth << " layer of reads");
			if (reads.hasPairs())
				picked += selector.pickBestCoveringSubsetPairs(depth,
						minDepth, Options::getMinReadLength(), FilterReadsOptions::getBothPairs());
			else
				picked += selector.pickBestCoveringSubsetReads(depth,
						minDepth, Options::getMinReadLength());
			LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());
		}

		if (picked > 0 && !outputFilename.empty()) {
			LOG_VERBOSE_OPTIONAL(1, true, "Writing " << picked << " reads to output file(s)");
			selector.writePicks(ofmap, oldPicked);
		}
		LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());
		oldPicked += picked;


	} else {

		int maxDepth = FilterReadsOptions::getPartitionByDepth();
		if (maxDepth < 0) {
			maxDepth = 1;
		}

		for (unsigned int depth = maxDepth; depth >= 1; depth /= 2) {

			string ofname = outputFilename;
			if (maxDepth > 1) {
				ofname += "-PartitionDepth" + boost::lexical_cast< string >( depth );
			}
			ofmap = OfstreamMap(ofname, ".fastq");
			float tmpMinDepth = std::max(minDepth, depth);
			if (Options::getKmerSize() == 0) {
				tmpMinDepth = 0;
				depth = 0;
			}
			LOG_VERBOSE_OPTIONAL(1, true, "Selecting reads over depth: " << depth << " (" << tmpMinDepth << ") ");

			if (reads.hasPairs()) {
				picked = selector.pickAllPassingPairs(tmpMinDepth,
						Options::getMinReadLength(),
						FilterReadsOptions::getBothPairs());
			} else {
				picked = selector.pickAllPassingReads(tmpMinDepth,
						Options::getMinReadLength());
			}
			LOG_VERBOSE_OPTIONAL(2, true, "At or above coverage: " << depth << " Picked " << picked
			<< " / " << reads.getSize() << " reads");
			LOG_DEBUG_OPTIONAL(1, true, MemoryUtils::getMemoryUsage());

			if (picked > 0 && !outputFilename.empty()) {
				LOG_VERBOSE_OPTIONAL(1, true, "Writing " << picked << " reads  to output files");
				selector.writePicks(ofmap, oldPicked);
			}
			oldPicked += picked;

			if (minDepth > depth) {
				break;
			}
		}
	}
	ofmap.clear();
	LOG_VERBOSE_OPTIONAL(1, true, "Done.  Cleaning up. " << MemoryUtils::getMemoryUsage());

	return oldPicked;
}


#endif /* FILTER_READS_TEMPLATE_H_ */
