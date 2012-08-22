/*
 * FilterAndSortSam.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: regan
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "sam_header.h"
#include "sam.h"
#include "bam.h"
#include "kstring.c"
#include "ksort.h"

int bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn,
                                        int flag, const char *reg);

int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: FilterAndSortSam [options] input.sam sortedOutput.bam\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-u\tPlace unaligned reads in a separate file: sortedOutput.bam.unaligned.fastq.gz\n");
	fprintf(stderr, "\t-U\tDo not include unaligned reads in the sorted BAM\n");
	fprintf(stderr, "\n\nThis program is optimized for BAMs with very large headers, so it is in memory only once\n\n");
	return -1;
}

int bamIsMapped(bam1_t *b) {
	return ((b->core.flag & BAM_FUNMAP) == BAM_FUNMAP) ? 0 : 1;;
}
int bamIsKept(bam1_t *b, int noUnaligned) {
	if (noUnaligned == 0)
		return 1;
	if ((b->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
		// paired, keep if at least one of pair is mapped
		return ((b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP)) == (BAM_FUNMAP | BAM_FMUNMAP)) ? 0 : 1;
	} else {
		// single, keep if mapped
		return bamIsMapped(b);
	}
}

void getPair(bam1_t *b, char *pairTag) {
	if ((b->core.flag & BAM_FPAIRED) == BAM_FPAIRED) {
		sprintf(pairTag, "/%d", (b->core.flag & BAM_FREAD1) ? 1 : 2);
	}
}

void printFastq(FILE *f, bam1_t *b) {
	int i;
	char pairTag[3];
	uint8_t *s = bam1_seq(b), *t = bam1_qual(b);
	pairTag[0] = '\0';
	getPair(b, pairTag);
	const bam1_core_t *c = &b->core;
        kstring_t str;
        str.l = str.m = 0; str.s = 0;
	kputc('@', &str);
        kputsn(bam1_qname(b), c->l_qname-1, &str);
        kputsn(pairTag, strlen(pairTag), &str);
	kputc('\n', &str);
	if (c->l_qseq) {
                for (i = 0; i < c->l_qseq; ++i) kputc(bam_nt16_rev_table[bam1_seqi(s, i)], &str);
                kputc('\n', &str);
                kputc('+', &str);
                kputc('\n', &str);
                if (t[0] == 0xff) for(i=0; i < c->l_qseq; ++i) kputc(93, &str);
                else for (i = 0; i < c->l_qseq; ++i) kputc(t[i] + 33, &str);
        } else kputsn("N\n+\nA", 6, &str);

	//fprintf(f, "@%s%s\n%s\n+\n%s\n", bam1_qname(b), pairTag, bam1_seq(b), bam1_qual(b));
	fprintf(f, "%s\n", str.s);
}

/* Derived from bam_sort.c ; modified modestly */

typedef bam1_t *bam1_p;
static inline int bam1_lt2(const bam1_p a, const bam1_p b)
{
        return (((uint64_t)a->core.tid<<32|(a->core.pos+1)) < ((uint64_t)b->core.tid<<32|(b->core.pos+1)));
}
KSORT_INIT(sort2, bam1_p, bam1_lt2)

static void sort_blocks(int n, int k, bam1_p *buf, const char *prefix, const bam_header_t *h, int is_stdout)
{
	char *name, mode[3];
	int i;
	bamFile fp;
	fprintf(stderr, "Sorting %d reads\n", k);
	ks_mergesort(sort2, k, buf, 0);
	name = (char*)calloc(strlen(prefix) + 30, 1);
	if (n >= 0) {
		sprintf(name, "%s.%.4d.tmpsort", prefix, n);
		strcpy(mode, "w1");
	} else {
		sprintf(name, "%s", prefix);
		strcpy(mode, "w");
	}
	fp = is_stdout? bam_dopen(fileno(stdout), mode) : bam_open(name, mode);
	if (fp == 0) {
		fprintf(stderr, "[sort_blocks] fail to create file %s.\n", name);
		free(name);
		// FIXME: possible memory leak
		return;
	}
	fprintf(stderr, "Writing %d sorted reads to %s\n", k, name);
	free(name);
	bam_header_write(fp, h);
	for (i = 0; i < k; ++i)
		bam_write1_core(fp, &buf[i]->core, buf[i]->data_len, buf[i]->data);
	bam_close(fp);
}

int main(int argc, char *argv[]) {
	samfile_t *in = 0;
	char *outBam, *cmd;
	int c, unalignedFastq = 0, noUnaligned = 0;
	int bytesRead;
	bam1_t **bamBuf;
	bam_header_t *header = 0;
	size_t mem = 0, max_mem = 512*1024*1024; // 512MB
	int numSortedFiles = 0, bamBufOffset = 0;
	unsigned long mapped = 0, unmapped = 0;

	FILE *unalignedFastqFile = 0;
	while ((c = getopt(argc, argv, "uU")) >= 0) {
		switch (c) {
		case 'u' : unalignedFastq = 1; break;
		case 'U' : noUnaligned = 1; break;
		default: return usage();
		}
	}
	if (argc != optind + 2) {
		return usage();
	}
	outBam = strdup(argv[optind+1]);

	if((in = samopen(argv[optind], "r", 0)) == 0) {
		fprintf(stderr, "Could not open %s as a SAM file!\n", argv[optind]);
		return -1;
	}
	if (in->header == 0) {
		fprintf(stderr, "[main_samview] fail to read the header!\n");
		return -1;
	}
	header = in->header;

	if (unalignedFastq) {
		cmd = (char *) malloc(strlen(outBam) + 128);
		fprintf(stderr, "Writing unaligned reads to %s.unaligned.fastq.gz\n", outBam);
		sprintf(cmd, "gzip -c > %s.unaligned.fastq.gz", outBam);

		if((unalignedFastqFile = popen(cmd, "w")) == 0) {
			fprintf(stderr, "Could not open '%s' for unaligned fastq output\n", cmd);
			return -1;
		}
		free(cmd); cmd = 0;
	}

	bamBuf = (bam1_t**)calloc((max_mem / BAM_CORE_SIZE) + 1, sizeof(bam1_t*));
	for (;;) {
		//while ((r = samread(in, b)) >= 0) {
		if(bamBuf[bamBufOffset] == 0)
			bamBuf[bamBufOffset] = (bam1_t*) calloc(1, sizeof(bam1_t));
		bam1_t *b = bamBuf[bamBufOffset];
		if ((bytesRead = samread(in, b)) < 0) break;
		int isMapped = bamIsMapped(b);
		if (isMapped == 1) {
			mapped++;
		} else {
			unmapped++;
		}
		if (unalignedFastq && isMapped == 0) {
			printFastq(unalignedFastqFile, b);
		}
		if (bamIsKept(b, noUnaligned) == 1) {
			mem += bytesRead;
			++bamBufOffset;
			if (mem >= max_mem) {
				sort_blocks(numSortedFiles++, bamBufOffset, bamBuf, outBam, header, 0);
				mem = 0;
				bamBufOffset = 0;
			}
		}
	}

	fprintf(stderr, "Finished reading SAM, finalizing sort\n");
	if (numSortedFiles == 0) {
		sort_blocks(-1, bamBufOffset, bamBuf, outBam, header, 0);
	} else { // then merge
		sort_blocks(numSortedFiles++, bamBufOffset, bamBuf, outBam, header, 0);
	}

	// free memory of headers, so merge can use it
	if (header != in->header)
		bam_header_destroy(header);
	samclose(in);

	// free memory
	for (bamBufOffset = 0; bamBufOffset < max_mem / BAM_CORE_SIZE; ++bamBufOffset) {
		if (bamBuf[bamBufOffset]) {
			free(bamBuf[bamBufOffset]->data);
			free(bamBuf[bamBufOffset]);
		}
	}
	free(bamBuf);

	// free memory, close pipe
	if (unalignedFastq) {
		if(pclose(unalignedFastqFile) != 0) {
			fprintf(stderr, "Could not close %s.unaligned.fastq.gz!", outBam);
			return -1;
		}
	}

	// merge the sorted files now
	if (numSortedFiles != 0) {
		fprintf(stderr, "Merging %d sorted files\n", numSortedFiles);
		char **fns;
		int i;
		fns = (char**)calloc(numSortedFiles, sizeof(char*));
		for (i = 0; i < numSortedFiles; ++i) {
			fns[i] = (char*)calloc(strlen(outBam) + 30, 1);
			sprintf(fns[i], "%s.%.4d.tmpsort", outBam, i);
		}

		char *reg = 0;
		bam_merge_core(0, outBam, 0, numSortedFiles, fns, 0, reg);
		for (i = 0; i < numSortedFiles; ++i) {
			unlink(fns[i]);
			free(fns[i]);
		}
		free(fns);
	}
	fprintf(stderr, "%ld mapped, %ld unmapped, %0.4lf%%\n", mapped, unmapped, 100.0 * (double) mapped / ((double) mapped + unmapped));

	return 0;
}
