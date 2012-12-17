# **Kmernator**

* An MPI Toolkit for large scale genomic analysis

Kmernator is a large, scalable toolkit designed to efficiently analyse and
reduce genome, metagenome, transcriptome and metatranscriptome raw data.  There
are many features to remove redundant or errant data and the code runs and they
scale well on hunderds of computers tackling tera-base sized datasets.   

*****************
*****************

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

*****************
*****************


See INSTALL for instruction on how to build and cmake-flags for platform
specific hints


*****************
*****************


### *The main features of the toolkit incude*:

* A re-usable MPI based C++ header library core for parallel and scalable genomic data
  IO and analysis.

#### *With implementations of the following applications*:

*FilterReads* and *FilterReads-P* (the MPI version):

* A memory efficient, parallel and scalable kmer counter and raw data filter which
  can simultaneuosly remove errant data below a minimum kmer count and normalize overly
  abundant or redundant data.  This combination dramatically reduces the computation
  costs of subsequent assembly and alignment steps and typically improves the asseblies
  as well.

*DistributedNucleatingAssembler*:

* A parallel and scalable local greedy assembler capable of accurately assembling
  multiple alleles of genes and operons which are represented within a metagenome by
  extending targeted seed sequences.

*TNFDistance*:

*  A tetra-nucleotide clustering algorithm designed to optimally collect short contigs
   and scaffolds into a bag of genes representing a nearly complete genome.

*SplitSequenceOnTheFly*:

*  A parallel and scalable tool to partition genomic datasets and run arbitrary software
   on each partition without consuming additional disk space.

*RandomlySample*:

*  A fast and efficient tool to randomly sample genomic data sets.

*BamSort-P*:

*  A parallel and scalable tool to efficiently sort, and optionally filter, many or large
   SAM and BAM files.


*****************
*****************

##### Example meta-genomic workflow:

*  Calculate the total size of your raw, uncompressed fastq files
*  Schedule a job in your batch scheduler (like SGE).  Size the job
   so that collective RAM >= 3 * (total_input_file_size).
*  Execute FilterReads-P in the job (this example assumes mpirun is
   integrated with the batch scheduler to automatically understand the
   job size and shape)

>    This will remove or trim all reads from the input set which contain
>    a 31-mer only observed once.  Additionally it will normalize
>    highly abundant and/or redundant reads to an average kmer depth
>    of 100.  Illumina primer-dimer and long homopolymer regions are 
>    also removed by default

    mpirun FilterReads-P --output-file FilteredData \  
      --max-kmer-depth 100 --min-kmer-depth 2 \  
      31 input*.fastq

*  Assemble with your favorite assembler on FilteredData*.fastq (not shown)

######  Align your original reads to your assembly 
*  (optional) Split your assembly into 4 parts (if >4G as needed by bowtie2)

>   This will parse assembly.fa, remove any contigs lt 350 bases
>   and partition them into 4 separate files "{Uniq}" is expanded by SplitSequenceOnTheFly

    mpirun -n 4 SplitSequenceOnTheFly \  
       --format-output 1 --min-read-length 350 \  
       --output-file "myassembly-{Uniq}" assembly.fa

*  (optional) run bowtie2-align on each of the 4 partitions of the assembly


>    This will start 4 bowtie2-build processes on the 4 parations of the assembly

    mpirun -n 4 SplitSequenceOnTheFly \
       --fork-command "bowtie2-build myassembly-{Uniq} myassembly-{Uniq}-bowtie2idx" /dev/null


* Align the raw data to the assembly in parallel using a large job


>   Assuming an interleaved fastq with pairs alternating, this command will
>   spawn multiple bowtie2 processes and provide the input data split
>   through unix fifo files

    f1=".{Uniq}.p{UniqSecond}.1"
    f2=".{Uniq}.p{UniqSecond}.2"
    s=".{Uniq}.p{UniqSecond}.sam"
    runbowtie="bowtie2 -p 2 --reorder -x myassembly-{Uniq}-bowtie2idx -1 $f1 -2 $f2 -S $s"

    mpirun SplitSequenceOnTheFly --trim-pair-in-name \
      --second-dim 4 \ # leave out if the index is *not* partitioned
      --output-fifo 1 --output-file  "$f1" --split-file "$f2" \
      --fork-command "$runbowtie" input*.fastq

* Sort the resulting sam files

>   This will read all the sam files into memory, remove and compress
>   the unmapped read-pairs (where both are unmapped)
>   make sure that the total ram of the jobs is >= 2x the total size of the sams

    mpirun BamSort-P \  
       --num-partitions 4 \  # leave out if the index is *not* partitioned
       --unmapped-read-pairs myassembly-unmapped.fastq.gz \
       myassembly.bam .*of*.p*of*.sam && rm -f .*of*.p*of*.sam


