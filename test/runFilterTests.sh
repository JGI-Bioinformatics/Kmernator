#!/bin/bash

FR=../apps/FilterReads
FRP=../apps/FilterReads-P
procs=$(($(lscpu -p | tail -1 | awk -F, '{print $2}')+1) || echo 2)

TMP=$(mktemp testXXXXXX)
export TMPDIR=/tmp
if [ ! -f "$TMP" ]
then
  exit 1
fi
rm $TMP
KEEP=${KEEP:=0}
clean()
{
  [ X$KEEP != X0 ] || rm -rf $TMP*
}
trap clean 0 1 2 3 15

IN=
GOOD=

check()
{
  opts=" --kmer-scoring-type MEDIAN --mask-simple-repeats 0 --artifact-edit-distance 1 --out $TMP 31 $IN"
  echo "Executing: $@ $opts"
  if $@ $opts
  then
    if ! diff -w -q $TMP-MinDepth2-$IN $GOOD
    then
       echo "FAILED $@ --out $TMP 31 $IN"
       wc $TMP-MinDepth2-$IN $GOOD
       diff -w $TMP-MinDepth2-$IN $GOOD | head -50
       exit 1
    fi
  else
    echo "FAILED with exit status $?: $@ --out $TMP 31 $IN"
    exit 1
  fi
}

// make sure base quality conversions work fine
IN=1000.fastq
GOOD=1000-Filtered-0.85.std.fastq
check $FR --fastq-output-base-quality 33 --min-read-length 0.85
GOOD=1000-Filtered-0.85.fastq
check $FR --fastq-output-base-quality 64 --min-read-length 0.85
IN=1000.std.fastq
GOOD=1000-Filtered-0.85.std.fastq
check $FR --fastq-output-base-quality 33 --min-read-length 0.85
GOOD=1000-Filtered-0.85.fastq
check $FR --fastq-output-base-quality 64 --min-read-length 0.85

IN=1000.fastq
GOOD=1000-Filtered-0.85.fastq
check $FR --fastq-output-base-quality 64 --min-read-length 0.85
GOOD=1000-Filtered-readlength.fastq
check $FR --fastq-output-base-quality 64 --min-read-length 1 
GOOD=1000-Filtered-readlength-both.fastq
check $FR --fastq-output-base-quality 64 --min-read-length 1 --min-passing-in-pair 2
rm -f $TMP*


IN=1000.fastq
GOOD=1000-Filtered.fastq

for thread in {1..3}
do
  check $FR --fastq-output-base-quality 64 --min-read-length 25 --thread $thread
  rm -f $TMP*
  check $FR --fastq-output-base-quality 64 --min-read-length 25 --thread $thread --save-kmer-mmap 1
  mv $TMP-mmap $TMP-mmap-saved
  check $FR --fastq-output-base-quality 64 --min-read-length 25 --thread $thread --load-kmer-mmap $TMP-mmap-saved
  rm -f $TMP*
done

MPI=""
MPI_OPTS=""

if mpirun /bin/true
then
  MPI="mpirun"
  MPI_OPTS="-np"
elif aprun -n 1 /bin/true
then
  MPI="aprun"
  MPI_OPTS="-n"
fi

if [ -n "$MPI" ]
then
  for mpi in 1 2 3 4 6 7 8 12 13 16 20 24 32
  do
    if [ $mpi -gt $procs ]
    then
      break
    fi
    export OMP_NUM_THREADS=$(((procs+mpi-1)/mpi))
    GOOD=1000-Filtered-0.85.fastq
    check $MPI $MPI_OPTS $mpi $FRP --fastq-output-base-quality 64 --min-read-length 0.85
    rm -f $TMP*
 
    GOOD=1000-Filtered-readlength.fastq
    check $MPI $MPI_OPTS $mpi $FRP --fastq-output-base-quality 64 --min-read-length 1     
    rm -f $TMP*
 
    GOOD=1000-Filtered-readlength-both.fastq
    check $MPI $MPI_OPTS $mpi $FRP --fastq-output-base-quality 64 --min-read-length 1 --min-passing-in-pair 2
    rm -f $TMP*
 
    GOOD=1000-Filtered.fastq
    check $MPI $MPI_OPTS $mpi $FRP --fastq-output-base-quality 64 --min-read-length 25 --thread 1
    rm -f $TMP*
    check $MPI $MPI_OPTS $mpi $FRP --fastq-output-base-quality 64 --min-read-length 25
    rm -f $TMP*
    
    # TODO restore save/load kmer map in MPI version...
    #check $MPI $MPI_OPTS $mpi $FRP --fastq-output-base-quality 64 --min-read-length 25 --thread 1 --save-kmer-mmap 1
    #mv $TMP-mmap $TMP-mmap-saved
    #check $MPI $MPI_OPTS $mpi $FRP --fastq-output-base-quality 64 --min-read-length 25 --thread 1 --load-kmer-mmap $TMP-mmap-saved
    #check $FR --min-read-length 25 --fastq-output-base-quality 64 --load-kmer-mmap $TMP-mmap-saved
    #rm -f $TMP*
    check $FR --min-read-length 25 --fastq-output-base-quality 64 --save-kmer-mmap 1 
    mv $TMP-mmap $TMP-mmap-saved
    check $MPI $MPI_OPTS $mpi $FRP --fastq-output-base-quality 64 --min-read-length 25 --load-kmer-mmap $TMP-mmap-saved
    rm -f $TMP*
  done
fi


