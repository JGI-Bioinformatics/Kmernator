#!/bin/bash

FR=../apps/FilterReads
FRP=../apps/FilterReads-P

TMP=$(mktemp)
rm $TMP

IN=1000.fastq
GOOD=1000-Filtered.fastq

check()
{
  if $@ --mask-simple-repeats 0 --artifact-edit-distance 1 --out $TMP 31 $IN
  then
    if ! diff -q $TMP-MinDepth2-1000.fastq $GOOD
    then
       echo "FAILED $@ --out $TMP 31 $IN"
       diff $TMP-MinDepth2-1000.fastq $GOOD | head -50
       rm -f $TMP*
       exit 1
    fi
  else
    echo "FAILED with exit status $?: $@ --out $TMP 31 $IN"
    rm -f $TMP*
    exit 1
  fi
}

for thread in {1..3}
do
  check $FR --thread $thread
  rm -f $TMP*
  check $FR --thread $thread --save-kmer-mmap 1
  check $FR --thread $thread --load-kmer-mmap 1
  rm -f $TMP*
done

if mpirun /bin/true
then
  for mpi in {1..5}
  do
    check mpirun -np $mpi $FRP --thread 1
    rm -f $TMP*
    check mpirun -np $mpi $FRP --thread 1 --save-kmer-mmap 1
    check mpirun -np $mpi $FRP --thread 1 --load-kmer-mmap 1
    check $FR --load-kmer-mmap 1
    rm -f $TMP*
    check $FR --thread $thread --save-kmer-mmap 1 
    check mpirun -np $mpi $FRP --thread 1 --load-kmer-mmap 1
    rm -f $TMP*
  done
fi

rm -f $TMP*

