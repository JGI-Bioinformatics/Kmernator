#!/bin/bash

FR=../apps/FilterReads
FRP=../apps/FilterReads-P

TMPDIR=${TMPDIR:=/tmp}
TMP=$(mktemp $TMPDIR/test$$)
if [ ! -f $TMP ]
then
  exit 1
fi
rm $TMP

IN=1000.fastq
GOOD=1000-Filtered.fastq

check()
{
  echo "Executing: $@ --mask-simple-repeats 0 --artifact-edit-distance 1 --out $TMP 31 $IN"
  if $@ --mask-simple-repeats 0 --artifact-edit-distance 1 --out $TMP 31 $IN
  then
    if ! diff -q $TMP-MinDepth2-1000.fastq $GOOD
    then
       echo "FAILED $@ --out $TMP 31 $IN"
       wc $TMP-MinDepth2-1000.fastq $GOOD
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
  mv $TMP-mmap $TMP-mmap-saved
  check $FR --thread $thread --load-kmer-mmap $TMP-mmap-saved
  rm -f $TMP*
done

MPI=""

if mpirun /bin/true
then
  MPI="mpirun -np"
elif aprun -n 1 /bin/true
then
  MPI="aprun -n"
fi

if [ -n "$MPI" ]
then
  for mpi in {1..5}
  do
    check $MPI $mpi $FRP --thread 1
    rm -f $TMP*
    check $MPI $mpi $FRP --thread 1 --save-kmer-mmap 1
    mv $TMP-mmap $TMP-mmap-saved
    check $MPI $mpi $FRP --thread 1 --load-kmer-mmap $TMP-mmap-saved
    check $FR --load-kmer-mmap $TMP-mmap-saved
    rm -f $TMP*
    check $FR --thread $thread --save-kmer-mmap 1 
    mv $TMP-mmap $TMP-mmap-saved
    check $MPI $mpi $FRP --thread 1 --load-kmer-mmap $TMP-mmap-saved
    rm -f $TMP*
  done
fi

rm -f $TMP*

