#!/bin/bash


FRP=../apps/MeraculousCounter

TMPDIR=${TMPDIR:=/tmp}
TMP=$(mktemp $TMPDIR/test$$)
if [ ! -f $TMP ]
then
  exit 1
fi
rm $TMP

IN=1000.fastq
TEST1=${TMP}.mercount.m21
GOOD1=phix.mercount.m21
TEST2=${TMP}.mergraph.m21.D2
GOOD2=phix.mergraph.m21.D2

check()
{
  opts=" --min-kmer-quality=0 --min-quality-score=2 --kmer-size 21 --out $TMP $IN"
  echo "Executing: $@ $opts"
  if $@ $opts
  then
    if ! sort $TEST1 | diff -q - $GOOD1
    then
       echo "FAILED $@ $opts"
       wc $TEST1 $GOOD1
       sort $TEST1 | diff - $GOOD1 | head -50
       rm -f $TMP*
       exit 1
    fi
    if ! sort $TEST2 | diff -q - $GOOD2
    then
       echo "FAILED $@ $opts"
       wc $TEST2 $GOOD2
       sort $TEST2 | diff - $GOOD2 | head -50
       rm -f $TMP*
       exit 1
    fi
    
  else
    echo "FAILED with exit status $?: $@ $opts"
    rm -f $TMP*
    exit 1
  fi
}

for thread in {1..4}
do
  check $FRP --thread $thread
  rm -f $TMP*
done

if mpirun /bin/true
then
  for thread in {1..3}
  do
  for mpi in {1..5}
  do
    check mpirun -np $mpi $FRP --thread $thread
    rm -f $TMP*
  done
  done
fi

rm -f $TMP*

