#!/bin/bash

SSOTF=../apps/SplitSequenceOnTheFly

TMPDIR=${TMPDIR:=/tmp}
TMP=$(mktemp $TMPDIR/testXXXXXX)
if [ ! -f "$TMP" ]
then
  exit 1
fi
rm $TMP

set -e
set -x

IN=1000.fastq
total=$(grep -c ^@ $IN)

for p in {1..10}
do
  rm -f ${TMP}*
  mpirun -np $p $SSOTF --output-file "$TMP-{Uniq}.fastq" $IN
  cat ${TMP}-*.fastq | diff -q - $IN

  rm -f ${TMP}*
  mpirun -np $p $SSOTF --output-file "$TMP-{Uniq}.fastq.1" --split-file "$TMP-{Uniq}.fastq.2" $IN
  c1=$(cat ${TMP}*.1 | grep -c ^@)
  c2=$(cat ${TMP}*.2 | grep -c ^@)
  if ((c1+c2 != total))
  then
    echo "Split counts differ $c1 + $c2 != $total"
    exit 1
  fi

  if ((p%2 == 0))
  then
    rm -f ${TMP}*
    mpirun -np $p $SSOTF --second-dim 2 --output-file "$TMP-{UniqFirst}x{UniqSecond}.fastq" $IN
    cat ${TMP}-*x000000of000002.fastq | diff -q - $IN
    cat ${TMP}-*x000001of000002.fastq | diff -q - $IN

    rm -f ${TMP}*
    mpirun -np $p $SSOTF --second-dim 2 --output-file "$TMP-{UniqFirst}x{UniqSecond}.fastq.1" --split-file "$TMP-{UniqFirst}x{UniqSecond}.fastq.2" $IN
    c1=$(cat ${TMP}*x000000of000002.fastq.1 | grep -c ^@)
    c2=$(cat ${TMP}*x000000of000002.fastq.2 | grep -c ^@)
    if ((c1+c2 != total))
    then
      echo "Split counts differ $c1 + $c2 != $total"
      exit 1
    fi

    c1=$(cat ${TMP}*x000001of000002.fastq.1 | grep -c ^@)
    c2=$(cat ${TMP}*x000001of000002.fastq.2 | grep -c ^@)
    if ((c1+c2 != total))
    then
      echo "Split2 counts differ $c1 + $c2 != $total"
      exit 1
    fi

  fi
done
