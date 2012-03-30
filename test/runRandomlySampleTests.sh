#!/bin/bash

RS=../apps/RandomlySample

TMPDIR=${TMPDIR:=/tmp}
TMP=$(mktemp $TMPDIR/test$$)
if [ ! -f $TMP ]
then
  exit 1
fi
rm $TMP

IN=1000.fastq

check()
{
  samp=$1
  echo "Executing: $RS --num-samples $samp --out $TMP $IN"
  if $RS --num-samples $samp --out $TMP $IN
  then
    lines=$(cat $TMP | wc -l)
    if [ $((lines/8)) -ne $samp ]
    then
       echo "FAILED: $RS --num-samples $samp --out $TMP $IN"
       wc $TMP
       rm -f $TMP*
       exit 1
    fi
  else
    echo "FAILED with exit status $?: $RS --num-samples $samp --out $TMP $IN"
    rm -f $TMP*
    exit 1
  fi
}

for s in 1 2 5 10 46 50 100 200 400 450 475 490 495 496 496 498 499 500
do
  check $s
done
rm -f $TMP*

