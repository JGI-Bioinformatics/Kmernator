#!/bin/bash

RS=../apps/RandomlySample

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

