#!/bin/bash

MPI=""
procs=$(grep -c ^processor /proc/cpuinfo)

true=$(which true)

mt=time
if memtime $true
then
  mt=memtime
fi

ismpi=0
if mpirun $true
then
  MPI="mpirun -bycore -bind-to-core -np"
  ismpi=1
elif aprun -n 1 $true
then
  MPI="aprun -n"
  procs=$(aprun -B uname -n | wc -l)
fi


SSOTF=../apps/SplitSequenceOnTheFly

TMPDIR=${TMPDIR:=/tmp}
TMP=$(mktemp $TMPDIR/testXXXXXX)
if [ ! -f "$TMP" ]
then
  exit 1
fi
rm $TMP

fail()
{
  echo $@
  /bin/false
}

set -e
set -x

IN=1000.fastq
total=$(grep -c ^@ $IN)

if [ -n "$MPI" ]
then

  for mpi in $(seq 1 ${procs})
  do
    export OMP_NUM_THREADS=1
    rm -f ${TMP}*
    $mt $MPI $mpi $SSOTF --output-file "$TMP-{Uniq}.fastq" $IN
    cat ${TMP}-*.fastq | diff -q - $IN

    rm -f ${TMP}*
    $mt $MPI $mpi $SSOTF --output-file "$TMP-{Uniq}.fastq.1" --split-file "$TMP-{Uniq}.fastq.2" $IN
    c1=$(cat ${TMP}*.1 | grep -c ^@)
    c2=$(cat ${TMP}*.2 | grep -c ^@)
    if ((c1+c2 != total))
    then
      fail "Split counts differ $c1 + $c2 != $total"
    fi

    if ((mpi%2 == 0))
    then
      rm -f ${TMP}*
      $mt $MPI $mpi  $SSOTF --second-dim 2 --output-file "$TMP-{UniqFirst}x{UniqSecond}.fastq" $IN
      cat ${TMP}-*x000000of000002.fastq | diff -q - $IN
      cat ${TMP}-*x000001of000002.fastq | diff -q - $IN

      rm -f ${TMP}*
      $mt $MPI $mpi $SSOTF --second-dim 2 --output-file "$TMP-{UniqFirst}x{UniqSecond}.fastq.1" --split-file "$TMP-{UniqFirst}x{UniqSecond}.fastq.2" $IN
      c1=$(cat ${TMP}*x000000of000002.fastq.1 | grep -c ^@)
      c2=$(cat ${TMP}*x000000of000002.fastq.2 | grep -c ^@)
      if ((c1+c2 != total))
      then
        fail "Split counts differ $c1 + $c2 != $total"
      fi

      c1=$(cat ${TMP}*x000001of000002.fastq.1 | grep -c ^@)
      c2=$(cat ${TMP}*x000001of000002.fastq.2 | grep -c ^@)
      if ((c1+c2 != total))
      then
        fail "Split2 counts differ $c1 + $c2 != $total"
      fi

    fi

  done
fi

