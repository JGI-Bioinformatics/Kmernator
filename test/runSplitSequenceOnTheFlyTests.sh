#!/bin/bash

MPI=""
procs=$(($(lscpu -p | tail -1 | awk -F, '{print $2}')+1))

true=$(which true)

ismpi=0
if mpirun $true
then
  MPI="mpirun -bycore -bind-to-core -np"
  ismpi=1
elif aprun -n 1 $true
then
  MPI="aprun -n"
  procs=$(aprun -B -q uname -n | wc -l)
fi


SSOTF=../apps/SplitSequenceOnTheFly

TMP=$(mktemp testXXXXXX)
export TMPDIR=/tmp
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
KEEP=${KEEP:=0}
clean()
{
  [ X$KEEP != X0 ] || rm -rf $TMP*
}
trap clean 0 1 2 3 15

set -e
set -x

IN=1000.fastq
total=$(grep -c ^@ $IN)

if [ -n "$MPI" ]
then

  for mpi in 1 2 3 4 6 7 8 12 13 16 20 24 32
  do
    if [ $mpi -gt $procs ]
    then
      break
    fi
    export OMP_NUM_THREADS=$(((procs+mpi-1) / mpi))
    rm -f ${TMP}*
    $MPI $mpi $SSOTF --output-file "$TMP-{Uniq}.fastq" $IN
    cat ${TMP}-*.fastq | diff -q - $IN

    rm -f ${TMP}*
    $MPI $mpi $SSOTF --output-file "$TMP-{Uniq}.fastq.1" --split-file "$TMP-{Uniq}.fastq.2" $IN
    c1=$(cat ${TMP}*.1 | grep -c ^@)
    c2=$(cat ${TMP}*.2 | grep -c ^@)
    if ((c1+c2 != total))
    then
      fail "Split counts differ $c1 + $c2 != $total"
    fi

    if ((mpi%2 == 0))
    then
      rm -f ${TMP}*
      $MPI $mpi  $SSOTF --second-dim 2 --output-file "$TMP-{UniqFirst}x{UniqSecond}.fastq" $IN
      cat ${TMP}-*x000000of000002.fastq | diff -q - $IN
      cat ${TMP}-*x000001of000002.fastq | diff -q - $IN

      rm -f ${TMP}*
      $MPI $mpi $SSOTF --second-dim 2 --output-file "$TMP-{UniqFirst}x{UniqSecond}.fastq.1" --split-file "$TMP-{UniqFirst}x{UniqSecond}.fastq.2" $IN
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

