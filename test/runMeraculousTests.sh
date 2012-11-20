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
  procs=$(aprun -B -q uname -n | wc -l)
fi


FRP=../apps/MeraculousCounter

TMPDIR=${TMPDIR:=/tmp}
TMP=$(mktemp $TMPDIR/testXXXXXX)
if [ ! -f "$TMP" ]
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

if ((isaprun==0))
then
  for thread in {1..4}
  do
    check $FRP --thread $thread
    rm -f $TMP*
  done
fi

if [ -n "$MPI" ]
then

  for mpi in $(seq 1 ${procs})
  do
    export OMP_NUM_THREADS=1
    check $mt $MPI $mpi $FRP
  done
  
  if ((ismpi)) && mpirun -n $procs -bysocket -bind-to-socket $true
  then
   MPI="mpirun -bysocket -np"
  fi
   
  echo "Running in hybrid mode"
 
  for threads in $(seq 1 $procs)
  do
    socket=
    if ((ismpi)) && ((threads<=procs/2))
    then
      socket='-bind-to-socket'
    fi

    mpi=$(((procs+threads-1)/threads))
    
    if ((ismpi==0))
    then
      socket="-d $threads"     
    fi

    export OMP_NUM_THREADS=$threads
    check $mt $MPI $mpi $socket $FRP
  done

fi

rm -f $TMP*

