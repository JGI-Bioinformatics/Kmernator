#!/bin/bash

MPI=""
MPI_OPTS=""
procs=$(($(lscpu -p | tail -1 | awk -F, '{print $2}')+1))

true=$(which true)

ismpi=0
isaprun=0
if mpirun $true
then
  MPI="mpirun"
  if $MPI -bycore -bind-to-core $true
  then
    MPI_OPTS="-bycore -bind-to-core -np"
  else
    MPI_OPTS="-np"
  fi
  ismpi=1
elif aprun -n 1 $true
then
  MPI="aprun"
  MPI_OPTS="-n"
  procs=$(aprun -B -q uname -n | wc -l)
  isaprun=1
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

  for mpi in 1 2 3 4 6 7 8 12 13 16 20 24 32
  do
    if [ $mpi -gt $procs ]
    then
      break
    fi
    export OMP_NUM_THREADS=1
     check $MPI $MPI_OPTS $mpi $FRP
  done
  
  if ((ismpi)) && $MPI -bysocket -bind-to-socket $true
  then
   MPI_OPTS="-bysocket -bind-to-socket -np"
  elif ((ismpi))
  then
   MPI_OPTS="-np"
  fi
   
  echo "Running in hybrid mode"
 
  for threads in $(seq 1 $procs)
  do
    socket=
    if ((ismpi)) && ((threads<=procs/2))
    then
      socket=''
    fi

    mpi=$(((procs+threads-1)/threads))
    
    if ((ismpi==0))
    then
      socket="-d $threads"     
    fi

    export OMP_NUM_THREADS=$threads
    check $MPI $MPI_OPTS $mpi $socket $FRP
  done

fi

rm -f $TMP*

