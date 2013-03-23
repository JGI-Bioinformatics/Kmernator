#!/bin/bash
test=TestMPI
MPI=""
procs=$(($(lscpu -p | tail -1 | awk -F, '{print $2}')+1))

true=$(which true)
if mpirun $true
then
  MPI="mpirun -np"
elif aprun -n 1 $true
then
  MPI="aprun -n"
fi

export TMPDIR=/tmp
set -e
set -x
if [ -n "$MPI" ]
then
  for mpi in 1 2 3 4 6 7 8 12 13 16 20 24 32
  do
    if [ $mpi -gt $procs ]
    then
      break
    fi
    export OMP_NUM_THREADS=$(((procs+mpi-1) / mpi))
    $MPI $mpi $test
  done
fi

