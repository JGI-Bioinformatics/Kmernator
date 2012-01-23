#!/bin/bash
test=TestMPI
MPI=""

true=$(which true)
if mpirun $true
then
  MPI="mpirun -np"
elif aprun -n 1 $true
then
  MPI="aprun -n"
fi

set -e
set -x
if [ -n "$MPI" ]
then
  for mpi in {1..8}
  do
    $MPI $mpi $test
  done
fi

