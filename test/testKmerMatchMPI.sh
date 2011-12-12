#!/bin/bash
test=TestKmerMatchMPI
MPI=""

if mpirun /bin/true
then
  MPI="mpirun -np"
elif aprun -n 1 /bin/true
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

