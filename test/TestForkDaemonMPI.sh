#!/bin/bash
test=TestForkDaemonMPI
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
  for mpi in ${procs}
  do
    $MPI $mpi $test
  done
fi

