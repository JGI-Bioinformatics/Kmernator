#!/bin/bash

set -e
set -x
procs=$(grep -c ^processor /proc/cpuinfo)

T=$(mktemp -d)
export TMPDIR=/tmp

KEEP=${KEEP:=0}
cleanup() 
{
  echo "Cleaning up" 1>&2
  [ X$KEEP != X0 ] || rm -rf ${T}
}
abort()
{
  cleanup
  exit 1
}
failed()
{
  echo "The last test failed: $@"
  exit 1
}
trap abort 1 2 3 15
trap cleanup 0

MPI=""
procs=$(grep -c ^processor /proc/cpuinfo)

true=$(which true)

mt=time
if memtime $true
then
  mt=memtime
fi

ismpi=0
isaprun=0
if mpirun $true
then
  MPI="mpirun -bycore -bind-to-core -np"
  ismpi=1
elif aprun -n 1 $true
then
  MPI="aprun -n"
  procs=$(aprun -B -q uname -n | wc -l)
  isaprun=1
fi

samtools=$(which samtools || /bin/true)
if [ -z "$samtools" ]
then
  samtools=$(echo ../samtools*/samtools)
  if [ ! -x $samtools ]
  then
    ( cd ../samtools* ; make samtools )
  fi
fi

test=${1}
test=${test:=10k.bam}
testout=${T}/${test##*/}-test.bam
testoutsorted=${testout}-sort.bam
testoutsam=${T}/${test##*/}-test.sam
testoutsam2=${testoutsam}2
testoutsam3=${testoutsam}3
$samtools view -h ${test} > ${testoutsam}

testsampart1=${testoutsam}.1.sam
testsampart2=${testoutsam}.2.sam
$samtools view -H ${test} > ${testsampart1}
$samtools view -H ${test} > ${testsampart2}
$samtools view ${test} | perl -ne 'if ($c++ < 10000) {print STDOUT;} else {print STDERR;}' >> ${testsampart1} 2>> ${testsampart2} 

testbampart1=${testsampart1}.bam
testbampart2=${testsampart2}.bam
$samtools view -Sb ${testsampart1} > ${testbampart1}
$samtools view -Sb ${testsampart2} > ${testbampart2}

sorted=${T}/${test##*/}-sort
$samtools sort ${test} ${sorted}
sortedsam=${sorted}.sam
sorted=${sorted}.bam
sortedorder=${sorted}.sam.order
$samtools view -h ${sorted} > ${sortedsam}
awk '{print $3" "$4}'  ${sortedsam} > ${sortedorder}

if ((isaprun==0))
then

  ./SamUtilsTest ${test} ${testout} ${testoutsorted}
  $samtools view -h $testout > ${testoutsam2}
  $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  diff -q ${testoutsam2} ${testoutsam} || failed SamUtilsTest bam copy failed
  diff -q ${testoutsam3} ${sortedorder} || failed SamUtilsTest bam sort-order failed

  ./SamUtilsTest ${testoutsam} ${testout} ${testoutsorted}
  $samtools view -h $testout > ${testoutsam2}
  $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  diff -q ${testoutsam2} ${testoutsam} || failed SamUtilsTest sam copy failed
  diff -q ${testoutsam3} ${sortedorder} || failed SamUtilsTest sam sort-order failed

  ./SamUtilsTest $sorted $testout ${testoutsorted}
  $samtools view -h $testout > ${testoutsam2}
  $samtools view -h ${testoutsorted} | awk '{print $3" "$4}'> ${testoutsam3}
  diff -q ${testoutsam2} ${sortedsam} || failed SamUtilsTest sorted-bam copy failed
  diff -q ${testoutsam3} ${sortedorder} || failed SamUtilsTest sorted-bam sort-order failed

fi

if [ -n "$MPI" ]
then
  for mpi in 1 2 3 4 6 7 8 12 13 16 20 24 32
  do
    if [ $mpi -gt $procs ]
    then
      break
    fi
    export OMP_NUM_THREADS=$(((procs+mpi-1)/mpi))
 
    $mt $MPI $mpi  ./SamUtilsTest ${test} ${testout} ${testoutsorted}
    $samtools view -h ${testout} > ${testoutsam2}
    $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
    #diff -q ${testoutsam2} ${testoutsam} || failed mpi $i SamUtilsTest bam copy failed
    diff -q ${testoutsam3} ${sortedorder} || failed mpi $i SamUtilsTest bam sort-order failed

    $mt $MPI $mpi ./SamUtilsTest ${testoutsam} ${testout} ${testoutsorted}
    $samtools view -h ${testout} > ${testoutsam2}
    $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
    #diff -q ${testoutsam2} ${testoutsam} || failed mpi $i SamUtilsTest sam copy failed
    diff -q ${testoutsam3} ${sortedorder} || failed mpi $i SamUtilsTest sam sort-order failed

    $mt $MPI $mpi  ./SamUtilsTest ${sorted} ${testout} ${testoutsorted}
    $samtools view -h $testout > ${testoutsam2}
    $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
    #diff -q ${testoutsam2} ${sortedsam} || failed mpi $i SamUtilsTest sorted-bam copy failed
    diff -q ${testoutsam3} ${sortedorder}  || failed mpi $i SamUtilsTest sorted-bam sort-order failed

    $mt $MPI $mpi  ../apps/BamSort-P ${testoutsorted} ${test}
    $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
    diff -q ${testoutsam3} ${sortedorder} || failed mpi $i BamSort-P bam sort-order failed

    $mt $MPI $mpi  ../apps/BamSort-P ${testoutsorted} ${testoutsam}
    $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
    diff -q ${testoutsam3} ${sortedorder} || failed mpi $i BamSort-P sam sort-order failed
 
    $mt $MPI $mpi  ../apps/BamSort-P ${testoutsorted} ${testbampart1} ${testbampart2}
    $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
    diff -q ${testoutsam3} ${sortedorder} || failed mpi $i BamSort-P multiple-bam sort-order failed

    $mt $MPI $mpi  ../apps/BamSort-P ${testoutsorted} ${testsampart1} ${testsampart2}
    $samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
    diff -q ${testoutsam3} ${sortedorder} || failed mpi $i BamSort-P multiple-sam sort-order failed

  done
  
fi



