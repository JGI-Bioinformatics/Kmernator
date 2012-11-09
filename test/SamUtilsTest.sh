#!/bin/bash

set -e
set -x

T=$(mktemp -d)
cleanup() 
{
  echo "Cleaning up" 1>&2
  rm -rf ${T}
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

test=${1}
test=${test:=10k.bam}
testout=${T}/${test##*/}-test.bam
testoutsorted=${testout}-sort.bam
testoutsam=${T}/${test##*/}-test.sam
testoutsam2=${testoutsam}2
testoutsam3=${testoutsam}3
samtools view -h ${test} > ${testoutsam}

testsampart1=${testoutsam}.1.sam
testsampart2=${testoutsam}.2.sam
samtools view -H ${test} > ${testsampart1}
samtools view -H ${test} > ${testsampart2}
samtools view ${test} | perl -ne 'if ($c++ < 10000) {print STDOUT;} else {print STDERR;}' >> ${testsampart1} 2>> ${testsampart2} 

testbampart1=${testsampart1}.bam
testbampart2=${testsampart2}.bam
samtools view -Sb ${testsampart1} > ${testbampart1}
samtools view -Sb ${testsampart2} > ${testbampart2}

sorted=${T}/${test##*/}-sort
samtools sort ${test} ${sorted}
sortedsam=${sorted}.sam
sorted=${sorted}.bam
sortedorder=${sorted}.sam.order
samtools view -h ${sorted} > ${sortedsam}
awk '{print $3" "$4}'  ${sortedsam} > ${sortedorder}


./SamUtilsTest ${test} ${testout} ${testoutsorted}
samtools view -h $testout > ${testoutsam2}
samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
diff -q ${testoutsam2} ${testoutsam} || failed SamUtilsTest bam copy failed
diff -q ${testoutsam3} ${sortedorder} || failed SamUtilsTest bam sort-order failed

./SamUtilsTest ${testoutsam} ${testout} ${testoutsorted}
samtools view -h $testout > ${testoutsam2}
samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
diff -q ${testoutsam2} ${testoutsam} || failed SamUtilsTest sam copy failed
diff -q ${testoutsam3} ${sortedorder} || failed SamUtilsTest sam sort-order failed

./SamUtilsTest $sorted $testout ${testoutsorted}
samtools view -h $testout > ${testoutsam2}
samtools view -h ${testoutsorted} | awk '{print $3" "$4}'> ${testoutsam3}
diff -q ${testoutsam2} ${sortedsam} || failed SamUtilsTest sorted-bam copy failed
diff -q ${testoutsam3} ${sortedorder} || failed SamUtilsTest sorted-bam sort-order failed


for i in {1..10}
do
  mpirun -n $i ./SamUtilsTest ${test} ${testout} ${testoutsorted}
  samtools view -h ${testout} > ${testoutsam2}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  #diff -q ${testoutsam2} ${testoutsam} || failed mpi $i SamUtilsTest bam copy failed
  diff -q ${testoutsam3} ${sortedorder} || failed mpi $i SamUtilsTest bam sort-order failed

  mpirun -n $i ./SamUtilsTest ${testoutsam} ${testout} ${testoutsorted}
  samtools view -h ${testout} > ${testoutsam2}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  #diff -q ${testoutsam2} ${testoutsam} || failed mpi $i SamUtilsTest sam copy failed
  diff -q ${testoutsam3} ${sortedorder} || failed mpi $i SamUtilsTest sam sort-order failed

  mpirun -n $i ./SamUtilsTest ${sorted} ${testout} ${testoutsorted}
  samtools view -h $testout > ${testoutsam2}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  #diff -q ${testoutsam2} ${sortedsam} || failed mpi $i SamUtilsTest sorted-bam copy failed
  diff -q ${testoutsam3} ${sortedorder}  || failed mpi $i SamUtilsTest sorted-bam sort-order failed

  mpirun -n $i ../apps/BamSort-P ${testoutsorted} ${test}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  diff -q ${testoutsam3} ${sortedorder} || failed mpi $i BamSort-P bam sort-order failed

  mpirun -n $i ../apps/BamSort-P ${testoutsorted} ${testoutsam}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  diff -q ${testoutsam3} ${sortedorder} || failed mpi $i BamSort-P sam sort-order failed

  mpirun -n $i ../apps/BamSort-P ${testoutsorted} ${testbampart1} ${testbampart2}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  diff -q ${testoutsam3} ${sortedorder} || failed mpi $i BamSort-P multiple-bam sort-order failed

  mpirun -n $i ../apps/BamSort-P ${testoutsorted} ${testsampart1} ${testsampart2}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  diff -q ${testoutsam3} ${sortedorder} || failed mpi $i BamSort-P multiple-sam sort-order failed


done

trap cleanup 0


