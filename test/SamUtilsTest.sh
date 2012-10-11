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
  echo The last test failed
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
diff -q ${testoutsam2} ${testoutsam}
diff -q ${testoutsam3} ${sortedorder}

./SamUtilsTest $sorted $testout ${testoutsorted}
samtools view -h $testout > ${testoutsam2}
samtools view -h ${testoutsorted} | awk '{print $3" "$4}'> ${testoutsam3}
diff -q ${testoutsam2} ${sortedsam}
diff -q ${testoutsam3} ${sortedorder}


for i in {1..10}
do
  mpirun -n $i ./SamUtilsTest ${test} ${testout} ${testoutsorted}
  samtools view -h ${testout} > ${testoutsam2}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  diff -q ${testoutsam2} ${testoutsam}
  diff -q ${testoutsam3} ${sortedorder}

  mpirun -n $i ./SamUtilsTest ${sorted} ${testout} ${testoutsorted}
  samtools view -h $testout > ${testoutsam2}
  samtools view -h ${testoutsorted} | awk '{print $3" "$4}' > ${testoutsam3}
  diff -q ${testoutsam2} ${sortedsam}
  diff -q ${testoutsam3} ${sortedorder}

done

trap cleanup 0


