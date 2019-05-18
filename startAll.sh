#!/bin/bash
echo "Hi"
numProcs=16
runs=10
rm logfile.log
touch logfile.log

for fileRun in "B2A2_SNDRCV.out" "B2A2_REDUCE.out" "B2A2_SNDRCV_Tree2Nodes.out" "B2A2_SNDRCV_Tree16Nodes.out"
do

for numProcs in 1 2 4 6 8 16 24 32 48 64
do
for i in {1..runs}
do
    mpirun -n $numProcs $fileRun >> logfile.log
done

done

done
