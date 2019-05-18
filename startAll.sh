#!/bin/bash
echo "Hello there!"


log=logfile.log
rm ${log}
touch ${log}

for fileRun in "B2A2_SNDRCV.out" "B2A2_REDUCE.out" "B2A2_SNDRCV_Tree2Nodes.out" "B2A2_SNDRCV_Tree16Nodes.out"
do
echo ${fileRun}
for numProcs in 1 2 4 6 8 16 24 32 48 64
do
for i in {1..100}
do
    mpirun -n ${numProcs} ${fileRun} >> ${log}
done

done

done
