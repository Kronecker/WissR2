#!/bin/bash
echo "Hello there!"
make B3A7scripted

log=logfile.log
rm ${log}
touch ${log}

for fileRun in "B3scripted.out"
do
echo ${fileRun}
for i in {1..5}
do
for numProcs in {1..16}
do
    mpirun -n ${numProcs} ${fileRun} >> ${log}
done

done

done
