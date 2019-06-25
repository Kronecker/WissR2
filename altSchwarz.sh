#!/bin/bash
echo "Hello there!"
make B5A10

log=slice.dat
rm ${log}
touch ${log}

echo ${fileRun}

for Nx in {51,52,53,60,90}
do
    ./B5_altSchwarz.out ${Nx} 50
done

gnuplot -p -e "plot for [IDX=0:4] 'slice.dat' i IDX u 1:2 w lines title columnheader
(1)"