#!/bin/bash
echo "Hi"

for var in "B2A2_REDUCE.out" "B2A2_SNDRCV_Tree2Nodes.out"
do
for i in {1..100}
do
    echo &var &i
done
done
