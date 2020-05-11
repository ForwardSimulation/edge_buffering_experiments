#!/bin/bash

echo "N method tsimplify time mem" > $(pwd)/benchmarks.txt

for N in 1000 5000 10000 25000
do
    runtime=`echo "5*$N"|bc -l`
    for tsimplify in 100 500 1000
    do
        SEED=$RANDOM
        /usr/bin/time -f "%e %M" -o classic.time ./wfbuffered --treefile classic.trees --N $N --simplify $tsimplify --seed $SEED --nsteps $runtime
        /usr/bin/time -f "%e %M" -o buffered.time ./wfbuffered --treefile buffered.trees --N $N --simplify $tsimplify --buffer --seed $SEED --nsteps $runtime
        python3 ../compare_treefiles.py $(pwd)/classic.trees $(pwd)/buffered.trees
        c=`cat classic.time`
        b=`cat buffered.time`
        echo $N "sort" $tsimplify $c >> benchmarks.txt
        echo $N "buffer" $tsimplify $b >> benchmarks.txt
    done
done
