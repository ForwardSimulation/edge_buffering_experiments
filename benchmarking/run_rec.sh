#!/bin/bash

# echo "simulator method N rho time mem" > neutrality_benchmark_rho.txt

# for N in 1000 
for N in 10000
do
    for rho in 100 1000 10000 100000
    do
    SEED=`shuf -i 0-42000000 -n 1`
    STEPS=`echo "5*$N"|bc -l`
    /usr/bin/time -f "%e %M" -o temp.txt ../build/wfbuffered --N $N --nsteps $STEPS --seed $SEED --rho $rho
    python3 numtrees.py
    RESULTS=`cat temp.txt`
    echo tskit tsk_sort" $N $rho $RESULTS" >> neutrality_benchmark_rho.txt
    /usr/bin/time -f "%e %M" -o temp.txt ../build/wfbuffered --N $N --nsteps $STEPS --seed $SEED --cppsort --rho $rho
    python3 numtrees.py
    RESULTS=`cat temp.txt`
    echo tskit cppsort" $N $rho $RESULTS" >> neutrality_benchmark_rho.txt
    /usr/bin/time -f "%e %M" -o temp.txt ../build/wfbuffered --N $N --nsteps $STEPS --seed $SEED --cppsort --parallel_sort --rho $rho
    python3 numtrees.py
    RESULTS=`cat temp.txt`
    echo tskit cppsort_par" $N $rho $RESULTS" >> neutrality_benchmark_rho.txt
    /usr/bin/time -f "%e %M" -o temp.txt ../build/wfbuffered --N $N --nsteps $STEPS --buffer --seed $SEED --rho $rho
    python3 numtrees.py
    RESULTS=`cat temp.txt`
    echo tskit buffer" $N $rho $RESULTS" >> neutrality_benchmark_rho.txt
    done
done

rm -f temp.txt
rm -f treefile.trees
