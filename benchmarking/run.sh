#!/bin/bash

echo "simulator method N time mem" > neutrality_benchmark.txt

for N in 1000 5000 10000 25000
do
    SEED=`shuf -i 0-42000000 -n 1`
    PYTHONPATH=$HOME/src/fwdpy11 /usr/bin/time -f "%e %M" -o temp.txt python3 run_fwdpy11.py $N $SEED
    RESULTS=`cat temp.txt`
    echo fwdpy11 sort" $N $RESULTS" >> neutrality_benchmark.txt
    bash slimulate.sh $N > slim.txt
    /usr/bin/time -f "%e %M" -o temp.txt ~/src/slim/build/slim slim.txt > /dev/null
    RESULTS=`cat temp.txt`
    echo slim sort" $N $RESULTS" >> neutrality_benchmark.txt
    STEPS=`echo "5*$N"|bc -l`
    /usr/bin/time -f "%e %M" -o temp.txt ../build/wfbuffered --N $N --nsteps $STEPS --seed $SEED
    RESULTS=`cat temp.txt`
    echo tskit sort" $N $RESULTS" >> neutrality_benchmark.txt
    /usr/bin/time -f "%e %M" -o temp.txt ../build/wfbuffered --N $N --nsteps $STEPS --buffer --seed $SEED
    RESULTS=`cat temp.txt`
    echo tskit buffer" $N $RESULTS" >> neutrality_benchmark.txt
done

rm -f slim.txt
rm -f temp.txt
rm -f slim.trees
rm -f fwdpy11.trees
