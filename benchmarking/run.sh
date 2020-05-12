#!/bin/bash

for N in 1000 5000 10000
do
    echo N = $N
    echo fwdpy11
    SEED=`shuf -i 0-42000000 -n 1`
    PYTHONPATH=$HOME/src/fwdpy11 /usr/bin/time -f "%e %M" python3 run_fwdpy11.py $N $SEED
    echo slim
    bash slimulate.sh $N > slim.txt
    /usr/bin/time -f "%e %M" ~/src/slim/build/slim slim.txt > /dev/null
    STEPS=`echo "5*$N"|bc -l`
    echo sorting
    /usr/bin/time -f "%e %M" ../build/wfbuffered --N $N --nsteps $STEPS --seed $SEED
    echo buffering
    /usr/bin/time -f "%e %M" ../build/wfbuffered --N $N --nsteps $STEPS --buffer --seed $SEED
done

rm -f slim.txt
