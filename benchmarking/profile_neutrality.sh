#!/usr/bin/bash

echo "simulator method N percent function" > neutrality_time_sorting.txt

for N in 1000 5000 10000 25000
do
    PYTHONPATH=$HOME/src/fwdpy11 LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=fwdpy11.prof python3 run_fwdpy11.py $N $RANDOM
    google-pprof --text `which python3` fwdpy11.prof | grep introsort | head -n 1 | sed 's/ \+ /\t/g' | cut -f 6,7 | sed 's/%//g' > temp.txt
    x=`cat temp.txt`
    echo fwpdy11 cppsort $N $x >> neutrality_time_sorting.txt
    google-pprof --pdf `which python3` fwdpy11.prof > fwdpy11_neutrality_$N.pdf
    rm -f fwdpy11.prof

    STEPS=`echo "5*$N"|bc -l`
    LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=tskit.prof ../build/wfbuffered --N $N --nsteps $STEPS --seed $RANDOM
    google-pprof --text ../build/wfbuffered tskit.prof | grep qsort | head -n 1 | sed 's/ \+ /\t/g' | cut -f 6,7 | sed 's/%//g' > temp.txt
    x=`cat temp.txt`
    echo tskit tsk_sort $N $x >> neutrality_time_sorting.txt
    google-pprof --pdf ../build/wfbuffered tskit.prof > tskit_tsk_sort_neutrality_$N.pdf
    rm -f tskit.prof

    STEPS=`echo "5*$N"|bc -l`
    LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libprofiler.so.0 CPUPROFILE=tskit.prof ../build/wfbuffered --N $N --nsteps $STEPS --seed $RANDOM --cppsort
    google-pprof --text ../build/wfbuffered tskit.prof | grep introsort | head -n 1 | sed 's/ \+ /\t/g' | cut -f 6,7 | sed 's/%//g' > temp.txt
    x=`cat temp.txt`
    echo tskit cppsort $N $x >> neutrality_time_sorting.txt
    google-pprof --pdf ../build/wfbuffered tskit.prof > tskit_cppsort_neutrality_$N.pdf
    rm -f tskit.prof

done
