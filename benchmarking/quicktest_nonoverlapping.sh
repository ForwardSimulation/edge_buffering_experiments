#!/bin/bash

time ./wfbuffered --treefile classic.trees --psurvival 0.5
time ./wfbuffered --buffer --treefile buffered.trees --psurvival 0.5
python3 ../compare_treefiles.py $(pwd)/classic.trees $(pwd)/buffered.trees

