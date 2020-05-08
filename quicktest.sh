#!/bin/bash

time ./wfbuffered --treefile classic.trees
time ./wfbuffered --buffer --treefile buffered.trees
python3 ../compare_treefiles.py $(pwd)/classic.trees $(pwd)/buffered.trees
