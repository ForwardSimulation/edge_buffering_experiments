#!/bin/bash

N=$1
STEPS=`echo "5*$N"|bc -l`

cat <<HERE
initialize() {
        initializeTreeSeq(simplificationInterval=100);
        initializeMutationType("m1", 0.5, "f", 0.0);

        initializeGenomicElementType("g1", m1, 1.0);
        initializeGenomicElement(g1, 0, 600);
        initializeMutationRate(0);
        initializeRecombinationRate(0);
}
1 early() {
    sim.addSubpop("p1", $N);
}
$STEPS late() {
    sim.treeSeqOutput("slim.trees");
}
HERE
