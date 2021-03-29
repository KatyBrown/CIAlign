#!/bin/bash
for i in {1..100};
do
    mkdir sim\_$i
    EvolveAGene4 -f human_gapdh.txt  -n 8 -b 0.62
    mv human*True* sim\_$i
    mv human_Trees* sim\_$i
    mv human*Unaligned* sim\_$i
done
