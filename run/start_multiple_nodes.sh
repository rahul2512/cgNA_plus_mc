#!/bin/bash 

mkdir $2 

for i in `seq 1 $1`; do 
    sbatch Parallel.run $i $2
done 


