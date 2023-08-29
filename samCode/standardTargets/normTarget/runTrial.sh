#!/bin/bash

export OMP_NUM_THREADS=1

# this script will take an input and return the number of rows
num_of_trials=$(Rscript ../../utilityFunctions/num_of_lines.R input/$1.rds)

echo $num_of_trials

for ((i = 1; i <= $num_of_trials; i++)); do
  #echo $i
  echo "Rscript norm_$1.R $i > log/$1/$1$i.txt 2>&1" 
done | parallel --jobs 20
wait

echo "finished"
