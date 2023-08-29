#!/bin/bash

export OMP_NUM_THREADS=1

for file in gprior*.R
do
  fileName="${file/.R/.txt}"
  fileClass="${file/.R/}"
  echo "${fileName}"
  for i in {1..100..1}
  do
    echo "Rscript $file $i > log/$fileClass/$fileClass$i.txt 2>&1"
  done | parallel --jobs 50 
  wait 
done
