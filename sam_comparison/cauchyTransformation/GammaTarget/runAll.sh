#!/bin/bash
for file in gamma_[glrst]*.R
do
    fileName="${file/.R/.txt}"
    Rscript $file > log/$fileName 2>&1 &
done

