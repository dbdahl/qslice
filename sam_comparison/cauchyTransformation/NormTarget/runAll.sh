#!/bin/bash
for file in norm_[glsrt]*.R
do
	fileName="${file/.R/.txt}"
	Rscript $file > log/$fileName 2>&1 &
done
