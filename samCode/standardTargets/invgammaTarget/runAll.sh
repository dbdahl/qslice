#!/bin/bash
#for file in norm_[glrst]*.R
#do
#  fileName="${file/.R/.txt}"
#  echo $file
#  Rscript $file > log/$fileName 2>&1 &
#done

bash runTrial.sh gess
wait
bash runTrial.sh steppingout
wait
bash runTrial.sh latent
wait
bash runTrial.sh randwalk
wait
bash runTrial.sh transform
echo 'Finished'
