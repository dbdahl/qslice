#!/bin/bash

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
