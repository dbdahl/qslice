#!/bin/bash
Rscript gamma_gess.R > log/logGess.txt 2>&1 &
Rscript gamma_latent.R > log/logLatent.txt 2>&1 &
Rscript gamma_rand_walk.R > log/logRandWalk.txt 2>&1 &
Rscript gamma_stepping_out.R > log/logSteppingOut.txt 2>&1 &
Rscript gamma_transform.R > log/logTransform.txt 2>&1 &
