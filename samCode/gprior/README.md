# Description of files and directories

## R files 

**setupGprior.R** is a script that sets up each of the samplers.

**gprior_.R** are scripts that run the gprior simulation for the _ method.

**exampleGprior.R** is the initial script to run the gprior sampler.

**experiment.R** is a script that confirms each method is coming to the same result.

**derivationOfLapalce.R** is a script to confirm the derivation of the Laplace function.

**formatingFunctions.R** is a script that contains functions that the format the results from the gprior simulations.

**summary.R** is a script that produces plots for the gprior results.

## bash files

**runAll.sh** is a script that runs all the samplers in parallel.

## directories

**output** is a directory that stores the results from teh simulations.

**log** is a directory that stores the log of each trial run.

**images** is a directory that stores the stores the plots from **summary.R** as PDF files.

**explorationPlots** is a directory that stores the plots from **experiment.R**. as PDF files.

