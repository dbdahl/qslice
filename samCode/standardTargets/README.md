# Content of each directory

Each of the Target directories contains similar files that have been adapted to the selected target. Each directory contains eight R files two bash scripts and four directories.

## R files

**_setup.R** is a script that needs to be run that specifies the target distribution.

**_trials.R** is a script where a user can specify the parameter values for each sampler and the number of times that want each set of parameter values run.

**_gess.R** is a script that runs each set of parameter values for the generalized elliptical slice sampler.

**_latent.R**  is a script that runs each set of parameter values for the latent slice sampler.

 **_randwalk.R** is a script that runs each set of parameter values for the Metropolis random walk sampler.

 **_steppingout.R** is a script that runs each set of parameter values for the stepping out slice sampler.

 **_transform.R** is a script that runs each set of parameter values for the transform slice sampler.

 **summary.R** is a script that creates plots for each sampler. 

 ## bash files

 **runTrial.sh** is a bash script that runs all the trials for a specific sampler. It is run using this command
 - bash runTrial.sh (sampling method)
 - eg: bash runTrial.sh steppingout.

 **runAll.sh** is a bash script that runs all sampling method trials.

 ## directories

 **input** is a directory that stores all the trials to be run. Is populated by running the **_trials.R** script

 **output** is a directory that stores the results from each trial in a separate CSV file. 

 **images** is a directory that stores all the plots from **summary.R** as PDF files.

 **log** is a directory that has all the logs for each trial run. 


The **images** directory found inside the **standardTargets** directory contains scripts used to summarize the information from each of the three targets.

