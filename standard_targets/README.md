# Standard Targets

Benchmarking and performance comparisons among samplers on three standard target 
distributions: normal, gamma, inverse gamma.

The R scripts in this folder are numbered by their order in the workflow sequence. 
The general workflow is to schedule runs for the desired targets and samplers by round 
(round 1 is for initial testing of settings within samplers; round 2 is for full
randomized comparison at selected settings).

A user should begin with script `1_setup_trials.R`, where they can produce an 
initial schedule for a desired target, selected samplers, and round. 
The schedule is saved in the `input` folder. 
All schedules can be combined and randomized using `1_schedule_all.R`, 
which will save a complete schedule in `input`. 

The script `2_run_trials_all.R` executes the sampler for each scheduled job. It 
can also be run in a single instance to explore sampler output. The shell script 
`runallTrials.sh` can be used to run all jobs in the schedule with GNU Parallel. 
Summaries for each run are saved in the `output` folder. 

Once all jobs in a round are completed, scripts `3_combine_output.R` and 
`4_summarize_combined.R` summarize sampler performance.
