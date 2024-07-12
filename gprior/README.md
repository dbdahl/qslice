---
editor_options: 
  markdown: 
    wrap: 72
---

# gprior Folder

Benchmarking and performance comparisons among samplers in a Bayesian
linear regression.

The R scripts in this folder are numbered by their order in the workflow
sequence. The general workflow is to schedule runs for the desired
targets and samplers. In contrast with the simulations in
`standard_targets`, adaptive tuning occurs prior to each timed run.

The script `1_schedule.R` creates a randomized simulation schedule that
is carried out by the shell script `runallTrials.sh`. Each run calls
`2_run_trials_all.R`, which outputs a summary to the `output`
folder.

Once all jobs in a round are completed, scripts `3_combine_output.R` and
`4_summarize_combined.R` summarize sampler performance.
