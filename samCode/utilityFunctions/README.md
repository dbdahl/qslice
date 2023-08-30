## Description of what functions are found in each script

**combine_best_func** is a function that takes the results output from each of the sampling methods and puts them in the same format. In addition to this, it also filters so that only the best set of parameters is kept for each method. 

**combine_func** is a function similar to **combine_best_func** but it returns every set of parameter values for the transform slice sampler.

**comparison_functions** contains functions for each sampling method where the inputs are the target function and and number of samples needed. The output for this function includes among other things the number of effective samples and the time it took to get those samples.

**effect_of_transform** is a function that takes a target and a pseudo target and returns a plot showing what the transformation will look like.

**fit_trunc_Cauchy** is a function that selects a pseudo target using samples by finding MLE estimates for the location and scale of a Cauchy distribution. 

**fn_AUC_hist** is a function that takes samples and returns the area under the curve for those samples using a histogram as an approximation. 

**fn_expectedSliceWidth** is a function that finds the expected slice width of a function.

**fn_invgamma** contains functions for the pdf, cdf, and inverse cdf of an inverse gamma distribution.

**fn_shrinkslice_utility** contains functions that find the utility metric (AUC - water). It also has functions that take a location, scale, and degrees of freedom and return the pdf, log pdf, and inverse cdf needed for a t distribution.

**fn_skew_hist** is a function that finds the skewness of a distribution using samples.

**fn_waterArea_hist** is a function that finds the amount of water a function will hold using a histogram as an approximation. 

**fn_waterArea_int** is a function that finds the amount of water a function will hold using numerical integration. 

**format_df** is a function that takes in a data frame of results returned from any of the sampling methods and returns a data frame with three columns: number of evaluations, effective samples per CPU second, and the parameter values of the method. 

**gprior_functions** is a script that contains functions that are used to run the gprior simulation.

**lapproxt** is a function that takes in a target and returns the Laplace approximation.

**num_of_lines** is a function that returns the number of lines in .rds file

**read_in_data** is a function that takes multiple CSV files and combines them into one data frame.

