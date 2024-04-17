n_samples <- 1e3

samples_train <- truth$q(runif(n_samples))

if (isTRUE(run_info$subtype == "MSW_samples")) {
  
  new_pseu <- opt_t(samples = samples_train,
                    type = "samples",
                    lb = truth$lb,
                    ub = truth$ub,
                    degf = ifelse(target == "igamma",
                                  c(1, 5), c(1, 5, 20)),
                    coeffs = trials[["Qslice"]][["MSW_samples"]]$pseudo$pseu$coeffs,
                    use_meanSliceWidth = TRUE,
                    plot = FALSE)
  
  trials[["Qslice"]][["MSW_samples"]][["pseudo"]] <- new_pseu

} else if (isTRUE(run_info$subtype == "AUC_samples")) {
  
  new_pseu <- opt_t(samples = samples_train,
                    type = "samples",
                    lb = truth$lb,
                    ub = truth$ub,
                    degf = ifelse(target == "igamma",
                                  c(1, 5), c(1, 5, 20)),
                    coeffs = trials[["Qslice"]][["AUC_samples"]]$pseudo$pseu$coeffs,
                    use_meanSliceWidth = FALSE,
                    plot = FALSE)
  
  trials[["Qslice"]][["AUC_samples"]][["pseudo"]] <- new_pseu
  
} else if (isTRUE(run_info$subtype == "MM_Cauchy")) {
  
  new_pseu <- list(pseu =  # for uniformity of structure
                     pseudo_t_list(loc = mean(samples_train),
                                   sc = sd(samples_train),
                                   degf = 1.0,
                                   lb = truth$lb,
                                   ub = truth$ub))
  
  trials[["Qslice"]][["MM_Cauchy"]][["pseudo"]] <- new_pseu
  
}

trials[[run_info$type]][[run_info$subtype]]$algo_descrip <- trials[[run_info$type]][[run_info$subtype]]$pseudo$pseu$t

rm(n_samples, samples_train, new_pseu)