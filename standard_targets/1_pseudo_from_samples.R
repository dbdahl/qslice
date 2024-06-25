n_samples <- 1e3

samples_train <- truth$q(runif(n_samples))

if (isTRUE(run_info$subtype == "MSW_samples")) {

  new_pseu <- pseudo_opt(samples = samples_train,
                         type = "samples",
                         family = "t",
                         degf = ifelse(target == "igamma",
                                  c(1, 5), c(1, 5, 20)),
                         lb = truth$lb,
                         ub = truth$ub,
                         utility_type = "MSW",
                         plot = FALSE)

  trials[["Qslice"]][["MSW_samples"]][["pseudo"]] <- new_pseu

} else if (isTRUE(run_info$subtype == "AUC_samples")) {

  new_pseu <- pseudo_opt(samples = samples_train,
                         type = "samples",
                         family = "t",
                         degf = ifelse(target == "igamma",
                                       c(1, 5), c(1, 5, 20)),
                         lb = truth$lb,
                         ub = truth$ub,
                         utility_type = "AUC",
                         plot = FALSE)

  trials[["Qslice"]][["AUC_samples"]][["pseudo"]] <- new_pseu

} else if (isTRUE(run_info$subtype == "MM_Cauchy")) {

  new_pseu <- list(pseu =  # for uniformity of structure
                     pseudo_list(family = "cauchy",
                                 params = list(loc = mean(samples_train),
                                               sc = sd(samples_train)),
                                 lb = truth$lb,
                                 ub = truth$ub))

  trials[["Qslice"]][["MM_Cauchy"]][["pseudo"]] <- new_pseu

}

trials[[run_info$type]][[run_info$subtype]]$algo_descrip <- trials[[run_info$type]][[run_info$subtype]]$pseudo$pseu$txt

rm(n_samples, samples_train, new_pseu)
