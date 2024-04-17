
## transform (Qslice)

Qtypes <- c("MSW", "MSW_samples", "AUC", "AUC_wide", "AUC_samples",
            "MM_Cauchy", "Laplace_Cauchy")

coeffs <- c(1.0, 0.0) # no water penalty

trials[["Qslice"]] <- list()

trials[["Qslice"]][["MSW"]] <- list(target = target,
                                    type = "Qslice",
                                    subtype = "MSW",
                                    n_iter = n_iter, x0 = x0,
                                    pseudo = opt_t(target = truth,
                                                   type = "function",
                                                   lb = truth$lb,
                                                   ub = truth$ub,
                                                   coeffs = coeffs,
                                                   use_meanSliceWidth = TRUE))

trials[["Qslice"]][["MSW_samples"]] <- list(target = target,
                                            type = "Qslice",
                                            subtype = "MSW_samples",
                                            n_iter = n_iter, x0 = x0,
                                            pseudo = list(pseu = list(t = "TBD",
                                                                      coeffs = coeffs)))

trials[["Qslice"]][["AUC"]] <- list(target = target,
                                    type = "Qslice",
                                    subtype = "AUC",
                                    n_iter = n_iter, x0 = x0,
                                    pseudo = opt_t(target = truth,
                                                   type = "function",
                                                   lb = truth$lb,
                                                   ub = truth$ub,
                                                   coeffs = coeffs,
                                                   use_meanSliceWidth = FALSE))

trials[["Qslice"]][["AUC_wide"]] <- list(target = target,
                                         type = "Qslice",
                                         subtype = "AUC_wide",
                                         n_iter = n_iter, x0 = x0,
                                         pseudo = list(pseu =  # for uniformity of structure
                                                         pseudo_t_list(loc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$loc,
                                                                       sc = wide_factor*trials[["Qslice"]][["AUC"]]$pseudo$pseu$sc,
                                                                       degf = trials[["Qslice"]][["AUC"]]$pseudo$pseu$degf,
                                                                       lb = truth$lb,
                                                                       ub = truth$ub)))

trials[["Qslice"]][["AUC_samples"]] <- list(target = target,
                                            type = "Qslice",
                                            subtype = "AUC_samples",
                                            n_iter = n_iter, x0 = x0,
                                            pseudo = list(pseu = list(t = "TBD",
                                                                      coeffs = coeffs)))

trials[["Qslice"]][["MM_Cauchy"]] <- list(target = target,
                                          type = "Qslice",
                                          subtype = "MM_Cauchy",
                                          n_iter = n_iter,
                                          x0 = x0,
                                          pseudo = list(pseu = list(t = "TBD")))

trials[["Qslice"]][["Laplace_Cauchy"]] <- list(target = target,
                                               type = "Qslice",
                                               subtype = "Laplace_Cauchy",
                                               n_iter = n_iter,
                                               x0 = x0,
                                               pseudo = list(pseu =  # for uniformity of structure
                                                               lapproxt(truth$d, init = x0,
                                                                        sc_adj = 1,
                                                                        lb = truth$lb,
                                                                        ub = truth$ub)))


stopifnot(all.equal(Qtypes, names(trials[["Qslice"]])))
for (xx in Qtypes) {
  trials[["Qslice"]][[xx]][["algo_descrip"]] <- trials[["Qslice"]][[xx]]$pseudo$pseu$t
}


## Independence M-H

trials[["imh"]] <- list()

trials[["imh"]][["AUC"]] <- list(target = target,
                                 type = "imh",
                                 subtype = "AUC",
                                 n_iter = n_iter, x0 = x0,
                                 pseudo = list(pseu =  # for uniformity of structure
                                                 pseudo_t_list(loc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$loc,
                                                               sc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$sc,
                                                               degf = trials[["Qslice"]][["AUC"]]$pseudo$pseu$degf,
                                                               lb = truth$lb,
                                                               ub = truth$ub)))

trials[["imh"]][["AUC_wide"]] <- list(target = target,
                                      type = "imh",
                                      subtype = "AUC_wide",
                                      n_iter = n_iter, x0 = x0,
                                      pseudo = list(pseu =  # for uniformity of structure
                                                      pseudo_t_list(loc = trials[["Qslice"]][["AUC"]]$pseudo$pseu$loc,
                                                                    sc = wide_factor*trials[["Qslice"]][["AUC"]]$pseudo$pseu$sc,
                                                                    degf = trials[["Qslice"]][["AUC"]]$pseudo$pseu$degf,
                                                                    lb = truth$lb,
                                                                    ub = truth$ub))
)

for (xx in names(trials[["imh"]])) {
  trials[["imh"]][[xx]][["algo_descrip"]] <- trials[["imh"]][[xx]]$pseudo$pseu$t
}



## create schedule
sched0[["Qslice"]] <- expand.grid(target = target,
                                  type = "Qslice", subtype = Qtypes,
                                  setting = 1,  # we're not varying the settings on Qslice
                                  rep = 1:n_rep,
                                  stringsAsFactors = FALSE)
sched0[["Qslice"]]$algo_descrip <- sapply(sched0[["Qslice"]]$subtype, function(xx) {
  trials[["Qslice"]][[xx]]$algo_descrip
}
)

head(sched0[["Qslice"]]); tail(sched0[["Qslice"]])

sched0[["imh"]] <- expand.grid(target = target,
                               type = "imh", subtype = names(trials[["imh"]]),
                               setting = 1,  # we're not varying the settings
                               rep = 1:n_rep,
                               stringsAsFactors = FALSE)
sched0[["imh"]]$algo_descrip <- sapply(sched0[["imh"]]$subtype, function(xx) {
  trials[["imh"]][[xx]]$algo_descrip
}
)

head(sched0[["imh"]]); tail(sched0[["imh"]])

