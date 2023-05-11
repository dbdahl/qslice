rm(list=ls())

library("cucumber")
library("coda")

sessionInfo()

lf <- function(x) {
  dunif(x, log=TRUE)
}

nsim <- 1e5

### Stepping out

set.seed(1)

counter_st <- 0
draws_st <- numeric(nsim) + 0.5
time_st <- system.time({
  for ( i in seq.int(2,length(draws_st)) ) {
    out <- slice_sampler_stepping_out(x = draws_st[i-1],
                                      target = lf,
                                      w=2.0,
                                      log = TRUE)
    draws_st[i] <- out$x
    counter_st <- counter_st + out$nEvaluations
  }
})
counter_st
coda::effectiveSize(draws_st)
coda::effectiveSize(draws_st) / time_st['user.self']

quit(save="no")
