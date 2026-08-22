# Focused sampler smoke tests used to reproduce the dynamic portion of the
# qslice package audit. These are short diagnostics, not study schedules and
# not substitutes for formal Markov-chain regression tests.

library("qslice")

qslice_ns <- asNamespace("qslice")

qfun <- function(name) get(name, envir = qslice_ns, inherits = FALSE)

slice_stepping_out <- qfun("slice_stepping_out")
slice_quantile <- qfun("slice_quantile")
slice_latent <- qfun("slice_latent")
slice_elliptical <- qfun("slice_elliptical")
slice_genelliptical <- qfun("slice_genelliptical")
slice_hyperrect <- qfun("slice_hyperrect")
slice_quantile_mv <- qfun("slice_quantile_mv")
slice_quantile_mv_seq <- qfun("slice_quantile_mv_seq")
slice_elliptical_mv <- qfun("slice_elliptical_mv")
slice_genelliptical_mv <- qfun("slice_genelliptical_mv")
imh_pseudo <- qfun("imh_pseudo")
pseudo_list <- qfun("pseudo_list")

n_iter <- 10e3
burn <- 2e3
keep <- seq.int(burn + 1L, n_iter, by = 10)

run_univariate <- function(initial, update) {
  draws <- numeric(n_iter)
  draws[1L] <- initial
  evaluations <- 0L
  for (i in 2:n_iter) {
    ans <- update(draws[i - 1L])
    draws[i] <- ans$x
    evaluations <- evaluations + ans$nEvaluations
  }
  list(
    draws = draws,
    mean = mean(draws[keep]),
    sd = sd(draws[keep]),
    evaluations_per_transition = evaluations / (n_iter - 1L)
  )
}

run_multivariate <- function(initial, update) {
  k <- length(initial)
  draws <- matrix(NA_real_, nrow = n_iter, ncol = k)
  draws[1L, ] <- initial
  evaluations <- 0L
  for (i in 2:n_iter) {
    ans <- update(draws[i - 1L, ])
    draws[i, ] <- ans$x
    evaluations <- evaluations + ans$nEvaluations
  }
  list(
    draws = draws,
    mean = colMeans(draws[keep, , drop = FALSE]),
    sd = apply(draws[keep, , drop = FALSE], 2L, sd),
    evaluations_per_transition = evaluations / (n_iter - 1L)
  )
}

rows <- list()
record <- function(name, expected, cdf_truth_list, result) {
  if (length(expected) == 1) {
    pvals <- ks.test(result$draws, cdf_truth_list[[1]])$p.value
    correlation <- 0.0
  } else if (length(expected) == 2) {
    pvals <- sapply(1:length(expected), function(k) {
      ks.test(result$draws[, k], cdf_truth_list[[k]])$p.value
    })
    correlation <- rep(
      cor(result$draws[, 1], result$draws[, 2]),
      length(expected)
    )
  }

  rows[[length(rows) + 1L]] <<- data.frame(
    sampler = name,
    component = seq_along(expected),
    expected_mean = expected,
    observed_mean = result$mean,
    observed_sd = result$sd,
    correlation = correlation,
    ks_pval = pvals,
    evaluations_per_transition = result$evaluations_per_transition,
    row.names = NULL
  )
}

# Univariate Beta(3, 4) target.
beta_log_target <- function(x) dbeta(x, shape1 = 3, shape2 = 4, log = TRUE)
beta_expected <- 3 / 7
beta_cdf <- list(function(x) pbeta(x, 3, 4))
beta_pseudo <- pseudo_list("beta", list(shape1 = 2, shape2 = 2))

set.seed(1)
record(
  "stepping-out: Beta(3,4)",
  beta_expected,
  beta_cdf,
  run_univariate(0.5, function(x) {
    slice_stepping_out(x, beta_log_target, w = 0.5, max = Inf)
  })
)

set.seed(2)
record(
  "QSlice: Beta(3,4)",
  beta_expected,
  beta_cdf,
  run_univariate(0.5, function(x) {
    slice_quantile(x, beta_log_target, beta_pseudo)
  })
)

set.seed(3)
latent_s <- 1
record(
  "latent slice: Beta(3,4)",
  beta_expected,
  beta_cdf,
  run_univariate(0.5, function(x) {
    ans <- slice_latent(x, latent_s, beta_log_target, rate = 1)
    latent_s <<- ans$s
    ans
  })
)

set.seed(4)
record(
  "GESS: Beta(3,4)",
  beta_expected,
  beta_cdf,
  run_univariate(0.5, function(x) {
    slice_genelliptical(x, beta_log_target, mu = 0.5, sigma = 1, df = 5)
  })
)

set.seed(5)
record(
  "IMH: Beta(3,4)",
  beta_expected,
  beta_cdf,
  run_univariate(0.5, function(x) imh_pseudo(x, beta_log_target, beta_pseudo))
)

# Univariate Gaussian ESS: N(0,1) prior times flat likelihood.
ess_expected_mean <- 0.0
ess_log_lik <- function(x) 1
ess_cdf <- list(function(x) pnorm(x))

set.seed(6)
record(
  "ESS: Gaussian posterior",
  ess_expected_mean,
  ess_cdf,
  run_univariate(0, function(x) {
    slice_elliptical(x, ess_log_lik, mu = 0, sigma = 1)
  })
)

# Independent bivariate beta target.
mv_beta_log_target <- function(x) {
  dbeta(x[1L], 3, 4, log = TRUE) + dbeta(x[2L], 5, 3, log = TRUE)
}
mv_beta_expected <- c(3 / 7, 5 / 8)
mv_beta_margcdf <- list(function(x) pbeta(x, 3, 4), function(x) pbeta(x, 5, 3))
mv_beta_pseudo <- list(
  pseudo_list("beta", list(shape1 = 1.5, shape2 = 1.5)),
  pseudo_list("beta", list(shape1 = 1.5, shape2 = 1.5))
)

set.seed(7)
record(
  "hyperrectangle: independent betas",
  mv_beta_expected,
  mv_beta_margcdf,
  run_multivariate(c(0.5, 0.5), function(x) {
    slice_hyperrect(x, mv_beta_log_target, w = c(0.5, 0.5))
  })
)

set.seed(8)
record(
  "multivariate QSlice: independent betas",
  mv_beta_expected,
  mv_beta_margcdf,
  run_multivariate(c(0.5, 0.5), function(x) {
    slice_quantile_mv(x, mv_beta_log_target, mv_beta_pseudo)
  })
)

set.seed(9)
record(
  "multivariate GESS: independent betas",
  mv_beta_expected,
  mv_beta_margcdf,
  run_multivariate(c(0.5, 0.5), function(x) {
    slice_genelliptical_mv(
      x,
      mv_beta_log_target,
      mu = c(0.5, 0.5),
      Sig = diag(2),
      df = 5
    )
  })
)

set.seed(10)
record(
  "multivariate IMH: independent betas",
  mv_beta_expected,
  mv_beta_margcdf,
  run_multivariate(c(0.5, 0.5), function(x) {
    imh_pseudo(x, mv_beta_log_target, mv_beta_pseudo)
  })
)

# Multivariate Gaussian ESS with an N(0,I) prior and flat likelihood.
mv_ess_expected <- c(0, 0)
mv_ess_log_lik <- function(x) 1

set.seed(11)
record(
  "multivariate ESS: Gaussian posterior",
  mv_ess_expected,
  list(function(x) pnorm(x), function(x) pnorm(x)),
  run_multivariate(c(0, 0), function(x) {
    slice_elliptical_mv(x, mv_ess_log_lik, mu = c(0, 0), Sig = diag(2))
  })
)

# Sequential QSlice with a pseudo-target exactly equal to a correlated
# bivariate standard normal target.
rho <- 0.7
target_cov <- matrix(c(1, rho, rho, 1), nrow = 2)
target_precision <- solve(target_cov)
seq_log_target <- function(x) {
  -0.5 * drop(crossprod(x, target_precision %*% x))
}
seq_control <- list(
  pseudo_init = pseudo_list("normal", list(loc = 0, sc = 1)),
  loc_fn = function(previous) rho * previous[length(previous)],
  sc_fn = function(previous) sqrt(1 - rho^2),
  lb = c(-Inf, -Inf),
  ub = c(Inf, Inf)
)

set.seed(12)
record(
  "sequential QSlice: exact Gaussian pseudo-target",
  c(0, 0),
  list(function(x) pnorm(x), function(x) pnorm(x)),
  run_multivariate(c(0, 0), function(x) {
    slice_quantile_mv_seq(x, seq_log_target, seq_control)
  })
)

results <- do.call(rbind, rows)
print(results, row.names = FALSE, digits = 6)

invisible(results)
