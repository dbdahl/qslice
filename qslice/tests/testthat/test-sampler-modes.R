run_chain <- function(x0, update, diagnostics, n_iter = 30L, seed = 1L) {
  set.seed(seed)
  x <- x0
  draws <- matrix(NA_real_, nrow = n_iter, ncol = length(x0))
  n_eval <- 0L
  for (i in seq_len(n_iter)) {
    out <- update(x, diagnostics)
    x <- out$x
    draws[i, ] <- x
    n_eval <- n_eval + out$nEvaluations
  }
  list(draws = draws, n_eval = n_eval)
}

expect_modes_identical <- function(x0, update) {
  fast <- run_chain(x0, update, diagnostics = FALSE)
  diag <- run_chain(x0, update, diagnostics = TRUE)
  expect_identical(fast$draws, diag$draws)
  expect_identical(fast$n_eval, diag$n_eval)
}

lt_gamma <- function(x) ifelse(x > 0.0, 1.5 * log(x) - x, -Inf)
lt_beta2 <- function(x) {
  dbeta(x[1], 3, 4, log = TRUE) + dbeta(x[2], 5, 3, log = TRUE)
}
pseudo_t <- pseudo_list("t", list(loc = 2, sc = 2, df = 5), lb = 0)
pseudo_beta2 <- list(
  pseudo_list("beta", list(shape1 = 2, shape2 = 2)),
  pseudo_list("beta", list(shape1 = 3, shape2 = 2))
)
Sig <- matrix(c(0.5, 0.25, 0.25, 0.5), nrow = 2)

test_that("univariate samplers give identical output in fast and diagnostic modes", {
  expect_modes_identical(1, function(x, d) {
    slice_stepping_out(x, lt_gamma, w = 1, diagnostics = d)
  })
  expect_modes_identical(1, function(x, d) {
    slice_stepping_out(x, lt_gamma, w = 0.3, max = 5, diagnostics = d)
  })
  expect_modes_identical(1, function(x, d) {
    slice_quantile(x, lt_gamma, pseudo_t, diagnostics = d)
  })
  latent_update <- function() {
    s <- 1
    function(x, d) {
      out <- slice_latent(x, s, lt_gamma, rate = 0.5, diagnostics = d)
      s <<- out$s
      out
    }
  }
  expect_identical(
    run_chain(1, latent_update(), diagnostics = FALSE),
    run_chain(1, latent_update(), diagnostics = TRUE)
  )
  expect_modes_identical(0.5, function(x, d) {
    slice_elliptical(x, function(z) dbeta(z, 8, 2, log = TRUE),
                     mu = 0.2, sigma = 0.3, diagnostics = d)
  })
  expect_modes_identical(1, function(x, d) {
    slice_genelliptical(x, lt_gamma, mu = 1.5, sigma = 1.5, df = 5,
                        diagnostics = d)
  })
  expect_modes_identical(1, function(x, d) {
    imh_pseudo(x, lt_gamma, pseudo_t, diagnostics = d)
  })
  expect_modes_identical(1, function(x, d) {
    rwm_norm(x, lt_gamma, cand_sd = 1.5, support = c(0, Inf), diagnostics = d)
  })
})

test_that("multivariate samplers give identical output in fast and diagnostic modes", {
  expect_modes_identical(c(0.5, 0.5), function(x, d) {
    slice_hyperrect(x, lt_beta2, w = c(0.5, 0.5), diagnostics = d)
  })
  expect_modes_identical(c(0.5, 0.5), function(x, d) {
    slice_hyperrect(x, lt_beta2, L = c(0, 0), R = c(1, 1), diagnostics = d)
  })
  expect_modes_identical(c(0.5, 0.5), function(x, d) {
    slice_quantile_mv(x, lt_beta2, pseudo_beta2, diagnostics = d)
  })
  expect_modes_identical(c(0.5, 0.5), function(x, d) {
    imh_pseudo(x, lt_beta2, pseudo_beta2, diagnostics = d)
  })
  expect_modes_identical(c(0.3, 0.3), function(x, d) {
    slice_elliptical_mv(x, lt_beta2, mu = c(0.5, 0.5), Sig = Sig,
                        diagnostics = d)
  })
  expect_modes_identical(c(0.3, 0.3), function(x, d) {
    slice_genelliptical_mv(x, lt_beta2, mu = c(0.5, 0.5), Sig = Sig, df = 5,
                           diagnostics = d)
  })

  rho <- 0.7
  precision <- solve(matrix(c(1, rho, rho, 1), nrow = 2))
  control <- list(
    pseudo_init = pseudo_list("t", list(loc = 0, sc = 1.2, df = 5)),
    loc_fn = function(previous) rho * previous[length(previous)],
    sc_fn = function(previous) 1.2 * sqrt(1 - rho^2),
    lb = c(-Inf, -Inf),
    ub = c(Inf, Inf)
  )
  expect_modes_identical(c(0, 0), function(x, d) {
    slice_quantile_mv_seq(x, function(z) -0.5 * drop(crossprod(z, precision %*% z)),
                          control, diagnostics = d)
  })
})

test_that("the qslice.diagnostics option selects diagnostic mode", {
  old <- options(qslice.diagnostics = TRUE)
  on.exit(options(old))
  expect_error(
    slice_stepping_out(1, function(x) NaN, w = 1),
    class = "qslice_target_error"
  )
})

test_that("deprecated auc() accepts a function for y", {
  value <- suppressWarnings(auc(y = function(x) dbeta(x, 2, 2)))
  expect_equal(value, 1 / 1.5, tolerance = 1e-4)
})

test_that("deprecated auc() passes nbins through", {
  set.seed(1)
  u <- rbeta(500, 2, 2)
  for (nbins in c(10, 30, 50)) {
    expect_identical(
      suppressWarnings(auc(u = u, nbins = nbins)),
      utility_shrinkslice(u = u, nbins = nbins, utility_type = "AUC")
    )
  }
})
