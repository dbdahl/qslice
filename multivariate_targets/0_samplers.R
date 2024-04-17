gess_sampler <- function(lt, n_iter, x0, loc, ScL, degf, is_chol) {
  draws <- matrix(rep(x0, n_iter + 1), nrow = n_iter + 1, byrow = TRUE)
  n_eval <- 0
  for (i in 2:(n_iter + 1)) {
    tmp <- slice_genelliptical_mv(draws[i-1,],
                                  target = lt,
                                  mu = loc, Sig = ScL, df = degf,
                                  is_chol = is_chol)
    draws[i,] <- tmp$x
    n_eval <- n_eval + tmp$nEvaluations
  }
  list(draws = draws[-1,], n_eval = n_eval)
}

Qslice_sampler <- function(lt, n_iter, x0, ps_lpdf, ps_icdf, ps_cdf) {
  draws <- matrix(rep(x0, n_iter + 1), nrow = n_iter + 1, byrow = TRUE)
  draws_u <- matrix(rep(x0, n_iter), nrow = n_iter, byrow = TRUE)
  n_eval <- 0
  for (i in 2:(n_iter + 1)) {
    tmp <- slice_mv_transform(draws[i-1,],
                              target = lt,
                              pseudo_log_pdf = ps_lpdf,
                              pseudo_inv_cdf = ps_icdf,
                              pseudo_cdf = ps_cdf)
    draws[i,] <- tmp$x
    draws_u[i-1,] <- tmp$u
    n_eval <- n_eval + tmp$nEvaluations
  }
  list(draws = draws[-1,], draws_u = draws_u, n_eval = n_eval)
}


Qslice_seq_sampler <- function(lt, n_iter, x0, pseudo_control) {
  draws <- matrix(rep(x0, n_iter + 1), nrow = n_iter + 1, byrow = TRUE)
  draws_u <- matrix(rep(x0, n_iter), nrow = n_iter, byrow = TRUE)
  n_eval <- 0
  for (i in 2:(n_iter + 1)) {
    tmp <- slice_mv_transform_seq(draws[i-1,],
                                  target = lt,
                                  pseudo_control = pseudo_control)
    draws[i,] <- tmp$x
    draws_u[i-1,] <- tmp$u
    n_eval <- n_eval + tmp$nEvaluations
  }
  list(draws = draws[-1,], draws_u = draws_u, n_eval = n_eval)
}

