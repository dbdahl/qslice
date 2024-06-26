% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hyperrectangle_slice.R
\name{slice_quantile_mv_seq}
\alias{slice_quantile_mv_seq}
\title{Multivariate Quantile Slice Sampler from a sequence of conditional pseudo-targets}
\usage{
slice_quantile_mv_seq(x, log_target, pseudo_control)
}
\arguments{
\item{x}{The current state (as a numeric vector).}

\item{log_target}{A function taking numeric vector that evaluates the log-target
density, returning a numeric scalar.}

\item{pseudo_control}{A list with

\code{pseudo_init}, a list output from
\link[qslice]{pseudo_list} describing the marginal pseudo-target for \code{x[1]}.
Attributes of \code{pseudo_init} will be used in subsequent pseudo-targets,
except for location and scale parameters.

\code{loc_fn}, a function that specifies the location of a conditional
pseudo-target given the elements in \code{x} that precede it.

\code{sc_fn}, a function that specifies the scale of a conditional
pseudo-target given the elements in \code{x} that precede it.

\code{lb}, a numeric vector (same length as \code{x}) specifying the lower
bound of support for each conditional pseudo-target.

\code{ub}, a numeric vector (same length as \code{x}) specifying the upper
bound of support for each conditional pseudo-target.}
}
\value{
A list containing three elements: "x" is the new state, "u" is a vector
of values of the sequence of conditional CDFs of the psuedo-targets associated
with the returned value, and "nEvaluations is the number of evaluations of the
target function used to obtain the new state.
}
\description{
Quantile slice sampler for a random vector (Heiner et al., 2024+). The pseudo-target is specified as
a sequence of growing conditional distributions.
}
\examples{
# Funnel distribution from Neal (2003).
K <- 10
n_iter <- 50 # MCMC iterations; set to 10e3 for more complete illustration
n <- 100 # number of iid samples from the target; set to 10e3 for more complete illustration
Y <- matrix(NA, nrow = n, ncol = K) # iid samples from the target
Y[,1] <- rnorm(n, 0.0, 3.0)
for (i in 1:n) {
  Y[i, 2:K] <- rnorm(K-1, 0.0, exp(0.5*Y[i,1]))
}
ltarget <- function(x) {
dnorm(x[1], 0.0, 3.0, log = TRUE) +
  sum(dnorm(x[2:K], 0.0, exp(0.5*x[1]), log = TRUE))
}
pseudo_control <- list(
  loc_fn = function(x) {
    0.0
  },
  sc_fn = function(x) {
    if (is.null(x)) {
      out <- 3.0
    } else {
      out <- exp(0.5*x[1])
    }
    out
  },
  pseudo_init = pseudo_list(family = "t",
                            params = list(loc = 0.0, sc = 3.0, degf = 20),
                            lb = -Inf, ub = Inf),
  lb = rep(-Inf, K),
  ub = rep(Inf, K)
)
x0 <- runif(K)
draws <- matrix(rep(x0, n_iter + 1), nrow = n_iter + 1, byrow = TRUE)
draws_u <- matrix(rep(x0, n_iter), nrow = n_iter, byrow = TRUE)
n_eval <- 0
for (i in 2:(n_iter + 1)) {
  tmp <- slice_quantile_mv_seq(draws[i-1,],
                                log_target = ltarget,
                                pseudo_control = pseudo_control)
  draws[i,] <- tmp$x
  draws_u[i-1,] <- tmp$u
  n_eval <- n_eval + tmp$nEvaluations
}
# (es <- coda::effectiveSize(coda::as.mcmc(draws)))
# mean(es)
n_eval / n_iter
sapply(1:K, function (k) auc(u = draws_u[,k]))
hist(draws_u[,1])
plot(draws[,1], draws[,2])
points(Y[,1], Y[,2], col = "blue", cex = 0.5)
}
\references{
Heiner, M. J., Johnson, S. B., Christensen, J. R., and Dahl, D. B. (2024+), "Quantile Slice Sampling," \emph{arXiv preprint arXiv:###}
}
