free_locscale <- function(family, params) {
  loc <- params$loc
  sc <- params$sc

  if (family == "t") {
    dens_std <- function(x) dt(x, df = params$df)
    ldens_std <- function(x) dt(x, df = params$df, log = TRUE)
    cdf_std <- function(x) pt(x, df = params$df)
    invcdf_std <- function(u) qt(u, df = params$df)
  } else if (family == "cauchy") {
    dens_std <- dcauchy
    ldens_std <- function(x) dcauchy(x, log = TRUE)
    cdf_std <- pcauchy
    invcdf_std <- qcauchy
  } else if (family == "normal") {
    dens_std <- dnorm
    ldens_std <- function(x) dnorm(x, log = TRUE)
    cdf_std <- pnorm
    invcdf_std <- qnorm
  } else if (family == "logistic") {
    dens_std <- dlogis
    ldens_std <- function(x) dlogis(x, log = TRUE)
    cdf_std <- plogis
    invcdf_std <- qlogis
  } else {
    stop("Unsupported location-scale family in test helper.", call. = FALSE)
  }

  list(
    d = function(x) dens_std((x - loc) / sc) / sc,
    ld = function(x) ldens_std((x - loc) / sc) - log(sc),
    p = function(x) cdf_std((x - loc) / sc),
    q = function(u) invcdf_std(u) * sc + loc
  )
}

location_scale_cases <- list(
  t = list(loc = 0.7, sc = 1.3, df = 4),
  cauchy = list(loc = 0.7, sc = 1.3),
  normal = list(loc = 0.7, sc = 1.3),
  logistic = list(loc = 0.7, sc = 1.3)
)

location_scale_bounds <- list(
  unbounded = c(-Inf, Inf),
  upper_only_left = c(-Inf, -2),
  both_left = c(-4, -2),
  straddling = c(-2, 3),
  both_right = c(2, 4),
  lower_only_right = c(2, Inf),
  opposite_extremes = c(-40, 40),
  extreme_right = c(8, 9),
  extreme_left = c(-9, -8)
)

location_scale_probabilities <- c(
  0,
  1e-12,
  1e-8,
  0.01,
  0.25,
  0.5,
  0.75,
  0.99,
  1 - 1e-8,
  1 - 1e-12,
  1
)

beta_shapes <- list(
  left_skewed = c(5, 2),
  right_skewed = c(2, 5),
  symmetric = c(2, 2),
  u_shaped = c(0.3, 0.4),
  concentrated_left = c(25, 80),
  concentrated_right = c(80, 25)
)

beta_bounds <- list(
  unbounded_request = c(-Inf, Inf),
  lower_natural = c(0, 0.15),
  lower_interior = c(0.05, 0.35),
  central = c(0.2, 0.8),
  upper_interior = c(0.65, 0.95),
  upper_natural = c(0.85, 1),
  clipped_left = c(-2, 0.7),
  clipped_right = c(0.3, 2)
)

beta_probabilities <- c(
  1e-12,
  1e-8,
  1e-5,
  0.01,
  0.25,
  0.5,
  0.75,
  0.99,
  1 - 1e-5,
  1 - 1e-8,
  1 - 1e-12
)
