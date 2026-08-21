test_that("untruncated beta pseudo-target matches stats functions", {
  pseudo <- pseudo_list("beta", list(shape1 = 2, shape2 = 5))
  x <- c(0, seq(0.01, 0.99, length.out = 101), 1)
  u <- c(0, seq(0.001, 0.999, length.out = 101), 1)

  expect_equal(pseudo$lb, 0)
  expect_equal(pseudo$ub, 1)
  expect_equal(pseudo$d(x), dbeta(x, 2, 5), tolerance = 1e-14)
  expect_equal(pseudo$ld(x), dbeta(x, 2, 5, log = TRUE),
               tolerance = 1e-14)
  expect_equal(pseudo$p(x), pbeta(x, 2, 5), tolerance = 1e-14)
  expect_equal(pseudo$q(u), qbeta(u, 2, 5), tolerance = 1e-14)
  expect_equal(
    pseudo$q(log(u), log.p = TRUE),
    qbeta(log(u), 2, 5, log.p = TRUE),
    tolerance = 1e-14
  )
})

test_that("beta truncation intersects the natural support", {
  params <- list(shape1 = 2, shape2 = 5)

  left_clipped <- pseudo_list("beta", params, lb = -2, ub = 0.7)
  right_clipped <- pseudo_list("beta", params, lb = 0.3, ub = 2)
  interior <- pseudo_list("beta", params, lb = 0.2, ub = 0.8)

  expect_equal(c(left_clipped$lb, left_clipped$ub), c(0, 0.7))
  expect_equal(c(right_clipped$lb, right_clipped$ub), c(0.3, 1))
  expect_equal(c(interior$lb, interior$ub), c(0.2, 0.8))
  expect_match(interior$txt, "I\\(0.2 < x < 0.8\\)")

  expect_error(
    pseudo_list("beta", params, lb = 1, ub = 2),
    "must overlap the beta support"
  )
  expect_error(
    pseudo_list("beta", params, lb = -2, ub = 0),
    "must overlap the beta support"
  )
})

test_that("skewed and symmetric truncated beta distributions invert and normalize", {
  for (shape_name in names(beta_shapes)) {
    shape <- beta_shapes[[shape_name]]
    params <- list(shape1 = shape[1], shape2 = shape[2])

    for (configuration in names(beta_bounds)) {
      bounds <- beta_bounds[[configuration]]
      pseudo <- pseudo_list(
        "beta",
        params,
        lb = bounds[1],
        ub = bounds[2]
      )

      expect_equal(pseudo$lb, max(bounds[1], 0),
                   info = paste(shape_name, configuration))
      expect_equal(pseudo$ub, min(bounds[2], 1),
                   info = paste(shape_name, configuration))

      q <- pseudo$q(beta_probabilities)
      p <- pseudo$p(q)
      expect_true(all(diff(q) >= 0),
                  info = paste(shape_name, configuration))
      expect_true(all(q >= pseudo$lb & q <= pseudo$ub),
                  info = paste(shape_name, configuration))

      ## U-shaped beta quantiles can round to a support endpoint. Check
      ## inversion wherever the returned x remains strictly interior.
      keep <- q > pseudo$lb & q < pseudo$ub
      if (any(keep)) {
        expect_equal(
          p[keep],
          beta_probabilities[keep],
          tolerance = 2e-8,
          info = paste(shape_name, configuration)
        )
      }

      expect_equal(
        pseudo$q(log(beta_probabilities), log.p = TRUE),
        q,
        tolerance = 2e-12,
        info = paste(shape_name, configuration)
      )
      expect_equal(
        pseudo$p(c(pseudo$lb, pseudo$ub)),
        c(0, 1),
        tolerance = 0,
        info = paste(shape_name, configuration)
      )
      expect_equal(
        pseudo$q(c(0, 1)),
        c(pseudo$lb, pseudo$ub),
        tolerance = 0,
        info = paste(shape_name, configuration)
      )

      area <- integrate(pseudo$d, pseudo$lb, pseudo$ub,
                        rel.tol = 1e-9)$value
      expect_equal(area, 1, tolerance = 2e-7,
                   info = paste(shape_name, configuration))
    }
  }
})

test_that("beta functions are vectorized and preserve names", {
  pseudo <- pseudo_list(
    "beta",
    list(shape1 = 2, shape2 = 5),
    lb = 0.1,
    ub = 0.9
  )
  x <- setNames(c(0, 0.1, 0.3, 0.7, 0.9, 1), letters[1:6])
  u <- setNames(c(0, 0.1, 0.5, 0.9, 1),
                c("lower", "low", "middle", "high", "upper"))

  expect_equal(pseudo$d(x), vapply(x, pseudo$d, numeric(1)), tolerance = 0)
  expect_equal(pseudo$ld(x), vapply(x, pseudo$ld, numeric(1)), tolerance = 0)
  expect_equal(pseudo$p(x), vapply(x, pseudo$p, numeric(1)), tolerance = 0)
  expect_equal(pseudo$q(u), vapply(u, pseudo$q, numeric(1)), tolerance = 0)
  expect_named(pseudo$d(x), names(x))
  expect_named(pseudo$ld(x), names(x))
  expect_named(pseudo$p(x), names(x))
  expect_named(pseudo$q(u), names(u))
  expect_length(pseudo$d(numeric()), 0)
  expect_length(pseudo$ld(numeric()), 0)
  expect_length(pseudo$p(numeric()), 0)
  expect_length(pseudo$q(numeric()), 0)
})

test_that("extreme beta truncation intervals remain numerically usable", {
  cases <- list(
    list(shape = c(2, 5), bounds = c(1e-100, 1e-80), tolerance = 1e-10),
    list(
      shape = c(2, 5),
      bounds = c(1 - 1e-10, 1 - 1e-12),
      tolerance = 1e-5
    ),
    list(shape = c(0.2, 8), bounds = c(1e-100, 1e-40), tolerance = 1e-10),
    list(
      shape = c(8, 0.2),
      bounds = c(1 - 1e-10, 1 - 1e-12),
      tolerance = 1e-5
    )
  )
  u <- c(1e-8, 1e-5, 0.01, 0.5, 0.99, 1 - 1e-5, 1 - 1e-8)

  for (case in cases) {
    pseudo <- pseudo_list(
      "beta",
      list(shape1 = case$shape[1], shape2 = case$shape[2]),
      lb = case$bounds[1],
      ub = case$bounds[2]
    )
    q <- pseudo$q(u)
    p <- pseudo$p(q)

    expect_true(all(diff(q) >= 0))
    expect_true(all(q >= pseudo$lb & q <= pseudo$ub))
    keep <- q > pseudo$lb & q < pseudo$ub
    if (any(keep)) {
      expect_equal(p[keep], u[keep], tolerance = case$tolerance)
    }
  }
})

test_that("beta quantiles reject invalid probabilities", {
  pseudo <- pseudo_list("beta", list(shape1 = 2, shape2 = 5))

  expect_error(pseudo$q(-0.1), "must be in \\[0, 1\\]")
  expect_error(pseudo$q(1.1), "must be in \\[0, 1\\]")
  expect_error(pseudo$q(NA_real_), "must not be NA or NaN")
  expect_error(pseudo$q(0.2, log.p = 1), "must be TRUE or FALSE")
  expect_error(pseudo$q(0.2, log.p = NA), "must be TRUE or FALSE")
})
