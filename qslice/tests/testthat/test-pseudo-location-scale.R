test_that("untruncated location-scale families match stats distributions", {
  for (family in names(location_scale_cases)) {
    params <- location_scale_cases[[family]]
    free <- free_locscale(family, params)
    pseudo <- pseudo_list(family = family, params = params)
    x <- seq(params$loc - 5 * params$sc,
             params$loc + 5 * params$sc,
             length.out = 101)
    u <- seq(0.001, 0.999, length.out = 101)

    expect_equal(pseudo$d(x), free$d(x), tolerance = 1e-13,
                 info = family)
    expect_equal(pseudo$ld(x), free$ld(x), tolerance = 1e-13,
                 info = family)
    expect_equal(pseudo$p(x), free$p(x), tolerance = 1e-13,
                 info = family)
    expect_equal(pseudo$q(u), free$q(u), tolerance = 1e-12,
                 info = family)
  }
})

test_that("routine truncated location-scale families match direct formulas", {
  for (family in names(location_scale_cases)) {
    params <- location_scale_cases[[family]]
    free <- free_locscale(family, params)
    loc <- params$loc
    sc <- params$sc
    lb <- loc - 1.7 * sc
    ub <- loc + 2.3 * sc
    pseudo <- pseudo_list(family, params, lb = lb, ub = ub)

    p_lb <- free$p(lb)
    interval_probability <- free$p(ub) - p_lb
    x <- seq(lb + 0.01 * sc, ub - 0.01 * sc, length.out = 101)
    u <- seq(0.001, 0.999, length.out = 101)

    expect_equal(
      pseudo$d(x),
      free$d(x) / interval_probability,
      tolerance = 1e-12,
      info = family
    )
    expect_equal(
      pseudo$p(x),
      (free$p(x) - p_lb) / interval_probability,
      tolerance = 1e-12,
      info = family
    )
    expect_equal(
      pseudo$q(u),
      free$q(p_lb + u * interval_probability),
      tolerance = 1e-11,
      info = family
    )
  }
})

test_that("all location-scale truncation configurations are normalized", {
  for (family in names(location_scale_cases)) {
    params <- location_scale_cases[[family]]
    loc <- params$loc
    sc <- params$sc

    for (configuration in names(location_scale_bounds)) {
      standardized_bounds <- location_scale_bounds[[configuration]]
      lb <- loc + sc * standardized_bounds[1]
      ub <- loc + sc * standardized_bounds[2]
      pseudo <- pseudo_list(family, params, lb = lb, ub = ub)

      expect_equal(
        pseudo$p(c(lb, ub)),
        c(0, 1),
        tolerance = 0,
        info = paste(family, configuration)
      )
      expect_equal(
        pseudo$q(c(0, 1)),
        c(lb, ub),
        tolerance = 0,
        info = paste(family, configuration)
      )

      q <- pseudo$q(location_scale_probabilities)
      keep <- is.finite(q) & q > lb & q < ub
      if (any(keep)) {
        expect_equal(
          pseudo$p(q[keep]),
          location_scale_probabilities[keep],
          tolerance = 2e-8,
          info = paste(family, configuration)
        )
      }

      area <- integrate(pseudo$d, lb, ub, rel.tol = 1e-9)$value
      expect_equal(area, 1, tolerance = 2e-7,
                   info = paste(family, configuration))
    }
  }
})

test_that("location-scale pseudo-target functions are vectorized", {
  for (family in names(location_scale_cases)) {
    params <- location_scale_cases[[family]]
    loc <- params$loc
    sc <- params$sc
    pseudo <- pseudo_list(
      family,
      params,
      lb = loc - 2 * sc,
      ub = loc + 3 * sc
    )

    x <- setNames(
      c(loc - 3 * sc, loc - sc, loc, loc + sc, loc + 4 * sc),
      c("outside_left", "left", "center", "right", "outside_right")
    )
    u <- setNames(c(0, 0.1, 0.5, 0.9, 1),
                  c("lower", "low", "middle", "high", "upper"))

    expect_equal(pseudo$d(x), vapply(x, pseudo$d, numeric(1)),
                 tolerance = 0, info = family)
    expect_equal(pseudo$ld(x), vapply(x, pseudo$ld, numeric(1)),
                 tolerance = 0, info = family)
    expect_equal(pseudo$p(x), vapply(x, pseudo$p, numeric(1)),
                 tolerance = 0, info = family)
    expect_equal(pseudo$q(u), vapply(u, pseudo$q, numeric(1)),
                 tolerance = 0, info = family)

    expect_named(pseudo$d(x), names(x))
    expect_named(pseudo$ld(x), names(x))
    expect_named(pseudo$p(x), names(x))
    expect_named(pseudo$q(u), names(u))

    expect_length(pseudo$d(numeric()), 0)
    expect_length(pseudo$ld(numeric()), 0)
    expect_length(pseudo$p(numeric()), 0)
    expect_length(pseudo$q(numeric()), 0)
    expect_length(pseudo$d(loc), 1)
    expect_length(pseudo$ld(loc), 1)
    expect_length(pseudo$p(loc), 1)
    expect_length(pseudo$q(0.5), 1)
  }
})

test_that("location-scale quantiles reject invalid probabilities", {
  pseudo <- pseudo_list("normal", list(loc = 0, sc = 1))

  expect_error(pseudo$q(c(0.2, 1.1)), "must be in \\[0, 1\\]")
  expect_error(pseudo$q(c(0.2, -0.1)), "must be in \\[0, 1\\]")
  expect_error(pseudo$q(c(0.2, NA_real_)), "must be in \\[0, 1\\]")
  expect_error(pseudo$q("0.2"), "must be numeric")
})
