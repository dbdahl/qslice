loc <- -1.0
sc <- 0.1
df <- 2

(lb_vec <- c(-Inf, loc - 2.0, loc + 1.0))
(ub_vec <- c(loc - 1.0, loc + 2.0, Inf))
(x_vec <- c(loc - 0.1, loc + 0.1))

for (lb in lb_vec) {
  for (ub in ub_vec) {
    for (x in x_vec) {
      txt0 <- paste0("lb = ", lb, "  ub = ", ub, "  x = ", x, "\n")
      cat(txt0)

      ps <- tryCatch(
        {
          value <- .pseudo_t_list(loc = loc, sc = sc, df = df, lb = lb, ub = ub)
        },
        error = function(e) {
          print(e$message)
        }
      )

      if (is.list(ps)) {
        dens0 <- ifelse(
          x > lb & x < ub,
          dt((x - loc) / sc, df = df) /
            sc /
            (pt((ub - loc) / sc, df = df) - pt((lb - loc) / sc, df = df)),
          0
        )
        dens1 <- ps$d(x)

        u0 <- ifelse(
          x > lb,
          ifelse(
            x < ub,
            (pt((x - loc) / sc, df = df) - pt((lb - loc) / sc, df = df)) /
              (pt((ub - loc) / sc, df = df) - pt((lb - loc) / sc, df = df)),
            1
          ),
          0
        )
        u1 <- ps$p(x)

        q0 <- qt(
          (pt((lb - loc) / sc, df = df) +
            u1 * (pt((ub - loc) / sc, df = df) - pt((lb - loc) / sc, df = df))),
          df = df
        ) *
          sc +
          loc
        q1 <- ps$q(u1)

        txt <- paste0(
          "density = ",
          dens0,
          " density (pseudo) = ",
          dens1,
          " diff = ",
          dens0 - dens1,
          "\ncdf = ",
          u0,
          " cdf (pseudo) = ",
          u1,
          " diff = ",
          u0 - u1,
          "\nquantile = ",
          q0,
          " quantile (pseudo) = ",
          q1,
          " diff = ",
          q0 - q1,
          "\n\n"
        )

        cat(txt)
      }
    }
  }
}
