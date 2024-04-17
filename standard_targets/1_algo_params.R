
if (rnd == 1) {  # initial round for tuning

  algo_params <- list(

    normal = list(
      rw = expand.grid(c = seq(0.5, 6.5, by = 0.25)),
      stepping = expand.grid(w = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
                                   5.0, 6.0, 10.0, 20.0, 50.0)),
      gess = expand.grid(loc = c(0.0), sc = seq(0.6, 2.2, by = 0.4), degf = c(1, 5, 20)),
      latent = expand.grid(rate = c(0.01, 0.02, 0.05, 0.10, 0.2, 0.5, 1.0),
                           s_init = c(0.5, 1.0, 3.0, 5.0))
    ),

    gamma = list(
      rw = expand.grid(c = seq(0.5, 6.5, by = 0.25)),
      stepping = expand.grid(w = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
                                   5.0, 6.0, 7.0, 10.0, 20.0, 50.0)),
      gess = expand.grid(loc = seq(0.0, 2.5, by = 0.5),
                         sc = seq(0.5, 3.0, by = 0.5),
                         degf = c(1, 5, 20)),
      latent = expand.grid(rate = c(0.01, 0.02, 0.05, 0.10, 0.2, 0.5, 1.0),
                           s_init = c(0.5, 1.0, 3.0, 5.0))
    ),

    gammalog = list(
      rw = expand.grid(c = seq(0.5, 6.5, by = 0.25)),
      stepping = expand.grid(w = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
                                   5.0, 6.0, 10.0, 20.0, 50.0)),
      gess = expand.grid(loc = seq(0.0, 1.25, by = 0.25),
                         sc = seq(0.5, 2.0, by = 0.5),
                         degf = c(1, 5, 20)),
      latent = expand.grid(rate = c(0.01, 0.02, 0.05, 0.10, 0.2, 0.5, 1.0),
                           s_init = c(0.5, 1.0, 3.0, 5.0))
    ),

    igamma = list(
      rw = expand.grid(c = c(0.5, 1.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 15.0)),
      stepping = expand.grid(w = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
                                   5.0, 6.0, 10.0, 20.0, 50.0)),
      gess = expand.grid(loc = seq(0.0, 1.0, by = 0.25),
                         sc = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 2.0),
                         degf = c(1, 5)),
      latent = expand.grid(rate = c(0.01, 0.02, 0.05, 0.10, 0.2, 0.5, 1.0),
                           s_init = c(1.0, 3.0))
    ),

    igammalog = list(
      rw = expand.grid(c = seq(0.5, 6.5, by = 0.25)),
      stepping = expand.grid(w = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
                                   5.0, 6.0, 10.0, 20.0, 50.0)),
      gess = expand.grid(loc = seq(-1.0, 0, by = 0.2),
                         sc = seq(0.5, 2.5, by = 0.5),
                         degf = c(1, 5, 20)),
      latent = expand.grid(rate = c(0.01, 0.02, 0.05, 0.10, 0.2, 0.5, 1.0),
                           s_init = c(0.5, 1.0, 3.0, 5.0))
    )

  )

} else if (rnd == 2) {  # using optimal settings

  algo_params <- list(

    normal = list(
      rw = expand.grid(c = 2.5), # was 2.5, then 3, then 2.5
      stepping = expand.grid(w = 2.5), # was 2.5, then 2.5, then anything in (2,4)
      gess = expand.grid(loc = 0.0, sc = 1.0, degf = 20), # was 0, 1, 20
      latent = expand.grid(rate = 0.05, s_init = 0.5) # was 0.05
    ),

    gamma = list(
      rw = expand.grid(c = 4.0), # was 3.5, then 3.75
      stepping = expand.grid(w = 6.0), # was 3.5, 7 about same
      gess = expand.grid(loc = 2.0, sc = 1.5, degf = 1), # was 2.5, 2.0, 5 (is second or third best)
      latent = expand.grid(rate = 0.05, s_init = 3.0) # was 0.05
    ),

    igamma = list(
      rw = expand.grid(c = 7), # was 5.5, then 5, then 7 (which is probably safer)
      stepping = expand.grid(w = 1.5), # was 1.5
      gess = expand.grid(loc = 0.5, sc = 0.4, degf = 1), # was 0.5 0.25 1, .5,.4,1 won again
      latent = expand.grid(rate = 0.02, s_init = 3.0) # was 0.01, uncertain (.02 looks good)
    ),

    gammalog = list(
      rw = expand.grid(c = 1.75),
      stepping = expand.grid(w = 2.0),
      gess = expand.grid(loc = 0.75, sc = 1.0, degf = 5),
      latent = expand.grid(rate = 0.1, s_init = 1.0)
    ),

    igammalog = list(
      rw = expand.grid(c = 2.0),
      stepping = expand.grid(w = 2.0),
      gess = expand.grid(loc = -0.4, sc = 1.0, degf = 5),
      latent = expand.grid(rate = 0.1, s_init = 0.5)
    )

  )

}
