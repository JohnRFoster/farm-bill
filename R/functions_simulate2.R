## functions to simulate pig population abundance

mu_p <- -2.2
sigma_l <- 0.2
sigma_p <- 0.6
sigma_b <- 0.1
size <- seq(0.5, 5, length.out = n_county)

simulate_swine <- function(mu_p,
                           sigma_l,
                           sigma_p,
                           sigma_b,
                           size,
                           effort,
                           method,
                           gamma,
                           sample_occ
){


  county_areas <- c(1451.16, 1891.04, 1688.15)
  property_area <- c(9.73, 73.84, 74.81, 81.48, 30.75, 28.72, 62.85, 22.27, 47.77,
                     86.45, 67.69, 69.71, 40.28, 13.01, 71.96)
  PA <- matrix(property_area, n_county, n_property_per_county, byrow = TRUE)
  PC <- matrix(1:n_property, n_county, n_property_per_county, byrow = TRUE)

  # constants
  n_county <- 3 # number of counties
  n_property_per_county <- 5 # number of properties within each county
  n_pp <- 10 # number of primary periods
  n_passes <- 5 # number of removals per sampled primary period
  n_method <- 3 # number of removal methods
  n_property <- n_property_per_county * n_county # total number of properties
  n_units <- n_property*n_pp # total number of spatio-temporal units
  n_removals <- n_units*n_passes - n_county*n_passes*n_pp # total number of removal events

  unit <- matrix(1:n_units, n_property, n_pp, byrow = TRUE)
  county <- rep(1:3, each = n_property_per_county)

  Xdat <- cbind(1, 1:n_pp/n_pp)
  X <- cbind(1, bs(Xdat[,2], 5))

  logit_p <- rnorm(n_method, mu_p, sigma_p)
  p <- ilogit(logit_p)
  beta <- rnorm(ncol(X), 0, sigma_b)
  alpha <- rnorm(n_property, 0, 1)

  # initial density of pigs in each county
  # density of pigs is 7 km^2
  M <- matrix(NA, n_county, n_pp)
  M[,1] <- round(county_areas*7)
  M_avail <- M
  N <- matrix(NA, n_property, n_pp)
  C <- matrix(NA, n_units, n_passes)      # counts
  removed <- matrix(NA, n_property, n_pp) # number removed

  for(i in 1:n_county){
    for(t in 1:n_pp){

      # property areas for county
      pi_n <- PA[i, ]

      # proportion of county sampled in PP
      phi <- sum(pi_n) / county_areas[i]

      # the number of pigs available to be captured
      M_avail[i, t] <- rbinom(1, M[i, t], phi)

      # place available pigs into properties based on property size
      N[PC[i,], t] <- rmultinom(1, M_avail[i, t], pi_n)

      # data model - determine take
      # determine which properties are sampled
      properties_sampled <- PC[i, which(!is.na(sample_occ[PC[i,], t]))]

      for(j in seq_along(properties_sampled)){
        prop <- properties_sampled[j]
        theta <- rep(NA, n_passes)
        for(l in seq_len(n_passes)){
          theta[l] <- 1 - (1 - p[method[unit[prop, t], l]]) ^ effort[unit[prop, t], l]
          if(l == 1){
            p_star <- gamma[unit[prop, t], l] * theta[l]
            N_avail <- N[prop, t]
          } else {
            p_star <- gamma[unit[prop, t], l] * theta[l] *
              exp(
                sum(
                  log(
                    (1 - gamma[unit[prop, t], 1:(l - 1)]) +
                      (gamma[unit[prop, t], 1:(l - 1)] * (1 - theta[1:(l - 1)]))
                  )
                )
              )
            N_avail <- N[prop, t] - sum(C[unit[prop, t], 1:(l-1)])
          }

          if(p_star >= 1) stop("Observation rate greater than 1!")

          # for the random draws use min so we don't remove more than available
          C[unit[prop, t], l] <- min(rpois(1, p_star * N_avail), N_avail)
          removed[prop, t] <- sum(C[unit[prop, t], ])
        }
      }

      if(t < n_pp){
        # determine growth rate
        log_lambda_mu <- inprod(X[t,], beta)
        log_lambda <- rnorm(1, log_lambda_mu, sigma_l)
        z <- M[i, t] - sum(removed[PC[i,], t], na.rm = TRUE)

        # predict
        M[i, t+1] <- rpois(1, exp(log(z) + log_lambda))
      }

    }
  }

  for(i in 1:n_county){
    plot(log(M[i,]+1), type = "l", ylim = c(0, max(log(M[i,]+1))))
    for(j in 1:n_property_per_county) lines(log(N[PC[i, j],]+1), col = j+1)
  }

  colnames(M) <- 1:n_pp
  colnames(N) <- 1:n_pp
  colnames(C) <- 1:n_passes

  M_long <- M |>
    as_tibble() |>
    mutate(County = 1:n()) |>
    pivot_longer(cols = -County,
                 values_to = "County_abundance",
                 names_to = "PPNum")

  N_long <- N |>
    as_tibble() |>
    mutate(County = rep(1:n_county,
                        each = n_property_per_county),
           Property = 1:n()) |>
    pivot_longer(cols = -c(County, Property),
                 values_to = "Property_abundance",
                 names_to = "PPNum")

  prop_area_df <- tibble(
    area_property = property_area,
    Property = 1:length(property_area)
  )

  county_area_df <- tibble(
    area_county = county_areas,
    County = 1:length(county_areas)
  )

  take <- C |>
    as_tibble() |>
    mutate(County = rep(1:n_county, each = n_property_per_county*n_pp),
           Property = rep(1:n_property, each = n_pp),
           Property_county = rep(rep(1:n_property_per_county, each = n_pp), n_county),
           PPNum = rep(rep(1:n_pp), n_county*n_property_per_county)) |>
    left_join(prop_area_df) |>
    left_join(county_area_df)

  # keep1 <- which(!is.na(y$`1`))

  y <- take |>
    filter(!is.na(`1`))

  county_pp_samples <- y |>
    select(County, PPNum) |>
    distinct() |>
    group_by(County) |>
    mutate(timestep = 1:n()) |>
    ungroup() |>
    pivot_wider(names_from = timestep,
                values_from = PPNum)

  county_pp_samples_mat <- county_pp_samples |>
    select(-County) |>
    as.matrix()

  start_pp <- county_pp_samples_mat[,1]
  end_pp <- apply(county_pp_samples_mat, 1, function(x) max(which(!is.na(x))))

  y |>
    select(County, Property, PPNum, area_property) |>
    pivot_wider(values_from = area_property,
                )


    n_timestep <- y |>
      slice(keep1) |>
      group_by(Property) |>
      mutate(timestep = 1:n()) |>
      summarise(n = max(timestep)) |>
      pull(n)

    obs_units <- y |>
      slice(keep1) |>
      mutate(obs_id = as.numeric(paste0(County, PPNum)))

    colnames(N) <- stringr::str_replace(1:n_pp, "\\d+",
                                        function(x) sprintf("%02d", as.integer(x)))

    p_obs_id <- y |>
      slice(keep1) |>
      mutate(obs_id = paste(Property, PPNum, sep = "_")) |>
      pull(obs_id)



    n_pc <- y |>
      group_by(County, PPNum) |>
      tally() |>
      pull(n)

    pp <- y |>
      slice(keep1) |>
      select(Property, PPNum) |>
      group_by(Property) |>
      mutate(timestep_idx = 1:n()) |>
      ungroup() |>
      pivot_wider(names_from = timestep_idx,
                  values_from = PPNum,
                  id_cols = Property) |>
      select(-Property) |>
      as.matrix()

    npc <- y |>
      slice(keep1) |>
      select(County, PPNum, Property_county) |>
      pivot_wider(names_from = Property_county,
                  values_from = Property_county,
                  id_cols = c(County, PPNum)) |>
      arrange(County, PPNum) |>
      select(-County, -PPNum) |>
      as.matrix()


    y_removed <- matrix(NA, nrow(removed), ncol(removed))
    for(i in 1:nrow(removed)){
      vec <- removed[i, which(!is.na(removed[i,]))]
      if(length(vec) == 0){
        y_removed[i,] <- NA
      } else {
        n0 <- length(vec)
        y_removed[i, 1:n0] <- vec
      }

    }

    out[[s]] <- list(
      y = y,                          # data with missing observations
      N = N_long,                     # latent property abundance
      M = N_sums,                          # latent county abundance
      y_removed = y_removed,          # the number of pigs removed at each timestep
      r_c = r_c,                          # mean growth rate
      alpha_p = alpha_p,
      beta = beta,
      x1 = x1,
      x2 = x2,
      phi = phi,                          # mean survival rate
      rho = rho,                          # mean reproduction rate
      size = size,
      sigma_proc = sigma_proc,        # process error precision
      sigma_p = sigma_p,        # process error precision
      p = p,                          # observation rate by method
      effort = effort[keep1,],   # removal effort
      gamma = gamma[keep1,],     # area impacted
      method = method[keep1,],   # trapping method
      n_method = n_method,            # number of methods
      n_property = n_property,        # number of distinct properties
      n_timestep = n_timestep,        # number of timesteps in each property
      n_pc = n_pc,                    # number of properties sampled in each county in each PP
      pp = pp,                        # timestep index within each property
      # n_county_units = nrow(M),       # number of distinct county by PP units
      n_prop_county = n_prop_county1, # properties within counties that are sampled in each PP
      n_prop = n_prop                 # how many properties are sampled within counties that are sampled in each PP
    )
  }

  names(out) <- proportion_surveyed

  return(
    out
  )
}


# build init function
make_inits <- function(ds, process){
  n_pp <- 10
  N_mat_na <- ds$N |>
    pivot_wider(names_from = PPNum,
                values_from = Property_abundance,
                values_fill = NA) |>
    select(sprintf("%02d", as.integer(1:n_pp))) |>
    as.matrix()
  n_mat <- matrix(NA, nrow(N_mat_na), ncol(N_mat_na))
  for(i in 1:nrow(N_mat_na)){
    vec <- N_mat_na[i, which(!is.na(N_mat_na[i,]))]
    n0 <- length(vec)
    n_mat[i, 1:n0] <- vec
  }

  y_init <- ds$y

  keep1 <- which(!is.na(y_init$V1))
  keep2 <- which(y_init$Property %in% c(6, 12, 18))
  keep_index <- sort(c(keep1, keep2))

  y_init <- y_init |>
    slice(keep_index) |>
    select(1:5) |> as.matrix()

  y_na <- is.na(y_init)
  y_init[y_na] <- rpois(1, 5000)
  y_init[!y_na] <- NA


  inits <- function(){
    n_init <- (n_mat+max(ds$y_removed, na.rm = TRUE))
    log_mu_proc <- log(n_init)
    z <- n_init - ds$y_removed
    z <- cbind(NA, z[,-ncol(z)])
    # n_init[c(6,12,18),] <- NA
    log_mu_proc <- log_mu_proc[-c(6,12,18),]

    l1 <- list(
      # y = y_init,
      n = n_init,
      z = z,
      log_mu_proc = log_mu_proc,
      M = ds$M$sumAbundance,
      z1 = log(n_init[,1]),
      tau_p = 1/runif(1, 2, 4)^2,
      mu_p = boot::logit(runif(3, 0.2, 0.6)),
      size = runif(3, 0.01, 5)
    )

    if(process == "dm"){
      z <- n_mat - ds$y_removed
      z <- cbind(NA, z[,-ncol(z)])
      s <- r <- dm <- array(NA, dim = c(ds$n_property, max(ds$n_timestep), n_pp))

      phi <- runif(1, 0.2, 0.7)
      rho <- runif(1)
      mu_phi <- boot::logit(phi)
      log_rho <- log(rho)

      # phi <- runif(1, 0.7, 0.9)
      # rho <- runif(1, 0, 0.4)

      for(i in 1:ds$n_property){
        for(t in 2:ds$n_timestep[i]){
          dm[i, t, ds$pp[i, t-1]] <- z[i, t]
          for(j in ds$pp[i, t-1]:(ds$pp[i, t]-1)){
            s[i, t, j] <- rbinom(1, dm[i, t, j], phi)
            r[i, t, j] <- rpois(1, dm[i, t, j] * rho)
            # r[i, t, j] <- rpois(1, rho)
            dm[i, t, j+1] <- s[i, t, j] + r[i, t, j]
          }
        }
      }

      l2 <- list(
        # mu_w = rnorm(1),
        # rho = rho,
        # phi = phi,
        s = s + 1,
        r = r,
        # ddm = dm
        log_rho = log_rho,
        mu_phi = mu_phi
      )

    } else if(process == "exponential") {
      sigma_p <- runif(1, 0, 0.2)
      log_mean_r <- runif(3, 0.8, 1)
      l2 <- list(
        log_mean_r = log_mean_r,
        alpha_p = rnorm(ds$n_property, 0, sigma_p),
        tau_lambda = 1/sigma_p^2,
        beta = rnorm(2)
      )
    }

    append(l1, l2)

  }
  return(inits)
}
