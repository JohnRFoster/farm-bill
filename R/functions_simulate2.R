############################################################################
###
### Simulate swine abundance with population dynamics at the county level
###
############################################################################

simulate_swine <- function(mu_p,
                           sigma_l,
                           sigma_p,
                           sigma_b,
                           size,
                           effort,
                           method,
                           sample_occ,
                           n_county,
                           n_property_per_county,
                           n_pp,
                           n_passes,
                           n_method
){


  # constants
  n_property <- n_property_per_county * n_county # total number of properties
  n_units <- n_property*n_pp # total number of spatio-temporal units
  n_removals <- n_units*n_passes - n_county*n_passes*n_pp # total number of removal events

  county_areas <- c(1451.16, 1891.04, 1688.15)
  property_area <- c(runif(n_property_per_county, 5, 180),
                     runif(n_property_per_county, 50, 200),
                     runif(n_property_per_county, 200, 300))

  PA <- matrix(property_area, n_county, n_property_per_county, byrow = TRUE)
  PC <- matrix(1:n_property, n_county, n_property_per_county, byrow = TRUE)

  unit <- matrix(1:n_units, n_property, n_pp, byrow = TRUE)
  county <- rep(1:3, each = n_property_per_county)

  # simulate the area impacted by each removal event
  # area impacted will differ by method within each property
  gamma <- matrix(NA, n_units, n_passes)
  g <- c(0.05, 0.13, 0.02)
  for(i in 1:n_property){
    # g <- (ga * property_area[i])
    for(j in 1:n_pp){
      for(k in 1:n_method){
        gamma[unit[i, j], which(method[unit[i, j],] == k)] <- g[k]
      }
    }
  }

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
    left_join(county_area_df) |>
    arrange(County, Property, PPNum)

  y <- take |>
    filter(!is.na(`1`))

    keep1 <- which(!is.na(y$`1`))

  county_pp_samples <- y |>
    select(County, PPNum) |>
    distinct() |>
    arrange(County, PPNum) |>
    group_by(County) |>
    mutate(timestep = 1:n()) |>
    ungroup() |>
    pivot_wider(names_from = timestep,
                values_from = PPNum)

  pp <- county_pp_samples |>
    select(-County) |>
    as.matrix()

  start_pp <- pp[,1]
  end_pp <- apply(pp, 1, function(x) max(which(!is.na(x))))
  n_prop <- apply(PC, 1, function(x) length(which(!is.na(x))))

  proportion_surveyed <- take |>
    select(County, Property, PPNum, area_property) |>
    group_by(County, PPNum) |>
    summarise(sum_area = sum(area_property)) |>
    ungroup() |>
    left_join(county_area_df) |>
    mutate(proportion_surveyed = sum_area / area_county) |>
    pivot_wider(id_cols = County,
                names_from = PPNum,
                values_from = proportion_surveyed) |>
    select(-County) |>
    as.matrix()

  y$take <- rowSums(y[,1:5])
  y_mat <- as.matrix(y[,1:5])
  n_passes_y <- apply(y_mat, 1, function(x) length(which(!is.na(x))))

  y_removed <- y |>
    select(County, PPNum, take) |>
    group_by(County, PPNum) |>
    summarise(total_removed = sum(take)) |>
    ungroup() |>
    pivot_wider(names_from = PPNum,
                values_from = total_removed) |>
    select(-County) |>
    as.matrix()

  out <- list(
    take = take,                         # all data
    y = y,                               # data with missing observations
    y_mat = y_mat,                       # take matrix for nimble
    N = N_long,                          # latent property abundance
    M = M_long,                          # latent county abundance
    y_removed = y_removed,               # the number of pigs removed at each PP at each county
    beta = beta,                         # basis function coefficients
    X = X,                               # basis expansions
    sigma_p = sigma_p,                   # error across observation rates
    sigma_b = sigma_b,                   # error across betas
    sigma_l = sigma_l,                   # error in growth (process error)
    mu_p = mu_p,                         # mean observation rate
    p = p,                               # observation rate by method
    effort = effort[keep1,],             # removal effort
    gamma = gamma[keep1,],               # area impacted
    method = method[keep1,],             # trapping method
    n_method = n_method,                 # number of methods
    n_property = n_property,             # number of distinct properties
    pp = pp,                             # PP index within each county,
    start_pp = start_pp,                 # first PP for each county
    end_pp = end_pp,                     # last PP for each county
    prop_surveyed = proportion_surveyed, # proportion of each county surveyed in each PP,
    PA = PA,                             # property areas
    PC = PC,                             # property indexes by county
    n_prop = n_prop,                     # the number of properties in each county
    n_passes_y = n_passes_y              # the number of passes for each property_pp unit
  )

  return(out)
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

  inits <- function(){
    n_init <- (n_mat+max(ds$y_removed, na.rm = TRUE))
    log_mu_proc <- log(n_init)
    z <- n_init - ds$y_removed
    z <- cbind(NA, z[,-ncol(z)])
    # n_init[c(6,12,18),] <- NA

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
