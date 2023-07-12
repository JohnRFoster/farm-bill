## functions to simulate pig population abundance

# n_county <- 3
# n_property_per_county <- 5 # number of properties in each county
# n_property <- n_property_per_county * n_county
# n_pp <- 10
# n_passes <- 5
#
# county_areas_low <- runif(n_county, 100, 999)
# county_areas_mid <- runif(n_county, 999, 4999)
# county_areas_high <- runif(n_county, 5000, 9999)
#
# get_property_areas <- function(county_areas){
#   property_areas <- matrix(NA, n_county, n_property_per_county)
#   for(i in 1:n_county){
#     pa <- runif(n_property_per_county, 20, county_areas[i]*0.8)
#     while(sum(pa) > county_areas[i]*0.8){
#       pa <- runif(n_property_per_county, 20, county_areas[i]*0.8)
#     }
#     property_areas[i, ] <- pa
#   }
#   property_areas
# }
#
# property_areas_low <- get_property_areas(county_areas_low)
# property_areas_mid <- get_property_areas(county_areas_mid)
# property_areas_high <- get_property_areas(county_areas_high)

simulate_swine <- function(process,
                           likelihood,
                           r_c,
                           phi,
                           rho,
                           beta,
                           x1,
                           x2,
                           alpha_p,
                           sigma_proc,
                           sigma_p,
                           p,
                           size,
                           effort,
                           method,
                           gamma,
                           sample_occ
                           ){


  county_areas_low <- c(451.16, 891.04, 688.15)
  county_areas_mid <- c(4643.96, 2645.81, 2759.02)
  county_areas_high <- c(6797.46, 7412.5, 9148.72)

  property_area <- c(9.73, 73.84, 74.81, 81.48, 30.75, 28.72, 62.85, 22.27, 47.77,
                     86.45, 67.69, 69.71, 40.28, 13.01, 71.96, 45.46, 47.51, 38.46)

  # the area not sampled in each county
  unsampled <- function(ca, pa){
    c(ca[1] - sum(pa[1:6]),
      ca[2] - sum(pa[7:12]),
      ca[3] - sum(pa[13:18]))
  }

  # the unsampled area is the area of the completely latent property
  unsampled_low <- unsampled(county_areas_low, property_area)
  unsampled_mid <- unsampled(county_areas_mid, property_area)
  unsampled_high <- unsampled(county_areas_high, property_area)
  U <- matrix(c(unsampled_low, unsampled_mid, unsampled_high), 3, 3)

  P <- matrix(property_area, 6, 3)

  # constants
  n_county <- 3 # number of counties
  n_property_per_county <- 10 # number of properties within each county
  n_pp <- 10 # number of primary periods
  n_passes <- 5 # number of removals per sampled primary period
  n_method <- 3 # number of removal methods
  n_property <- n_property_per_county * n_county # total number of properties
  n_units <- n_property*n_pp # total number of spatio-temporal units
  n_removals <- n_units*n_passes - n_county*n_passes*n_pp # total number of removal events

  unit <- matrix(1:n_units, n_property, n_pp, byrow = TRUE)
  county <- rep(1:3, each = n_property_per_county)

  # initial density of pigs in each county
  # max density of pigs is 10 km^2
  M_low <- round(county_areas_low*10)
  M_mid <- round(county_areas_mid*10)
  M_high <- round(county_areas_high*10)
  M_vec <- c(M_low, M_mid, M_high)
  M <- matrix(M_vec, 3, 3)

  initial_n <- function(M, pa, u){
    area_vec <- c(pa, u)
    N <- rmultinom(1, M, area_vec / sum(area_vec))
    N <- as.vector(N)
    l <- round(N[length(N)] / 4)
    c(N[-length(N)], rep(l, 4))
  }

  proportion_surveyed <- c("low", "med", "high")

  # loop over area scenarios
  out <- list()
  for(s in 1:3){
    # initial log abundance for each property
    n <- vector()
    for(m in 1:3){
      N1 <- initial_n(M[m, s], P[,m], U[m, s])
      n <- c(n, N1)
    }

    # storage
    XN <- N <- matrix(NA, n_property, n_pp)       # latent abundance [property, primary period] w/ stochasticity
    x <- matrix(NA, n_property, n_pp)       # expected log abundance
    M <- matrix(NA, n_county, n_pp)
    C <- matrix(NA, n_units, n_passes)      # counts
    removed <- matrix(NA, n_property, n_pp) # number removed

    M[1,1] <- sum(n[1:10])
    M[2,1] <- sum(n[11:20])
    M[3,1] <- sum(n[21:30])
    N[,1] <- n
    XN[,1] <- n
    x[,1] <- log(n)


    # simulate dynamics
    for(i in seq_len(n_property)){
      for(j in seq_len(n_pp)){

        if(j %in% sample_occ[i,]){ # if a removal occasion, determine how many pigs removed
          theta <- rep(NA, n_passes)
          for(l in 1:n_passes){
            theta[l] <- 1 - (1 - p[method[unit[i, j], l]]) ^ effort[unit[i, j], l]
            if(l == 1){
              p_star <- gamma[unit[i, j], l] * theta[l]
              N_avail <- N[i, j]
            } else {
              p_star <- gamma[unit[i, j], l] * theta[l] *
                exp(
                  sum(
                    log(
                      (1 - gamma[unit[i, j], 1:(l - 1)]) +
                        (gamma[unit[i, j], 1:(l - 1)] * (1 - theta[1:(l - 1)]))
                    )
                  )
                )
              N_avail <- N[i, j] - sum(C[unit[i, j], 1:(l-1)])
            }

            if(p_star >= 1) stop("Observation rate greater than 1!")

            ex_removed <- p_star * N_avail
            if(likelihood == "deterministic"){
              C[unit[i, j], l] <- round(ex_removed)
            } else if(likelihood == "poisson"){
              # for the random draws use min so we don't remove more than available
              C[unit[i, j], l] <- min(rpois(1, ex_removed), N_avail)
            } else if(likelihood == "nb"){
              C[unit[i, j], l] <- min(rnbinom(1, mu = ex_removed, size = size[county[i]]), N_avail)
            }

          }
          removed[i, j] <- sum(C[unit[i, j], ])
        } else {
          removed[i, j] <- NA
        }

        z <- N[i, j] - if_else(is.na(removed[i, j]), 0, removed[i, j])
        if(z < 0) stop("Number removed more than population!")

        if(j < n_pp){

          if(process == "exponential"){



            growth <- exp(log(r_c[county[i]]) +
                            beta[1]*x1[county[i], j] +
                            beta[2]*x2[county[i], j] +
                            alpha_p[i])
            log_mu <- log(growth) + log(z)
            x[i, j+1] <- rnorm(1, log_mu, sigma_proc) # add process error
            N[i, j+1] <- rpois(1, exp(x[i, j+1]))     # add stochasticity
          } else if(process == "dm"){
            S <- rbinom(1, z, phi) # survival process
            R <- rpois(1, z * (rho[county[i]] + alpha_p[i])) # recruitment process
            # R <- rZIP(1, rho, w) # recruitment process
            x[i, j+1] <- S + R
            N[i, j+1] <- rpois(1, x[i, j+1])     # add stochasticity
          }

        }
      }
    }

    N_sums <- tibble(
      county_1 = colSums(N[1:6,]),
      county_2 = colSums(N[7:12,]),
      county_3 = colSums(N[13:18,])
    ) |>
      mutate(PPNum = 1:n()) |>
      pivot_longer(cols = -PPNum,
                   names_to = "County",
                   values_to = "sumAbundance") |>
      mutate(proportion_surveyed = proportion_surveyed[s])

    colnames(C) <- paste0("V", 1:n_passes)

    y <- C |>
      as_tibble() |>
      mutate(County = rep(1:n_county, each = n_property_per_county*n_pp),
             Property = rep(1:n_property, each = n_pp),
             Property_county = rep(rep(1:n_property_per_county, each = n_pp), n_county),
             PPNum = rep(rep(1:n_pp), n_county*n_property_per_county))

    keep1 <- which(!is.na(y$V1))
    keep2 <- which(y$Property %in% c(6, 12, 18))
    keep_index <- sort(c(keep1, keep2))

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

    N_long <- N |>
      as_tibble() |>
      mutate(County = rep(1:n_county,
                          each = n_property_per_county),
             Property = 1:n()) |>
      pivot_longer(cols = -c(County, Property),
                   values_to = "Property_abundance",
                   names_to = "PPNum") |>
      mutate(obs_id = paste(Property, as.numeric(PPNum), sep = "_")) |>
      # filter(obs_id %in% p_obs_id) |>
      select(-obs_id)

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

    n_prop_county1 <- matrix(NA, nrow(npc), ncol(npc))
    n_prop <- rep(NA, nrow(npc))
    for(i in 1:nrow(npc)){
      vec <- npc[i, which(!is.na(npc[i,]))]
      n_prop[i] <- length(vec)
      n_prop_county1[i, 1:n_prop[i]] <- vec
    }

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
