## functions to simulate pig population abundance

# n_county <- 3
# n_property_per_county <- 5 # number of properties in each county
# n_pp <- 10
# n_passes <- 5

simulate_swine <- function(n_county, n_property_per_county, n_pp, n_passes, likelihood){

  # constants
  n_method <- 3
  n_property <- n_property_per_county * n_county
  n_units <- n_property*n_pp
  property <- seq_len(n_property)
  county <- rep(1:n_county, each = n_property)
  county_area <- runif(n_county, 1000, 3000)
  property_areas <- runif(n_property, 10, 150)

  # ecological process parameters
  mean_r <- runif(1, 0.8, 2) # mean growth rate
  sigma_proc <- runif(1, 0.1, 1) # process error
  size <- runif(n_county, 0, 5) # overdispersion

  # observation process parameters
  p <- runif(n_method, 0.2, 0.6)

  # observation "data"
  effort <- matrix(sample(1:3, n_units*n_passes, c(0.9, 0.07, 0.03), replace = TRUE), n_units, n_passes)
  method <- matrix(round(runif(n_units*n_passes, 0.5, 2.8)), n_units, n_passes)

  unit <- matrix(seq_len(n_units), n_property, n_pp, byrow = TRUE)

  gamma <- matrix(NA, n_units, n_passes)
  for(i in 1:n_property){
    g <- runif(3, 0, 0.8)
    for(j in 1:n_pp){
      for(k in 1:n_method){
        gamma[unit[i, j], which(method[unit[i, j],] == k)] <- g[k]
      }
    }
  }

  # storage
  N <- matrix(NA, n_property, n_pp)  # latent abundance [property, primary period] w/ stochasticity
  x <- matrix(NA, n_property, n_pp)  # expected log abundance
  C <- matrix(NA, n_units, n_passes) # counts
  removed <- matrix(NA, n_property, n_pp)


  # process model
  for(i in seq_len(n_property)){

    # initial abundance for each property
    x[i, 1] <- log(max(round(rnorm(1, 100, 50)), 10))
    N[i, 1] <- rpois(1, exp(x[i, 1]))

    # when we sample (between 3 and 8 random PP are sampled)
    n_samps <- round(runif(1, 2.6, 8.4))
    samps <- sample(n_pp, n_samps)

    for(j in seq_len(n_pp)){

      if(j %in% samps){ # if a removal occasion, determine how many pigs removed
        theta <- rep(NA, n_passes)
        for(l in 1:n_passes){
          theta[l] <- 1 - (1 - p[method[unit[i, j], l]]) ^ effort[unit[i, j], l]
          if(l == 1){
            p_star <- gamma[unit[i, j], l] *
              theta[l] *
              (1 - gamma[unit[i, j], l]) +
              (gamma[unit[i, j], l] * (1 - theta[l]))
            N_avail <- N[i, j]
          } else {
            p_star <- gamma[unit[i, j], l] *
              theta[l] *
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
          # N_avail <- N[i, j]
          if(p_star >= 1) stop("Observation rate greater than 1!")

          if(likelihood == "poisson"){
            C[unit[i, j], l] <- min(rpois(1, p_star * N_avail), N_avail)
          } else if(likelihood == "nb"){
            C[unit[i, j], l] <- min(rnbinom(1, mu = p_star * N_avail, size = size[county[i]]), N_avail)
          }

        }
        removed[i, j] <- sum(C[unit[i, j], ])
      } else {
        removed[i, j] <- NA
      }

      z <- N[i, j] - if_else(is.na(removed[i, j]), 0, removed[i, j])
      if(z < 0) stop("Number removed more than population!")
      mu <- log(mean_r) + log(z)

      if(j < n_pp){
        x[i, j+1] <- rnorm(1, mu, sigma_proc) # add process error
        N[i, j+1] <- rpois(1, exp(x[i, j+1])) # add stochasticity
        # N[i, j+1] <- max(N[i, j+1], 1)
      }

    }
  }

  y <- C |>
    as_tibble() |>
    mutate(County = rep(1:n_county, each = n_property_per_county*n_pp),
           Property = rep(1:n_property, each = n_pp),
           Property_county = rep(rep(1:n_property_per_county, each = n_pp), n_county),
           PPNum = rep(rep(1:n_pp), n_county*n_property_per_county),
           prop_area = rep(property_areas, each = n_pp))

  n_timestep <- y |>
    filter(!is.na(V1)) |>
    group_by(Property) |>
    mutate(timestep = 1:n()) |>
    summarise(n = max(timestep)) |>
    pull(n)

  obs_units <- y |>
    filter(!is.na(V1)) |>
    group_by(County, PPNum) |>
    summarise(total_area_sampled = sum(prop_area)) |>
    ungroup() |>
    mutate(obs_id = as.numeric(paste0(County, PPNum)))

  colnames(N) <- stringr::str_replace(1:n_pp, "\\d+",
                              function(x) sprintf("%02d", as.integer(x)))

  p_obs_id <- y |>
    filter(!is.na(V1)) |>
    mutate(obs_id = as.numeric(paste0(Property, PPNum))) |>
    pull(obs_id)

  N_long <- N |>
    as_tibble() |>
    mutate(County = rep(1:n_county,
                        each = n_property_per_county),
           Property = 1:n(),
           county_area = rep(county_area, each = n_property_per_county),
           property_area = property_areas) |>
    pivot_longer(cols = -c(County, county_area, Property, property_area),
                 values_to = "Property_abundance",
                 names_to = "PPNum") |>
    mutate(obs_id = as.numeric(paste0(Property, as.numeric(PPNum)))) |>
    filter(obs_id %in% p_obs_id) |>
    select(-obs_id)

  # county level abundance is the sum across properties in each timestep
  Nsum <- N_long |>
    group_by(County, PPNum, county_area) |>
    summarise(Nsum = sum(Property_abundance)) |>
    mutate(PPNum = as.numeric(PPNum),
           obs_id = as.numeric(paste0(County, PPNum))) |>
    ungroup()

  size_tb <- tibble(
    County = 1:n_county,
    size = size
  )

  # simulate county abundance
  M <- left_join(obs_units, Nsum) |>
    left_join(size_tb) |>
    mutate(mu = Nsum / total_area_sampled * county_area,
           County_abundnace = rnbinom(length(1:n()), mu = mu, size = size),
           County_abundnace = pmax(County_abundnace, Nsum))

  n_pc <- y |>
    group_by(County, PPNum) |>
    tally() |>
    pull(n)

  pp <- y |>
    filter(!is.na(V1)) |>
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
    filter(!is.na(V1)) |>
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
    n0 <- length(vec)
    y_removed[i, 1:n0] <- vec
  }

  return(
    list(
      y = y,                          # data with missing observations
      N = N_long,                     # latent property abundance
      M = M,                          # latent county abundance
      y_removed = y_removed,          # the number of pigs removed at each timestep
      size = size,                    # overdispersion parameter
      mean_r = mean_r,                # mean growth rate
      sigma_proc = sigma_proc,        # process error
      p = p,                          # observation rate by method,
      effort = effort,                # removal effort
      gamma = gamma,                  # area impacted
      method = method,                # trapping method
      county_area = county_area,      # the area of each county
      n_method = n_method,            # the number of methods
      n_property = n_property,        # the number of distinct properties
      n_timestep = n_timestep,        # the number of timesteps in each property
      n_pc = n_pc,                    # the number of properties sampled in each county in each PP
      pp = pp,                        # timestep index within each property
      n_county_units = nrow(M),       # the number of distinct county by PP units
      n_prop_county = n_prop_county1, # properties within counties that are sampled in each PP
      n_prop = n_prop                 # how many properties are sampled within counties that are sampled in each PP
    )
  )
}



