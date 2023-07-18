library(config)
library(tidyverse)
library(nimble)
library(lubridate)
library(targets)
library(coda)
library(reshape2)
library(splines)

fit_simulated_data <- function(data_simulate, inits, lambda_config, likelihood, process){

  latent_property <- c(7:10, 17:20, 27:30)

  counts <- data_simulate$y

  keep1 <- which(!is.na(counts$V1))
  keep2 <- which(counts$Property %in% latent_property)

  C <- counts |>
    slice(sort(c(keep1, keep2))) |>
    group_by(Property) |>
    mutate(timestep = 1:n()) |>
    ungroup()

  y_attributes <- C |> filter(!is.na(V1))

  y <- y_attributes |>
    select(-County, -Property, -Property_county, -PPNum, -timestep) |>
    as.matrix()

  y_removed <- data_simulate$y_removed
  y_removed <- y_removed[-latent_property,]

  data <- list(
    y = y,
    y_removed = y_removed,
    effort = data_simulate$effort,
    gamma = data_simulate$gamma,
    x1 = data_simulate$x1,
    x2 = data_simulate$x2
  )

  M_groups <- C |>
    select(County, PPNum) |>
    distinct() |>
    mutate(m_id = 1:n())

  M_join <- C |>
    mutate(n_id = 1:n()) |>
    select(County, PPNum, n_id, Property_county) |>
    left_join(M_groups)

  M_wide <- M_join |>
    pivot_wider(id_cols = m_id,
                values_from = n_id,
                names_from = Property_county) |>
    select(-m_id) |>
    as.matrix()

  M_lookup <- matrix(NA, nrow(M_wide), ncol(M_wide))
  for(i in 1:nrow(M_wide)){
    vec <- which(!is.na(M_wide[i,]))
    M_lookup[i,1:length(vec)] <- M_wide[i, vec]
  }

  n_prop <- apply(M_lookup, 1, function(x) max(which(!is.na(x))))

  constants <- list(
    n_units = nrow(y),
    n_units_all = nrow(counts),
    n_passes = rep(5, nrow(y)),
    n_property = data_simulate$n_property,
    n_latent_property = length(latent_property),
    n_latent_timestep = 10,
    n_county = 3,
    n_timestep = data_simulate$n_timestep,
    n_methods = data_simulate$n_method,
    # n_beta = n_beta,
    n_county_units = nrow(data_simulate$M),
    # n_pc = data_simulate$n_pc,
    n_monitor = nrow(C),
    property = y_attributes$Property,
    propertyN = y_attributes$Property |> unique(),
    n_propertyN = length(y_attributes$Property |> unique()),
    latent_property = latent_property,
    # latent_timestep = 1:10,
    county = data_simulate$M$County |> as.factor() |> as.numeric(),
    countyN = y_attributes$County,
    timestep = y_attributes$timestep,
    property_x = counts$Property,
    timestep_x = counts$PPNum,
    county_x = counts$County,
    method = data_simulate$method,
    pp = data_simulate$pp,
    # nt = data_simulate$n_pp,
    # pp_c = data_simulate$M$PPNum,
    # maxp = max(data_simulate$M$PPNum),
    n_prop = n_prop,
    M_lookup = M_lookup,
    # cnty_property = C$Property_county,
    # PPNum = C$PPNum,
    property_m = C$Property,
    timestep_m = C$timestep
    # sum_prop_area = data_simulate$M$total_area_sampled,
    # county_area = data_simulate$M$county_areas |> unique()
    # k_mu = rep(0, 4),
    # k_cov = diag(1, 4, 4)
  )


  lambda_static <- if_else(lambda_config == "static", TRUE, FALSE)
  lambda_data <- if_else(lambda_config == "data", TRUE, FALSE)
  lambda_beta <- if_else(lambda_config == "beta", TRUE, FALSE)
  lambda_county <- if_else(lambda_config == "county", TRUE, FALSE)

  include_all_pp <- TRUE
  # include_beta <- if_else(intercept_only, FALSE, TRUE)
  nb_likelihood <- if_else(likelihood == "nb", TRUE, FALSE)
  pois_likelihood <- if_else(likelihood == "poisson", TRUE, FALSE)
  exponential_process <- if_else(process == "exponential", TRUE, FALSE)
  dm_process <- if_else(process == "dm", TRUE, FALSE)

  # nimbleOptions(MCMCusePredictiveDependenciesInCalculations = FALSE)

  Rmodel <- nimbleModel(
    code = modelCode,
    constants = constants,
    inits = inits(),
    data = data
  )

  Rmodel$initializeInfo()

  # default MCMC configuration
  mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)
  # addHMC(mcmcConf)
  mcmcConf$addMonitors(c("nm", "M", "N_mu", "alpha_p"))
  # if(process == "dm"){
  #   mcmcConf$removeSampler(c("log_rho", "mu_phi"))
  #   mcmcConf$addSampler(c("log_rho", "mu_phi"), "AF_slice")
  #   # mcmcConf$addSampler("mu_phi", "slice")
  # }
  # mcmcConf$removeSampler("n")

  Rmcmc <- buildMCMC(mcmcConf)
  Cmodel <- compileNimble(Rmodel)
  # Cmodel$calculate()
  Cmcmc <- compileNimble(Rmcmc)

  n_iter <- 75000
  n_chains <- 3

  samples <- runMCMC(
    Cmcmc,
    nburnin = n_iter/2,
    niter = n_iter,
    nchains = n_chains,
    # thin = 5,
    samplesAsCodaMCMC = TRUE
  )
  # plot(samples)
  return(as.matrix(samples))
}



