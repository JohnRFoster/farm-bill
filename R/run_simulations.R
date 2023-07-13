library(config)
library(tidyverse)
library(nimble)
library(nimbleHMC)
library(lubridate)
library(targets)
library(coda)
source("R/functions_simulate.R")
source("R/nimble_amy.R")
source("R/functions_nimble.R")
source("R/fit_simulated_data.R")

out_dir <- file.path("out/simulationMultipleLatentPropertiesCountyRTimeVarying")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
run_dm <- FALSE
run_exponential <- TRUE

ex_diffs <- c("low", "high")

exponential_experiments <- expand_grid(
  r = ex_diffs,
  prob = "low"
)

dail_madsen_experiments <- expand_grid(
  phi = ex_diffs,
  rho = ex_diffs,
  prob = "low"
)

n_county <- 3 # number of counties
n_property_per_county <- 5 # number of properties within each county
n_pp <- 10 # number of primary periods
n_passes <- 5 # number of removals per sampled primary period
n_method <- 3 # number of removal methods
n_property <- n_property_per_county * n_county # total number of properties
n_units <- n_property*n_pp # total number of spatio-temporal units
n_removals <- n_units*n_passes # total number of removal events

# simulate effort
effort_vec <- sample(1:3, n_removals, c(0.7, 0.2, 0.1), replace = TRUE)
effort <- matrix(effort_vec, n_units, n_passes)

# simulate which methods are used when
method_vec <- sample(1:n_method, n_removals, c(0.7, 0.2, 0.1), replace = TRUE)
method <- matrix(method_vec, n_units, n_passes)

# simulate the area impacted by each removal event
gamma <- matrix(NA, n_units, n_passes)

# row index
unit <- matrix(seq_len(n_units), n_property, n_pp, byrow = TRUE)

# area impacted will differ by method within each property
for(i in 1:n_property){
  g <- runif(n_method, 0, 0.7)
  for(j in 1:n_pp){
    for(k in 1:n_method){
      gamma[unit[i, j], which(method[unit[i, j],] == k)] <- g[k]
    }
  }
}

# determine which primary periods are sampled for each property
sample_occ <- matrix(NA, n_property, n_pp)
for(i in seq_len(n_property)){
  n_samps <- round(runif(1, 4.6, n_pp - 1.6))
  s <- sample(n_pp, n_samps)
  sample_occ[i, sort(s)] <- sort(s)
}
# sample_occ[c(7:10, 17:20, 27:30),] <- NA


#### need to loop over overdispersion, area group, and likelihood for mcmc fits
sim_likelihood <- "nb" #c("deterministic", "poisson", "nb")
fit_likelihood <- "nb" #c("poisson", "nb")
area_group <- c("low", "med", "high")


r_c <- c(1, 1.05, 1.1)
phi <- 0.9
rho <- c(0.2, 0.4, 0.6)
sigma_p <- 0.1
sigma_proc <- 0.5
p <- c(0.2, 0.3, 0.4)
size <- c(0.5, 1, 2)
alpha_p <- rnorm(n_property, 0, sigma_p)

x1 <- matrix(rnorm(30), 3, 10)
x2 <- matrix(rnorm(30), 3, 10)

beta <- c(-0.2, 0.3)

exp(log(r_c[1]) + beta[1]*x1[1,] + beta[2]*x2[1,])
exp(log(r_c[2]) + beta[1]*x1[2,] + beta[2]*x2[2,])
exp(log(r_c[3]) + beta[1]*x1[3,] + beta[2]*x2[3,])

## simulated data sets (one fore each process; exponential and DM)
# 1. deterministic catch
# 2. poisson catch
# 3. negative binomial catch

i=j=k=1

if(run_exponential){

  for(i in seq_along(sim_likelihood)){
    sim_data <- simulate_swine(
      process = "exponential",
      likelihood = sim_likelihood[i],
      r_c = r_c,
      phi = NA,
      rho = NA,
      beta = beta,
      x1 = x1,
      x2 = x2,
      alpha_p = alpha_p,
      sigma_proc = sigma_proc,
      sigma_p = sigma_p,
      p = p,
      size = size,
      effort = effort,
      method = method,
      gamma = gamma,
      sample_occ = sample_occ
    )

    for(j in seq_along(fit_likelihood)){
      for(k in seq_along(area_group)){
        message("============================")
        message("simulted data: ", sim_likelihood[i])
        message("fit using: ", fit_likelihood[j])
        message("area surveyed: ", area_group[k], "\n")

        ag <- area_group[k]
        fit_data <- sim_data[[ag]]

        inits <- make_inits(fit_data, "exponential")


        out_name <- paste(
          "exponential",
          paste0("sim_", sim_likelihood[i]),
          paste0("fit_", fit_likelihood[j]),
          paste0("area", area_group[k]),
          sep = "_"
        )

        dest <- paste0(file.path(out_dir, out_name), ".rds")
        if(file.exists(dest)) next

        source("R/nimble_amy.R")
        samps <- fit_simulated_data(
          data_simulate = fit_data,
          inits = inits,
          lambda_config = "county",
          likelihood = fit_likelihood[j],
          process = "exponential"
        )

        write_rds(
          list(
            sim_likelihood = sim_likelihood[i],
            fit_likelihood = fit_likelihood[j],
            area_group = area_group[k],
            sim_data = fit_data,
            samps = samps
          ),
          file = dest
        )

      } # k; area_group
    } # j; fit_likelihood

  } # i; sim_likelihood

}



i=j=1


if(run_dm){
  for(i in seq_len(nrow(dail_madsen_experiments))){
    simulation_params <- dail_madsen_experiments |>
      slice(i)

    message("============================")
    message(i, " of ", nrow(dail_madsen_experiments))
    print(simulation_params)
    message("============================")

    sim_data <- simulate_swine(
      process = "dm",
      phi_ex = pull(simulation_params, phi),
      rho_ex = pull(simulation_params, rho),
      prob = pull(simulation_params, prob),
      effort = effort,
      method = method,
      gamma = gamma,
      sample_occ = sample_occ
    )

    inits <- make_inits(sim_data, "dm")

    for(j in seq_len(nrow(fits))){

      message("   =========================")
      message("   ", j, " of ", nrow(fits))
      print(slice(fits, j))
      message("   =========================")

      out_name <- paste(
        "dm",
        paste0("phi_", pull(simulation_params, phi)),
        paste0("rho_", pull(simulation_params, rho)),
        paste0("prob_", pull(simulation_params, prob)),
        paste0(fits$area_group[j]),
        paste0("overdispersion_", fits$overdispersion[j]),
        paste0("likelihood_", fits$likelihood[j]),
        sep = "_"
      )

      dest <- paste0(file.path(out_dir, out_name), ".rds")

      if(file.exists(dest)) next

      fit_data <- sim_data

      fit_data$M <- fit_data$M |>
        filter(area_group == fits$area_group[j],
               overdispersion == fits$overdispersion[j])

      samps <- fit_simulated_data(
        data_simulate = fit_data,
        inits = inits,
        intercept_only = TRUE,
        likelihood = fits$likelihood[j],
        process = "dm"
      )

      write_rds(
        list(
          simulation_params = simulation_params,
          fit_params = slice(fits, j),
          sim_data = fit_data,
          samps = samps
        ),
        file = dest
      )
    }

  }
}




