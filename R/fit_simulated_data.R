library(config)
library(tidyverse)
library(nimble)
library(lubridate)
library(targets)
library(coda)

n_county <- 3
n_property_per_county <- 5
n_pp <- 10
n_passes <- 5
likelihood <- "nb"

source("R/functions_simulate.R")
data_simulate <- simulate_swine(
  n_county = n_county,
  n_property_per_county = n_property_per_county,
  n_pp = n_pp,
  n_passes = n_passes,
  likelihood = likelihood
)

counts <- data_simulate$y
keep_index <- which(!is.na(counts$V1))

C <- counts |>
  slice(keep_index) |>
  group_by(Property) |>
  mutate(timestep = 1:n()) |>
  ungroup()

y <- C |>
  select(-County, -Property, -Property_county, -PPNum, -prop_area, -timestep) |>
  as.matrix()


data <- list(
  y = y,
  y_removed = data_simulate$y_removed,
  effort = data_simulate$effort,
  gamma = data_simulate$gamma
)

constants <- list(
  n_units = nrow(y),
  n_units_all = nrow(counts),
  n_passes = rep(n_passes, nrow(y)),
  n_property = data_simulate$n_property,
  n_county = n_county,
  n_timestep = data_simulate$n_timestep,
  n_methods = data_simulate$n_method,
  # n_beta = n_beta,
  n_county_units = data_simulate$n_county_units,
  n_pc = data_simulate$n_pc,
  property = C$Property,
  county = data_simulate$M$County,
  countyN = C$County,
  timestep = C$timestep,
  property_x = counts$Property,
  timestep_x = counts$PPNum,
  method = data_simulate$method[keep_index, ],
  pp = data_simulate$pp,
  # nt = data_simulate$n_pp,
  pp_c = data_simulate$M$PPNum,
  # maxp = max(data_simulate$M$PPNum),
  n_prop_cnty = data_simulate$n_prop_county,
  n_prop = data_simulate$n_prop,
  cnty_property = C$Property_county,
  PPNum = C$PPNum,
  sum_prop_area = data_simulate$M$total_area_sampled,
  county_area = data_simulate$county_area
  # k_mu = rep(0, 4),
  # k_cov = diag(1, 4, 4)
)

include_county <- TRUE
include_proc <- TRUE
include_all_pp <- TRUE
intercept_only <- TRUE
include_beta <- if_else(intercept_only, FALSE, TRUE)
nb_likelihood <- if_else(likelihood == "nb", TRUE, FALSE)
pois_likelihood <- if_else(likelihood == "poisson", TRUE, FALSE)

N_mat_na <- data_simulate$N |>
  pivot_wider(names_from = PPNum,
              values_from = Property_abundance) |>
  select(sprintf("%02d", as.integer(1:n_pp))) |>
  as.matrix()
n_mat <- matrix(NA, nrow(N_mat_na), ncol(N_mat_na))
for(i in 1:nrow(N_mat_na)){
  vec <- N_mat_na[i, which(!is.na(N_mat_na[i,]))]
  n0 <- length(vec)
  n_mat[i, 1:n0] <- vec
}

inits <- function(){
  list(
    n = n_mat + 1,
    log_mu_proc = log(n_mat + 1),
    M = data_simulate$M$County_abundnace,
    log_mu1 = log(n_mat[,1] + 1),
    # log_lambda = log(matrix(runif(n_property * max(n_timestep), 0.01, 1), n_property, max(n_timestep))),
    mean_r = log(runif(1, 0.8, 2)),
    tau_p = 1/runif(1, 0.1, 1)^2,
    # beta = matrix(rnorm(n_property * n_beta), n_property, n_beta),
    mu_p = boot::logit(runif(3, 0.2, 0.6)),
    size = runif(3, 0.01, 5)
    # beta_k = rnorm(4)
  )
}

inits_test <- inits()
str(data)
str(constants)
# str(inits_test)

# options(error = recover)
source("R/nimble_amy.R")
source("R/functions_nimble.R")
Rmodel <- nimbleModel(
  code = modelCode,
  constants = constants,
  inits = inits_test,
  data = data
)

Rmodel$initializeInfo()

# default MCMC configuration
mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)
mcmcConf$addMonitors(c("n", "M"))
# mcmcConf$removeSampler("n")
Rmcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc)

n_iter <- 2e4
n_chains <- 3

samples <- runMCMC(
  Cmcmc,
  nburnin = n_iter/2,
  niter = n_iter,
  nchains = n_chains,
  samplesAsCodaMCMC = TRUE
)

samps <- as.matrix(samples)

x_lab <- function(obs, pred){
  paste("Observed:", round(obs, 2),
        "\nPredicted:", round(mean(pred), 2))
}

hist_check <- function(obs, pred, main){
  hist(pred, xlab = x_lab(obs, pred), main = main)
  abline(v = obs, col = "red", lwd = 2)
}


par(mfrow = c(1, 1))
hist_check(data_simulate$mean_r, exp(samps[,"mean_r"]), "Mean r")

par(mfrow = c(1, 1))
hist_check(data_simulate$sigma_proc, 1/sqrt(samps[,"tau_p"]), "Process error")

par(mfrow = c(2, 3))
hist_check(data_simulate$p[1], boot::inv.logit(samps[,"mu_p[1]"]), "p[1]")
hist_check(data_simulate$p[2], boot::inv.logit(samps[,"mu_p[2]"]), "p[2]")
hist_check(data_simulate$p[3], boot::inv.logit(samps[,"mu_p[3]"]), "p[3]")
hist_check(data_simulate$size[1], samps[,"size[1]"], "size[1]")
hist_check(data_simulate$size[2], samps[,"size[2]"], "size[2]")
hist_check(data_simulate$size[3], samps[,"size[3]"], "size[3]")

par(mfrow = c(5, 5))
N <- data_simulate$N |>
  group_by(Property) |>
  mutate(timestep = 1:n()) |>
  ungroup()
mu_all <- tibble()
for(i in seq_len(nrow(N))){
  n_node <- paste0("n[", N$Property[i], ", ", N$timestep[i], "]")
  hist_check(N$Property_abundance[i], samps[,n_node],
             paste("Property", N$Property[i], "PP", N$PPNum[i]))
  mu_i <- tibble(
    lower = quantile(samps[,n_node], 0.025),
    med = quantile(samps[,n_node], 0.5),
    upper = quantile(samps[,n_node], 0.975),
    Property = N$Property[i],
    timestep = N$timestep[i],
    PPNum = N$PPNum[i]
  )
  mu_all <- bind_rows(mu_all, mu_i)
}

N_samps <- left_join(N, mu_all)

N_samps |>
  # mutate(property_area = as.character(property_area)) |>
  ggplot() +
  aes(x = Property_abundance, y = med) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed",
       y = "Median predicted") +
  theme_bw()



par(mfrow = c(5, 6))
M <- data_simulate$M
for(i in seq_len(nrow(M))){
  m_node <- paste0("M[", i, "]")
  hist_check(M$County_abundnace[i], samps[, m_node],
             paste("County", M$County[i], "PP", M$PPNum[i]))
}
