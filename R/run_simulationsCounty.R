library(config)
library(tidyverse)
library(nimble)
library(nimbleHMC)
library(lubridate)
library(targets)
library(coda)
library(splines)
source("R/functions_simulate2.R")
source("R/nimble_DynamicMultiMethodCounty.R")
source("R/functions_nimble.R")

n_county <- 3 # number of counties
n_property_per_county <- 5 # number of properties within each county
n_pp <- 10 # number of primary periods
n_passes <- 5 # number of removals per sampled primary period
n_method <- 3 # number of removal methods
n_property <- n_property_per_county * n_county # total number of properties
n_units <- n_property*n_pp # total number of spatio-temporal units
n_removals <- n_units*n_passes # total number of removal events

# simulate effort
effort_vec <- sample(1:3, n_removals, c(0.7, 0.25, 0.05), replace = TRUE)
effort <- matrix(effort_vec, n_units, n_passes)

# simulate which methods are used when
method_vec <- sample(1:n_method, n_removals, c(0.7, 0.2, 0.1), replace = TRUE)
method <- matrix(method_vec, n_units, n_passes)

# row index
unit <- matrix(seq_len(n_units), n_property, n_pp, byrow = TRUE)

# determine which primary periods are sampled for each property
sample_occ <- matrix(NA, n_property, n_pp)
for(i in seq_len(n_property)){
  n_samps <- round(runif(1, 4, n_pp - 2))
  s <- sample(n_pp, n_samps)
  sample_occ[i, sort(s)] <- sort(s)
}

mu_p <- -2.2
sigma_l <- 0.3
sigma_p <- 0.6
sigma_b <- 0.1

simulated_data <- simulate_swine(
  mu_p = mu_p,
  sigma_l = sigma_l,
  sigma_p = sigma_p,
  sigma_b = sigma_b,
  size = size,
  effort = effort,
  method = method,
  sample_occ = sample_occ,
  n_county = n_county,
  n_property_per_county = n_property_per_county,
  n_pp = n_pp,
  n_passes = n_passes,
  n_method = n_method
)

start_pp <- simulated_data$start_pp
end_pp <- simulated_data$end_pp
n_beta <- length(simulated_data$beta)

X <- matrix(NA, 1, n_beta)
for(i in seq_len(n_county)){
  Xdat <- cbind(1, start_pp[i]:end_pp[i]/end_pp[i])
  X <- rbind(X, cbind(1, bs(Xdat[,2], 5)))
}
X <- X[-1,]

# nimble data
data <- list(
  y = simulated_data$y_mat,
  y_removed = simulated_data$y_removed,
  effort = simulated_data$effort,
  gamma = simulated_data$effort,
  X = X
)

# nimble constants
constants <- list(
  n_units = nrow(simulated_data$y),
  n_units_all = nrow(X),
  n_passes = simulated_data$n_passes_y,
  n_prop = simulated_data$n_prop,
  n_county = n_county,
  n_methods = n_method,
  n_beta = n_beta,
  property = simulated_data$y$Property,
  PPNum = simulated_data$y$PPNum,
  method = simulated_data$method,
  pp = simulated_data$pp,
  county_x = simulated_data$M$County,
  pp_x = as.numeric(simulated_data$M$PPNum),
  PA = simulated_data$PA,
  PC = simulated_data$PC,
  start_pp = start_pp,
  end_pp = end_pp,
  proportion_surveyed = simulated_data$prop_surveyed,
  alpha_p = 1,
  beta_p = 1,
  alpha_b = 1,
  beta_b = 1,
  alpha_l = 1,
  beta_l = 1
)

county_areas <- simulated_data$take$area_county |> unique()
property_areas <- simulated_data$take$area_property |> unique()

inits <- function(){
  mu_p = rnorm(1, -2, 0.01)
  sigma_p = runif(1)
  logit_p = rnorm(n_method, mu_p, sigma_p)
  sigma_b = runif(1)
  beta = rnorm(constants$n_beta, 0, sigma_b)
  M = matrix(round(county_areas * runif(n_pp * n_county, 9, 10)), n_county, n_pp, byrow = TRUE)
  M_avail = round(M * simulated_data$prop_surveyed)
  list(
    z1 = log(M[,1]),
    mu_p = mu_p,
    logit_p = logit_p,
    sigma_p = sigma_p,
    sigma_l = runif(1),
    sigma_b = sigma_b,
    beta = beta,
    M = M,
    M_avail = M_avail
    # n = matrix(round(property_areas * runif(n_pp * n_property, 8, 9)), n_property, n_pp, byrow = TRUE)
  )
}

Rmodel <- nimbleModel(
  code = modelCode,
  constants = constants,
  inits = inits(),
  data = data
)

Rmodel$initializeInfo()

mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)
# addHMC(mcmcConf)
mcmcConf$addMonitors(c("xn", "M", "beta", "logit_p"))
Rmcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(Rmodel)
# Cmodel$calculate()
Cmcmc <- compileNimble(Rmcmc)

n_iter <- 100000
n_chains <- 3

samples <- runMCMC(
  Cmcmc,
  nburnin = n_iter/2,
  niter = n_iter,
  nchains = n_chains,
  # thin = 5,
  samplesAsCodaMCMC = TRUE
)

write_rds(
  list(
    samples = samples,
    simulated_data = simulated_data,
    constants = constants,
    data = data
  ),
  file = "out/CountySimulation.rds"
)

str(simulated_data)

ls <- read_rds("out/CountySimulation.rds")
attach(ls)

samples <- ls$samples
simulated_data <- ls[[2]]
ss <- as.matrix(samples)


hist(ss[,"sigma_l"])
abline(v = simulated_data$sigma_l, col = "red", lwd = 2)

hist(ss[,"sigma_p"])
abline(v = simulated_data$sigma_p, col = "red", lwd = 2)

hist(ss[,"mu_p"])
abline(v = simulated_data$mu_p, col = "red", lwd = 2)

par(mfrow = c(1, 3))
for(i in 1:simulated_data$n_method){
  hist(ilogit(ss[,paste0("logit_p[", i, "]")]), main = paste0("p[", i, "]"))
  abline(v = simulated_data$p[i], col = "red", lwd = 2)
}

hist(ss[,"sigma_b"])
abline(v = simulated_data$sigma_b, col = "red", lwd = 2)

for(i in 1:6){
  hist(ss[,paste0("beta[", i, "]")], main = paste0("beta[", i, "]"))
  abline(v = simulated_data$beta[i], col = "red", lwd = 2)
}

M <- simulated_data$M
par(mfrow = c(2, 5))
for(i in 1:30){
  node <- paste0("M[", M$County[i], ", ", M$PPNum[i], "]")
  hist(ss[,node], main = node)
  abline(v = M$County_abundance[i], col = "red", lwd = 2)
}

dev.off()
hist(ss[,"xn[40]"])




