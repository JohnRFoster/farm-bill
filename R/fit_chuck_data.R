library(tidyverse)
library(nimble)
library(lubridate)
library(targets)
library(coda)

dir_data <- "data"
take_csv <- "all_chuck_data.csv"

dir_out <- "out"
dir_model <- "basic"

chuck_data <- read_csv(file.path(dir_data, take_csv))
take_df <- chuck_data |>
  arrange(Property, timestep) |>
  mutate(method_idx = as.numeric(as.factor(Methods)),
         property_idx = as.numeric(as.factor(Property)),
         county_idx = as.numeric(as.factor(County)))

area_csv <- "PropertyAreas.csv"
property_areas <- read_csv(file.path(dir_data, area_csv)) |>
  select(Property, `Area (km2)`) |>
  rename(area_km2 = `Area (km2)`)

passes <- take_df |>
  left_join(property_areas) |>
  # select(property_idx, timestep, Take, Effort, method_idx, area_km2) |>
  group_by(property_idx, timestep) |>
  mutate(pass = 1:n()) |>
  filter(max(pass) > 1) |>
  arrange(property_idx, timestep, pass) |>
  group_by(property_idx) |>
  mutate(timestep = timestep - min(timestep) + 1) |>
  filter(max(timestep) > 1) |>
  ungroup() |>
  arrange(property_idx, timestep, pass)


source("R/functions_data_processing.R")

data_ls <- get_site_data(passes)

ymat <- data_ls$ymat       # take, each row is a spatio-temporal unit (property x primary period), passes in columns
Areaper <- data_ls$Areaper # area of impact
Xd <- data_ls$Xd           # basis function coefficients
gbe <- data_ls$gbe         # effort





# the number of passes within each property x timestep
n_passes <- passes |>
  group_by(property_idx, timestep) |>
  summarise(n_passes = max(pass)) |>
  group_by(property_idx) |>
  mutate(timestep_idx = 1:n()) |>
  ungroup()

make_wide <- function(df, v){
  df |>
    select(property_idx, timestep, pass, .data[[v]]) |>
    group_by(property_idx, timestep) |>
    pivot_wider(names_from = pass,
                values_from = .data[[v]]) |>
    ungroup()
}

# take <- passes |> make_wide("Take")
take <- ymat |>
  select(-Site) |>
  mutate(property_idx = as.numeric(as.factor(Property))) |>
  rename(timestep = ts)
methods <- passes |> make_wide("method_idx") |>
  select(-property_idx, -timestep) |>
  as.matrix()
area <- passes |> make_wide("area_km2") |>
  select(-property_idx, -timestep) |>
  as.matrix()

y_sums <- take |>
  select(-property_idx, -timestep, -Property, -PPNum) |>
  rowSums(na.rm = TRUE)

ysum <- take |>
  select(property_idx, timestep) |>
  mutate(y_sums = y_sums) |>
  group_by(property_idx) |>
  mutate(y_cumsum = cumsum(y_sums)) |>
  mutate(pass = 1:n())

y_removed <- ysum|>
  select(-y_cumsum) |>
  pivot_wider(names_from = pass,
              values_from = y_sums,
              id_cols = property_idx) |>
  ungroup() |>
  select(-property_idx) |>
  as.matrix()


data <- list(
  y = take |> select(-timestep, -Property, -PPNum, -property_idx) |> as.matrix(),
  y_removed = y_removed,
  effort = gbe |> select(-ts, -Property) |> as.matrix(),
  gamma = Areaper |> select(-ts, -Property) |> as.matrix(),
  X = Xd |> select(-ts, -Property) |> as.matrix()
)

n_property <- n_passes |> pull(property_idx) |> unique() |> length()
n_timestep <- n_passes |> group_by(property_idx) |> tally() |> pull(n)
n_methods <- max(methods, na.rm = TRUE)

constants <- list(
  n_units = nrow(ymat),
  n_passes = n_passes |> pull(n_passes),
  n_property = n_property,
  n_timestep = n_timestep,
  n_methods = n_methods,
  property = n_passes |> pull(property_idx) |> as.factor() |> as.numeric(),
  timestep = n_passes |> pull(timestep_idx),
  method = methods
)

lambda_monitor <- n_monitor <- vector()
for(i in 1:n_property){
  for(t in 1:(n_timestep[i]-1)){
    ll <- paste0("log_lambda[", i, ", ", t, "]")
    lambda_monitor <- c(lambda_monitor, ll)
    nn <- paste0("n[", i, ", ", t, "]")
    n_monitor <- c(n_monitor, nn)
  }
  nn <- paste0("n[", i, ", ", n_timestep[i], "]")
  n_monitor <- c(n_monitor, nn)
}


inits <- function(){
  # z_init <- matrix(rpois(n_property * max(n_timestep), 50), n_property, max(n_timestep))
  list(
    n = y_removed + 10,
    log_mu1 = log(rowSums(data$y_removed, na.rm = TRUE)),
    # z = z_init,
    # mu_log = log(n_init + 1),
    log_lambda = log(matrix(runif(n_property * max(n_timestep), 0.01, 1), n_property, max(n_timestep))),
    mu_p = runif(constants$n_methods, 1, 2)
    # tau_p = runif(1, 0, 2)
  )
}

inits_test <- inits()
str(data)
str(constants)
str(inits_test)


source("R/nimble_amy.R")
Rmodel <- nimbleModel(
  code = modelCode,
  constants = constants,
  inits = inits_test,
  data = data
)
# warnings()
# check initialization
Rmodel$initializeInfo()

# default MCMC configuration
mcmcConf_conj <- configureMCMC(Rmodel, useConjugacy = TRUE)
# mcmcConf_no_conj <- configureMCMC(Rmodel, useConjugacy = FALSE)
mcmcConf <- mcmcConf_conj
mcmcConf$addMonitors(c("n"))
# these print statements will not display when running in parallel
mcmcConf$printMonitors()
# mcmcConf$printSamplers(byType = TRUE)

Rmcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc)

n_iter <- 2e5
n_chains <- 3

samples <- runMCMC(
  Cmcmc,
  nburnin = n_iter/2,
  niter = n_iter,
  nchains = n_chains,
  samplesAsCodaMCMC = TRUE
)


subset_out <- function(mcmc_out, state.col){

  mat2mcmc_list <- function(w) {
    temp <- list()
    chain.col <- which(colnames(w) == "CHAIN")
    for (i in unique(w[, "CHAIN"])) {
      temp[[i]] <- coda:::as.mcmc(w[w[, "CHAIN"] == i, -chain.col])
    }
    return(as.mcmc.list(temp))
  }

  mfit <- as.matrix(mcmc_out, chains = TRUE)
  pred.cols <- grep(state.col, colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out <- mat2mcmc_list(mfit[, c(chain.col, pred.cols)])
  return(out)
}

n <- subset_out(samples, "n")
mu_p <- subset_out(samples, "mu_p")
log_mu1 <- subset_out(samples, "log_mu1")
log_lambda <- subset_out(samples, "log_lambda")

# check convergence and effective sample size on specified node
subset_check_mcmc <- function(node){
  s <- mcmc[,grep(node, colnames(mcmc[[1]]), value = TRUE, fixed = TRUE)]
  df <- data.frame(
    psrf = gelman.diag(s, multivariate = FALSE)$psrf[,2], # checking the upper CI
    effective_samples = effectiveSize(s)
  )
  if(is.null(dim(s$chain1))) rownames(df) <- node
  return(df)
}

subset_check_burnin <- function(node, plot = FALSE){
  s <- mcmc[,grep(node, colnames(mcmc[[1]]), value = TRUE, fixed = TRUE)]
  if(plot){
    GBR <- gelman.plot(s)
  } else {
    ff <- tempfile()
    png(filename = ff)
    GBR <- gelman.plot(s)
    dev.off()
    unlink(ff)
  }
  shrink <- GBR$shrink[, , 2]
  if(is.null(dim(s$chain1))){
    burnin <- GBR$last.iter[tail(which(shrink > 1.1, 1), 1) + 1]
  } else {
    burnin <- GBR$last.iter[tail(which(apply(shrink > 1.1, 1, any)), 1) + 1]
  }
  return(burnin)
}

nodes1 <- c("mu_p", "log_mu1")
nodes2 <- c(lambda_monitor, n_monitor)
nodes_check <- c(nodes1, lambda_monitor)

mcmc <- samples
checks <- map_dfr(lapply(nodes_check, subset_check_mcmc), as.data.frame)
checks

burnin <- map_dbl(lapply(nodes_check, subset_check_burnin), as.vector)
burnin
max(burnin)
samples_burn <- window(samples, start = max(burnin))

method_lookup <- passes |>
  select(Methods, method_idx) |>
  distinct()

property_idxs <- passes |>
  select(Property, property_idx) |>
  distinct() |>
  mutate(property_idx = 1:n())

property_lookup <- passes |>
  select(Property, PPNum, pp_start_date, pp_end_date) |>
  distinct() |>
  group_by(Property) |>
  arrange(Property, PPNum) |>
  mutate(timestep = 1:n()) |>
  ungroup() |>
  left_join(property_idxs)


path <- file.path(dir_out, dir_model)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

write_rds(
  list(
    samples = samples_burn,
    nimble_data = data,
    nimble_constants = constants,
    method_lookup = method_lookup,
    property_lookup = property_lookup,
    nodes1 = nodes1,
    nodes2 = nodes2
  ),
  file.path(path, "samples.rds")
)








