library(config)
library(tidyverse)
library(nimble)
library(lubridate)
library(targets)
library(coda)

config_name <- "county_secondary_all_pp_pkt"

Sys.setenv(R_CONFIG_ACTIVE = config_name)
config <- config::get()
dir_out <- config$dir_out
dir_model <- config$dir_model
include_beta <- config$include_beta
intercept_only <- config$intercept_only
include_proc <- config$include_proc
include_county <- config$include_county
include_all_pp <- config$include_all_pp
not_include_all_pp <- !include_all_pp

dir_data <- "data"

take_csv <- "all_chuck_data.csv"
chuck_data <- read_csv(file.path(dir_data, take_csv))
take_df <- chuck_data |>
  mutate(method_idx = as.numeric(as.factor(Methods)),
         property_idx = as.numeric(as.factor(Property)),
         county_idx = as.numeric(as.factor(County)),
         Date = ymd(Date)) |>
  arrange(Property, Date, timestep)

area_csv <- "PropertyAreas.csv"
property_areas <- read_csv(file.path(dir_data, area_csv)) |>
  select(Property, `Area (km2)`) |>
  rename(area_km2 = `Area (km2)`)

dist_cum <- function(var){
  sapply(seq_along(var), function(x) length(unique(head(var, x))))
}

source("R/functions_data_processing.R")

# take_df$SPNum <- secondaryperiod(1, take_df$Date, take_df)

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
  arrange(property_idx, Date, timestep, pass)

all_county_units <- passes |>
  left_join(property_areas) |>
  select(Property, County, county_idx, PPNum, area_km2) |>
  distinct() |>
  group_by(Property) |>
  mutate(timestep = PPNum - min(PPNum) + 1) |>
  ungroup() |>
  group_by(County) |>
  mutate(cnty_property = dist_cum(Property)) |>
  ungroup()

sum_prop_area <- all_county_units |>
  group_by(County, PPNum) |>
  summarise(sum_area = sum(area_km2)) |>
  ungroup()

n_prop_county <- all_county_units |>
  select(County, PPNum, cnty_property) |>
  pivot_wider(names_from = cnty_property,
              values_from = cnty_property,
              id_cols = c(County, PPNum)) |>
  arrange(County, PPNum) |>
  select(-County, -PPNum) |>
  as.matrix()

n_prop_county1 <- matrix(NA, nrow(n_prop_county), ncol(n_prop_county))
n_prop <- rep(NA, nrow(n_prop_county))
for(i in 1:nrow(n_prop_county)){
  vec <- n_prop_county[i, which(!is.na(n_prop_county[i,]))]
  n_prop[i] <- length(vec)
  n_prop_county1[i, 1:n_prop[i]] <- vec
}

pp_c <- sum_prop_area |> pull(PPNum)

n_pc <- all_county_units |>
  group_by(County) |>
  summarise(n = max(cnty_property)) |>
  ungroup() |>
  pull(n)

county_names <- all_county_units |>
  pull(County) |>
  unique()

county_vec <- paste0(county_names, ", SC")

county_areas_xl <- readxl::read_xls("../pigs-statistical/data/counties/area.xls")
county_areas <- county_areas_xl |>
  select(Areaname, LND110210D) |>
  filter(Areaname %in% county_vec) |>
  mutate(area_km2 = LND110210D * 2.589)


pp <- all_county_units |>
  group_by(Property) |>
  mutate(timestep_idx = 1:n()) |>
  ungroup() |>
  select(Property, timestep, timestep_idx) |>
  pivot_wider(names_from = timestep_idx,
              values_from = timestep)

data_ls <- get_site_data(passes)

ymat <- data_ls$ymat         # take, each row is a spatio-temporal unit (property x primary period), passes in columns
#ymat[is.na(ymat)] <- -1
Areaper <- data_ls$Areaper   # area of impact
Xd <- data_ls$Xd             # basis function coefficients
Xd_all <- data_ls$Xd_all     # basis function coefficients for all primary periods
n_pp_all <- data_ls$n_pp_all # the number of primary periods for each property from first to last removal
gbe <- data_ls$gbe           # effort
methods <- data_ls$method |>    # methods
  pivot_longer(cols = -c(ts, Property)) |>
  mutate(value = if_else(value == 0, NA, value),
         method_fac = as.numeric(as.factor(value))) |>
  select(-value) |>
  pivot_wider(names_from = name,
              values_from = method_fac)


make_wide <- function(df, v){
  df |>
    select(property_idx, timestep, pass, all_of(v)) |>
    group_by(property_idx, timestep) |>
    pivot_wider(names_from = pass,
                values_from = all_of(v)) |>
    ungroup()
}

# take <- passes |> make_wide("Take")
take <- ymat |>
  select(-Site) |>
  mutate(property_idx = as.numeric(as.factor(Property))) |>
  rename(timestep = ts)

y <- take |>
  select(-timestep, -Property, -PPNum, -property_idx) |>
  as.matrix()

# the number of passes within each property x timestep
n_passes <- apply(y, 1, function(x) max(which(x != 0)))

area <- take |>
  left_join(property_areas) |>
  select(Property, area_km2) |>
  distinct() |>
  arrange(Property)

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

n_beta <- if_else(intercept_only, 1, 6)

data <- list(
  y = y,
  y_removed = y_removed,
  effort = gbe |> select(-ts, -Property) |> as.matrix(),
  gamma = Areaper |> select(-ts, -Property) |> as.matrix(),
  X = Xd |> select(-ts, -Property) |> as.matrix(),
  Xall = Xd_all |> select(-ts, -Property) |> as.matrix(),
  Xk = cbind(rep(1, 3), matrix(rnorm(9), 3, 3))
)

property <- take |> pull(Property) |> as.factor() |> as.numeric()
n_property <- max(property)
county <- sum_prop_area |> pull(County) |> as.factor() |> as.numeric()
n_county <- max(county)
n_timestep <- take |> group_by(Property) |> tally() |> pull(n)

methods <- methods |>
  select(-ts, -Property) |>
  as.matrix()
n_methods <- max(methods, na.rm = TRUE)

constants <- list(
  n_units = nrow(ymat),
  n_units_all = nrow(Xd_all),
  n_passes = n_passes,
  n_property = n_property,
  n_county = n_county,
  n_timestep = n_timestep,
  n_methods = n_methods,
  n_beta = n_beta,
  n_county_units = nrow(sum_prop_area),
  n_pc = n_pc,
  property = property,
  county = county,
  countyN = all_county_units |> pull(county_idx),
  timestep = take$timestep,
  property_x = Xd_all |> pull(Property) |> as.factor() |> as.numeric(),
  timestep_x = Xd_all |> pull(ts),
  method = methods,
  pp = pp |> select(-Property) |> as.matrix(),
  nt = max(n_timestep),
  pp_c = pp_c,
  maxp = max(pp_c),
  n_prop_cnty = n_prop_county1,
  n_prop = n_prop,
  cnty_property = all_county_units |> pull(cnty_property),
  PPNum = all_county_units |> pull(PPNum),
  sum_prop_area = sum_prop_area |> pull(sum_area),
  county_area = county_areas$area_km2,
  k_mu = rep(0, 4),
  k_cov = diag(1, 4, 4)
)

lambda_monitor <- n_monitor <- m_monitor <- vector()
for(i in 1:n_property){
  for(t in 1:n_timestep[i]){
    ll <- paste0("log_lambda[", i, ", ", t, "]")
    lambda_monitor <- c(lambda_monitor, ll)
    nn <- paste0("n[", i, ", ", t, "]")
    n_monitor <- c(n_monitor, nn)
  }
}


inits <- function(){
  list(
    n = y_removed + rpois(1, 10),
    # M = (y_removed + rpois(1, 10))*rpois(1, 3),
    log_mu1 = log(rowSums(data$y_removed, na.rm = TRUE)),
    log_mu_proc = log(y_removed + 10),
    log_lambda = log(matrix(runif(n_property * max(n_timestep), 0.01, 1), n_property, max(n_timestep))),
    tau_p = runif(1, 0, 2),
    beta = matrix(rnorm(n_property * n_beta), n_property, n_beta),
    mu_p = matrix(runif(n_methods*max(n_timestep), 0, 1), n_methods, max(n_timestep)),
    # size = runif(3, 0.01, 2),
    beta_k = rnorm(4)
  )
}

inits_test <- inits()
str(data)
str(constants)
str(inits_test)

# options(error = recover)
source("R/nimble_amy.R")
source("R/function_nimble.R")
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
mcmcConf$addMonitors(c("n", "M", "size", "N_disp", "size"))
# these print statements will not display when running in parallel
# mcmcConf$printMonitors()
# mcmcConf$printSamplers(byType = TRUE)

Rmcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc)

n_iter <- 1e5
n_chains <- 3

samples <- runMCMC(
  Cmcmc,
  nburnin = n_iter/2,
  niter = n_iter,
  nchains = n_chains,
  samplesAsCodaMCMC = TRUE
)

ss <- as.matrix(samples)
hist(ss[,"N_disp[3]"])
hist(ss[,"M[3]"])
hist(ss[,"beta_k[1]"])
hist(ss[,"beta_k[2]"])
hist(ss[,"size[1]"])
hist(ss[,"size[2]"])
hist(ss[,"size[3]"])


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
  if(all(shrink < 1.1)){
    burnin <- 1
  } else {
    if(is.null(dim(s$chain1))){
      burnin <- GBR$last.iter[tail(which(shrink > 1.1), 1) + 1]
    } else {
      burnin <- GBR$last.iter[tail(which(apply(shrink > 1.1, 1, any)), 1) + 1]
    }
  }
  return(burnin)
}

nodes1 <- config$nodes1
nodes2 <- n_monitor
if(config_name %in% c("default", "basic_process")){
  nodes2 <- c(nodes2, lambda_monitor)
}

nodes_check <- c(nodes1, nodes2)

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
print(path)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

write_rds(
  list(
    samples = samples_burn,
    nimble_data = data,
    nimble_constants = constants,
    method_lookup = method_lookup,
    property_lookup = property_lookup,
    all_county_units = all_county_units,
    county_level_lookup = sum_prop_area |> mutate(m_idx = 1:n()),
    nodes1 = nodes1,
    nodes2 = nodes2
  ),
  file.path(path, "samples.rds")
)







