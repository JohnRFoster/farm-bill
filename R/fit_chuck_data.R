library(config)
library(tidyverse)
library(nimble)
library(lubridate)
library(targets)
library(coda)

config_name <- "dynamic_multi_method_poisson"

Sys.setenv(R_CONFIG_ACTIVE = config_name)
config <- config::get()
dir_out <- config$dir_out
phi_config <- config$phi_config
rho_config <- config$rho_config
lambda_config <- config$lambda_config
dir_model <- config_name
nb_likelihood <- if_else(config$nb_likelihood, TRUE, FALSE)
pois_likelihood <- !nb_likelihood
phi_static <- config$phi_static
rho_static <- config$rho_static
lambda_static <- config$lambda_static

phi_data <- if_else(phi_config == "data", TRUE, FALSE)
rho_data <- if_else(rho_config == "data", TRUE, FALSE)
phi_beta <- if_else(phi_config %in% c("basis", "covs"), TRUE, FALSE)
rho_beta <- if_else(rho_config %in% c("basis", "covs"), TRUE, FALSE)

phi_hb <- if_else(phi_config == "hb", TRUE, FALSE)

lambda_data <- if_else(lambda_config == "data", TRUE, FALSE)
lambda_beta <- if_else(lambda_config %in% c("basis", "covs"), TRUE, FALSE)

dir_data <- "data"

take_csv <- "all_chuck_data.csv"
chuck_data <- read_csv(file.path(dir_data, take_csv))
take_df <- chuck_data |>
  mutate(method_idx = as.numeric(as.factor(Methods)),
         property_idx = as.numeric(as.factor(Property)),
         county_idx = as.numeric(as.factor(County)),
         Date = ymd(Date)) |>
  arrange(Property, Date, timestep)

area_csv <- "PropertyAreas_2kmBuffer.csv"
property_areas <- read_csv(file.path(dir_data, area_csv)) |>
  select(Property, `With 2km Buffer`) |>
  rename(area_km2 = `With 2km Buffer`)

dist_cum <- function(var){
  sapply(seq_along(var), function(x) length(unique(head(var, x))))
}

source("R/functions_data_processing.R")

# take_df$SPNum <- secondaryperiod(1, take_df$Date, take_df)

passes <- take_df |>
  left_join(property_areas) |>
  # select(property_idx, timestep) |>
  group_by(property_idx, timestep) |>
  mutate(pass = 1:n()) |>
  filter(max(pass) > 1) |>
  arrange(property_idx, timestep, pass) |>
  ungroup() |>
  group_by(property_idx) |>
  mutate(timestep = timestep - min(timestep) + 1,
         Effort = if_else(is.na(Effort), 1, Effort)) |>
  filter(max(timestep) > 1) |>
  ungroup() |>
  mutate(property_idx = as.numeric(as.factor(Property))) |>
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
  mutate(n_idx = 1:n()) |>
  select(County, PPNum, cnty_property, n_idx) |>
  pivot_wider(names_from = cnty_property,
              values_from = n_idx,
              id_cols = c(County, PPNum)) |>
  arrange(County, PPNum) |>
  select(-County, -PPNum) |>
  as.matrix()

n_prop <- apply(n_prop_county, 1, function(x) length(which(!is.na(x))))

n_prop_county1 <- matrix(NA, nrow(n_prop_county), max(n_prop))
for(i in 1:nrow(n_prop_county)){
  vec <- n_prop_county[i, which(!is.na(n_prop_county[i,]))]
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
              values_from = timestep,
              id_cols = Property)

data_ls <- suppressWarnings(get_site_data(passes, 1))

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


n_beta <- 6

data <- list(
  y = y,
  y_removed = y_removed,
  effort = gbe |> select(-ts, -Property) |> as.matrix(),
  gamma = Areaper |> select(-ts, -Property) |> as.matrix(),
  # X = Xd |> select(-ts, -Property) |> as.matrix(),
  Xall = Xd_all |> select(-ts, -Property) |> as.matrix()
)

property <- take |> pull(Property) |> as.factor() |> as.numeric()
n_property <- max(property)
county <- sum_prop_area |> pull(County) |> as.factor() |> as.numeric()
n_county <- max(county)
n_timestep <- take |> group_by(Property) |> tally() |> pull(n)

methods_mat <- methods |>
  select(-ts, -Property) |>
  as.matrix()
n_methods <- max(methods_mat, na.rm = TRUE)

phi_prior <- c(3.294963, 0.4520253)

# par(mfrow = c(2,2))
# n <- 10000
# xlim = c(0.5, 1)
# var(ilogit(rnorm(n, phi_prior[1], phi_prior[2])))
# var(ilogit(rnorm(n, phi_prior[1], phi_prior[2]*2)))
# var(ilogit(rnorm(n, phi_prior[1], phi_prior[2]*5)))
# var(ilogit(rnorm(n, phi_prior[1], phi_prior[2]*12)))


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
  n_monitor = nrow(take),
  property = property,
  county = county,
  countyN = all_county_units |> pull(county_idx),
  timestep = take$timestep,
  property_x = Xd_all |> pull(Property) |> as.factor() |> as.numeric(),
  timestep_x = Xd_all |> pull(ts),
  method = methods_mat,
  pp = pp |> select(-Property) |> as.matrix(), # will need to reparamaterize to fit nimble code
  nt = max(n_timestep),
  pp_c = pp_c,
  maxp = max(pp_c),
  M_lookup = n_prop_county1,
  n_prop = n_prop,
  cnty_property = all_county_units |> pull(cnty_property),
  PPNum = all_county_units |> pull(PPNum),
  sum_prop_area = sum_prop_area |> pull(sum_area),
  county_area = county_areas$area_km2,
  phi_prior = phi_prior,
  property_m = take$property_idx,
  timestep_m = take$timestep,
  alpha_p = 1,
  beta_p = 1,
  alpha_b = 1,
  beta_b = 1,
  alpha_l = 1,
  beta_l = 1
)

inits <- function(){
  mu_p <- rnorm(1)
  sigma_p <- runif(1)
  sigma_b <- runif(1)
  sigma_l <- runif(1)
  logit_p <- rnorm(n_methods, mu_p, sigma_p)
  n_init <- round(y_removed / min(boot::inv.logit(logit_p)))
  list(
    n = n_init,
    z1 = n_init[,1],
    mu_p = mu_p,
    logit_p = logit_p,
    sigma_p = sigma_p,
    sigma_l = sigma_l,
    sigma_b = runif(1),
    beta = rnorm(n_beta, 0, sigma_b),
    size = runif(n_county, 0, 5)
  )
}

inits_test <- inits()
str(data)
str(constants)
str(inits_test)

# options(error = recover)
source("R/functions_nimble.R")
source("R/nimble_DynamicMultiMethod.R")
Rmodel <- nimbleModel(
  code = modelCode,
  constants = constants,
  # inits = inits_test,
  inits = inits(),
  data = data
)
# Rmodel$n
# Rmodel$n - y_removed

# warnings()
# check initialization
Rmodel$initializeInfo()

Cmodel <- compileNimble(Rmodel)
Cmodel$calculate()

# default MCMC configuration
mcmcConf <- configureMCMC(Cmodel, useConjugacy = TRUE)

mcmcConf$addMonitors(c("xn", "M", "beta", "logit_p", "N_mu"))
# if(phi_hb) mcmcConf$addMonitors(c("mu_phi_prop"))

# mcmcConf$removeSampler("n")

Rmcmc <- buildMCMC(mcmcConf)
Cmcmc <- compileNimble(Rmcmc)

n_iter <- 120000
n_chains <- 3

samples <- runMCMC(
  Cmcmc,
  nburnin = 0,
  thin = 1,
  niter = n_iter,
  nchains = n_chains,
  samplesAsCodaMCMC = TRUE
)

nodes_check <- config$nodes1
all_nodes <- colnames(samples[[1]])
j <- unlist(lapply(nodes_check, function(x) grep(x, all_nodes)))

params <- samples[,j]

dest <- file.path(dir_out, dir_model, if_else(config$nb_likelihood, "negbinom", "poissonGamma"))
print(dest)
if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

message("Creating traceplots...")
png(filename = file.path(dest, "mcmcTimeseries%03d.png"))
plot(params)
dev.off()
message("  done")

message("Calculating PSRF...")
psrf <- gelman.diag(params, multivariate = FALSE)
print(psrf)

message("Calculating effective samples...")
effective_samples <- effectiveSize(params)
print(effective_samples)

message("Calculating burnin...")
ff <- tempfile()
png(filename = ff)
GBR <- gelman.plot(params)
dev.off()
unlink(ff)

burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
print(burnin)

if(is.na(burnin)) burnin <- 85000


samples_burn <- window(samples, start = max(burnin, na.rm = TRUE))

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

print(warnings())

write_rds(
  list(
    samples = samples_burn,
    nimble_data = data,
    nimble_constants = constants,
    method_lookup = method_lookup,
    property_lookup = property_lookup,
    all_county_units = all_county_units,
    county_level_lookup = sum_prop_area |> mutate(m_idx = 1:n()),
    nodes1 = config$nodes1
    # nodes2 = nodes2
  ),
  file.path(dest, "samples.rds")
)






