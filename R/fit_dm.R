library(config)
library(tidyverse)
library(nimble)
library(lubridate)
library(targets)
library(coda)

config_name <- "mis_multi_method_poisson"

Sys.setenv(R_CONFIG_ACTIVE = config_name)
config <- config::get()
dir_model <- config_name
likelihood <- config$likelihood

dir_out <- "out"
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

# mean litter size year from VerCauteren et al. 2019 pg 63
data_litter_size <- c(5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
                      5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
                      4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4)

data_litter_size <- round(data_litter_size)

passes <- take_df |>
  left_join(property_areas) |>
  # select(-timestep) |>
  group_by(property_idx, PPNum) |>
  mutate(pass = 1:n()) |>
  filter(max(pass) > 1) |>
  ungroup() |>
  mutate(PPNum = PPNum - min(PPNum) + 1)

good_properties <- passes |>
  select(Property, PPNum) |>
  distinct() |>
  group_by(Property) |>
  mutate(timestep = 1:n()) |>
  filter(max(timestep) > 1) |>
  ungroup()

passes <- passes |>
  select(-timestep) |>
  filter(Property %in% good_properties$Property) |>
  arrange(Property, PPNum, pass) |>
  mutate(property_idx = as.numeric(as.factor(Property))) |>
  left_join(good_properties) |>
  arrange(property_idx, Date, timestep, pass)

reps <- passes |>
  select(Property, PPNum) |>
  group_by(Property, PPNum) |>
  mutate(rep = 1:n()) |>
  ungroup() |>
  left_join(good_properties)

data_ls <- suppressWarnings(get_site_data(passes, 0))
Areaper <- data_ls$Areaper
gamma <- Areaper |>
  pivot_longer(cols = -c(ts, Property),
               names_to = "rep") |>
  mutate(rep = as.numeric(str_extract(rep, "\\d*$"))) |>
  rename(timestep = ts) |>
  filter(!is.na(value)) |>
  right_join(reps)

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

county_names <- all_county_units |>
  pull(County) |>
  unique()

county_vec <- paste0(county_names, ", SC")

county_areas_xl <- readxl::read_xls("../pigs-statistical/data/counties/area.xls")
county_areas <- county_areas_xl |>
  select(Areaname, LND110210D) |>
  filter(Areaname %in% county_vec) |>
  mutate(area_km2 = LND110210D * 2.589,
         County = county_names)

survey_obs <- read_csv("../pigs-statistical/data/covariates/FINAL.Process.Model.Observation.Covariates.12Jan2018.csv")

# generate centered and scaled versions of these numeric variables
center_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

c_survey_obs <- survey_obs |>
  filter(STATE_NAME == "South Carolina",
         NAME %in% county_names) |>
  select(NAME, rural.road.density, mean.ruggedness, mean.canopy.density) |>
  mutate(c_road_den = center_scale(rural.road.density),
         c_rugged = center_scale(mean.ruggedness),
         c_canopy = center_scale(mean.canopy.density))

X <- c_survey_obs |>
  select(starts_with("c_")) |>
  as.matrix()
X <- cbind(1, X)

pp <- all_county_units |>
  group_by(Property) |>
  mutate(timestep_idx = 1:n()) |>
  ungroup() |>
  select(Property, timestep, timestep_idx) |>
  pivot_wider(names_from = timestep_idx,
              values_from = timestep,
              id_cols = Property)

n_property <- max(passes$property_idx)
all_pp <- tibble()
for(i in 1:n_property){
  sub <- filter(passes, property_idx == i)
  reps <- tibble(
    property_idx = i,
    PPNum = min(sub$PPNum):max(sub$PPNum),
    timestep = PPNum
  ) |>
    mutate(timestep = timestep - min(timestep) + 1)
  all_pp <- bind_rows(all_pp, reps)
}

all_pp_wide_prop <- all_pp |>
  pivot_wider(names_from = timestep,
              values_from = PPNum)

pop_growth_lookup <- all_pp_wide_prop |>
  pivot_longer(cols = -property_idx,
               names_to = "timestep",
               values_to = "PPNum") |>
  filter(!is.na(PPNum)) |>
  group_by(property_idx) |>
  filter(PPNum < max(PPNum)) |>
  ungroup() |>
  select(-PPNum) |>
  mutate(H = 1:n()) |>
  pivot_wider(values_from = H,
              names_from = timestep) |>
  select(-property_idx)

all_pp_wide <- all_pp_wide_prop |>
  select(-property_idx)
n_pp_include <- apply(all_pp_wide, 1, function(x) max(which(!is.na(x))))


# Generate start and end indices for previous surveys ---------------------
take <- passes |>
  mutate(order = pass) |>
  mutate(method = as.numeric(as.factor(Methods)))

take$start <- 0
take$end <- 0

pb <- txtProgressBar(max = nrow(take), style = 3)
for (i in 1:nrow(take)) {
  if (take$order[i] > 1) {
    idx <- which(take$County == take$County[i] &
                   take$Property == take$Property[i] &
                   take$timestep == take$timestep[i] &
                   take$order < take$order[i])
    take$start[i] <- idx[1]
    take$end[i] <- idx[length(idx)]
    tar_assert_identical(idx, take$start[i]:take$end[i])
  }
  setTxtProgressBar(pb, i)
}
close(pb)

sampled_units <- take |>
  select(County, Property, property_idx, PPNum, timestep) |>
  distinct() |>
  mutate(n_id = 1:n())

county_sampled_units <- sampled_units |>
  select(County, PPNum) |>
  distinct() |>
  mutate(m_id = 1:n()) |>
  left_join(county_areas) |>
  select(-Areaname, -LND110210D)

sampled_units <- left_join(sampled_units, county_sampled_units)

n_timesteps <- sampled_units |>
  group_by(Property) |>
  count() |>
  pull(n)

timestep <- sampled_units |>
  select(-n_id, -m_id, -area_km2, -property_idx) |>
  pivot_wider(names_from = timestep,
              values_from = PPNum) |>
  select(-County, -Property)

y_rem <- take |>
  group_by(property_idx, PPNum) |>
  summarise(ysum = sum(Take)) |>
  ungroup()

y_sum_wide <- left_join(all_pp, y_rem) |>
  mutate(ysum = if_else(is.na(ysum), 0, ysum)) |>
  select(-PPNum) |>
  pivot_wider(names_from = timestep,
              values_from = ysum) |>
  select(-property_idx)

sum_area_surveyed <- take |>
  group_by(County, PPNum) |>
  distinct() |>
  summarise(sum_area_surveyed = sum(area_km2))

constants <- list(
  n_timesteps = n_timesteps,
  n_county = length(unique(take$County)),
  n_survey = nrow(take),
  n_ls = length(data_litter_size),
  n_property = n_property,
  n_county_units = nrow(county_sampled_units),
  n_prop = n_prop,
  n_first_survey = length(which(take$order == 1)),
  n_not_first_survey = length(which(take$order != 1)),
  n_units = nrow(sampled_units),
  n_method = 3,
  n_pp = max(n_pp_include),
  n_pp_prop = n_pp_include,
  all_pp = as.matrix(all_pp_wide),
  m_p = ncol(X),
  first_survey = which(take$order == 1),
  not_first_survey = which(take$order != 1),
  p_property_idx = take$property_idx,
  p_pp_idx = take$PPNum,
  start = take$start,
  end = take$end,
  PPNum = as.matrix(timestep),
  method = take$method,
  property_x = sampled_units$property_idx,
  pp_x = sampled_units$PPNum,
  pp_len = 90,
  phi_prior_mean = 3.3,
  phi_prior_tau = 1,
  M_lookup = n_prop_county1,
  county = as.numeric(as.factor(take$County)),
  pH = as.matrix(pop_growth_lookup)
)

data <- list(
  y = take$Take,
  J = data_litter_size,
  rem = as.matrix(y_sum_wide),
  sum_prop_area = sum_area_surveyed$sum_area_surveyed,
  county_area = county_sampled_units$area_km2,
  X_p = X,
  log_gamma = log(gamma$value)
)

inits <- function(){
  sigma_phi <- runif(1, 0.2, 0.5)
  logit_mean_phi <- rnorm(1, 3, sigma_phi)
  mean_lpy <- 1
  mean_ls <- round(runif(1, 3.5, 6.4))
  zeta <- rep(mean_lpy / 365 * constants$pp_len * mean_ls, constants$n_pp)
  dm <- matrix(NA, n_property, max(constants$all_pp, na.rm = TRUE))
  S <- R <- dm
  for(i in 1:constants$n_property){
    dm[i, constants$all_pp[i, 1]] <- round(runif(1, 100, 500))
    for(t in 2:constants$n_pp_prop[i]){
      phi <- nimble::ilogit(rnorm(1, logit_mean_phi, sigma_phi))
      n_avail <- max(0, dm[i, constants$all_pp[i, t-1]] - data$rem[i, t-1])
      if(is.na(n_avail)) print(t)
      S[i, t-1] <- rbinom(1, n_avail, phi)
      R[i, t-1] <- rpois(1, zeta[constants$all_pp[i, t-1]] * n_avail/2)
      dm[i, constants$all_pp[i, t]] <- S[i, t-1] + R[i, t-1]
    }
  }

  N <- matrix(NA, constants$n_property, max(take$PPNum))
  for(i in 1:constants$n_property){
    for(t in 1:constants$n_timesteps[i]){ # loop through sampled PP only
      N[i, constants$PPNum[i, t]] <- dm[i, constants$PPNum[i, t]]
    }
  }
  list(
    z1 = apply(N, 1, function(x) x[min(which(!is.na(x)))]),
    beta_p = rnorm(ncol(X)),
    logit_mean_phi = logit_mean_phi,
    sigma_phi = sigma_phi,
    S = S,
    R = R,
    N = N,
    mean_ls = mean_ls
  )
}

likelihood_nb <- if_else(likelihood == "nb", TRUE, FALSE)
likelihood_binom <- ifelse(likelihood == "binomial", TRUE, FALSE)
likelihood_poisson <- ifelse(likelihood == "poisson", TRUE, FALSE)
spatial <- FALSE

source("R/nimble_DynamicMultiMethodDM.R")
source("R/functions_nimble.R")
Rmodel <- nimbleModel(
  code = modelCode,
  constants = constants,
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
mcmcConf$removeSamplers("beta_p")
mcmcConf$addSampler("beta_p", "AF_slice")

mcmcConf$addMonitors(c("xn", "M", "lambda"))

# mcmcConf$removeSampler("n")

Rmcmc <- buildMCMC(mcmcConf)
Cmcmc <- compileNimble(Rmcmc)

n_iter <- 50000
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

dest <- file.path(dir_out, dir_model)
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

if(is.na(burnin)) burnin <- 40000

samps_burn <- window(samples, start = burnin)
samps_mat <- as.matrix(samples)
samps_quants <- apply(samps_mat, 2, quantile, c(0.025, 0.5, 0.975))

subset_quants <- function(df, x){
  df[,grep(x, colnames(df), fixed = TRUE)] |>
    as_tibble() |>
    mutate(quant = c("low", "median", "high")) |>
    pivot_longer(cols = -quant,
                 names_to = "node") |>
    pivot_wider(names_from = quant,
                values_from = value)
}

point_size <- 4
line_width <- 2

subset_quants(samps_quants, "beta_p") |>
  mutate(Variable = factor(c("Intercept", "Road density", "Ruggedness", "Canopy density"),
    levels = c("Intercept", "Road density", "Ruggedness", "Canopy density"))) |>
  ggplot() +
  aes(y = Variable, x = median, xmin = low, xmax = high) +
  geom_point(size = point_size) +
  geom_linerange(linewidth = line_width) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect",
       y = "Covariate") +
  theme_bw()

samps_mat[,grep("logit_mean_phi", colnames(samps_mat), fixed = TRUE)] |>
  ilogit() |>
  quantile(c(0.025, 0.5, 0.975)) |>
  as_tibble() |>
  mutate(quant = c("low", "median", "high")) |>
  pivot_longer(cols = -quant,
               names_to = "node") |>
  pivot_wider(names_from = quant,
              values_from = value) |>
  mutate(node = "Mean survival") |>
  ggplot() +
  aes(y = node, x = median, xmin = low, xmax = high) +
  geom_point(size = point_size) +
  geom_linerange(linewidth = line_width) +
  labs(x = "90 day survival rate",
       y = "",
       title = "Mean survival") +
  coord_flip() +
  theme_bw()

samps_mat[,grep("sigma_phi", colnames(samps_mat), fixed = TRUE)] |>
  quantile(c(0.025, 0.5, 0.975)) |>
  as_tibble() |>
  mutate(quant = c("low", "median", "high")) |>
  pivot_longer(cols = -quant,
               names_to = "node") |>
  pivot_wider(names_from = quant,
              values_from = value) |>
  mutate(node = "Error in survival") |>
  ggplot() +
  aes(y = node, x = median, xmin = low, xmax = high) +
  geom_point(size = point_size) +
  geom_linerange(linewidth = line_width) +
  labs(x = "Standard dev.",
       y = "",
       title = "SD survival") +
  coord_flip() +
  theme_bw()

samps_mat[,grep("mean_ls", colnames(samps_mat), fixed = TRUE)] |>
  quantile(c(0.025, 0.5, 0.975)) |>
  as_tibble() |>
  mutate(quant = c("low", "median", "high")) |>
  pivot_longer(cols = -quant,
               names_to = "node") |>
  pivot_wider(names_from = quant,
              values_from = value) |>
  mutate(node = "Litter size") |>
  ggplot() +
  aes(y = node, x = median, xmin = low, xmax = high) +
  geom_point(size = point_size) +
  geom_linerange(linewidth = line_width) +
  labs(x = "Piglets per litter",
       y = "",
       title = "litter size") +
  coord_flip() +
  theme_bw()


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
  left_join(property_idxs) |>
  mutate(n_idx = 1:n())

property_addition <- all_county_units |>
  select(-timestep, -cnty_property)

p_lookup <- left_join(property_lookup, property_addition)

take_sums <- take |>
  group_by(Property, PPNum) |>
  summarise(Take = sum(Take)) |>
  ungroup() |>
  mutate(n_idx = 1:n())

draws <- sample.int(nrow(samps_mat), 5000, replace = TRUE)
N <- samps_mat[draws, grep("xn", colnames(samps_mat))]
n_long <- N |>
  as_tibble() |>
  mutate(iter = 1:n()) |>
  pivot_longer(cols = -iter,
               names_to = "node",
               values_to = "abundance") |>
  mutate(n_idx = as.numeric(str_extract(node, "(?<=\\[)\\d*"))) |>
  left_join(p_lookup) |>
  left_join(take_sums) |>
  mutate(density = abundance / area_km2,
         take_density = Take / area_km2)

N <- n_long |>
  group_by(node, n_idx, Property, PPNum, pp_start_date, pp_end_date,
           timestep, property_idx, County, county_idx, area_km2, Take, take_density) |>
  summarise(low = as.numeric(quantile(density, 0.025)),
            median = as.numeric(quantile(density, 0.5)),
            high = as.numeric(quantile(density, 0.975)),
            mean = mean(density),
            sd = sd(density),
            cv = sd / mean,
            ci_range = high-low) |>
  ungroup() |>
  mutate(metric = "density")

N |>
  filter(Property == "Hayes") |>
  mutate(time = mdy(pp_end_date)) |>
  ggplot() +
  aes(x = time, y = median, ymin = low, ymax = high, color = "Abundance") +
  geom_point(size = point_size) +
  geom_linerange(linewidth = line_width) +
  geom_point(aes(y = take_density, color = "Take")) +
  facet_wrap(~ Property) +
  labs(x = "Date",
       y = "Density") +
  theme_bw()











