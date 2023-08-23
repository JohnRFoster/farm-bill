library(config)
library(tidyverse)
library(nimble)
library(lubridate)
library(targets)
library(coda)

config_name <- "dynamic_multi_method_poisson"

Sys.setenv(R_CONFIG_ACTIVE = config_name)
config <- config::get()
dir_model <- config_name
nb_likelihood <- if_else(config$nb_likelihood, TRUE, FALSE)
pois_likelihood <- !nb_likelihood

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
  ungroup()

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

data_ls <- suppressWarnings(get_site_data(passes))
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
  group_by(Property, PPNum) |>
  mutate(order = 1:n()) |>
  ungroup() |>
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
  n_beta_p = 4,
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
  p_property_idx = take$Property,
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
  county = as.numeric(as.factor(county_sampled_units$County)),
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
