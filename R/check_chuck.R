library(coda)
library(nimble)
library(tidyverse)
library(stringr)

dir_data <- "data"
take_csv <- "all_chuck_data.csv"

dir_out <- "out"
dir_model <- "basic"

out <- read_rds(file.path(dir_out, dir_model, "samples.rds"))
samples <- out$samples                    # mcmc samples
data <- out$nimble_data                   # data fed to nimble
constants <- out$nimble_constants         # constants fed to nimble
method_lookup <- out$method_lookup        # method look up table
property_lookup <- out$property_lookup    # property look up table
nodes1 <- out$nodes1                      # nodes with consistent indexes
nodes2 <- out$nodes2                      # nodes with inconsistent indexes


samples_mat <- as.matrix(samples)

draws <- sample.int(nrow(samples_mat), 5000, replace = TRUE)

mcmc <- samples_mat |>
  as_tibble() |>
  slice(draws) |>
  select(contains("mu_p"), contains("log_mu1"), all_of(nodes2)) |>
  mutate(iter = 1:n()) |>
  pivot_longer(cols = -iter,
               names_to = "node") |>
  mutate(name = if_else(grepl("mu_p", node), "Capture Prob", "x"),
         name = if_else(grepl("log_mu1", node), "Initial Condition", name),
         name = if_else(grepl("lambda", node), "Growth Rate", name),
         name = if_else(grepl("n", node), "Abundance", name),
         value = if_else(grepl("mu_p", node), boot::inv.logit(value), value),
         value = if_else(grepl("log_mu1", node), exp(value), value),
         value = if_else(grepl("lambda", node), exp(value), value),
         property_idx = str_extract(node, "(?<=\\[)\\d"),
         property_idx = as.numeric(property_idx),
         property_idx = if_else(grepl("mu_p", node), -1, property_idx),
         timestep = str_extract(node, "\\d*(?=\\])"),
         timestep = as.numeric(timestep),
         method = if_else("mu_p[1]" == node, "Aerial", "x"),
         method = if_else("mu_p[2]" == node, "Shooting", method),
         method = if_else("mu_p[3]" == node, "Trap", method)) |>
  left_join(property_lookup)

param_stats <- mcmc |>
  group_by(node, name, property_idx, timestep, method, Property) |>
  summarise(lower = quantile(value, 0.025),
            median = quantile(value, 0.5),
            upper = quantile(value, 0.975)) |>
  ungroup() |>
  left_join(property_lookup)


param_stats |>
  filter(name == "Capture Prob") |>
  ggplot() +
  aes(x = method, ymin = lower, ymax = upper, y = median) +
  geom_pointrange(size = 1.2, linewidth = 1.2) +
  labs(title = "Capture Prob by method",
       x = "Method",
       y = "Probability") +
  theme_bw()


param_stats |>
  filter(name == "Growth Rate") |>
  ggplot() +
  aes(x = PPNum, ymin = lower, ymax = upper, y = median) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~ Property) +
  coord_cartesian(ylim = c(0, 15)) +
  labs(title = "Population growth rate",
       x = "Primary period index",
       y = "Growth rate") +
  theme_bw()

obs <- data$y_removed |>
  as_tibble() |>
  mutate(property_idx = 1:n()) |>
  pivot_longer(cols = -property_idx,
               names_to = "timestep",
               values_to = "Observed") |>
  mutate(timestep = as.numeric(timestep)) |>
  left_join(property_lookup)

param_stats |>
  filter(name == "Abundance") |>
  left_join(obs) |>
  ggplot() +
  aes(x = pp_end_date, ymin = lower, ymax = upper, y = median) +
  geom_pointrange() +
  geom_point(aes(y = Observed), color = "purple") +
  facet_wrap(~ Property) +
  # scale_y_log10() +
  labs(title = "Abundnace",
       x = "Timestep",
       y = "Number of pigs") +
  theme_bw()



