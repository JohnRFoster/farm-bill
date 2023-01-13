library(coda)
library(nimble)
library(tidyverse)
library(stringr)

############
# read and load model output
############
read_output <- function(filedest) {
  out <- read_rds(filedest)
  list2env(out, envir = .GlobalEnv)
  # out$samples            # mcmc samples
  # out$nimble_data        # data fed to nimble
  # out$nimble_constants   # constants fed to nimble
  # out$method_lookup      # method look up table
  # out$property_lookup    # property look up table
  # out$nodes1             # nodes with consistent indexes
  # out$nodes2             # nodes with inconsistent indexes
}


############
# modify output with correct names, time steps, and properties
############
name_nodes <- function(mcmc_samples, model_name, n_mcmc = 5000){

  samples_mat <- as.matrix(mcmc_samples)
  draws <- sample.int(nrow(samples_mat), n_mcmc, replace = TRUE)

  mcmc <- samples_mat |>
    as_tibble() |>
    slice(draws) |>
    # select(contains("mu_p"), contains("log_mu1"), all_of(nodes2)) |>
    mutate(iter = 1:n()) |>
    pivot_longer(cols = -iter,
                 names_to = "node") |>
    mutate(
      model = model_name,
      name = if_else(grepl("mu_p", node), "Capture Prob", "x"),
      name = if_else(grepl("log_mu1", node), "Initial Condition", name),
      name = if_else(grepl("lambda", node), "Growth Rate", name),
      name = if_else(grepl("n", node), "Abundance", name),
      name = if_else(grepl("tau", node), "Process Error", name),
      name = if_else(grepl("beta", node), "Basis Function Coefficients", name),
      value = if_else(grepl("mu_p", node), boot::inv.logit(value), value),
      value = if_else(grepl("log_mu1", node), exp(value), value),
      value = if_else(grepl("lambda", node), exp(value), value),
      value = if_else(grepl("tau", node), 1/sqrt(value), value),
      property_idx = str_extract(node, "(?<=\\[)\\d"),
      property_idx = as.numeric(property_idx),
      property_idx = if_else(grepl("mu_p", node), -1, property_idx),
      timestep = str_extract(node, "\\d*(?=\\])"),
      timestep = as.numeric(timestep),
      method = if_else("mu_p[1]" == node, "Aerial", "x"),
      method = if_else("mu_p[2]" == node, "Shooting", method),
      method = if_else("mu_p[3]" == node, "Trap", method)
      ) |>
    left_join(property_lookup)
  return(mcmc)
}


############
# create quantiles
############
get_quantiles <- function(mcmc_names){
  mcmc_names |>
    group_by(node, name, property_idx, timestep, method, Property, model) |>
    summarise(lower = quantile(value, 0.025),
              median = quantile(value, 0.5),
              upper = quantile(value, 0.975)) |>
    ungroup() |>
    left_join(property_lookup)
}


############
# plot capture rate
############
gg_capture_prob_pointrange <- function(mcmc_quants) {
  mcmc_quants |>
    filter(name == "Capture Prob") |>
    ggplot() +
    aes(x = method, ymin = lower, ymax = upper, y = median, color = model) +
    geom_pointrange(size = 1.2, linewidth = 1.2, position = position_dodge(width=0.5)) +
    labs(title = "Capture Prob by method",
         x = "Method",
         y = "Probability") +
    theme_bw()
}

############
# plot growth rate
############
gg_lambda_pointrange <- function(mcmc_quants){
  mcmc_quants |>
    filter(name == "Growth Rate") |>
    ggplot() +
    aes(x = pp_end_date, ymin = lower, ymax = upper, y = median, color = model) +
    geom_pointrange(position = position_dodge(width=45)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~ Property) +
    # coord_cartesian(ylim = c(0, 15)) +
    labs(title = "Population growth rate",
         x = "PP end date",
         y = "Growth rate") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

############
# plot abundance
############
gg_abundance_pointrange <- function(mcmc_quants, removed){
  obs <- removed |>
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
    aes(x = pp_end_date, ymin = lower, ymax = upper, y = median, color = model) +
    geom_pointrange(position = position_dodge(width=45)) +
    geom_point(aes(y = Observed), color = "black") +
    facet_wrap(~ Property) +
    # scale_y_log10() +
    labs(title = "Abundnace",
         x = "PP end date",
         y = "Number of pigs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


############
# calculate growth rate
############

calc_lambda <- function(beta_mcmc, X){

  l <- with(nimble_constants, {
    lambda_all <- tibble()
    for(i in 1:n_units){

      beta_node <- paste0("beta[", property[i])

      beta <- beta_mcmc |>
        filter(grepl(beta_node, node, fixed = TRUE)) |>
        select(iter, value, timestep) |>
        pivot_wider(names_from = timestep,
                    values_from = value) |>
        select(-iter) |>
        as.matrix()

      lambda <- vector()
      for(m in seq_len(nrow(beta))){
        lambda[m] <- exp(X[i,] %*% beta[m,])
      }

      lambda_tb <- tibble(
        property_idx = property[i],
        timestep = timestep[i],
        lower = quantile(lambda, 0.025),
        median = quantile(lambda, 0.5),
        upper = quantile(lambda, 0.975),
      )
    lambda_all <- bind_rows(lambda_all, lambda_tb)
    }
    lambda_all
  })
  left_join(l, property_lookup)
}
