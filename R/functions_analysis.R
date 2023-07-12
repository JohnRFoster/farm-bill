library(coda)
library(nimble)
library(tidyverse)
library(stringr)

############
# read and load model output
############
read_output <- function(filedest) {
  out <- read_rds(filedest)
  list2env(out, envir = globalenv()) |>
    invisible()
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
name_nodes <- function(mcmc_samples, n_nodes, m_nodes, model_name, dm_lookup, n_mcmc = 5000){

  method_lookup <- method_lookup |>
    rename(method = Methods)

  samples_mat <- as.matrix(mcmc_samples)
  draws <- sample.int(nrow(samples_mat), n_mcmc, replace = TRUE)

  mcmc <- samples_mat |>
    as_tibble() |>
    slice(draws) |>
    select(contains(nodes1), all_of(n_nodes), all_of(m_nodes)) |>
    mutate(iter = 1:n()) |>
    pivot_longer(cols = -iter,
                 names_to = "node") |> #filter(grepl("tau_p", node))
    mutate(
      model = model_name,
      name = if_else(grepl("mu_p[", node, fixed = TRUE), "Capture Prob", "x"),
      name = if_else(grepl("z1", node), "Initial Condition", name),
      name = if_else(grepl("mu_phi", node), "Survival Rate", name),
      name = if_else(grepl("log_rho", node), "Per Capita Recruitment Rate", name),
      name = if_else(grepl("n", node), "Abundance", name),
      name = if_else(grepl("M", node), "County Abundance", name),
      name = if_else(grepl("log_mean_r", node), "Growth Rate", name),
      # name = if_else(grepl("N_disp", node), "Property Abundance with dispersion", name),
      name = if_else(grepl("tau", node), "Process Error", name),
      name = if_else(grepl("beta", node), "Basis Function Coefficients", name),
      value = if_else(grepl("mu_p[", node, fixed = TRUE), boot::inv.logit(value), value),
      value = if_else(grepl("mu_phi", node), boot::inv.logit(value), value),
      value = if_else(grepl("log_rho", node), exp(value), value),
      value = if_else(grepl("log_mean_r", node), exp(value), value))

  my_summary <- function(df){
   df |> summarise(lower = quantile(value, 0.025),
                   median = quantile(value, 0.5),
                   upper = quantile(value, 0.975))
  }

  n_mcmc <- mcmc |>
    filter(grepl("n[", node, fixed = TRUE)) |>
    mutate(property_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
           timestep = as.numeric(str_extract(node, "\\d*(?=\\])"))) |>
    left_join(property_lookup) |>
    group_by(node, model, name, property_idx, timestep, Property, PPNum, pp_start_date, pp_end_date) |>
    my_summary() |>
    ungroup()

  m_mcmc <- mcmc |>
    filter(grepl("M[", node, fixed = TRUE)) |>
    mutate(m_idx = as.numeric(str_extract(node, "(?<=\\[)\\d*")),
           value = log(value + 1)) |>
    left_join(county_level_lookup, by = "m_idx") |>
    group_by(node, model, name, m_idx, County, PPNum, sum_area) |>
    my_summary() |>
    ungroup()

  p_mcmc <- mcmc |>
    filter(grepl("mu_p[", node, fixed = TRUE)) |>
    mutate(method_idx = as.numeric(str_extract(node, "(?<=mu_p\\[)\\d"))) |>
    left_join(method_lookup) |>
    group_by(node, model, name, method_idx, method) |>
    my_summary() |>
    ungroup()

  dm_mcmc <- mcmc |>
    filter(grepl("mu_phi", node, fixed = TRUE) |
             grepl("log_rho", node, fixed = TRUE)) |>
    mutate(property_timestep_idx = as.numeric(str_extract(node, "\\d*(?=\\])"))) |>
    left_join(dm_lookup) |>
    group_by(node, model, name, property_timestep_idx, Property, property_idx, PPNum, pp_start_date, pp_end_date, timestep) |>
    my_summary() |>
    ungroup()

  z1_mcmc <- mcmc |> filter(node == "z1")  |>
    group_by(node, model, name) |>
    my_summary() |>
    ungroup()

  return(bind_rows(n_mcmc, m_mcmc, p_mcmc, dm_mcmc, z1_mcmc))
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
# plot survival
############
gg_lambda_pointrange <- function(mcmc_quants){
  mcmc_quants |>
    filter(name == "Survival Rate") |>
    ggplot() +
    aes(x = pp_end_date, ymin = lower, ymax = upper, y = median, color = model) +
    geom_pointrange(position = position_dodge(width=45)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~ Property) +
    # coord_cartesian(ylim = c(0, 15)) +
    labs(title = "Survival Rate",
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

  mcmc_quants |>
    filter(name == "Abundance") |>
    left_join(obs) |>
    ggplot() +
    aes(x = pp_end_date, ymin = lower, ymax = upper, y = median) +
    geom_pointrange(position = position_dodge(width=45), color = "blue") +
    geom_line(position = position_dodge(width=45), color = "blue") +
    geom_line(aes(y = Observed), color = "black") +
    geom_point(aes(y = Observed), color = "black") +
    facet_wrap( ~ model) +
    # scale_y_log10() +
    labs(title = "Abundnace",
         x = "PP end date",
         y = "Number of pigs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


############
# plot county abundance
############
gg_county_abundance_pointrange <- function(mcmc_quants, county){
  obs <- removed |>
    as_tibble() |>
    mutate(property_idx = 1:n()) |>
    pivot_longer(cols = -property_idx,
                 names_to = "timestep",
                 values_to = "Observed") |>
    mutate(timestep = as.numeric(timestep)) |>
    left_join(property_lookup) |>
    mutate(County = if_else(Property %in% c("Chilton", "Stratos", "Youman 1", "Corbin"), "Hampton", "x"),
           County = if_else(Property %in% c("H&B", "McKinney"), "Newberry", County),
           County = if_else(Property %in% c("Mixon", "Okeetee"), "Jasper", County)) |>
    filter(!is.na(PPNum)) |>
    group_by(County, PPNum) |>
    summarise(observed = log(sum(Observed) + 1)) |>
    filter(County == county)

  mcmc_quants |>
    filter(name == "County Abundance",
           County == county) |>
    left_join(obs) |>
    ggplot() +
    aes(x = PPNum, ymin = lower, ymax = upper, y = median) +
    geom_pointrange(position = position_dodge(width=0.5), color = "blue") +
    geom_line(aes(y = observed), color = "black") +
    geom_point(aes(y = observed), color = "black") +
    facet_wrap(~ model) +
    # scale_y_log10() +
    labs(y = "log(County Abundnace + 1)",
         x = "PP") +
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

############
# plot process error
############
gg_process_error_pointrange <- function(mcmc_quants) {
  mcmc_quants |>
    filter(name == "Process Error") |>
    ggplot() +
    aes(x = model, ymin = lower, ymax = upper, y = median, color = model) +
    geom_pointrange(size = 1.2, linewidth = 1.2, position = position_dodge(width=0.5)) +
    labs(title = "Model structural error",
         x = "Model",
         y = "Standard deviation (swine/PP/property)") +
    theme_bw()
}

############
# plot betas - basis function coefficients
############
gg_beta_basis_pointrange <- function(mcmc_quants) {
  prop <- property_lookup |>
    select(Property, property_idx)
  mcmc_quants |>
    filter(name == "Basis Function Coefficients") |>
    rename(beta_n = timestep) |>
    select(-Property, -method, -PPNum, -pp_start_date, -pp_end_date) |>
    left_join(prop) |>
    ggplot() +
    aes(x = beta_n, ymin = lower, ymax = upper, y = median, color = model) +
    geom_pointrange(size = 1.2, linewidth = 1.2, position = position_dodge(width=0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~ Property) +
    labs(title = "Basis Function Coefficients",
         x = "Beta number (1 = intercept, 2:5 = basis functions)",
         y = "Effect") +
    theme_bw()
}

# Rcpp::cppFunction(
#   'int fibonacci(const int x) {
#               if (x == 0) return(0);
#               if (x == 1) return(1);
#               return (fibonacci(x - 1)) + fibonacci(x - 2);
#           }')
#
# fibonacci(7)



