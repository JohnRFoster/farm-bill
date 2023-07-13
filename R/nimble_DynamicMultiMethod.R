modelCode <- nimbleCode({

  ### Priors ###
  # hyperpriors for capture probability by method
  mu_p ~ dnorm(0, tau = 1)
  sigma_p ~ dinvgamma(alpha_p, beta_p)
  tau_p <- 1/sigma_p^2
  for(k in 1:n_methods){
    # prior for each capture probability by method
    logit_p[k] ~ dnorm(mu_p, tau = tau_p)
    logit(p[k]) <- logit_p[k]
  }

  # hyperpriors for beta
  sigma_b ~ dinvgamma(alpha_b, beta_b)
  tau_b <- 1/sigma_b^2
  for(i in 1:n_beta){
    # prior for coefficient of basis functions
    beta[i] ~ dnorm(0, tau = tau_b)
  }

  # hyperpriors for lambda
  sigma_l ~ dinvgamma(alpha_l, beta_l)
  tau_l <- 1/sigma_l^2

  # county dispersion
  for(i in 1:n_county){
    size[i] ~ dunif(0, 20)
  }

  # property level process model
  for(i in 1:n_property){

    # prior on first latent state
    n[i, 1] ~ dpois(mu[i, 1])
    log(mu[i, 1]) <- z1[i]
    z1[i] ~ dunif(0, 10)

    for(t in 2:n_timestep[i]){

      n[i, t] ~ dpois(mu[i, t])
      log(mu[i, t]) <- log_lambda[i, t] + log(z[i, t])
      z[i, t] <- n[i, t-1] - y_removed[i, t-1]
      log_lambda[i, t] <- sum(log_lambda_all[i, pp[i, t-1]:(pp[i, t]-1)])

    }
  }

  # calculate log growth rate across timesteps
  for(i in 1:n_units_all){
    log_lambda_mu[property_x[i], timestep_x[i]] <- inprod(Xall[i, 1:n_beta], beta[1:n_beta])
    log_lambda_all[property_x[i], timestep_x[i]] ~ dnorm(log_lambda_mu[property_x[i], timestep_x[i]], tau = tau_l)
  }

  for(i in 1:n_units){

    # likelihood - first pass
    if(nb_likelihood){
      y[i, 1] ~ dnegbin(py[i, 1], size[countyN[i]])
      py[i, 1] <- size[countyN[i]] / (y_mu[i, 1] + size[countyN[i]])
      y_mu[i, 1] <- pi[i, 1] * n[property[i], timestep[i]]
    }

    if(pois_likelihood){
      y[i, 1] ~ dpois(pi[i, 1] * n[property[i], timestep[i]])
    }

    pi[i, 1] <- gamma[i, 1] * theta[i, 1]
    theta[i, 1] <- 1 - pow((1 - p[method[i, 1]]), effort[i, 1])

    for(j in 2:n_passes[i]){

      # likelihood - subsequent passes
      if(nb_likelihood){
        y[i, j] ~ dnegbin(py[i, j], size[countyN[i]])
        py[i, j] <- size[countyN[i]] / (y_mu[i, j] + size[countyN[i]])
        y_mu[i, j] <- pi[i, j] * n[property[i], timestep[i]]
      }

      if(pois_likelihood){
        lam[i, j] <- pi[i, j] * n[property[i], timestep[i]]
        y[i, j] ~ dpois(lam[i, j])
      }

      pi[i, j] <- gamma[i, j] * theta[i, j] *
        exp(sum(log((1 - gamma[i, 1:(j-1)]) + (gamma[i, 1:(j-1)] * (1 - theta[i, 1:(j-1)])))))
      theta[i, j] <- 1 - pow((1 - p[method[i, j]]), effort[i, j])

    }

  }

  # county population extrapolation
  for(i in 1:n_county_units){ # county unit is a county x PP combination

    ### county dispersion ###
    # county level abundance
    M[i] ~ T(dnegbin(pM[i], size[county[i]]), N_mu[i], Inf)

    # re-parametrization with mean abundance
    pM[i] <- size[county[i]] / (M_mu[i] + size[county[i]])


    # convert abundance to density (across properties), scale to county abundance
    M_mu[i] <- N_mu[i] / sum_prop_area[i] * county_area[county[i]]

    # all pigs in county_timestep i
    N_mu[i] <- sum_properties(property = M_lookup[i, 1:n_prop[i]],
                              z = xn[1:n_monitor])


  }

  # long format for easier monitoring and use in county model
  for(i in 1:n_monitor){
    xn[i] <- n[property_m[i], timestep_m[i]]
  }


})

