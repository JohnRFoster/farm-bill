modelCode <- nimbleCode({

  ### Priors
  # capture probability by method
  for(k in 1:n_methods){
    # for(t in 1:nt){
    #   mu_p[k, t] ~ dnorm(0, tau = 1)
    #   logit(p[k, t]) <- mu_p[k, t]
    # }
      mu_p[k] ~ dnorm(0, tau = 1)
      logit(p[k]) <- mu_p[k]
  }

  if(lambda_static){
    log_mean_r ~ dnorm(0, tau = 1)
  } else if(lambda_county){

    for(i in 1:n_county){
      log_mean_r[i] ~ dnorm(0, tau = 1)
      # log(mean_r[i]) <- log_mean_r[i]
    }

    tau_lambda ~ dgamma(0.01, 0.01)

    for(i in 1:n_property){
      alpha_p[i] ~ dnorm(0, tau = tau_lambda)
    }

  }

  for(i in 1:2){
    beta[i] ~ dnorm(0, 0.1)
  }

  # process error
  tau_p ~ dgamma(0.01, 0.01)

  for(i in 1:n_property){
    n[i, 1] ~ dpois(mu[i, 1])
    log(mu[i, 1]) <- z1[i]
    z1[i] ~ dunif(0, 10)
  }

  # property level process model
  for(i in 1:(n_propertyN)){

    # growth rate coefficient priors
    if(lambda_beta){
      for(j in 1:n_beta){
        beta[i, j] ~ dnorm(0, tau = 1)
      }
    }

    for(t in 2:n_timestep[i]){

      n[propertyN[i], t] ~ dpois(mu[propertyN[i], t])

      log(mu[propertyN[i], t]) <- log_mu_proc[i, t]
      log_mu_proc[i, t] ~ dnorm(log_mu[i, t], tau = tau_p)
      log_mu[i, t] <- log_lambda[propertyN[i], t] + log(z[i, t])

      z[i, t] <- n[propertyN[i], t-1] - y_removed[i, t-1]

      if(include_all_pp){
        log_lambda[propertyN[i], t] <- sum(log_lambda_all[propertyN[i], pp[i, t-1]:(pp[i, t]-1)])
      }
    }
  }

  for(i in 1:n_latent_property){
    for(t in 2:n_latent_timestep){
      n[latent_property[i], t] ~ dpois(mu[latent_property[i], t])
      mu[latent_property[i], t] <- exp(log_lambda_all[latent_property[i], t]) * n[latent_property[i], t-1]
    }
  }

  for(i in 1:n_monitor){
    nm[i] <- n[property_m[i], timestep_m[i]]
  }


  if(include_all_pp){
    for(i in 1:n_units_all){

      if(lambda_static){
        log_lambda_all[property_x[i], timestep_x[i]] <- log_mean_r
      } else if(lambda_data) {
        log_mean_r[i] ~ dnorm(0, tau = 1)
        log_lambda_all[property_x[i], timestep_x[i]] <- log_mean_r[i]
      } else if(lambda_beta){
        log_lambda_all[property_x[i], timestep_x[i]] <- inprod(Xall[i, 1:n_beta], beta[property_x[i], 1:n_beta])
      } else if(lambda_county){
        log_lambda_all[property_x[i], timestep_x[i]] <- log_mean_r[county_x[i]] +
          alpha_p[property_x[i]] +
          beta[1]*x1[county_x[i], timestep_x[i]] +
          beta[2]*x2[county_x[i], timestep_x[i]]
      }
    }
  }

  for(i in 1:n_units){

    # lambda is growth rate (exponential population growth)
    if(!include_all_pp){
      if(lambda_beta){
        # currently using basis functions that vary by property and time
        log_lambda[property[i], timestep[i]] <- inprod(X[i, 1:n_beta], beta[property[i], 1:n_beta])
      } else {
        log_lambda[property[i], timestep[i]] ~ dnorm(0, tau = 1)
      }
    }

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

  # beta_k[1:4] ~ dmnorm(k_mu[1:4], k_cov[1:4, 1:4])

  for(i in 1:n_county){

    size[i] ~ dunif(0, 20) # property dispersion
    # size_m[i] <- dexp(1) # county dispersion
    # log(size[i]) <- inprod(beta_k[1:4], Xk[i, 1:4])
  }


  for(i in 1:n_county_units){ # county unit is a county x PP combination

    ### county dispersion ###
    # county level abundance
    M[i] ~ dnegbin(pM[i], size[county[i]])

    # re-parametrization with mean abundance
    pM[i] <- size[county[i]] / (N_mu[i] + size[county[i]])

    # convert abundance to density (across properties), scale to county abundance
    # M_mu[i] <- N_mu[i] / sum_prop_area[i] * county_area[county[i]]

    # all pigs in county_timestep i
    N_mu[i] <- sum_properties(property = M_lookup[i, 1:n_prop[i]],
                              z = nm[1:n_monitor])


  }



})

