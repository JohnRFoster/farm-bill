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

  # if(include_county){
  #   for(i in 1:n_county){
  #     mu_r[i] ~ dnorm(0, tau = 1)
  #     logit(r[i]) <- mu_k[i]
  #   }
  # }

  if(!include_beta){
    mean_r ~ dnorm(0, tau = 1)
  }

  # process error
  if(include_proc){
    tau_p ~ dgamma(0.01, 0.01)
  }

  # property level process model
  for(i in 1:n_property){

    # growth rate coefficient priors
    if(include_beta){
      for(j in 1:n_beta){
        beta[i, j] ~ dnorm(0, tau = 1)
      }
    }

    n[i, 1] ~ dpois(mu[i, 1])
    log(mu[i, 1]) <- log_mu_proc[i, 1]
    log_mu_proc[i, 1] <- log_mu1[i]
    log_mu1[i] ~ dnorm(0, tau = 0.1)

    for(t in 2:n_timestep[i]){

      if(include_proc){
        log(mu[i, t]) <- log_mu_proc[i, t]
        log_mu_proc[i, t] ~ dnorm(log_mu[i, t], tau = tau_p)
      } else {
        log(mu[i, t]) <- log_mu[i, t]
      }

      n[i, t] ~ dpois(mu[i, t])
      log_mu[i, t] <- log_lambda[i, t-1] + log(z[i, t-1])
      z[i, t-1] <- n[i, t-1] - y_removed[i, t-1]

      if(include_all_pp){
        log_lambda[i, t-1] <- sum(log_lambda_all[i, pp[i, t-1]:pp[i, t]])
      }
    }
  }

  if(include_all_pp){
    for(i in 1:n_units_all){
      if(intercept_only){
        log_lambda_all[property_x[i], timestep_x[i]] <- mean_r
      } else {
        log_lambda_all[property_x[i], timestep_x[i]] <- inprod(Xall[i, 1:n_beta], beta[property_x[i], 1:n_beta])
      }
    }
  }

  for(i in 1:n_units){

    # lambda is growth rate (exponential population growth)
    if(!include_all_pp){
      if(include_beta){
        # currently using basis functions that vary by property and time
        log_lambda[property[i], timestep[i]] <- inprod(X[i, 1:n_beta], beta[property[i], 1:n_beta])
      } else {
        log_lambda[property[i], timestep[i]] ~ dnorm(0, tau = 1)
      }
    }

    Nc[countyN[i], cnty_property[i], PPNum[i]] <- n[property[i], timestep[i]]

    # likelihood - first pass
    if(nb_likelihood){
      y[i, 1] ~ dnegbin(py[i, 1], size[countyN[i]])
      py[i, 1] <- size[countyN[i]] / (y_mu[i, 1] + size[countyN[i]])
      y_mu[i, 1] <- pi[i, 1] * n[property[i], timestep[i]]
    }

    if(pois_likelihood){
      y[i, 1] ~ dpois(pi[i, 1] * n[property[i], timestep[i]])
    }

    pi[i, 1] <- gamma[i, 1] * theta[i, 1] * (1 - gamma[i, 1]) + (gamma[i, 1] * (1 - theta[i, 1]))
    theta[i, 1] <- 1 - (1 - p[method[i, 1]]) ^ effort[i, 1]

    for(j in 2:n_passes[i]){

      # likelihood - subsequent passes
      if(nb_likelihood){
        y[i, j] ~ dnegbin(py[i, j], size[countyN[i]])
        py[i, j] <- size[countyN[i]] / (y_mu[i, j] + size[countyN[i]])
        y_mu[i, j] <- pi[i, j] * n[property[i], timestep[i]]
      }

      if(pois_likelihood){
        y[i, j] ~ dpois(pi[i, j] * n[property[i], timestep[i]])
      }

      pi[i, j] <- gamma[i, j] * theta[i, j] *
        exp(sum(log((1 - gamma[i, 1:(j-1)]) + (gamma[i, 1:(j-1)] * (1 - theta[i, 1:(j-1)])))))
      theta[i, j] <- 1 - (1 - p[method[i, j]]) ^ effort[i, j]

    }

  }

  # beta_k[1:4] ~ dmnorm(k_mu[1:4], k_cov[1:4, 1:4])

  for(i in 1:n_county){

    size[i] ~ dunif(0, 50) # property dispersion
    # size_m[i] <- dexp(1) # county dispersion
    # log(size[i]) <- inprod(beta_k[1:4], Xk[i, 1:4])
  }


  for(i in 1:n_county_units){ # county unit is a county x PP combination

    ### property dispersion ###
    # total abundance in sampled properties within county i
    # N_disp[i] ~ dnegbin(pN[i], size[county[i]])

    # re-parametrization with mean abundance
    # pN[i] <- size[county[i]] / (N_mu[i] + size[county[i]])

    # all pigs in county i at time j
    # N_mu[i, j] <- sum(Nc[i, 1:n_prop_cnty[i], pp_c[i, j]])
    N_mu[i] <- sum_properties(property = n_prop_cnty[i, 1:n_prop[i]],
                              # p = pp_c[i],
                              z = Nc[county[i], 1:n_pc[county[i]], pp_c[i]])

    ### county dispersion ###
    # county level abundance
    M[i] ~ dnegbin(pM[i], size[county[i]])

    # re-parametrization with mean abundance
    pM[i] <- size[county[i]] / (M_mu[i] + size[county[i]])

    # convert abundance to density (across properties), scale to county abundance
    M_mu[i] <- N_mu[i] / sum_prop_area[i] * county_area[county[i]]


  }



})

