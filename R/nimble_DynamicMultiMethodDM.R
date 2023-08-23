modelCode <- nimbleCode({

  # priors

  # non time varying coefficients - observation model
  for(i in 1:m_p){
    beta_p[i] ~ dnorm(0, tau = 1)
  }

  # estimate apparent survival
  logit_mean_phi ~ dnorm(phi_prior_mean, tau = phi_prior_tau)
  sigma_phi ~ dunif(0, 100)
  # log(sigma_phi) <- log_sigma_phi
  tau_phi <- 1/sigma_phi^2
  for(i in 1:n_property){
    for(t in 1:(n_pp_prop[i]-1)){
      logit_phi[pH[i, t]] ~ dnorm(logit_mean_phi, tau = tau_phi)
      logit(phi[pH[i, t]]) <- logit_phi[pH[i, t]]
    }
  }

  mean_ls ~ dinvgamma(3, 20)

  for(i in 1:n_ls){
    J[i] ~ dpois(mean_ls)
  }

  # convert to expected number of pigs per primary period
  zeta_mu <- 1 / 365 * pp_len * mean_ls
  for(i in 1:n_pp){
    zeta[i] <- zeta_mu
  }

  # property dispersion
  if(likelihood_nb){
    for(i in 1:n_county){
      size[i] ~ dunif(0, 20)
    }
  }


  for(i in 1:n_survey){

    # probability of capture, given that an individual is in the surveyed area
    log_theta[i] <- log(ilogit(inprod(X_p[county[i], 1:m_p], beta_p[1:m_p]))) +
      log_gamma[i]

    # data model
    if(likelihood_binom){
      y[i] ~ dbinom(p[i], z[i])
    }

    if(likelihood_nb){
      y[i] ~ dnegbin(py[i], size[p_county_idx[i]])
      py[i] <- size[p_county_idx[i]] / (y_mu[i] + size[p_county_idx[i]])
      y_mu[i] <- p[i] * z[i]
    }

    if(likelihood_poisson){
      y[i] ~ dpois(p[i] * N[p_property_idx[i], p_pp_idx[i]])
    }

  }

  # the probability an individual is captured on the first survey
  for(i in 1:n_first_survey){
    log(p[first_survey[i]]) <- log_theta[first_survey[i]]
  }

  # the probability an individual is captured after the first survey
  for(i in 1:n_not_first_survey){
    log(p[not_first_survey[i]]) <- log_theta[start[not_first_survey[i]]] +
      sum(log(1 - exp(log_theta[start[not_first_survey[i]]:end[not_first_survey[i]]])))
  }

  for(i in 1:n_property){

    # eps_property_pR[i] ~ dnorm(0, tau = 1) # property effect in observation model

    N[i, PPNum[i, 1]] ~ dpois(z1[i])
    z1[i] ~ dinvgamma(1, 1)

    for(t in 2:n_timesteps[i]){ # loop through sampled PP only
      N[i, PPNum[i, t]] ~ dpois(dm[i, PPNum[i, t]])
    }

    # population growth across time steps
    dm[i, all_pp[i, 1]] <- N[i, PPNum[i, 1]]
    for(j in 2:n_pp_prop[i]){ # loop through every PP, including missing ones

      Z[i, j-1] <- dm[i, all_pp[i, j-1]] - rem[i, j-1]
      S[i, j-1] ~ dbinom(phi[pH[i, j-1]], Z[i, j-1])
      R[i, j-1] ~ dpois(zeta[all_pp[i, j-1]] * Z[i, j-1] / 2)
      dm[i, all_pp[i, j]] <- S[i, j-1] + R[i, j-1]
      lambda[pH[i, j-1]] <- dm[i, all_pp[i, j]] / dm[i, all_pp[i, j-1]]
    }

  }

  # for easier monitoring of abundance - long format
  for(i in 1:n_units){
    xn[i] <- N[property_x[i], pp_x[i]]
  }

  for(i in 1:n_county_units){ # county unit is a county x PP combination

    ### county dispersion ###
    # county level abundance
    # a priori we know it cannot me less than total abundance across properties
    if(likelihood_nb){
      M[i] ~ T(dnegbin(pM[i], size[county[i]]), N_mu[i], 0)

      # re-parametrization with mean abundance
      pM[i] <- size[county[i]] / (M_mu[i] + size[county[i]])
    } else {
      M[i] ~ dpois(M_mu[i])
    }

    # convert abundance to density (across properties), scale to county abundance
    M_mu[i] <- N_mu[i] / sum_prop_area[i] * county_area[i]

    # all pigs in county_timestep i
    N_mu[i] <- sum_properties(property = M_lookup[i, 1:n_prop[i]],
                              z = xn[1:n_units])

  }

})
