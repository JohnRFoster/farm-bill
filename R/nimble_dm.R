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

  if(phi_static | phi_hb){
    mu_phi ~ dnorm(phi_prior[1], tau = 1/phi_prior[2]^2)
    if(phi_static){
      logit(phi) <- mu_phi
    }
    if(phi_hb){
      tau_surv ~ dgamma(0.001, 0.001)
    }
  }

  if(rho_static){
    log_rho ~ dnorm(1, tau = 1)
    log(rho) <- log_rho
  }

  for(i in 1:n_county){
    size[i] ~ dunif(0, 20) # property dispersion
  }

  z1 ~ dgamma(1, 0.001)

  # property level process model
  for(i in 1:n_property){

    # survival coefficient priors
    if(phi_beta){
      for(j in 1:n_beta){
        beta_p[i, j] ~ dnorm(0, tau = 1)
      }
    }

    # if(phi_hb){
    #   mu_phi_prop[i] ~ dnorm(mu_phi, tau = tau_surv)
    #   logit(phi[i]) <- mu_phi_prop[i]
    # }

    # recruitment coefficient priors
    if(rho_beta){
      for(j in 1:n_beta){
        beta_r[i, j] ~ dnorm(0, tau = 1)
      }
    }

    n[i, 1] ~ dpois(z1)
    # z1[i] ~ dunif(0, 10000)

    for(t in 2:n_timestep[i]){

      n[i, t] ~ dpois(z[i, t, pp[i, t]])

      z[i, t, pp[i, t-1]] <- n[i, t-1] - y_removed[i, t-1]
      for(j in pp[i, t-1]:(pp[i, t]-1)){

        s[i, t, j] ~ dbinom(phi_all[i, j], z[i, t, j])
        r[i, t, j] ~ dpois(rho_all[i, j] * z[i, t, j])
        z[i, t, j+1] <- s[i, t, j] + r[i, t, j]

      } # j
    } # t
  } # i

  for(i in 1:n_units_all){

    if(phi_static){
      phi_all[property_x[i], timestep_x[i]] <- phi
    } else if(phi_hb){
      mu_phi_prop[i] ~ dnorm(mu_phi, tau = tau_surv)
      logit(phi_all[property_x[i], timestep_x[i]]) <- mu_phi_prop[i]
    } else if(phi_data){
      mu_phi[i] ~ dnorm(phi_prior[1], tau = 1/phi_prior[2]^2)
      logit(phi_all[property_x[i], timestep_x[i]]) <- mu_phi[i]
    } else {
      logit(phi_all[property_x[i], timestep_x[i]]) <- inprod(Xall[i, 1:n_beta], beta_p[property_x[i], 1:n_beta])
    }

    if(rho_static){
      rho_all[property_x[i], timestep_x[i]] <- rho
    } else if(rho_data){
      log_rho[i] ~ dnorm(1, tau = 1)
      log(rho_all[property_x[i], timestep_x[i]]) <- log_rho[i]
    } else {
      log(rho_all[property_x[i], timestep_x[i]]) <- inprod(Xall[i, 1:n_beta], beta_r[property_x[i], 1:n_beta])
    }

  }

  for(i in 1:n_units){

    Nc[countyN[i], cnty_property[i], PPNum[i]] <- n[property[i], timestep[i]]

    # likelihood - first pass
    if(nb_likelihood){
      y[i, 1] ~ dnegbin(py[i, 1], size[countyN[i]])
      py[i, 1] <- size[countyN[i]] / (y_mu[i, 1] + size[countyN[i]])
      y_mu[i, 1] <- pi[i, 1] * n[property[i], timestep[i]]
    } else {
      y[i, 1] ~ dpois(pi[i, 1] * n[property[i], timestep[i]])
    }

    pi[i, 1] <- gamma[i, 1] * theta[i, 1]
    theta[i, 1] <- 1 - (1 - p[method[i, 1]]) ^ effort[i, 1]

    for(j in 2:n_passes[i]){

      # likelihood - subsequent passes
      if(nb_likelihood){
        y[i, j] ~ dnegbin(py[i, j], size[countyN[i]])
        py[i, j] <- size[countyN[i]] / (y_mu[i, j] + size[countyN[i]])
        y_mu[i, j] <- pi[i, j] * n[property[i], timestep[i]]
      } else {
        lam[i, j] <- pi[i, j] * n[property[i], timestep[i]]
        y[i, j] ~ dpois(lam[i, j])
      }

      pi[i, j] <- gamma[i, j] * theta[i, j] *
        exp(sum(log((1 - gamma[i, 1:(j-1)]) + (gamma[i, 1:(j-1)] * (1 - theta[i, 1:(j-1)])))))
      theta[i, j] <- 1 - (1 - p[method[i, j]]) ^ effort[i, j]

    }

  }

  for(i in 1:n_county_units){ # county unit is a county x PP combination


    ### county dispersion ###
    # county level abundance
    M[i] ~ dnegbin(pM[i], size[county[i]])

    # re-parametrization with mean abundance
    pM[i] <- size[county[i]] / (M_mu[i] + size[county[i]])

    # convert abundance to density (across properties), scale to county abundance
    M_mu[i] <- N_mu[i] / sum_prop_area[i] * county_area[county[i]]

    # all pigs in county_timestep i
    N_mu[i] <- sum_properties(property = n_prop_cnty[i, 1:n_prop[i]],
                              z = Nc[county[i], 1:n_pc[county[i]], pp_c[i]])

  }



})

