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

  # county level process model
  for(i in 1:n_county){

    # prior on first latent state
    M[i, start_pp[i]] ~ dpois(mu[i, start_pp[i]])
    log(mu[i, start_pp[i]]) <- z1[i]
    z1[i] ~ dunif(0, 15)

    M_avail[i, start_pp[i]] ~ dbinom(proportion_surveyed[i, start_pp[i]], M[i, start_pp[i]])
    n[PC[i, start_pp[i]]:PC[i, n_prop[i]], start_pp[i]] ~ dmulti(PA[i, 1:n_prop[i]], M_avail[i, start_pp[i]])

    for(t in (start_pp[i]+1):end_pp[i]){

      M[i, t] ~ dpois(mu[i, t])
      log(mu[i, t]) <- log_lambda[i, t] + log(z[i, t])
      z[i, t] <- M[i, t-1] - y_removed[i, t-1]
      log_lambda[i, t] <- sum(log_lambda_all[i, pp[i, t-1]:(pp[i, t]-1)])

      M_avail[i, t] ~ dbinom(proportion_surveyed[i, t], M[i, t])
      n[PC[i, 1]:PC[i, n_prop[i]], t] ~ dmulti(PA[i, 1:n_prop[i]], M_avail[i, t])

    }
  }

  # calculate log growth rate across timesteps
  for(i in 1:n_units_all){
    log_lambda_mu[county_x[i], pp_x[i]] <- inprod(X[i, 1:n_beta], beta[1:n_beta])
    log_lambda_all[county_x[i], pp_x[i]] ~ dnorm(log_lambda_mu[county_x[i], pp_x[i]], tau = tau_l)
  }

  for(i in 1:n_units){

    # long format for easier monitoring
    xn[i] <- n[property[i], PPNum[i]]

    # likelihood - first pass
    lam[i, 1] <- pi[i, 1] * n[property[i], PPNum[i]]
    y[i, 1] ~ dpois(lam[i, 1])

    pi[i, 1] <- gamma[i, 1] * theta[i, 1]
    theta[i, 1] <- 1 - pow((1 - p[method[i, 1]]), effort[i, 1])

    for(j in 2:n_passes[i]){

      # likelihood - subsequent passes
      lam[i, j] <- pi[i, j] * n[property[i], PPNum[i]]
      y[i, j] ~ dpois(lam[i, j])

      pi[i, j] <- gamma[i, j] * theta[i, j] *
        exp(sum(log((1 - gamma[i, 1:(j-1)]) + (gamma[i, 1:(j-1)] * (1 - theta[i, 1:(j-1)])))))
      theta[i, j] <- 1 - pow((1 - p[method[i, j]]), effort[i, j])

    }

  }

})

