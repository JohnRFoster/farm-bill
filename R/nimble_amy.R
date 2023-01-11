modelCode <- nimbleCode({

  ### Priors
  # capture probability by method
  for(k in 1:n_methods){
    mu_p[k] ~ dnorm(0, tau = 1)
    logit(p[k]) <- mu_p[k]
  }

  # process error
  if(include_proc){
    tau_p ~ dexp(1)
  }

  for(i in 1:n_property){

    # growth rate coefficient priors
    if(include_beta){
      for(j in 1:n_beta){
        beta[i, j] ~ dnorm(0, tau = 1)
      }
    }

    n[i, 1] ~ dpois(mu1[i])
    log(mu1[i]) <- log_mu1[i]
    log_mu1[i] ~ dnorm(0, tau = 0.1)

    for(t in 2:n_timestep[i]){

      n[i, t] ~ dpois(mu[i, t])

      if(include_proc){
        log(mu[i, t]) <- log_mu_proc[i, t]
        log_mu_proc[i, t] ~ dnorm(log_mu[i, t], tau = tau_p)
      } else {
        log(mu[i, t]) <- log_mu[i, t]
      }

      log_mu[i, t] <- log_lambda[i, t-1] + log(z[i, t-1])
      z[i, t-1] <- n[i, t-1] - y_removed[i, t-1]

    }
  }

  for(i in 1:n_units){

    # lambda is growth rate (exponential population growth)
    if(include_beta){
      # currently using basis functions that vary by property and time
      log_lambda[property[i], timestep[i]] <- inprod(X[i, 1:n_beta], beta[property[i], 1:n_beta])
    } else {
      log_lambda[property[i], timestep[i]] ~ dnorm(0, tau = 1)
    }

    # poisson likelihood - first pass
    y[i, 1] ~ dpois(pi[i, 1] * n[property[i], timestep[i]])
    pi[i, 1] <- gamma[i, 1] * theta[i, 1] *
      (1 - gamma[i, 1]) + (gamma[i, 1] * (1 - theta[i, 1]))
    theta[i, 1] <- 1 - (1 - p[method[i, 1]]) ^ effort[i, 1]

    for(j in 2:n_passes[i]){

      # poisson likelihood - subsequent passes
      y[i, j] ~ dpois(pi[i, j] * n[property[i], timestep[i]])
      pi[i, j] <- gamma[i, j] * theta[i, j] *
        exp(sum(log((1 - gamma[i, 1:(j-1)]) + (gamma[i, 1:(j-1)] * (1 - theta[i, 1:(j-1)])))))
      theta[i, j] <- 1 - (1 - p[method[i, j]]) ^ effort[i, j]

    }
  }



})

