modelCode <- nimbleCode({

  for(i in 1:n_property){

    n[i, 1] ~ dpois(mu1[i])
    log(mu1[i]) <- log_mu1[i]
    log_mu1[i] ~ dnorm(0, tau = 0.1)

    for(t in 2:n_timestep[i]){

      # lambda is growth rate (exponential population growth)
      # log(lambda[i]) ~ dnorm(X %*% beta, sigma_l)
      # log(lambda[i, t]) <- lambda_mu[i, t]
      log_lambda[i, t-1] ~ dnorm(0, tau = 0.1)

      n[i, t] ~ dpois(mu[i, t])
      log(mu[i, t]) <- log_mu_proc[i, t]
      log_mu_proc[i, t] ~ dnorm(log_mu[i, t], tau = tau_p)
      log_mu[i, t] <- log_lambda[i, t-1] + log(z[i, t-1])
      z[i, t-1] <- n[i, t-1] - y_removed[i, t-1]
    }
  }

  tau_p ~ dexp(1)

  for(i in 1:n_units){

    # poisson likelihood - first pass
    y[i, 1] ~ dpois(pi[i, 1] * n[property[i], timestep[i]])
    pi[i, 1] <- gamma[i, 1] * theta[i, 1]
    theta[i, 1] <- 1 - (1 - p[method[i, 1]]) ^ effort[i, 1]


    for(j in 2:n_passes[i]){

      # poisson likelihood - subsequent passes
      y[i, j] ~ dpois(pi[i, j] * n[property[i], timestep[i]])
      pi[i, j] <- gamma[i, j] * theta[i, j] *
        exp(sum(log((1 - gamma[i, 1:(j-1)]) + (gamma[i, 1:(j-1)] * (1 - theta[i, 1:(j-1)])))))
      theta[i, j] <- 1 - (1 - p[method[i, j]]) ^ effort[i, j]

    }
  }

  for(k in 1:n_methods){
    mu_p[k] ~ dnorm(0, tau = 1)
    logit(p[k]) <- mu_p[k]
  }

  # for(i in 1:n_beta){
  #   beta[i] ~ dnorm(0, 1)
  # }

})

