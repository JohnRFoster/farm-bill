tau_2_sigma <- function(tau){
  1/sqrt(tau)
}

sigma_2_tau <- function(sigma){
  1/sigam^2
}

sigma_2_var <- function(sigma){
  sigma^2
}

var_2_sigma <- function(var){
  sqrt(var)
}

var_2_tau <- function(var){
  1/var
}

tau_2_var <- function(tau){
  1/tau
}

sum_properties <- nimbleFunction(

  run = function(
    property = double(1),  # property vector index
    z = double(1)          # latent abundance vector
  ){
    N <- sum(z[property])
    return(N)
    returnType(double(0))
   }
)

# C_sum_properties <- compileNimble(sum_properties)
# z <- matrix(rpois(12, 10), 4, 3)
# p <- c(1, 3)
# pp <- 2
# z
# C_sum_properties(property = p, PPNum = pp, z = z)
# z
# sum(A[p, pp])


rDM <- nimbleFunction(
  run = function(
    n = integer(0),
    x1 = double(0),
    phi = double(0),
    lambda = double(0),
    zeroProb = double(0)
  ){

    returnType(double(0))
    if(n != 1) print("rDM only allows n = 1; using n = 1.")

    s <- rbinom(1, x1, phi)

    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(s)
    return(s + rpois(1, lambda))

  }
)

# phi <- rep(0.8, 20)
# rho <- rep(0.3, 20)
# x1 <- 12
# rDM(n = 1, x1 = x1, phi = phi, rho = rho)
# C_rDM <- compileNimble(rDM)
# C_rDM(n = 1, x1 = x1, phi = phi, rho = rho)


# dDM <- nimbleFunction(
#   run = function(x = double(0), x1 = double(0), phi = double(0), rho = double(0),
#                  log = logical(0, default = 0)){
#
#     returnType(double(0))
#
#     logProb <- lfactorial(x1) -
#       (lfactorial(x) + lfactorial(x1 - x)) +
#       log(phi) * x +
#       log(1 - phi) * (x1 - x) +
#       log(rho) * x -
#       rho -
#       lfactorial(x)
#
#     if(log){
#       return(logProb)
#     } else {
#       return(exp(logProb))
#     }
#
#   }
# )


dDM <- nimbleFunction(
  run = function(
    x = double(0),
    x1 = double(0),
    phi = double(0),
    lambda = double(0),
    zeroProb = double(0),
    log = logical(0, default = 0)
  ) {

    returnType(double(0))

    log_binom <- dbinom(x, x1, phi, log = TRUE)

    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb) + log_binom)
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dpois(x, lambda, log = FALSE) * exp(log_binom))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero) + log_binom)
    return(totalProbZero * exp(log_binom))

  }
)

# x <- rDM(n = 1, x = x1, phi = phi[1:5], rho = rho[1:5])
# dDM(x = x, x1 = x1, phi = phi[1:5], rho = rho[1:5], log = 0)
# dDM(x = x, x1 = x1, phi = phi[1:5], rho = rho[1:5], log = 1)

registerDistributions(list(
  dDM = list(
    BUGSdist = "dDM(x1, phi, lambda, zeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c("value = double(0)",
              "x1 = double(0)",
              "phi = double(0)",
              "lambda = double(0)",
              "zeroProb = double(0)")
  )))

# code <- nimbleCode({
#   y ~ dDM(20, phi[1:n], rho[1:n])
#   for(i in 1:n){
#     phi[i] ~ dunif(0, 1)
#     rho[i] ~ dunif(0, 2)
#   }
# })
#
# model <- nimbleModel(code, constants = list(n = 2))
# model <- nimbleModel(code, constants = list(n = 3))
# model <- nimbleModel(code, constants = list(n = 4))

# p <- 0.4
# q <- 1 - p
# n <- 20
# lambda <- n
# x <- 9
#
# # binomial pmf
# factorial(n)/(factorial(x)*factorial(n-x))*(p^x)*(q^(n-x))
# dbinom(x, n, p)
#
# # binomial pmf log scale
# lfactorial(n) - (lfactorial(x) + lfactorial(n-x)) + log(p)*x + log(q)*(n-x)
#
# # poisson pmf
# (lambda^x)*exp(-lambda)/factorial(x)
# dpois(x, lambda)
#
# # poisson pmf log scale
# log(lambda)*x - lambda - lfactorial(x)

rZIP <- nimbleFunction(
  run = function(n = integer(0), lambda = double(0), zeroProb = double(0)) {
    returnType(double())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })


dZIP <- nimbleFunction(
  run = function(x = double(0), lambda = double(0),
                 zeroProb = double(0), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })

registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, zeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = double(0)', 'lambda = double(0)', 'zeroProb = double(0)')
  )))



inits_dm <- function(){

  tau_surv <- NA
  mu_phi_prop <- NA
  if(phi_beta){
    beta_p <- matrix(rnorm(n_property * n_beta, ilogit(phi_prior[1]), phi_prior[2]), n_property, n_beta)
    phi_all <- matrix(NA, constants$n_property, max(constants$timestep_x))
    for(i in 1:constants$n_units_all){
      s <- inprod(data$Xall[i,], beta_p[constants$property_x[i],])
      phi_all[constants$property_x[i], constants$timestep_x[i]] <- ilogit(s)
    }
    mu_phi <- NA
  } else {
    beta_p <- NA
    if(phi_data){
      mu_phi <- rnorm(constants$n_units_all, phi_prior[1], phi_prior[2])
      phi_all <- matrix(NA, constants$n_property, max(constants$timestep_x))
      for(i in 1:constants$n_units_all){
        phi_all[constants$property_x[i], constants$timestep_x[i]] <- ilogit(mu_phi[i])
      }
    } else if(phi_static){
      mu_phi <- rnorm(1, phi_prior[1], phi_prior[2])
      phi_all <- matrix(ilogit(mu_phi), constants$n_property, max(constants$timestep_x))
    } else if(phi_hb){
      tau_surv <- rgamma(1, 1, 1)
      mu_phi <- rnorm(1, phi_prior[1], phi_prior[2])
      mu_phi_prop <- rnorm(constants$n_units_all, mu_phi, 1/sqrt(tau_surv))
      phi_all <- matrix(NA, constants$n_property, max(constants$timestep_x))
      for(i in 1:constants$n_units_all){
        phi_all[constants$property_x[i], constants$timestep_x[i]] <- ilogit(mu_phi_prop[i])
      }
    }
  }

  if(rho_beta){
    beta_r <- matrix(rnorm(n_property * n_beta, -1, 1), n_property, n_beta)
    rho_all <- matrix(NA, constants$n_property, max(constants$timestep_x))
    for(i in 1:constants$n_units_all){
      r <- inprod(data$Xall[i,], beta_r[constants$property_x[i],])
      rho_all[constants$property_x[i], constants$timestep_x[i]] <- exp(r)
    }
    log_rho <- NA
  } else {
    if(rho_data){
      rho <- runif(constants$n_units_all, 0, 1)
      rho_all <- matrix(NA, constants$n_property, max(constants$timestep_x))
      for(i in 1:constants$n_units_all){
        rho_all[constants$property_x[i], constants$timestep_x[i]] <- rho[i]
      }
    } else if(rho_static){
      rho <- runif(1)
      rho_all <- matrix(rho, constants$n_property, max(constants$timestep_x))
    }
    log_rho <- log(rho)
    beta_r <- NA
  }

  n_init <- matrix(NA, constants$n_property, max(constants$n_timestep))
  n_init[,1] <- y_removed[,1]
  dm <- y_removed
  dm <- cbind(NA, dm[,-ncol(dm)])
  s <- r <- z <- array(NA, dim = c(constants$n_property, max(constants$n_timestep), max(n_pp_all)+1))

  for(i in 1:constants$n_property){
    for(t in 2:constants$n_timestep[i]){
      z[i, t, constants$pp[i, t-1]] <- dm[i, t]
      for(j in constants$pp[i, t-1]:(constants$pp[i, t]-1)){
        # if(is.na(phi_all[i, j]) | is.na(rho_all[i, j])) stop("NA")
        s[i, t, j] <- rbinom(1, z[i, t, j], phi_all[i, j])
        # s[i, t, j] <- min(s[i, t, j], z[i, t, j]-1)
        r[i, t, j] <- rpois(1, z[i, t, j] * rho_all[i, j])
        z[i, t, j+1] <- s[i, t, j] + r[i, t, j]
      }
      n_init[i, t] <- z[i, t, j]
    }
  }

  a <- 500
  list(
    n = n_init + a,
    # z = z, # don't include z!
    # M = n_init*rpois(1, 3),
    z1 = mean(rowSums(data$y_removed, na.rm = TRUE)) + a,
    # z1 = rowSums(data$y_removed, na.rm = TRUE),
    beta_p = beta_p,
    beta_r = beta_r,
    mu_p = runif(n_methods),
    # mu_p = matrix(runif(n_methods*max(n_timestep), 0, 1), n_methods, max(n_timestep)),
    size = runif(3, 0.01, 5),
    s = s+1,
    r = r,
    log_rho = log_rho,
    mu_phi = mu_phi,
    mu_phi_prop = mu_phi_prop,
    tau_surv = tau_surv
  )
}

inits_lambda <- function(){

}
