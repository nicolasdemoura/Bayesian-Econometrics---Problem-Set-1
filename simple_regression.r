###############################################################################
# Topic: Bayesian Econometrics - Problem Set 1
# Instructor: Eduardo Mendes
# Course: Bayesian Econometrics
# Autor: Nícolas de Moura
# Goal: Implement a Bayesian model to estimate the parameters of a non-linear model
###############################################################################
# Organize the working environment
###############################################################################

# Clean the working environment
rm(list = ls())
load.lib <- c("dplyr", "ipumsr", "ggplot2", "splines", "stargazer", "Hmisc", "AER","readxl", "tidyverse", "data.table", "stargazer", "lubridate", "fixest", "ggplot2", "pracma", "dplyr", "remotes", "tidyr", "mvProbit", "ipw", "MASS", "xtable", "quantreg", "nprobust", "chron", "WDI")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib, require, character=TRUE)
remove(install.lib, lib, load.lib)

# Set the seed
set.seed(20250324)

###############################################################################
# Understanding the metropolis algorithm for a simple Bayesian regression model
###############################################################################
set.seed(42)

# Simulate Data
n <- 1000
p <- 2
X <- cbind(1, rnorm(n))  # Design matrix with intercept
beta_true <- c(7, 2)
sigma <- 10
y <- X %*% beta_true + rnorm(n, sd=sigma)

# Prior settings
tau2 <- 100  # Prior variance for beta

# Burn-in
burn_in <- 1000

# Metropolis-Hastings Algorithm
metropolis_hastings <- function(X, y, n_iter = 10000, sigma2 = 100, proposal_sd = 1) {
  p <- ncol(X)
  beta_current <- rep(0, p)  # Initialize
  samples <- matrix(NA, nrow = n_iter, ncol = p)
  
  log_posterior <- function(beta) {
    ll <- sum(dnorm(y, mean = X %*% beta, sd = sqrt(sigma2), log = TRUE))
    lp <- sum(dnorm(beta, mean = 0, sd = sqrt(tau2), log = TRUE))
    return(ll + lp)  # Posterior up to constant
  }
  
  for (i in 1:n_iter) {
    beta_proposal <- beta_current + rnorm(p, mean = 0, sd = proposal_sd)
    log_acceptance_ratio <- log_posterior(beta_proposal) - log_posterior(beta_current)
    
    if (log(runif(1)) < log_acceptance_ratio) {
      beta_current <- beta_proposal
    }
    samples[i, ] <- beta_current
  }
  
  return(samples)
}

# Run MCMC
samples <- metropolis_hastings(X, y, n_iter = 50000)

# Plot posterior distributions
par(mfrow = c(1,2))
hist(samples[burn_in:50000,1], breaks = 30, main = "Posterior of beta0", xlab = "beta0")
hist(samples[burn_in:50000,2], breaks = 30, main = "Posterior of beta1", xlab = "beta1")