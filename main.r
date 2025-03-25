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
# Set the function to simulate the data
###############################################################################

nelson_siegel <- function(tau, gamma){
    # Tau : Time to maturity
    # gamma[1] : Beta_0 - Level
    # gamma[2] : Beta_1 - Slope
    # gamma[3] : Beta_2 - Curvature
    # gamma[4] : lambda - Factor of decay

    # Compute the Nelson-Siegel model
    y <- gamma[1] + gamma[2] * (1 - exp(-tau / gamma[4])) / (tau / gamma[4]) + gamma[3] * ((1 - exp(-tau / gamma[4])) / (tau / gamma[4]) - exp(-tau / gamma[4]))
}


###############################################################################
# Simulate the data
###############################################################################

# Set the parameters
n <- 1000
tau <- c(1:120)
beta <- c(0, 10, 0)
lambda <- 1

# Simulate the data
y <- nelson_siegel(tau, c(beta, lambda)) + rnorm(n, sd = 0.1)

# Metropolis-Hastings Algorithm
metropolis_hastings <- function(y, tau, n_iter = 100, sigma2 = 1, proposal_sd = 1) {
    p <- 4
    gamma_current <- rep(0, p)  # Initialize
    samples <- matrix(NA, nrow = n_iter, ncol = p)
    
    log_posterior <- function(gamma) {
        ll <- sum(dnorm(y, mean = nelson_siegel(tau, gamma), sd = sqrt(sigma2), log = TRUE))
        lp <- sum(dnorm(gamma, mean = 0, sd = sqrt(sigma2), log = TRUE))
        return(ll + lp)  # Posterior up to constant
    }
    
    for (i in 1:n_iter) {
        gamma_proposal <- gamma_current + rnorm(p, mean = 0, sd = proposal_sd)
        log_acceptance_ratio <- log_posterior(gamma_proposal) - log_posterior(gamma_current)
            
        if (is.nan(log_acceptance_ratio)) {
            log_acceptance_ratio <- - Inf 
        }
        
        if(log(runif(1)) < log_acceptance_ratio) {
            gamma_current <- gamma_proposal
        }
        samples[i, ] <- gamma_current
    }
    
    return(samples)
}

# 
n_iter <- 10000

# Run MCMC
samples <- metropolis_hastings(y, tau, n_iter = n_iter, proposal_sd = 0.1)

# Discard the burn-in
burn_in <- 0.1
samples <- samples[round(burn_in * n_iter):n_iter, ]

# Plot the results
par(mfrow = c(2, 4))
for(i in 1:4){
    hist(samples[,i], main = paste("gamma", i), xlab = paste("gamma", i))
    plot (samples[,i], type = "l", main = paste("gamma", i), xlab = "Iteration", ylab = paste("gamma", i))
}



# Plot the Nelson-Siegel curve
#plot(tau, y, pch = 20, col = "blue", xlab = "Time to Maturity", ylab = "Yield")

