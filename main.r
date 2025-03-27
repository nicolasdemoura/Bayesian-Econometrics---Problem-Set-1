#####################################################################
# Topic: Bayesian Econometrics - Problem Set 1
# Instructor: Eduardo Mendes
# Course: Bayesian Econometrics
# Autor: Nícolas de Moura
# Goal: Implement a Bayesian model to estimate the parameters of a non-linear model
#####################################################################
# Organize the working environment
#####################################################################

# Clean the working environment
rm(list = ls())
load.lib <- c("dplyr", "ipumsr", "ggplot2", "splines", "stargazer", "Hmisc", "AER","readxl", "tidyverse", "data.table", "stargazer", "lubridate", "fixest", "ggplot2", "pracma", "dplyr", "remotes", "tidyr", "mvProbit", "ipw", "MASS", "xtable", "quantreg", "nprobust", "chron", "WDI")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib, require, character=TRUE)
remove(install.lib, lib, load.lib)

# Set the seed
set.seed(20250324)

# Import the functions
source("functions.R")

#####################################################################
# Simulate the data
#####################################################################

# Set the parameters
n <- 100000
taus <- c(6,12,24,30,36,48,60,72,84,96,108,180,240,360)
tau <- sample(taus, n, replace = TRUE)
beta <- c(5, -3, 2)
lambda <- 1
sigma2 <- 0.1

# Simulate the data
y <- nelson_siegel(tau, beta, lambda) + rnorm(n, sd = sigma2)

#####################################################################
# Run the MCMC
#####################################################################

# Set the simulation parameters
n_iter <- 10000
burn_in <- 0.25
c_proposal <- 0.0001
sigma_proposal <- 1 * diag(c(1,1,1,1)) * c_proposal

# Run MCMC
metropolis <- metropolis_hastings(y, tau, n_iter = n_iter, sigma_proposal=sigma_proposal)
samples <- metropolis[,1:4]
acceptance <- metropolis[,5]

samples <- as.data.frame(samples)
colnames(samples) <- c("Beta_0", "Beta_1", "Beta_2", "Lambda")

# Discard the burn-in
samples <- samples[round(burn_in * n_iter):n_iter, ]

#####################################################################
# Estimate the parameters
#####################################################################

beta_hat <- colMeans(samples[,1:3])
lambda_hat <- mean(samples$Lambda)

#####################################################################
# Plot the results
#####################################################################

# Plot and save the histograms of the parameters
plot_hist(samples, " ")
plot_time_series(samples, " ")

# Plot the Nelson-Siegel curve 
plot_nelson_siegel(beta_hat,lambda_hat, y, tau, taus, " ")

print(beta_hat)
print(lambda_hat)
print(mean(acceptance))

# Print the correlation matrix
samples_corr <- samples
samples_corr$Lambda <- log(samples$Lambda)
samples_corr <- cov(samples)
print(samples_corr)
