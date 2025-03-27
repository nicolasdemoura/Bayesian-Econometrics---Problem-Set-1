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
n <- 1000
taus <- c(6,12,24,30,36,48,60,72,84,96,108,180,240,360)
tau <- sample(taus, n, replace = TRUE)
beta <- c(3, 10, 7)
lambda <- 1

# Simulate the data
y <- nelson_siegel(tau, c(beta, lambda)) + rnorm(n, sd = 0)

# Set the number of iterations
n_iter <- 10000
burn_in <- 0.1

# Run MCMC
samples <- metropolis_hastings(y, tau, n_iter = n_iter, proposal_sd = 0.2
samples <- as.data.frame(samples)
colnames(samples) <- c("Beta_0", "Beta_1", "Beta_2", "Lambda")

# Discard the burn-in
samples <- samples[round(burn_in * n_iter):n_iter, ]

#####################################################################
# Estimate the parameters
#####################################################################

gamma_hat <- colMeans(samples)

#####################################################################
# Plot the results
#####################################################################

# Plot and save the histograms of the parameters
plot_hist(samples, " ")
plot_time_series(samples, " ")

# Plot the Nelson-Siegel curve 
plot_nelson_siegel(gamma_hat, y, tau, taus, " ")

print(gamma_hat)
