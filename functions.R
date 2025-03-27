###############################################################################
# Functions 
###############################################################################

###############################################################################
# Prior, Likelihood, and Posterior
###############################################################################

# Prior
log_prior <- function(beta, lambda) {
    lp <- sum(dnorm(beta, mean = c(0, 0, 0), sd = sqrt(sigma2*100), log = TRUE))
    lp <- lp + dnorm(log(lambda), mean = 0, sd = sqrt(sigma2*100), log = TRUE)
    return(lp)  # Prior up to constant
}

# Log-likelihood
log_likelihood <- function(beta, lambda) {
    ll <- sum(dnorm(y, mean = nelson_siegel(tau, beta, lambda), sd = sqrt(sigma2), log = TRUE))
    return(ll)  # Likelihood up to constant
}

# Log-posterior
log_posterior <- function(beta, lambda) {
    return(log_prior(beta, lambda) + log_likelihood(beta, lambda))
}

###############################################################################
# Simulate the data
###############################################################################

# Nelson-Siegel Model
nelson_siegel <- function(tau, beta, lambda){
    # Tau : Time to maturity
    # beta[1] : Beta_0 - Level
    # beta[2] : Beta_1 - Slope
    # beta[3] : Beta_2 - Curvature
    # lambda : lambda - Factor of decay

    # Compute the Nelson-Siegel model
    y <- beta[1] + beta[2] * (1 - exp(-tau / lambda)) / (tau / lambda) + beta[3] * ((1 - exp(-tau / lambda)) / (tau / lambda) - exp(-tau / lambda))
}

# Metropolis-Hastings Algorithm
metropolis_hastings <- function(y, tau, n_iter = 100, sigma_proposal = diag(c(1,1,1,1))) {
    p <- 4

    beta_current <- rnorm(p - 1, mean = 0, sd = 1)
    lambda_current <- exp(rnorm(1, mean = 0, sd = 1))

    samples <- matrix(NA, nrow = n_iter, ncol = p+1)
    
    for (i in 1:n_iter) {
      
        step <- mvrnorm(1, mu = rep(0, p), Sigma = sigma_proposal)

        beta_proposal <- beta_current + step[1:(p-1)]
        lambda_proposal <- exp(log(lambda_current) + step[p])
        accepted_proposal <- 0

        log_acceptance_ratio <- log_posterior(beta_proposal, lambda_proposal) - log_posterior(beta_current, lambda_current)
            
        if (is.nan(log_acceptance_ratio)) {
            print(paste("Iteration", i, ":", "NaN"))
            print(paste("     beta_proposal, lambda_proposal : ", beta_proposal[1], beta_proposal[2], beta_proposal[3], lambda_proposal))
            print(paste("     beta_current, lambda_current : ", beta_current[1], beta_current[2], beta_current[3], lambda_current))
            print(paste("     log_posterior(beta_proposal, lambda_proposal) : ", log_posterior(beta_proposal, lambda_proposal)))
            print(paste("           log_prior(beta_proposal, lambda_proposal) : ", log_prior(beta_proposal, lambda_proposal)))
            print(paste("           log_likelihood(beta_proposal, lambda_proposal) : ", log_likelihood(beta_proposal, lambda_proposal)))
            print(paste("     log_posterior(beta_current, lambda_current) : ", log_posterior(beta_current, lambda_current)))
            print(paste("           log_prior(beta_current, lambda_current) : ", log_prior(beta_current, lambda_current)))
            print(paste("           log_likelihood(beta_current, lambda_current) : ", log_likelihood(beta_current, lambda_current)))
            log_acceptance_ratio <- - Inf 
        }
        
        if(log(runif(1)) < log_acceptance_ratio) {
            beta_current <- beta_proposal
            lambda_current <- lambda_proposal
            accepted_proposal <- 1
        }
        samples[i, 1:3] <- beta_current
        samples[i, 4] <- lambda_current
        samples[i, 5] <- accepted_proposal
    }
    
    return(samples)
}

###############################################################################
# Plot the results
###############################################################################

# Plot the histogram of the data using ggplot2
plot_hist <- function(data, title) {
    for(i in 1:length(data)){
        # Plot the histogram

        if(title != " "){
            gg <- ggplot(data, aes(x = data[, i])) +
                geom_histogram(bins = 50, fill = "#1E88E5", color = "black", alpha = 0.5) +
                labs(x = colnames(data)[i],
                    y = "Frequency",
                    title = paste("Histogram of", colnames(data)[i])) +
                theme_bw(base_size = 25) +
                theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
                    legend.position = "none",
                    legend.text = element_text(size = 15),
                    legend.title = element_text(size = 16),
                    legend.key.size = unit(1, "cm"),
                    legend.background = element_rect(color = "black", size = 0.5)) +
                coord_cartesian()
        } else {
            gg <- ggplot(data, aes(x = data[, i])) +
                geom_histogram(bins = 50, fill = "#1E88E5", color = "black", alpha = 0.5) +
                labs(x = colnames(data)[i],
                    y = "Frequency") +
                theme_bw(base_size = 25) +
                theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
                    legend.position = "none",
                    legend.text = element_text(size = 15),
                    legend.title = element_text(size = 16),
                    legend.key.size = unit(1, "cm"),
                    legend.background = element_rect(color = "black", size = 0.5)) +
                coord_cartesian()
        }
        ggsave(paste("figures/hist_", colnames(data)[i], ".png", sep = ""), gg, width = 10, height = 7)
    }
}

# Plot the time-series of the data using ggplot2
plot_time_series <- function(data, title) {
    for(i in 1:length(data)){
        # Plot the time-series
        if(title != " "){
            gg <- ggplot(data, aes(x = 1:nrow(data), y = data[, i])) +
                geom_line(size = 1.5, color = "#1E88E5") +
                labs(x = "Iteration",
                    y = colnames(data)[i],
                    title = paste("Time-series of", colnames(data)[i])) +
                theme_bw(base_size = 25) +
                theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
                    legend.position = "none",
                    legend.text = element_text(size = 15),
                    legend.title = element_text(size = 16),
                    legend.key.size = unit(1, "cm"),
                    legend.background = element_rect(color = "black", size = 0.5)) +
                coord_cartesian()
        } else {
            gg <- ggplot(data, aes(x = 1:nrow(data), y = data[, i])) +
                geom_line(size = 1.5, color = "#1E88E5") +
                labs(x = "Iteration",
                    y = colnames(data)[i]) +
                theme_bw(base_size = 25) +
                theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
                    legend.position = "none",
                    legend.text = element_text(size = 15),
                    legend.title = element_text(size = 16),
                    legend.key.size = unit(1, "cm"),
                    legend.background = element_rect(color = "black", size = 0.5)) +
                coord_cartesian()
        }
        ggsave(paste("figures/ts_", colnames(data)[i], ".png", sep = ""), gg, width = 10, height = 7)
    }
}

# Plot the Nelson-Siegel curve and the data points
plot_nelson_siegel <- function(beta_hat,lambda_hat, y, tau, taus, title) {
    # Plot the Nelson-Siegel curve and the data points
    if(title != " "){            
        gg <- ggplot(data.frame(tau = taus, y = nelson_siegel(taus, beta_hat,lambda_hat)), aes(x = tau, y = y)) +
            geom_line(size = 1.5, color = "#1E88E5") +
            geom_point(aes(y = y), data = data.frame(tau = tau, y = y), color = "red", size = 3) +
            theme_bw(base_size = 25) +
            theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
                legend.position = "none",
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 16),
                legend.key.size = unit(1, "cm"),
                legend.background = element_rect(color = "black", size = 0.5)) +
            labs(title = "Nelson-Siegel Curve", x = "Time to Maturity", y = "Yield") +
            coord_cartesian()
        print(gg)
        ggsave("figures/nelson_siegel_curve.png", gg, width = 10, height = 10)
    } else {            
        gg <- ggplot(data.frame(tau = taus, y = nelson_siegel(taus,beta_hat,lambda_hat)), aes(x = tau, y = y)) +
            geom_line(size = 1.5, color = "#1E88E5") +
            geom_point(aes(y = y), data = data.frame(tau = tau, y = y), color = "red", size = 3) +
            theme_bw(base_size = 25) +
            theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
                legend.position = "none",
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 16),
                legend.key.size = unit(1, "cm"),
                legend.background = element_rect(color = "black", size = 0.5)) +
            labs(x = "Time to Maturity", y = "Yield") +
            coord_cartesian()
        print(gg)
        ggsave("figures/nelson_siegel_curve.png", gg, width = 10, height = 10)
    }
}
