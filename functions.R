###############################################################################
# Functions 
###############################################################################

###############################################################################
# Simulate the data
###############################################################################

# Nelson-Siegel Model
nelson_siegel <- function(tau, gamma){
    # Tau : Time to maturity
    # gamma[1] : Beta_0 - Level
    # gamma[2] : Beta_1 - Slope
    # gamma[3] : Beta_2 - Curvature
    # gamma[4] : lambda - Factor of decay

    # Compute the Nelson-Siegel model
    y <- gamma[1] + gamma[2] * (1 - exp(-tau / gamma[4])) / (tau / gamma[4]) + gamma[3] * ((1 - exp(-tau / gamma[4])) / (tau / gamma[4]) - exp(-tau / gamma[4]))
}


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
        beta_proposal <- gamma_current[1:p-1] + rnorm(p-1, mean = 0, sd = proposal_sd)
        lambda_proposal <- (sqrt(gamma_current[p]) + rnorm(1, mean = 0, sd = proposal_sd))**2

        gamma_proposal <- c(beta_proposal, lambda_proposal)
        log_acceptance_ratio <- log_posterior(gamma_proposal) - log_posterior(gamma_current)
            
        if (is.nan(log_acceptance_ratio)) {
            print(paste("Iteration", i, ":", "NaN"))
            print(paste("     log_posterior(gamma_proposal) : ", log_posterior(gamma_proposal)))
            print(paste("     log_posterior(gamma_current) : ", log_posterior(gamma_current)))
            print(paste("     gamma_proposal : ", gamma_proposal))
            print(paste("     gamma_current : ", gamma_current))
            log_acceptance_ratio <- - Inf 
        }
        
        if(log(runif(1)) < log_acceptance_ratio) {
            gamma_current <- gamma_proposal
        }
        samples[i, ] <- gamma_current
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
plot_nelson_siegel <- function(gamma_hat, y, tau, taus, title) {
    # Plot the Nelson-Siegel curve and the data points
    if(title != " "){            
        gg <- ggplot(data.frame(tau = taus, y = nelson_siegel(taus, gamma_hat)), aes(x = tau, y = y)) +
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
        gg <- ggplot(data.frame(tau = taus, y = nelson_siegel(taus,gamma_hat)), aes(x = tau, y = y)) +
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
