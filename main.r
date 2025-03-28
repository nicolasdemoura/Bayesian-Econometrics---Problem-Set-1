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
load.lib <- c("dplyr", "ipumsr", "ggplot2", "splines", "stargazer", "Hmisc", "AER","readxl", "tidyverse", "data.table", "stargazer", "lubridate", "fixest", "ggplot2", "pracma", "dplyr", "remotes", "tidyr", "mvProbit", "ipw", "MASS", "xtable", "quantreg", "nprobust", "chron", "WDI", "progress")
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

# Import the data
raw_data <- read_excel("data/Yields_Bloomberg.xlsx", sheet = "weekly yield")

# Transform the data to long format
temp_data <- raw_data %>%
  pivot_longer(cols = -DATA, names_to = "Maturity", values_to = "Yield") %>% # Replave the name to remove the TREASURY at the beginning of the column, change the , to . and read as.numeric
  mutate(Maturity = gsub("TREASURY", "", Maturity)) %>%
  mutate(Maturity = 12*as.numeric(gsub(",", ".", Maturity)))

# For each group of 26 weeks among the unique(temp_data$DATA), run the algorithm with the data
# and plot the Nielsen-Siegel curve
for (i in seq(1, length(unique(temp_data$DATA)), by = 8)) {
    # Select the data for the current group of 8 weeks
    data <- temp_data %>% filter(DATA %in% unique(temp_data$DATA)[i:(i+7)])

    # Run the algorithm with the data
    compute_nelson_siegel(data, i)

}
