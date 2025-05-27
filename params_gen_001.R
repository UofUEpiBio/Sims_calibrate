# Load required libraries
library(epiworldR)
library(data.table)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(cowplot)

# --------------------------
# Global Simulation Settings
# --------------------------
model_ndays <- 60   # simulation duration (days)
model_seed  <- 122  # seed for reproducibility
global_n    <- 5000  # population size (used in calibration)
N_SIMS      <- 600 # number of simulations to run

# --------------------------
# Generate Parameter Sets using Theta
# --------------------------
set.seed(model_seed)  # Ensure reproducibility

  # Make sure this matches your intended number of simulations

n_values <- sample(5000:10000, N_SIMS, replace = TRUE)  # population size constant

theta <- data.table(
  n      = n_values,
  preval = sample(100:2000, N_SIMS, replace = TRUE) / n_values,  # prevalence as a fraction
  crate  = runif(N_SIMS, 1, 5),       # contact rate
  recov  = runif(N_SIMS, 0.071, 0.25),# recovery rate (from your stats)
  R0     = runif(N_SIMS, 1.1, 5)      # basic reproduction number
)
# Calculate transmission probability (ptran)
theta[, ptran := (R0 * recov / crate)]

# Final dataset with needed columns
theta_use <- theta[, .(n, preval, crate, recov, R0, ptran)]


# Use only the needed columns
theta_use <- theta[, .(n, preval, crate, recov, ptran)]

# Print the true parameter sets for reference
# cat("True parameter sets:\n")
# print(theta_use)
saveRDS(theta_use, "theta_use.rds")

