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
N_SIMS      <- 1000 # number of simulations to run

# --------------------------
# Generate Parameter Sets using Theta
# --------------------------
set.seed(model_seed)  # Ensure reproducibility

  # Make sure this matches your intended number of simulations

n_values <- sample(5000:10000, N_SIMS, replace = TRUE)  # population size constant

theta <- data.table(
  n      = n_values,
  preval = runif(N_SIMS, 0.007, 0.02),
  crate  = runif(N_SIMS, 1, 5),
  recov  = runif(N_SIMS, 0.071, 0.25),
  R0     = runif(N_SIMS, 1.1, 5)
)
theta[, ptran := (R0 * recov / crate)]

# Calculate transmission probability (ptran)

# Final dataset with needed columns
theta_use <- theta[, .(n, preval, crate, recov, ptran)]


# Use only the needed columns
summary(theta_use)
# Print the true parameter sets for reference
# cat("True parameter sets:\n")
# print(theta_use)
saveRDS(theta_use, "theta_use.rds")

