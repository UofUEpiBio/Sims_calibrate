# Load libraries
library(epiworldR)
library(data.table)
library(ggplot2)
library(dplyr)

# Get task ID from Slurm
task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cat("Running simulation", task_id, "\n")

# Global Simulation Settings (load from file or define here)
model_ndays <- 30    # simulation duration (days)
model_seed  <- 122   # seed for reproducibility
global_n    <- 10000 # population size (used in calibration)

# Load parameters
theta_use <- readRDS("theta_use.rds")

# --------------------------
# Define Simulation Functions
# --------------------------

# Function to simulate the "observed" epidemic
simulate_epidemic_observed <- function(params, ndays = model_ndays, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize the model with the true parameters
  sim_model <- ModelSIRCONN(
    name              = "sim",
    n                 = as.integer(params[1]),
    prevalence        = params[2],
    contact_rate      = params[3],
    transmission_rate = params[5],
    recovery_rate     = params[4]
  )
  
  verbose_off(sim_model)
  
  # Run the simulation
  run(sim_model, ndays = ndays)
  
  # Get only the infected counts
  counts <- get_hist_total(sim_model)
  infected_counts <- counts[counts$state == "Infected", "counts"]
  
  return(infected_counts)
}

# Function for simulation during calibration
simulate_epidemic_calib <- function(params, ndays = model_ndays, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize model with calibration parameters
  sim_model <- ModelSIRCONN(
    name              = "sim",
    n                 = as.integer(global_n),
    prevalence        = params[1],
    contact_rate      = params[2],
    transmission_rate = params[4],
    recovery_rate     = params[3]
  )
  
  verbose_off(sim_model)
  
  # Run the simulation
  run(sim_model, ndays = ndays)
  
  # Get only the infected counts
  counts <- get_hist_total(sim_model)
  infected_counts <- counts[counts$state == "Infected", "counts"]
  
  return(infected_counts)
}

# --------------------------
# ABC (LFMCMC) Calibration Functions
# --------------------------

# This function simulates using the proposed parameters and returns all daily infected counts
simulation_fun <- function(params, lfmcmc_obj) {
  # Here we use the entire time series of infected counts for comparison
  infected_counts <- simulate_epidemic_calib(params, ndays = model_ndays)
  return(as.numeric(infected_counts))
}

# This returns the observed data (all daily infected counts)
summary_fun <- function(data, lfmcmc_obj) {
  return(as.numeric(data))
}

# Generate new parameter proposals
proposal_fun <- function(old_params, lfmcmc_obj) {
  # Proposals with appropriate step sizes for each parameter
  new_preval <- plogis(qlogis(old_params[1]) + rnorm(1, sd = 0.1))
  new_crate  <- old_params[2] * exp(rnorm(1, sd = 0.1))  # Log-normal proposal
  new_recov  <- old_params[3] * exp(rnorm(1, sd = 0.1))  # Log-normal proposal
  new_ptran  <- plogis(qlogis(old_params[4]) + rnorm(1, sd = 0.1))
  return(c(new_preval, new_crate, new_recov, new_ptran))
}

# Kernel function using sum of squared differences across all days
kernel_fun <- function(simulated_stat, observed_stat, epsilon, lfmcmc_obj) {
  # Calculate the sum of squared differences across all days
  diff <- sum((simulated_stat - observed_stat)^2)
  return(exp(-diff / (2 * epsilon^2)))
}

# Function to calculate confidence interval
calc_CI <- function(values, conf_level = 0.95) {
  # Calculate the mean
  mean_val <- mean(values)
  
  # Calculate standard error (still useful to return)
  n <- length(values)
  se <- sd(values) / sqrt(n)
  
  # Calculate quantile-based confidence intervals
  lower_ci <- quantile(values, 0.025)
  upper_ci <- quantile(values, 0.975)
  
  return(c(mean = mean_val, lower_ci = lower_ci, upper_ci = upper_ci, se = se))
}

# --------------------------
# Function: Simplified Simulation–Calibration Study
# --------------------------
simulate_and_calibrate <- function(true_params, sim_id) {
  cat("Running simulation", sim_id, "\n")
  
  # Step 1: Simulate the observed epidemic using the true parameters
  observed_infected <- simulate_epidemic_observed(true_params, ndays = model_ndays, 
                                                  seed = model_seed + sim_id)
  
  # Create a dummy model for LFMCMC
  dummy_model <- ModelSIRCONN(
    name              = "dummy",
    n                 = as.integer(global_n),
    prevalence        = 0.1,
    contact_rate      = 5.0,
    transmission_rate = 0.1,
    recovery_rate     = 0.1
  )
  
  # Setup and run LFMCMC with the full time series
  lfmcmc_obj <- LFMCMC(dummy_model)
  lfmcmc_obj <- set_simulation_fun(lfmcmc_obj, simulation_fun)
  lfmcmc_obj <- set_summary_fun(lfmcmc_obj, summary_fun)
  lfmcmc_obj <- set_proposal_fun(lfmcmc_obj, proposal_fun)
  lfmcmc_obj <- set_kernel_fun(lfmcmc_obj, kernel_fun)
  lfmcmc_obj <- set_observed_data(lfmcmc_obj, observed_infected)
  
  # Use the true parameters as the initial guess for faster convergence
  init_params <- as.numeric(true_params[2:5])
  
  # Run the LFMCMC
  n_samples_calib <- 500  # Increased from 100 for better convergence
  
  # Adaptive epsilon based on the magnitude of the observed data
  epsilon <- sqrt(sum(observed_infected^2)) * 0.05
  
  run_lfmcmc(
    lfmcmc = lfmcmc_obj,
    params_init = init_params,
    n_samples = n_samples_calib,
    epsilon = epsilon,
    seed = model_seed + sim_id + 100
  )
  
  # Get accepted parameters from the MCMC chain
  accepted <- get_all_accepted_params(lfmcmc_obj)
  
  # Calculate median parameters as our calibrated estimate
  if (!is.null(accepted) && nrow(accepted) > 0) {
    calibrated_params <- apply(accepted, 2, median)
  } else {
    calibrated_params <- true_params[2:5]
    cat("Warning: Using true parameters as calibration failed for simulation", sim_id, "\n")
  }
  
  # Step 3: Simulate predictions using the calibrated parameters
  predicted_infected <- simulate_epidemic_calib(calibrated_params, ndays = model_ndays, 
                                                seed = model_seed + sim_id + 200)
  
  # Create results dataframe for all days
  days <- 0:model_ndays
  daily_results <- data.frame(
    sim_id = sim_id,
    day = days,
    # True parameters (same for all days in this simulation)
    true_preval = true_params[2],
    true_crate = true_params[3],
    true_recov = true_params[4],
    true_ptran = true_params[5],
    # Calibrated parameters (same for all days in this simulation)
    calib_preval = calibrated_params[1],
    calib_crate = calibrated_params[2],
    calib_recov = calibrated_params[3],
    calib_ptran = calibrated_params[4],
    # Observed and predicted counts for each day
    observed_infected = observed_infected,
    predicted_infected = predicted_infected,
    # Bias calculation
    bias = predicted_infected - observed_infected,
    # Relative bias (to account for scale differences)
    rel_bias = ifelse(observed_infected > 0, 
                      (predicted_infected - observed_infected) / observed_infected, 
                      NA)
  )
  
  # Also return a parameters comparison
  param_comparison <- data.frame(
    sim_id = sim_id,
    parameter = c("prevalence", "contact_rate", "recovery_rate", "transmission_prob"),
    true_value = true_params[2:5],
    calibrated_value = calibrated_params,
    error_pct = (calibrated_params - true_params[2:5]) / true_params[2:5] * 100
  )
  
  result <- list(
    daily_results = daily_results,
    param_comparison = param_comparison
  )
  
  # Save the result immediately as a safeguard
  saveRDS(result, paste0("result_", sim_id, ".rds"))
  
  return(result)
}

# Run the simulation for this specific task ID
result <- simulate_and_calibrate(as.numeric(theta_use[task_id]), task_id)

# Result is already saved in the simulate_and_calibrate function
cat("Simulation", task_id, "completed\n")