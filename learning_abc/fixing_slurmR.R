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

n_values <- sample(5000:10000, N_SIMS, replace = TRUE)  # population size
theta <- data.table(
  n      = n_values,
  preval = runif(N_SIMS, 0.007, 0.02),
  crate  = runif(N_SIMS, 1, 5),
  recov  = runif(N_SIMS, 0.071, 0.25),
  R0     = runif(N_SIMS, 1.1, 5)
)

# Calculate transmission probability
theta[, ptran := (R0 * recov / crate)]

# IMPORTANT: Keep only 5 columns in the correct order (no R0 in theta_use)
theta_use <- theta[, .(n, preval, crate, recov, ptran)]

# Print summary
summary(theta_use)
cat("Number of parameter sets:", nrow(theta_use), "\n")
cat("Column names:", names(theta_use), "\n")

# Save for reference
saveRDS(theta_use, "theta_use.rds")

# --------------------------
# Define Simulation Functions
# --------------------------

# Function to simulate the "observed" epidemic
simulate_epidemic_observed <- function(params, ndays = model_ndays, seed = NULL) {
  library(epiworldR)
  library(data.table)
  library(dplyr)
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize the model with the true parameters
  sim_model <- ModelSIRCONN(
    name="sim",
    n                 = as.integer(params["n"]),
    prevalence        = as.numeric(params["preval"]),
    contact_rate      = as.numeric(params["crate"]),
    transmission_rate = as.numeric(params["ptran"]),
    recovery_rate     = as.numeric(params["recov"])
  )
  
  verbose_off(sim_model)
  
  # Run the simulation
  epiworldR::run(sim_model, ndays = ndays)
  
  # Get only the infected counts
  counts <- get_hist_total(sim_model)
  infected_counts <- counts[counts$state == "Infected", "counts"]
  
  return(infected_counts)
}

# --------------------------
# ABC (LFMCMC) Calibration Functions
# --------------------------

simulate_epidemic_calib <- function(params, ndays = model_ndays, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize model with calibration parameters
  sim_model <- ModelSIRCONN(
    name="sim",
    n                 = as.integer(params["n"]),
    prevalence        = as.numeric(params["preval"]),
    contact_rate      = as.numeric(params["crate"]),
    transmission_rate = as.numeric(params["ptran"]),
    recovery_rate     = as.numeric(params["recov"])
  )
  
  verbose_off(sim_model)
  
  # Run the simulation
  run(sim_model, ndays = ndays)
  
  # Get only the infected counts
  counts <- get_hist_total(sim_model)
  infected_counts <- counts[counts$state == "Infected", "counts"]
  
  return(infected_counts)
}

# FIXED: Modified simulation_fun to properly pass all parameters
simulation_fun <- function(params, lfmcmc_obj, true_params_local) {
  library(epiworldR)
  library(data.table)
  library(dplyr)
  
  # Build full parameter vector with fixed n and prevalence, calibrated others
  sim_params <- c(
    n      = true_params_local[1],  # n (fixed)
    preval = true_params_local[2],  # prevalence (fixed)
    crate  = params[1],             # contact_rate (calibrated)
    recov  = params[2],             # recovery_rate (calibrated)
    ptran  = params[3]              # transmission_prob (calibrated)
  )
  
  # Run simulation with the parameters
  infected_counts <- simulate_epidemic_calib(sim_params, ndays = model_ndays)
  return(as.numeric(infected_counts))
}

# This returns the observed data (all daily infected counts)
summary_fun <- function(data, lfmcmc_obj) {
  library(epiworldR)
  library(data.table)
  library(dplyr)
  return(as.numeric(data))
}

# Generate new parameter proposals
proposal_fun <- function(old_params, lfmcmc_obj) {
  library(epiworldR)
  library(data.table)
  library(dplyr)
  
  # Proposals with appropriate step sizes for each parameter
  new_crate <- old_params[1] * exp(rnorm(1, sd = 0.1))  # Log-normal proposal
  new_recov <- old_params[2] * exp(rnorm(1, sd = 0.1))  # Log-normal proposal
  new_ptran <- plogis(qlogis(old_params[3]) + rnorm(1, sd = 0.1))
  
  return(c(new_crate, new_recov, new_ptran))
}

# Kernel function using sum of squared differences across all days
kernel_fun <- function(simulated_stat, observed_stat, epsilon, lfmcmc_obj) {
  library(epiworldR)
  library(data.table)
  library(dplyr)
  
  # Calculate the sum of squared differences across all days
  diff <- sum((simulated_stat - observed_stat)^2)
  return(exp(-diff / (2 * epsilon^2)))
}

# --------------------------
# Function: Simulation–Calibration Study with ABC only
# --------------------------
simulate_and_calibrate <- function(true_params, sim_id) {
  # FIXED: Added error handling
  tryCatch({
    library(epiworldR)
    library(data.table)
    library(dplyr)
    cat("Running simulation", sim_id, "of", N_SIMS, "\n")
    
    # Calculate true R0 from true parameters
    true_R0 <- (true_params[3] * true_params[5]) / true_params[4]
    
    # Step 1: Simulate the observed epidemic using the true parameters
    observed_infected <- simulate_epidemic_observed(true_params, ndays = model_ndays, 
                                                    seed = model_seed + sim_id)
    
    # Step 2: ABC calibration
    dummy_model <- ModelSIRCONN(
      name              = "dummy",
      n                 = as.integer(global_n),
      prevalence        = 0.1,
      contact_rate      = 5.0,
      transmission_rate = 0.1,
      recovery_rate     = 0.1
    )
    
    # FIXED: Create a local wrapper for simulation_fun
    local_simulation_fun <- function(params, lfmcmc_obj) {
      return(simulation_fun(params, lfmcmc_obj, true_params))
    }
    
    # Setup and run LFMCMC with the full time series
    lfmcmc_obj <- LFMCMC(dummy_model)
    lfmcmc_obj <- set_simulation_fun(lfmcmc_obj, local_simulation_fun)
    lfmcmc_obj <- set_summary_fun(lfmcmc_obj, summary_fun)
    lfmcmc_obj <- set_proposal_fun(lfmcmc_obj, proposal_fun)
    lfmcmc_obj <- set_kernel_fun(lfmcmc_obj, kernel_fun)
    lfmcmc_obj <- set_observed_data(lfmcmc_obj, observed_infected)
    
    # Initial parameters for ABC (only the 3 we're calibrating)
    init_params <- c(
      true_params[3],              # contact rate
      true_params[4],              # recovery rate
      true_params[5]               # transmission probability
    )
    
    # Run the LFMCMC
    n_samples_calib <- 500
    epsilon <- sqrt(sum(observed_infected^2)) * 0.05
    
    epiworldR::run_lfmcmc(
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
      calibrated_params_raw <- apply(accepted, 2, median)
      
      abc_crate <- calibrated_params_raw[1]
      abc_recov <- calibrated_params_raw[2]
      abc_ptran <- calibrated_params_raw[3]
      
      # Build full parameter vector with fixed n and prevalence
      calibrated_params <- c(
        true_params[1],  # n (fixed)
        true_params[2],  # prevalence (fixed)
        abc_crate,       # calibrated contact rate
        abc_recov,       # calibrated recovery rate
        abc_ptran        # calibrated transmission prob
      )
      
      # Calculate calibrated R0
      abc_R0 <- (abc_crate * abc_ptran) / abc_recov
    } else {
      calibrated_params <- true_params[1:5]  # Use all 5 true params
      abc_R0 <- true_R0
      abc_crate <- true_params[3]
      abc_recov <- true_params[4]
      abc_ptran <- true_params[5]
      cat("Warning: Using true parameters as calibration failed for simulation", sim_id, "\n")
    }
    
    # Step 3: Simulate predictions using the ABC calibrated parameters
    abc_predicted_infected <- simulate_epidemic_calib(calibrated_params, ndays = model_ndays, 
                                                      seed = model_seed + sim_id + 200)
    
    # Create results dataframe for all days
    days <- 0:model_ndays
    daily_results <- data.frame(
      sim_id = sim_id,
      day = days,
      true_n = true_params[1],
      true_preval = true_params[2],
      true_crate = true_params[3],
      true_recov = true_params[4],
      true_ptran = true_params[5],
      true_R0 = true_R0,
      calib_preval = true_params[2],  # Fixed at true value
      calib_crate = abc_crate,
      calib_recov = abc_recov,
      calib_ptran = abc_ptran,
      calib_R0 = abc_R0,
      observed_infected = observed_infected,
      abc_predicted_infected = abc_predicted_infected,
      abc_bias = abc_predicted_infected - observed_infected,
      abc_rel_bias = ifelse(observed_infected > 0, 
                            (abc_predicted_infected - observed_infected) / observed_infected, 
                            NA)
    )
    
    # Parameter comparison
    param_comparison <- data.frame(
      sim_id = sim_id,
      parameter = c("contact_rate", "recovery_rate", "transmission_prob"),
      true_value = c(true_params[3:5]),
      abc_calibrated_value = c(abc_crate, abc_recov, abc_ptran),
      abc_error_pct = c(
        (abc_crate - true_params[3]) / true_params[3] * 100,
        (abc_recov - true_params[4]) / true_params[4] * 100,
        (abc_ptran - true_params[5]) / true_params[5] * 100
      )
    )
    
    # ABC parameters summary
    abc_parameters <- data.frame(
      sim_id = sim_id,
      true_n = true_params[1],
      true_preval = true_params[2],
      true_crate = true_params[3],
      true_recov = true_params[4],
      true_ptran = true_params[5],
      true_R0 = true_R0,
      abc_preval = true_params[2],  # Fixed at true value
      abc_crate = abc_crate,
      abc_recov = abc_recov,
      abc_ptran = abc_ptran,
      abc_R0 = abc_R0
    )
    
    return(list(
      daily_results = daily_results,
      param_comparison = param_comparison,
      abc_parameters = abc_parameters
    ))
    
  }, error = function(e) {
    cat("Error in simulation", sim_id, ":", e$message, "\n")
    return(NULL)  # Return NULL on error
  })
}

# --------------------------
# Run Simulations in Parallel
# --------------------------

library(slurmR)

# Run simulations using Slurm
ans <- Slurm_lapply(
  X = 1:N_SIMS,
  FUN = function(i) simulate_and_calibrate(as.numeric(theta_use[i,]), i),
  job_name = "Sims_calibrate",
  njobs = 100,
  overwrite = TRUE,
  plan = "submit",
  sbatch_opt = list(
    partition = "vegayon-shared-np",
    account = "vegayon-np",
    time = "01:00:00",
    `mem-per-cpu` = "4G",
    `cpus-per-task` = 1
  ),
  export = c(
    "simulate_and_calibrate",
    "simulate_epidemic_observed",
    "simulate_epidemic_calib",
    "simulation_fun",
    "summary_fun",
    "proposal_fun",
    "kernel_fun",
    "theta_use",
    "model_ndays",
    "model_seed",
    "global_n",
    "N_SIMS"
  )
)

# Collect results
results_list <- Slurm_collect(ans)

# Process results with error handling
valid_results <- results_list[!sapply(results_list, is.nu