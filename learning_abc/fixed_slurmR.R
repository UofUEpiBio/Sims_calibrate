# Load libraries
library(epiworldR)
library(data.table)
library(dplyr)
library(ggplot2)
library(slurmR)

# Source parameter generator (assumes theta_use is loaded)
source("~/Desktop/Sims_calibrate/params_gen_001.R")
theta_use

# Global settings
model_ndays <- 60
model_seed <- 122
N_SIMS <- nrow(theta_use)  # Make sure N_SIMS is defined

# Function to simulate the observed epidemic
simulate_epidemic_observed <- function(params, ndays = model_ndays, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sim_model <- ModelSIRCONN(
    name             = "sim",
    n                = as.double(params[[1]]),    # ◀ enforce double
    prevalence       = as.double(params[[2]]),    # ◀ enforce double
    contact_rate     = as.double(params[[3]]),    # ◀ enforce double
    transmission_rate= as.double(params[[5]]),    # ◀ enforce double
    recovery_rate    = as.double(params[[4]])     # ◀ enforce double
  )
  verbose_off(sim_model)
  run(sim_model, ndays = ndays)
  
  counts <- get_hist_total(sim_model)
  # Convert integer counts → double
  infected_counts <- as.numeric(counts[counts$state == "Infected", "counts"])
  return(infected_counts)
}

# Function for calibration-time simulation - FIXED to enforce double for all numeric inputs
simulate_epidemic_calib <- function(params, true_params_local, ndays = model_ndays, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  sim_model <- ModelSIRCONN(
    name             = "sim",
    n                = as.double(true_params_local[[1]]),  # ◀ enforce double
    prevalence       = as.double(true_params_local[[2]]),  # ◀ enforce double
    contact_rate     = as.double(params[[1]]),             # ◀ enforce double
    transmission_rate= as.double(params[[2]]),             # ◀ enforce double
    recovery_rate    = as.double(params[[3]])              # ◀ enforce double
  )
  verbose_off(sim_model)
  run(sim_model, ndays = ndays)
  
  counts <- get_hist_total(sim_model)
  # Convert integer counts → double
  infected_counts <- as.numeric(counts[counts$state == "Infected", "counts"])
  return(infected_counts)
}

# ABC Calibration Functions - unchanged except they rely on simulate_epidemic_calib
simulation_fun <- function(params, lfmcmc_obj, true_params_local) {
  simulate_epidemic_calib(params, true_params_local, ndays = model_ndays)
}

summary_fun <- function(data, lfmcmc_obj) {
  return(as.numeric(data))
}

proposal_fun <- function(old_params, lfmcmc_obj) {
  new_crate <- old_params[1] * exp(rnorm(1, sd = 0.1))
  new_ptran <- plogis(qlogis(old_params[2]) + rnorm(1, sd = 0.1))
  new_recov <- old_params[3] * exp(rnorm(1, sd = 0.1))
  return(c(new_crate, new_ptran, new_recov))
}

kernel_fun <- function(simulated_stat, observed_stat, epsilon, lfmcmc_obj) {
  diff <- sum((simulated_stat - observed_stat)^2)
  return(exp(-diff / (2 * epsilon^2)))
}

# Calibration wrapper - COMPLETELY FIXED to enforce double for dummy model and observed data
simulate_and_calibrate <- function(true_params, sim_id) {
  # Coerce entire true_params vector to double
  true_params <- as.double(true_params)
  
  tryCatch({
    cat("Running simulation", sim_id, "of", N_SIMS, "\n")
    
    true_R0 <- (true_params[3] * true_params[5]) / true_params[4]
    observed_infected <- simulate_epidemic_observed(true_params, seed = model_seed + sim_id)
    # Ensure observed_infected is double
    observed_infected <- as.numeric(observed_infected)
    
    # Build a dummy model with all numeric inputs as double
    dummy_model <- ModelSIRCONN(
      name             = "dummy",
      n                = as.double(true_params[1]),    # ◀ enforce double
      prevalence       = as.double(0.1),                # ◀ enforce double
      contact_rate     = as.double(5.0),                # ◀ enforce double
      transmission_rate= as.double(0.1),                # ◀ enforce double
      recovery_rate    = as.double(0.1)                 # ◀ enforce double
    )
    
    # Create a wrapper function that includes true_params
    local_simulation_fun <- function(params, lfmcmc_obj) {
      simulation_fun(params, lfmcmc_obj, true_params)
    }
    
    # Set up LFMCMC object
    lfmcmc_obj <- LFMCMC(dummy_model)
    lfmcmc_obj <- set_simulation_fun(lfmcmc_obj, local_simulation_fun)
    lfmcmc_obj <- set_summary_fun(lfmcmc_obj, summary_fun)
    lfmcmc_obj <- set_proposal_fun(lfmcmc_obj, proposal_fun)
    lfmcmc_obj <- set_kernel_fun(lfmcmc_obj, kernel_fun)
    lfmcmc_obj <- set_observed_data(lfmcmc_obj, observed_infected)
    
    # Initial parameters for ABC (as double)
    init_params <- c(
      true_params[3],  # contact_rate
      true_params[5],  # transmission_prob
      true_params[4]   # recovery_rate
    )
    
    epsilon <- sqrt(sum(observed_infected^2)) * 0.05
    
    run_lfmcmc(
      lfmcmc     = lfmcmc_obj,
      params_init= init_params,
      n_samples  = 500,
      epsilon    = epsilon,
      seed       = model_seed + sim_id + 100
    )
    
    accepted <- get_all_accepted_params(lfmcmc_obj)
    
    if (!is.null(accepted) && nrow(accepted) > 0) {
      calibrated_params_raw <- apply(accepted, 2, median)
      abc_crate    <- calibrated_params_raw[1]
      abc_ptran    <- calibrated_params_raw[2]
      abc_recov    <- calibrated_params_raw[3]
      abc_R0       <- (abc_crate * abc_ptran) / abc_recov
      calibrated_params <- c(abc_crate, abc_ptran, abc_recov)
    } else {
      abc_crate    <- true_params[3]
      abc_recov    <- true_params[4]
      abc_ptran    <- true_params[5]
      abc_R0       <- true_R0
      calibrated_params <- c(abc_crate, abc_ptran, abc_recov)
      cat("Warning: calibration failed on", sim_id, "\n")
    }
    
    # Use the calibrated parameters to predict
    abc_predicted_infected <- simulate_epidemic_calib(calibrated_params, true_params, seed = model_seed + sim_id + 200)
    # Ensure the predicted counts are double
    abc_predicted_infected <- as.numeric(abc_predicted_infected)
    
    days <- 0:model_ndays
    
    daily_results <- data.frame(
      sim_id                 = sim_id,
      day                    = days,
      true_n                 = true_params[1],
      true_preval            = true_params[2],
      true_crate             = true_params[3],
      true_recov             = true_params[4],
      true_ptran             = true_params[5],
      true_R0                = true_R0,
      abc_crate              = abc_crate,
      abc_ptran              = abc_ptran,
      abc_recov              = abc_recov,
      abc_R0                 = abc_R0,
      observed_infected      = observed_infected,
      abc_predicted_infected = abc_predicted_infected,
      abc_bias               = abc_predicted_infected - observed_infected,
      abc_rel_bias           = ifelse(
        observed_infected > 0,
        (abc_predicted_infected - observed_infected) / observed_infected,
        NA
      )
    )
    
    param_comparison <- data.frame(
      sim_id              = sim_id,
      parameter           = c("contact_rate", "transmission_prob", "recovery_rate"),
      true_value          = c(true_params[3], true_params[5], true_params[4]),
      abc_calibrated_value= c(abc_crate, abc_ptran, abc_recov),
      abc_error_pct       = c(
        (abc_crate - true_params[3]) / true_params[3] * 100,
        (abc_ptran - true_params[5]) / true_params[5] * 100,
        (abc_recov - true_params[4]) / true_params[4] * 100
      )
    )
    
    abc_parameters <- data.frame(
      sim_id      = sim_id,
      true_n      = true_params[1],
      true_preval = true_params[2],
      true_crate  = true_params[3],
      true_recov  = true_params[4],
      true_ptran  = true_params[5],
      true_R0     = true_R0,
      abc_crate   = abc_crate,
      abc_recov   = abc_recov,
      abc_ptran   = abc_ptran,
      abc_R0      = abc_R0
    )
    
    return(list(
      daily_results    = daily_results,
      param_comparison = param_comparison,
      abc_parameters   = abc_parameters
    ))
    
  }, error = function(e) {
    cat("Error in simulation", sim_id, ":", e$message, "\n")
    return(NULL)
  })
}

# Run Simulations with Slurm
ans <- Slurm_lapply(
  X       = 1:N_SIMS,
  FUN     = function(i) simulate_and_calibrate(as.numeric(theta_use[i, ]), i),
  job_name= "Sims_calibrate",
  njobs   = 100,
  overwrite = TRUE,
  plan    = "submit",
  sbatch_opt = list(
    partition     = "vegayon-shared-np",
    account       = "vegayon-np",
    time          = "01:00:00",
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
    "N_SIMS"
  )
)

# Collect results
results_list <- Slurm_collect(ans)

# Process results with error handling
valid_results <- results_list[!sapply(results_list, is.null)]
cat("Successfully completed", length(valid_results), "out of", N_SIMS, "simulations\n")


results_list <- Slurm_collect(ans)

# FIXED: Add result processing with error handling
valid_results <- results_list[!sapply(results_list, is.null)]
cat("Successfully completed", length(valid_results), "out of", N_SIMS, "simulations\n")

# Combine results if there are valid ones
if (length(valid_results) > 0) {
  # Combine daily results
  all_daily_results <- do.call(rbind, lapply(valid_results, function(x) x$daily_results))
  
  # Combine parameter comparisons
  all_param_comparisons <- do.call(rbind, lapply(valid_results, function(x) x$param_comparison))
  
  # Combine ABC parameters
  all_abc_parameters <- do.call(rbind, lapply(valid_results, function(x) x$abc_parameters))
  
  # Save results
  write.csv(all_daily_results, "daily_results.csv", row.names = FALSE)
  write.csv(all_param_comparisons, "param_comparisons.csv", row.names = FALSE)
  write.csv(all_abc_parameters, "abc_parameters.csv", row.names = FALSE)
  
  cat("Results saved to CSV files\n")
} else {
  cat("No valid results to process\n")
}