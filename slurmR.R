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
N_SIMS      <- 10  # number of simulations to run

# --------------------------
# Generate Parameter Sets using Theta
# --------------------------
set.seed(model_seed)  # Ensure reproducibility
n_values <- rep(5000, N_SIMS)  # population size (constant at 10000)
theta <- data.table(
  n      = n_values,
  preval = sample(100:2000, N_SIMS, replace = TRUE) / n_values,
  crate  = runif(N_SIMS, 5, 15),
  recov  = 0.1, #1 / runif(N_SIMS, 4, 14),
  R0     = runif(N_SIMS, 1, 5)
)
# Calculate transmission rate (ptran)
theta[, ptran := plogis(qlogis(R0 * recov / crate) + rnorm(.N))]

# Use only the needed columns
theta_use <- theta[, .(n, preval, crate, recov, ptran)]

# Print the true parameter sets for reference
cat("True parameter sets:\n")
print(theta_use)
saveRDS(theta_use, "theta_use.rds")

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
  # Extract parameters for simulation
  # We now have: contact_rate, recovery_rate, transmission_prob (3 parameters)
  # prevalence and R0 are not estimated/predicted
  
  sim_params <- c(
    true_params[2],  # prevalence (fixed at true value)
    params[1],       # contact_rate
    params[2],       # recovery_rate
    params[3]        # transmission_prob
  )
  
  # Run simulation with the parameters
  infected_counts <- simulate_epidemic_calib(sim_params, ndays = model_ndays)
  return(as.numeric(infected_counts))
}

# This returns the observed data (all daily infected counts)
summary_fun <- function(data, lfmcmc_obj) {
  return(as.numeric(data))
}

# Generate new parameter proposals
proposal_fun <- function(old_params, lfmcmc_obj) {
  # Proposals with appropriate step sizes for each parameter
  # old_params contains: contact_rate, recovery_rate, transmission_prob
  
  new_crate <- old_params[1] * exp(rnorm(1, sd = 0.1))  # Log-normal proposal
  new_recov <- old_params[2] * exp(rnorm(1, sd = 0.1))  # Log-normal proposal
  new_ptran <- plogis(qlogis(old_params[3]) + rnorm(1, sd = 0.1))
  
  return(c(new_crate, new_recov, new_ptran))
}

# Kernel function using sum of squared differences across all days
kernel_fun <- function(simulated_stat, observed_stat, epsilon, lfmcmc_obj) {
  # Calculate the sum of squared differences across all days
  diff <- sum((simulated_stat - observed_stat)^2)
  return(exp(-diff / (2 * epsilon^2)))
}

# --------------------------
# Function: Simulation–Calibration Study with ABC only (for Slurm)
# --------------------------
simulate_and_calibrate_slurm <- function(sim_id) {
  # Load the theta_use data
  theta_use <- readRDS("theta_use.rds")
  true_params <- as.numeric(theta_use[sim_id])
  
  cat("Running simulation", sim_id, "of", nrow(theta_use), "\n")
  
  # Calculate true R0 from true parameters
  # R0 = (contact_rate * transmission_rate) / recovery_rate
  true_R0 <- (true_params[3] * true_params[5]) / true_params[4]
  
  # Step 1: Simulate the observed epidemic using the true parameters
  observed_infected <- simulate_epidemic_observed(true_params, ndays = model_ndays, 
                                                  seed = model_seed + sim_id)
  
  # Step 2: ABC calibration (without prevalence and R0)
  # Create a dummy model for LFMCMC
  dummy_model <- ModelSIRCONN(
    name              = "dummy",
    n                 = as.integer(global_n),
    prevalence        = 0.1,
    contact_rate      = 5.0,
    transmission_rate = 0.1,
    recovery_rate     = 0.1
  )
  
  # Make true_params available to simulation_fun
  true_params <<- true_params
  
  # Setup and run LFMCMC with the full time series
  lfmcmc_obj <- LFMCMC(dummy_model)
  lfmcmc_obj <- set_simulation_fun(lfmcmc_obj, simulation_fun)
  lfmcmc_obj <- set_summary_fun(lfmcmc_obj, summary_fun)
  lfmcmc_obj <- set_proposal_fun(lfmcmc_obj, proposal_fun)
  lfmcmc_obj <- set_kernel_fun(lfmcmc_obj, kernel_fun)
  lfmcmc_obj <- set_observed_data(lfmcmc_obj, observed_infected)
  
  # Initial parameters for ABC: 
  # contact_rate, recovery rate, transmission probability (no prevalence, no R0)
  init_params <- c(
    true_params[3],              # contact rate
    true_params[4],              # recovery rate
    true_params[5]               # transmission probability
  )
  
  # Run the LFMCMC
  n_samples_calib <- 500  # Increased for better convergence
  
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
    calibrated_params_raw <- apply(accepted, 2, median)
    
    # Extract the calibrated parameters
    abc_crate <- calibrated_params_raw[1]
    abc_recov <- calibrated_params_raw[2]
    abc_ptran <- calibrated_params_raw[3]
    
    # Keep prevalence and R0 fixed at true values (not predicted)
    abc_preval <- true_params[2]
    abc_R0 <- true_R0
    
    # Create final calibrated parameter vector for simulation
    calibrated_params <- c(abc_preval, abc_crate, abc_recov, abc_ptran)
  } else {
    calibrated_params <- true_params[2:5]
    abc_R0 <- true_R0
    abc_crate <- true_params[3]
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
    # True parameters (same for all days in this simulation)
    true_preval = true_params[2],
    true_crate = true_params[3],
    true_recov = true_params[4],
    true_ptran = true_params[5],
    true_R0 = true_R0,
    # ABC calibrated parameters
    calib_preval = calibrated_params[1],
    calib_crate = calibrated_params[2],
    calib_recov = calibrated_params[3],
    calib_ptran = calibrated_params[4],
    calib_R0 = abc_R0,
    # Observed and predicted counts for each day
    observed_infected = observed_infected,
    abc_predicted_infected = abc_predicted_infected,
    # Bias calculation
    abc_bias = abc_predicted_infected - observed_infected,
    # Relative bias (to account for scale differences)
    abc_rel_bias = ifelse(observed_infected > 0, 
                          (abc_predicted_infected - observed_infected) / observed_infected, 
                          NA)
  )
  
  # Parameter comparison - Only the 3 calibrated parameters
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
  
  # Create a separate data frame for just the parameter values for this simulation
  abc_parameters <- data.frame(
    sim_id = sim_id,
    true_preval = true_params[2],
    true_crate = true_params[3],
    true_recov = true_params[4],
    true_ptran = true_params[5],
    true_R0 = true_R0,
    abc_preval = calibrated_params[1],
    abc_crate = calibrated_params[2],
    abc_recov = calibrated_params[3],
    abc_ptran = calibrated_params[4],
    abc_R0 = abc_R0
  )
  
  return(list(
    daily_results = daily_results,
    param_comparison = param_comparison,
    abc_parameters = abc_parameters
  ))
}

# --------------------------
# Manual Slurm Script Approach (workaround for slurmR bug)
# --------------------------

# Create chunks of simulation IDs
n_jobs <- 100  # Number of jobs to submit
sims_per_job <- ceiling(N_SIMS / n_jobs)
sim_chunks <- split(1:N_SIMS, ceiling(seq_along(1:N_SIMS) / sims_per_job))

# Create directories for Slurm scripts and output
dir.create("slurm_scripts", showWarnings = FALSE)
dir.create("slurm_output", showWarnings = FALSE)

# Function to create R script for each job
create_r_script <- function(chunk_id, sim_ids) {
  script_content <- paste0("
# Load required libraries
library(epiworldR)
library(data.table)
library(parallel)

# Source the main simulation function
source('simulation_functions.R')

# Get the simulation IDs for this job
sim_ids <- c(", paste(sim_ids, collapse = ", "), ")

# Run simulations for this chunk
results_list <- list()
for (i in seq_along(sim_ids)) {
  results_list[[i]] <- simulate_and_calibrate_slurm(sim_ids[i])
}

# Save results
saveRDS(results_list, paste0('chunk_', ", chunk_id, ", '_results.rds'))
")
  
  writeLines(script_content, paste0("slurm_scripts/job_", chunk_id, ".R"))
}

# Create bash script for each job
create_bash_script <- function(chunk_id) {
  bash_content <- paste0("#!/bin/bash
#SBATCH --job-name=epidemic_sim_", chunk_id, "
#SBATCH --partition=notchpeak-freecycle
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=slurm_output/job_", chunk_id, ".out
#SBATCH --error=slurm_output/job_", chunk_id, ".err

# Load R module if needed (uncomment and modify as needed)
# module load R

# Run the R script
cd $SLURM_SUBMIT_DIR
Rscript slurm_scripts/job_", chunk_id, ".R
")
  
  writeLines(bash_content, paste0("slurm_scripts/job_", chunk_id, ".sh"))
}

# Create a separate file with all the functions
functions_content <- "
# --------------------------
# All simulation functions from above
# --------------------------

# Global variables
model_ndays <- 60
model_seed <- 122
global_n <- 5000

# Function to simulate the \"observed\" epidemic
simulate_epidemic_observed <- function(params, ndays = model_ndays, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize the model with the true parameters
  sim_model <- ModelSIRCONN(
    name              = \"sim\",
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
  infected_counts <- counts[counts\$state == \"Infected\", \"counts\"]
  
  return(infected_counts)
}

# Function for simulation during calibration
simulate_epidemic_calib <- function(params, ndays = model_ndays, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize model with calibration parameters
  sim_model <- ModelSIRCONN(
    name              = \"sim\",
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
  infected_counts <- counts[counts\$state == \"Infected\", \"counts\"]
  
  return(infected_counts)
}

# ABC functions
simulation_fun <- function(params, lfmcmc_obj) {
  sim_params <- c(
    true_params[2],  # prevalence (fixed at true value)
    params[1],       # contact_rate
    params[2],       # recovery_rate
    params[3]        # transmission_prob
  )
  
  infected_counts <- simulate_epidemic_calib(sim_params, ndays = model_ndays)
  return(as.numeric(infected_counts))
}

summary_fun <- function(data, lfmcmc_obj) {
  return(as.numeric(data))
}

proposal_fun <- function(old_params, lfmcmc_obj) {
  new_crate <- old_params[1] * exp(rnorm(1, sd = 0.1))
  new_recov <- old_params[2] * exp(rnorm(1, sd = 0.1))
  new_ptran <- plogis(qlogis(old_params[3]) + rnorm(1, sd = 0.1))
  
  return(c(new_crate, new_recov, new_ptran))
}

kernel_fun <- function(simulated_stat, observed_stat, epsilon, lfmcmc_obj) {
  diff <- sum((simulated_stat - observed_stat)^2)
  return(exp(-diff / (2 * epsilon^2)))
}

# Main simulation function
simulate_and_calibrate_slurm <- function(sim_id) {
  theta_use <- readRDS(\"theta_use.rds\")
  true_params <- as.numeric(theta_use[sim_id])
  
  cat(\"Running simulation\", sim_id, \"of\", nrow(theta_use), \"\\n\")
  
  true_R0 <- (true_params[3] * true_params[5]) / true_params[4]
  
  observed_infected <- simulate_epidemic_observed(true_params, ndays = model_ndays, 
                                                  seed = model_seed + sim_id)
  
  dummy_model <- ModelSIRCONN(
    name              = \"dummy\",
    n                 = as.integer(global_n),
    prevalence        = 0.1,
    contact_rate      = 5.0,
    transmission_rate = 0.1,
    recovery_rate     = 0.1
  )
  
  true_params <<- true_params
  
  lfmcmc_obj <- LFMCMC(dummy_model)
  lfmcmc_obj <- set_simulation_fun(lfmcmc_obj, simulation_fun)
  lfmcmc_obj <- set_summary_fun(lfmcmc_obj, summary_fun)
  lfmcmc_obj <- set_proposal_fun(lfmcmc_obj, proposal_fun)
  lfmcmc_obj <- set_kernel_fun(lfmcmc_obj, kernel_fun)
  lfmcmc_obj <- set_observed_data(lfmcmc_obj, observed_infected)
  
  init_params <- c(
    true_params[3],
    true_params[4],
    true_params[5]
  )
  
  n_samples_calib <- 500
  epsilon <- sqrt(sum(observed_infected^2)) * 0.05
  
  run_lfmcmc(
    lfmcmc = lfmcmc_obj,
    params_init = init_params,
    n_samples = n_samples_calib,
    epsilon = epsilon,
    seed = model_seed + sim_id + 100
  )
  
  accepted <- get_all_accepted_params(lfmcmc_obj)
  
  if (!is.null(accepted) && nrow(accepted) > 0) {
    calibrated_params_raw <- apply(accepted, 2, median)
    abc_crate <- calibrated_params_raw[1]
    abc_recov <- calibrated_params_raw[2]
    abc_ptran <- calibrated_params_raw[3]
    abc_preval <- true_params[2]
    abc_R0 <- true_R0
    calibrated_params <- c(abc_preval, abc_crate, abc_recov, abc_ptran)
  } else {
    calibrated_params <- true_params[2:5]
    abc_R0 <- true_R0
    abc_crate <- true_params[3]
    cat(\"Warning: Using true parameters as calibration failed for simulation\", sim_id, \"\\n\")
  }
  
  abc_predicted_infected <- simulate_epidemic_calib(calibrated_params, ndays = model_ndays, 
                                                    seed = model_seed + sim_id + 200)
  
  days <- 0:model_ndays
  daily_results <- data.frame(
    sim_id = sim_id,
    day = days,
    true_preval = true_params[2],
    true_crate = true_params[3],
    true_recov = true_params[4],
    true_ptran = true_params[5],
    true_R0 = true_R0,
    calib_preval = calibrated_params[1],
    calib_crate = calibrated_params[2],
    calib_recov = calibrated_params[3],
    calib_ptran = calibrated_params[4],
    calib_R0 = abc_R0,
    observed_infected = observed_infected,
    abc_predicted_infected = abc_predicted_infected,
    abc_bias = abc_predicted_infected - observed_infected,
    abc_rel_bias = ifelse(observed_infected > 0, 
                          (abc_predicted_infected - observed_infected) / observed_infected, 
                          NA)
  )
  
  param_comparison <- data.frame(
    sim_id = sim_id,
    parameter = c(\"contact_rate\", \"recovery_rate\", \"transmission_prob\"),
    true_value = c(true_params[3:5]),
    abc_calibrated_value = c(abc_crate, abc_recov, abc_ptran),
    abc_error_pct = c(
      (abc_crate - true_params[3]) / true_params[3] * 100,
      (abc_recov - true_params[4]) / true_params[4] * 100,
      (abc_ptran - true_params[5]) / true_params[5] * 100
    )
  )
  
  abc_parameters <- data.frame(
    sim_id = sim_id,
    true_preval = true_params[2],
    true_crate = true_params[3],
    true_recov = true_params[4],
    true_ptran = true_params[5],
    true_R0 = true_R0,
    abc_preval = calibrated_params[1],
    abc_crate = calibrated_params[2],
    abc_recov = calibrated_params[3],
    abc_ptran = calibrated_params[4],
    abc_R0 = abc_R0
  )
  
  return(list(
    daily_results = daily_results,
    param_comparison = param_comparison,
    abc_parameters = abc_parameters
  ))
}
"

writeLines(functions_content, "simulation_functions.R")

# Create R and bash scripts for each chunk
for (i in seq_along(sim_chunks)) {
  create_r_script(i, sim_chunks[[i]])
  create_bash_script(i)
}

# Create a master script to submit all jobs
master_script <- paste0("#!/bin/bash

echo \"Submitting ", length(sim_chunks), " jobs to Slurm...\"

# Submit all jobs
job_ids=()
for i in {1..", length(sim_chunks), "}; do
  job_id=$(sbatch --parsable slurm_scripts/job_\${i}.sh)
  job_ids+=(\$job_id)
  echo \"Submitted job \$i with ID \$job_id\"
done

echo \"All jobs submitted. Job IDs: \${job_ids[@]}\"

# Optional: Create a script to check job status
cat > check_jobs.sh << 'EOF'
#!/bin/bash
squeue -u $USER -o \"%.8i %.8T %.4D %.8j %R\"
EOF

chmod +x check_jobs.sh

echo \"Use './check_jobs.sh' to check job status\"
echo \"Use 'squeue -u $USER' for detailed status\"
")

writeLines(master_script, "submit_all_jobs.sh")
system("chmod +x submit_all_jobs.sh")

# Create a collection script
collection_script <- "#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)

# Collect all results
results_list <- list()
result_files <- list.files(pattern = \"chunk_.*_results.rds\", full.names = TRUE)

cat(\"Found\", length(result_files), \"result files\\n\")

for (file in result_files) {
  chunk_results <- readRDS(file)
  results_list <- c(results_list, chunk_results)
}

cat(\"Total simulations collected:\", length(results_list), \"\\n\")

# Extract and combine results
daily_results <- bind_rows(lapply(results_list, function(x) x\$daily_results))
param_comparisons <- bind_rows(lapply(results_list, function(x) x\$param_comparison))
abc_parameters <- bind_rows(lapply(results_list, function(x) x\$abc_parameters))

# Save final results
saveRDS(abc_parameters, \"abc_predicted_parameters_1000_fixedrecn.rds\")
write.csv(abc_parameters, \"abc_predicted_parameters_1000_fixedrecn.csv\", row.names = FALSE)

cat(\"Results saved to:\\n\")
cat(\"  - abc_predicted_parameters_1000_fixedrecn.rds\\n\")
cat(\"  - abc_predicted_parameters_1000_fixedrecn.csv\\n\")
"

writeLines(collection_script, "collect_results.R")
system("chmod +x collect_results.R")

cat("Setup complete! To run the simulations:\n")
cat("1. Submit all jobs:     ./submit_all_jobs.sh\n")
cat("2. Check job status:    ./check_jobs.sh\n")
cat("3. Collect results:     Rscript collect_results.R\n")
cat("\nNumber of jobs created:", length(sim_chunks), "\n")
cat("Simulations per job:", sims_per_job, "\n")