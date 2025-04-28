# Load required libraries
library(epiworldR)
library(data.table)
library(parallel)
library(ggplot2)
library(dplyr)
library(slurmR)  # Add slurmR library

# --------------------------
# Global Simulation Settings
# --------------------------
model_ndays <- 30   # simulation duration (days)
model_seed  <- 122    # seed for reproducibility
global_n    <- 10000  # population size (used in calibration)
N_SIMS      <- 100  # number of simulations to run

# --------------------------
# Generate Parameter Sets using Theta
# --------------------------
set.seed(model_seed)  # Ensure reproducibility
n_values <- rep(10000, N_SIMS)  # population size (constant at 10000)
theta <- data.table(
  n      = n_values,
  preval = sample(100:2000, N_SIMS, replace = TRUE) / n_values,
  crate  = runif(N_SIMS, 5, 15),
  recov  = 1 / runif(N_SIMS, 4, 14),
  R0     = runif(N_SIMS, 1, 5)
)
# Calculate transmission rate (ptran)
theta[, ptran := plogis(qlogis(R0 * recov / crate) + rnorm(.N))]

# Use only the needed columns
theta_use <- theta[, .(n, preval, crate, recov, ptran)]

# Print the true parameter sets for reference
cat("True parameter sets:\n")
print(theta_use)

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

# --------------------------
# Function: Simplified Simulation–Calibration Study
# --------------------------
simulate_and_calibrate <- function(true_params, sim_id) {
  cat("Running simulation", sim_id, "of", N_SIMS, "\n")
  
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
  
  return(list(
    daily_results = daily_results,
    param_comparison = param_comparison
  ))
}

# --------------------------
# Run Simulations using slurmR instead of mclapply
# --------------------------
# Load required libraries
library(slurmR)

# Configure Slurm options correctly
opts <- list(
  account     = "my-account",     # Replace with your account
  partition   = "my-partition",   # Replace with your partition
  "cpus-per-task" = 1,
  "mem-per-cpu"   = "4G",
  time        = "01:00:00"
)

# Create a Slurm job
sjob <- Slurm_lapply(
  X        = 1:N_SIMS,
  FUN      = function(i) simulate_and_calibrate(as.numeric(theta_use[i]), i),
  njobs    = 20,           # Split into 20 jobs (5 simulations per job)
  mc.cores = 1,            # Cores per job
  tmp_path = getwd(),      # Use current directory for temp files
  sbatch_opt = opts,       # Pass Slurm options here
  job_name = "epidemic_sim"
)

# Submit the job
sjob <- sbatch(sjob)

# Later, collect results
results_list <- Slurm_collect(sjob)

# Option 2: Using makeSlurmCluster (alternative approach)
#Comment out Option 1 above and uncomment this section to use this approach
cl <- makeSlurmCluster(
  n_cores,
  sbatch_opt = list(
    time = "12:00:00",
    mem_per_cpu = "4G",
    job_name = "epi_sim"
  )
)

# Export any needed variables to the cluster
clusterExport(cl, c("model_ndays", "model_seed", "global_n", "theta_use", "N_SIMS"))

# Make sure all needed packages are loaded on the cluster nodes
clusterEvalQ(cl, {
  library(epiworldR)
  library(data.table)
  library(ggplot2)
  library(dplyr)
})

# Run the simulations
results_list <- parLapply(cl, 1:N_SIMS, function(i) {
  simulate_and_calibrate(as.numeric(theta_use[i]), i)
})

# Stop the cluster when done
stopCluster(cl)

# --------------------------
# Combine and Analyze Results
# --------------------------
# Extract and combine all daily results
daily_results <- bind_rows(lapply(results_list, function(x) x$daily_results))

# Extract and combine all parameter comparisons
param_comparisons <- bind_rows(lapply(results_list, function(x) x$param_comparison))

# --------------------------
# Calculate Coverage and Confidence Intervals (using all simulations)
# --------------------------
# For each day, calculate statistics across simulations
coverage_stats <- data.frame(day = 0:model_ndays)

# Function to calculate confidence interval
# Using t-distribution for small sample sizes
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
# For each day, calculate stats across simulations
for (day_num in 0:model_ndays) {
  day_data <- daily_results %>% filter(day == day_num)
  
  # Extract observed and predicted values for this day
  observed_values <- day_data$observed_infected
  predicted_values <- day_data$predicted_infected
  
  # Calculate confidence intervals for observed values
  obs_ci <- calc_CI(observed_values)
  coverage_stats$obs_mean[coverage_stats$day == day_num] <- obs_ci["mean"]
  coverage_stats$obs_lower_ci[coverage_stats$day == day_num] <- obs_ci["lower_ci"]
  coverage_stats$obs_upper_ci[coverage_stats$day == day_num] <- obs_ci["upper_ci"]
  coverage_stats$obs_se[coverage_stats$day == day_num] <- obs_ci["se"]
  
  # Calculate confidence intervals for predicted values
  pred_ci <- calc_CI(predicted_values)
  coverage_stats$pred_mean[coverage_stats$day == day_num] <- pred_ci["mean"]
  coverage_stats$pred_lower_ci[coverage_stats$day == day_num] <- pred_ci["lower_ci"]
  coverage_stats$pred_upper_ci[coverage_stats$day == day_num] <- pred_ci["upper_ci"]
  coverage_stats$pred_se[coverage_stats$day == day_num] <- pred_ci["se"]
  
  # Calculate prediction intervals (for coverage)
  # 95% prediction intervals based on the distribution of observed values
  coverage_stats$obs_lower_pi[coverage_stats$day == day_num] <- quantile(observed_values, 0.025)
  coverage_stats$obs_upper_pi[coverage_stats$day == day_num] <- quantile(observed_values, 0.975)
  
  # Count how many predicted values fall within the PI of observed values
  coverage <- mean(predicted_values >= coverage_stats$obs_lower_pi[coverage_stats$day == day_num] & 
                     predicted_values <= coverage_stats$obs_upper_pi[coverage_stats$day == day_num]) * 100
  coverage_stats$coverage[coverage_stats$day == day_num] <- coverage
  
  # Also calculate the reverse - how many observed values fall within the PI of predicted values
  coverage_stats$pred_lower_pi[coverage_stats$day == day_num] <- quantile(predicted_values, 0.025)
  coverage_stats$pred_upper_pi[coverage_stats$day == day_num] <- quantile(predicted_values, 0.975)
  
  reverse_coverage <- mean(observed_values >= coverage_stats$pred_lower_pi[coverage_stats$day == day_num] & 
                             observed_values <= coverage_stats$pred_upper_pi[coverage_stats$day == day_num]) * 100
  coverage_stats$reverse_coverage[coverage_stats$day == day_num] <- reverse_coverage
}

# Calculate daily statistics across simulations
daily_stats <- daily_results %>%
  group_by(day) %>%
  summarize(
    mean_observed = mean(observed_infected),
    mean_predicted = mean(predicted_infected),
    mean_bias = mean(bias),
    mean_abs_bias = mean(abs(bias)),
    median_bias = median(bias),
    mean_rel_bias = mean(rel_bias, na.rm = TRUE) * 100, # As percentage
    mean_preval_error = mean((calib_preval - true_preval)/true_preval) * 100,
    mean_crate_error = mean((calib_crate - true_crate)/true_crate) * 100,
    mean_recov_error = mean((calib_recov - true_recov)/true_recov) * 100,
    mean_ptran_error = mean((calib_ptran - true_ptran)/true_ptran) * 100,
    sd_observed = sd(observed_infected),
    sd_predicted = sd(predicted_infected),
    n_sims = n()
  )

# Join with coverage statistics
daily_stats <- left_join(daily_stats, coverage_stats, by = "day")

# --------------------------
# Display Results
# --------------------------
# Print mean bias by day
cat("\n=== Mean Bias by Day ===\n")
bias_by_day <- daily_stats %>% 
  select(day, mean_bias, mean_abs_bias, median_bias, mean_rel_bias)
print(bias_by_day)

# Print confidence intervals by day (selected days for readability)
cat("\n=== Confidence Intervals by Day ===\n")
ci_by_day <- daily_stats %>% 
  filter(day %in% c(0, 5, 10, 15, 20, 25, 30)) %>%
  select(day, obs_mean, obs_lower_ci, obs_upper_ci, pred_mean, pred_lower_ci, pred_upper_ci)
print(ci_by_day)

# Print coverage by day
cat("\n=== Coverage Statistics by Day ===\n")
coverage_by_day <- daily_stats %>% 
  filter(day %in% c(0, 5, 10, 15, 20, 25, 30)) %>%
  select(day, coverage, reverse_coverage)
print(coverage_by_day)

# Print parameter comparison table (wide format for easier reading)
param_wide <- param_comparisons %>%
  tidyr::pivot_wider(
    id_cols = sim_id,
    names_from = parameter,
    values_from = c(true_value, calibrated_value, error_pct)
  )

cat("\n=== True vs Calibrated Parameters ===\n")
print(param_wide)

# --------------------------
# Create Summary Tables
# --------------------------
# Parameter recovery by parameter type
param_summary <- param_comparisons %>%
  group_by(parameter) %>%
  summarize(
    median_error_pct = median(error_pct),
    mean_abs_error_pct = mean(abs(error_pct))
  )

cat("\n=== Parameter Recovery Summary ===\n")
print(param_summary)

# Coverage summary
coverage_summary <- data.frame(
  overall_coverage = mean(daily_stats$coverage, na.rm = TRUE),
  early_coverage = mean(daily_stats$coverage[daily_stats$day <= 10], na.rm = TRUE),
  middle_coverage = mean(daily_stats$coverage[daily_stats$day > 10 & daily_stats$day <= 20], na.rm = TRUE),
  late_coverage = mean(daily_stats$coverage[daily_stats$day > 20], na.rm = TRUE)
)

cat("\n=== Coverage Summary ===\n")
print(coverage_summary)

# --------------------------
# Generate Key Visualizations
# --------------------------

# 1. Mean bias over time (MAIN PLOT REQUESTED)
mean_bias_plot <- ggplot(daily_stats, aes(x = day)) +
  geom_line(aes(y = mean_bias), color = "blue", size = 1) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Mean Bias Across All Simulations By Day",
       subtitle = paste(N_SIMS, "simulations"),
       x = "Day", y = "Mean Bias (Predicted - Observed)") +
  theme_minimal()

print(mean_bias_plot)

epi_curve_pred_ci_plot <- ggplot(daily_stats, aes(x = day)) +
  # Only predicted confidence intervals (no observed CI)
  geom_ribbon(aes(ymin = pred_lower_ci, ymax = pred_upper_ci), fill = "red", alpha = 0.2) +
  # Both observed and predicted lines
  geom_line(aes(y = obs_mean), color = "blue", size = 1) +
  geom_line(aes(y = pred_mean), color = "red", size = 1, linetype = "dashed") +
  scale_color_manual(values = c("Observed" = "blue", "Predicted" = "red")) +
  labs(title = "Epidemic Curves with 95% Predicted Confidence Intervals",
       subtitle = paste("Based on", N_SIMS, "simulations"),
       x = "Day", y = "Infected Count") +
  theme_minimal() +
  # Add a legend to distinguish the lines
  annotate("text", x = max(daily_stats$day)*0.8, y = max(daily_stats$pred_mean)*0.9, 
           label = "Blue: Observed\nRed: Predicted", hjust = 0)

print(epi_curve_pred_ci_plot)

# 2. Modified coverage plot with just ONE line (predicted in observed PI)
single_coverage_plot <- ggplot(daily_stats, aes(x = day)) +
  geom_line(aes(y = coverage), color = "blue", size = 1) +
  geom_hline(yintercept = 95, color = "red", linetype = "dotted") +
  labs(title = "Coverage Rate Across Simulations",
       subtitle = "Percentage of predicted values that fall within the observed prediction interval",
       x = "Day", y = "Coverage Rate (%)") +
  ylim(0, 100) +
  theme_minimal()

print(single_coverage_plot)


simple_prevalence_plot <- ggplot(param_comparisons[param_comparisons$parameter == "prevalence", ], 
                                 aes(x = true_value, y = calibrated_value)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Initial Prevalence",
       x = "True Value", 
       y = "Calibrated Value") +
  theme_minimal()

# Contact rate plot - simplified
simple_contact_rate_plot <- ggplot(param_comparisons[param_comparisons$parameter == "contact_rate", ], 
                                   aes(x = true_value, y = calibrated_value)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Contact Rate",
       x = "True Value", 
       y = "Calibrated Value") +
  theme_minimal()

# Recovery rate plot - simplified
simple_recovery_rate_plot <- ggplot(param_comparisons[param_comparisons$parameter == "recovery_rate", ], 
                                    aes(x = true_value, y = calibrated_value)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Recovery Rate",
       x = "True Value", 
       y = "Calibrated Value") +
  theme_minimal()

# Transmission probability plot - simplified
simple_transmission_plot <- ggplot(param_comparisons[param_comparisons$parameter == "transmission_prob", ], 
                                   aes(x = true_value, y = calibrated_value)) +
  geom_point(color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(title = "Transmission Probability",
       x = "True Value", 
       y = "Calibrated Value") +
  theme_minimal()

# Print each plot
print(simple_prevalence_plot)
print(simple_contact_rate_plot)
print(simple_recovery_rate_plot)
print(simple_transmission_plot)

# Create a grid of all plots if gridExtra is available
if (requireNamespace("gridExtra", quietly = TRUE)) {
  simple_grid <- gridExtra::grid.arrange(
    simple_prevalence_plot, simple_contact_rate_plot,
    simple_recovery_rate_plot, simple_transmission_plot,
    ncol = 2
  )
  print(simple_grid)
} else {
  cat("Package 'gridExtra' not available for grid layout. Individual plots shown instead.\n")
}
# --------------------------
# Save Results to Files
# --------------------------
# Daily statistics across all simulations
write.csv(daily_stats, "epidemic_daily_stats_across_sims.csv", row.names = FALSE)

# Daily results for all days and all simulations
write.csv(daily_results, "epidemic_all_days_all_sims.csv", row.names = FALSE)

# Parameter comparisons
write.csv(param_comparisons, "epidemic_param_comparisons.csv", row.names = FALSE)

# Parameter summary
write.csv(param_summary, "epidemic_param_summary.csv", row.names = FALSE)

cat("\nSimulation and analysis complete!\n")