#Load required libraries
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
N_SIMS      <- 100   # number of simulations to run
N_CORES     <- 10    # number of cores to use

# --------------------------
# Generate Parameter Sets using Theta
# --------------------------
set.seed(model_seed)  # Ensure reproducibility
n_values <- rep(5000, N_SIMS)  # population size (constant at 5000)
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
# Function: Simulation–Calibration Study with ABC only
# --------------------------
simulate_and_calibrate <- function(true_params, sim_id) {
  # Load required libraries for parallel execution
  library(epiworldR)
  library(data.table)
  
  # Ensure true_params is a numeric vector
  true_params <- as.numeric(true_params)
  
  cat("Running simulation", sim_id, "of", N_SIMS, "\n")
  
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
  
  # Define ABC functions within the worker function to avoid scoping issues
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
# Run Simulations in Parallel
# --------------------------

# Setup cluster for parallel processing
cat("Setting up parallel cluster with", N_CORES, "cores...\n")
cl <- makeCluster(N_CORES)

# Export necessary objects and functions to the cluster
clusterExport(cl, c(
  "theta_use", "model_ndays", "model_seed", "global_n", "N_SIMS",
  "simulate_epidemic_observed", "simulate_epidemic_calib"
))

# Load required libraries on each worker
clusterEvalQ(cl, {
  library(epiworldR)
  library(data.table)
})

# Convert theta_use to a list of numeric vectors for parallel processing
param_list <- lapply(1:nrow(theta_use), function(i) as.numeric(theta_use[i, ]))

# Run simulations in parallel
cat("Running", N_SIMS, "simulations in parallel...\n")
start_time <- Sys.time()

results_list <- parLapply(cl, 1:N_SIMS, function(i) {
  simulate_and_calibrate(param_list[[i]], i)
})

end_time <- Sys.time()
cat("Parallel execution completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# Stop the cluster
stopCluster(cl)

# --------------------------
# Extract and Save Parameter Results
# --------------------------
# Extract and combine all daily results
cat("Processing results...\n")
daily_results <- bind_rows(lapply(results_list, function(x) x$daily_results))

# Extract and combine all parameter comparisons
param_comparisons <- bind_rows(lapply(results_list, function(x) x$param_comparison))

# Extract and combine all ABC parameter values
abc_parameters <- bind_rows(lapply(results_list, function(x) x$abc_parameters))

# Save the parameter results (this is what you specifically requested)
saveRDS(abc_parameters, "abc_predicted_parameters_100_fixedrecn.rds")
write.csv(abc_parameters, "abc_predicted_parameters_100_fixedrecn.csv", row.names = FALSE)

# Also save the other results for completeness
write.csv(daily_results, "abc_daily_results.csv", row.names = FALSE)
write.csv(param_comparisons, "abc_param_comparisons.csv", row.names = FALSE)

# --------------------------
# Calculate Coverage and Confidence Intervals
# --------------------------
# Function to calculate confidence interval
calc_CI <- function(values, conf_level = 0.95) {
  # Calculate the mean
  mean_val <- mean(values, na.rm = TRUE)
  
  # Calculate standard error
  n <- length(values[!is.na(values)])
  se <- sd(values, na.rm = TRUE) / sqrt(n)
  
  # Calculate quantile-based confidence intervals
  lower_ci <- quantile(values, 0.025, na.rm = TRUE)
  upper_ci <- quantile(values, 0.975, na.rm = TRUE)
  
  return(c(mean = mean_val, lower_ci = lower_ci, upper_ci = upper_ci, se = se))
}

# For each day, calculate stats across simulations
coverage_stats <- data.frame(day = 0:model_ndays)

for (day_num in 0:model_ndays) {
  day_data <- daily_results %>% filter(day == day_num)
  
  # Extract observed and predicted values for this day
  observed_values <- day_data$observed_infected
  abc_predicted_values <- day_data$abc_predicted_infected
  
  # Calculate confidence intervals for observed values
  obs_ci <- calc_CI(observed_values)
  coverage_stats$obs_mean[coverage_stats$day == day_num] <- obs_ci["mean"]
  coverage_stats$obs_lower_ci[coverage_stats$day == day_num] <- obs_ci["lower_ci"]
  coverage_stats$obs_upper_ci[coverage_stats$day == day_num] <- obs_ci["upper_ci"]
  
  # Calculate confidence intervals for ABC predicted values
  abc_pred_ci <- calc_CI(abc_predicted_values)
  coverage_stats$abc_pred_mean[coverage_stats$day == day_num] <- abc_pred_ci["mean"]
  coverage_stats$abc_pred_lower_ci[coverage_stats$day == day_num] <- abc_pred_ci["lower_ci"]
  coverage_stats$abc_pred_upper_ci[coverage_stats$day == day_num] <- abc_pred_ci["upper_ci"]
  
  # Calculate coverage for ABC
  abc_coverage <- mean(abc_predicted_values >= obs_ci["lower_ci"] & 
                         abc_predicted_values <= obs_ci["upper_ci"], na.rm = TRUE) * 100
  coverage_stats$abc_coverage[coverage_stats$day == day_num] <- abc_coverage
}

# Calculate daily statistics across simulations
daily_stats <- daily_results %>%
  group_by(day) %>%
  summarize(
    mean_observed = mean(observed_infected),
    abc_mean_predicted = mean(abc_predicted_infected),
    abc_mean_bias = mean(abc_bias),
    abc_mean_abs_bias = mean(abs(abc_bias)),
    n_sims = n()
  )

# Join with coverage statistics
daily_stats <- left_join(daily_stats, coverage_stats, by = "day")

# Save the aggregated statistics
write.csv(daily_stats, "abc_daily_stats.csv", row.names = FALSE)

# --------------------------
# Basic Visualization - ABC Only
# --------------------------

# 1. Bias plot
bias_plot <- ggplot(daily_stats, aes(x = day)) +
  geom_line(aes(y = abc_mean_bias), size = 1, color = "blue") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(title = "Mean Bias in ABC Calibration",
       subtitle = paste(N_SIMS, "simulations"),
       x = "Day", y = "Mean Bias (Predicted - Observed)") +
  theme_minimal()

# 2. Epidemic curves 
epidemic_curve_plot <- ggplot(daily_stats, aes(x = day)) +
  geom_line(aes(y = obs_mean, color = "Observed"), size = 1) +
  geom_line(aes(y = abc_pred_mean, color = "ABC"), size = 1, linetype = "dashed") +
  geom_ribbon(aes(ymin = abc_pred_lower_ci, ymax = abc_pred_upper_ci), fill = "blue", alpha = 0.1) +
  scale_color_manual(values = c("Observed" = "black", "ABC" = "blue")) +
  labs(title = "Epidemic Curves: Observed vs ABC",
       subtitle = paste("Based on", N_SIMS, "simulations"),
       x = "Day", y = "Infected Count") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 3. Parameter recovery comparison plots
param_plots <- list()
params <- c("contact_rate", "recovery_rate", "transmission_prob")

for (param in params) {
  param_data <- param_comparisons[param_comparisons$parameter == param, ]
  
  p <- ggplot(param_data, aes(x = true_value)) +
    geom_point(aes(y = abc_calibrated_value), alpha = 0.6, color = "blue") +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    labs(title = str_to_title(gsub("_", " ", param)),
         x = "True Value", 
         y = "ABC Estimated Value") +
    theme_minimal()
  
  param_plots[[param]] <- p
}

# Calculate aggregate statistics for parameter estimation
param_stats <- param_comparisons %>%
  group_by(parameter) %>%
  summarize(
    mean_true = mean(true_value, na.rm = TRUE),
    mean_abc = mean(abc_calibrated_value, na.rm = TRUE),
    mean_abc_error_pct = mean(abc_error_pct, na.rm = TRUE),
    rmse_abc = sqrt(mean((abc_calibrated_value - true_value)^2, na.rm = TRUE))
  )

# Save parameter statistics
write.csv(param_stats, "abc_parameter_statistics.csv", row.names = FALSE)

# Save plots
ggsave("abc_bias.png", plot = bias_plot, width = 10, height = 6)
ggsave("abc_epidemic_curve.png", plot = epidemic_curve_plot, width = 10, height = 6)

# Save individual parameter plots
for (param in params) {
  ggsave(paste0("abc_", param, "_recovery.png"), plot = param_plots[[param]], width = 8, height = 6)
}

cat("\nSimulation and calibration complete with ABC only. Predicted parameters saved to 'abc_predicted_parameters_100_fixedrecn.csv' and 'abc_predicted_parameters_100_fixedrecn.rds'!\n")