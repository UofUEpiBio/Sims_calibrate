source("~/Desktop/Sims_calibrate/params_gen_001.R")
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
# Function: Simulation–Calibration Study with ABC only
# --------------------------
simulate_and_calibrate <- function(true_params, sim_id) {
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
# Run Simulations in Parallel
# --------------------------
# Determine how many cores to use

# Sequential processing
results_list <- list()
for (i in 1:N_SIMS) {
  results_list[[i]] <- simulate_and_calibrate(as.numeric(theta_use[i]), i)
}

# --------------------------
# Extract and Save Parameter Results
# --------------------------
# Extract and combine all daily results
daily_results <- bind_rows(lapply(results_list, function(x) x$daily_results))

# Extract and combine all parameter comparisons
param_comparisons <- bind_rows(lapply(results_list, function(x) x$param_comparison))

# Extract and combine all ABC parameter values
abc_parameters <- bind_rows(lapply(results_list, function(x) x$abc_parameters))

# Save the parameter results (this is what you specifically requested)
saveRDS(abc_parameters, "abc_predicted_parameters_1000_fixedrecn.rds")
write.csv(abc_parameters, "abc_predicted_parameters_1000_fixedrecn.csv", row.names = FALSE)




#############
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
params <- c("prevalence", "contact_rate", "recovery_rate", "transmission_prob", "R0")

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

cat("\nSimulation and calibration complete with ABC only. Predicted parameters saved to 'abc_predicted_parameters.csv' and 'abc_predicted_parameters.rds'!\n")



############################
# Load required libraries
library(epiworldR)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)

# --------------------------
# Set global parameters
# --------------------------
n_population <- 5000     # Population size
recovery_rate <- 0.1     # Fixed recovery rate
model_ndays <- 50        # Simulation duration (days)
n_simulations <- 50      # Number of simulations per parameter set

# --------------------------
# Main function to run simulations for both true and ABC parameters
# --------------------------
run_and_compare <- function() {
  # Initialize lists to store results
  true_results_list <- list()
  abc_results_list <- list()
  
  # Process each row in abc_parameters
  for (i in 1:nrow(abc_parameters)) {
    cat("Processing parameter set", i, "of", nrow(abc_parameters), "\n")
    
    # Extract parameter values
    true_preval <- abc_parameters$true_preval[i]
    true_crate <- abc_parameters$true_crate[i]
    true_ptran <- abc_parameters$true_ptran[i]
    
    abc_preval <- abc_parameters$abc_preval[i]
    abc_crate <- abc_parameters$abc_crate[i]
    abc_ptran <- abc_parameters$abc_ptran[i]
    
    # Run simulations with true parameters
    cat("  Running simulations with true parameters...\n")
    model_true <- ModelSIRCONN(
      name = paste0("True_", i),
      n = n_population,
      prevalence = true_preval,
      contact_rate = true_crate,
      transmission_rate = true_ptran,
      recovery_rate = recovery_rate
    )
    
    # Create saver for true parameters
    saver_true <- make_saver("total_hist", "reproductive")
    
    # Run the simulations with true parameters
    run_multiple(model_true, ndays = model_ndays, nsims = n_simulations, 
                 saver = saver_true, nthread = 2)
    
    # Retrieve the results for true parameters
    true_results <- run_multiple_get_results(model_true)
    
    # Run simulations with ABC parameters
    cat("  Running simulations with ABC parameters...\n")
    model_abc <- ModelSIRCONN(
      name = paste0("ABC_", i),
      n = n_population,
      prevalence = abc_preval,
      contact_rate = abc_crate,
      transmission_rate = abc_ptran,
      recovery_rate = recovery_rate
    )
    
    # Create saver for ABC parameters
    saver_abc <- make_saver("total_hist", "reproductive")
    
    # Run the simulations with ABC parameters
    run_multiple(model_abc, ndays = model_ndays, nsims = n_simulations, 
                 saver = saver_abc, nthread = 2)
    
    # Retrieve the results for ABC parameters
    abc_results <- run_multiple_get_results(model_abc)
    
    # Process results for this parameter set
    true_data <- process_simulation_results(true_results, "true", i)
    abc_data <- process_simulation_results(abc_results, "abc", i)
    
    # Store results
    true_results_list[[i]] <- true_data
    abc_results_list[[i]] <- abc_data
  }
  
  # Combine all results
  all_true_results <- do.call(rbind, true_results_list)
  all_abc_results <- do.call(rbind, abc_results_list)
  
  # Combine true and ABC results
  all_results <- rbind(all_true_results, all_abc_results)
  
  return(all_results)
}

# Function to process simulation results
process_simulation_results <- function(results, param_type, param_set_id) {
  # Create a list to store results for each simulation
  processed_data <- list()
  
  # Process each simulation
  for (sim_idx in 1:length(results$total_hist)) {
    # Get data for this simulation
    sim_data <- results$total_hist[[sim_idx]]
    
    # Filter for Infected state
    infected_data <- sim_data[sim_data$state == "Infected", ]
    
    # Create data frame with results
    processed_data[[sim_idx]] <- data.frame(
      param_set = param_set_id,
      param_type = param_type,
      sim_id = sim_idx,
      day = infected_data$step,
      infected = infected_data$counts
    )
  }
  
  # Combine all simulations for this parameter set
  combined_data <- do.call(rbind, processed_data)
  return(combined_data)
}

# --------------------------
# Run simulations and process results
# --------------------------
cat("Running simulations for", nrow(abc_parameters), "parameter sets, with", 
    n_simulations, "replications each...\n")

# Run all simulations
all_results <- run_and_compare()

# Save the full results
saveRDS(all_results, "simulation_comparison_results.rds")

# --------------------------
# Analyze results
# --------------------------
# Calculate summary statistics by parameter set, day, and type
summary_stats <- all_results %>%
  group_by(param_set, day, param_type) %>%
  summarize(
    mean_infected = mean(infected),
    median_infected = median(infected),
    lower_ci = quantile(infected, 0.025),
    upper_ci = quantile(infected, 0.975),
    sd_infected = sd(infected)
  ) %>%
  ungroup()

# Calculate overall statistics by day and type across all parameter sets
overall_stats <- all_results %>%
  group_by(day, param_type) %>%
  summarize(
    mean_infected = mean(infected),
    median_infected = median(infected),
    lower_ci = quantile(infected, 0.025),
    upper_ci = quantile(infected, 0.975),
    sd_infected = sd(infected)
  ) %>%
  ungroup()

# --------------------------
# Create visualizations
# --------------------------

# 1. Overall average infected counts with confidence intervals
p1 <- ggplot(overall_stats, aes(x = day, y = mean_infected, color = param_type)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = param_type), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("true" = "blue", "abc" = "red"),
                     labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
  scale_fill_manual(values = c("true" = "blue", "abc" = "red"),
                    labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
  labs(title = "Average Infected Counts Over Time",
       subtitle = paste("Based on", nrow(abc_parameters), "parameter sets with", n_simulations, "replications each"),
       x = "Day", y = "Number of Infected",
       color = "Parameter Type", fill = "Parameter Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save overall plot
ggsave("overall_comparison.png", plot = p1, width = 10, height = 6)

# 2. Sample individual parameter sets
selected_param_sets <- sample(unique(summary_stats$param_set), min(6, length(unique(summary_stats$param_set))))

for (ps in selected_param_sets) {
  # Filter data for this parameter set
  ps_data <- summary_stats %>% filter(param_set == ps)
  
  # Create plot
  p <- ggplot(ps_data, aes(x = day, y = mean_infected, color = param_type)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = param_type), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("true" = "blue", "abc" = "red"),
                       labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
    scale_fill_manual(values = c("true" = "blue", "abc" = "red"),
                      labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
    labs(title = paste("Parameter Set", ps, "- Infected Counts Over Time"),
         x = "Day", y = "Number of Infected",
         color = "Parameter Type", fill = "Parameter Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Save individual plot
  ggsave(paste0("param_set_", ps, "_comparison.png"), plot = p, width = 8, height = 5)
}

# 3. Create data for error analysis
error_data <- summary_stats %>%
  select(param_set, day, param_type, mean_infected) %>%
  pivot_wider(names_from = param_type, values_from = mean_infected) %>%
  mutate(
    absolute_error = abc - true,
    relative_error = ifelse(true > 0, (abc - true) / true * 100, NA)
  )

# Create error plot
error_plot <- ggplot(error_data, aes(x = day, y = absolute_error, group = param_set)) +
  geom_line(alpha = 0.2, color = "gray") +
  geom_smooth(aes(group = 1), color = "red", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Absolute Error in Infected Counts: ABC vs True Parameters",
       x = "Day", y = "Error (ABC - True)") +
  theme_minimal()

# Save error plot
ggsave("absolute_error.png", plot = error_plot, width = 10, height = 6)

cat("\nSimulation and visualization complete!\n")
cat("Results saved to simulation_comparison_results.rds\n")
cat("Plots saved to current directory\n")



# Load required libraries
library(epiworldR)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)

# --------------------------
# Set global parameters
# --------------------------
n_population <- 5000     # Population size
recovery_rate <- 0.1     # Fixed recovery rate
model_ndays <- 50        # Simulation duration (days)
n_simulations <- 50      # Number of simulations per parameter set

# --------------------------
# Main function to run simulations for both true and ABC parameters
# --------------------------
run_and_compare <- function() {
  # Initialize lists to store results
  true_results_list <- list()
  abc_results_list <- list()
  
  # Process each row in abc_parameters
  for (i in 1:nrow(abc_parameters)) {
    cat("Processing parameter set", i, "of", nrow(abc_parameters), "\n")
    
    # Extract parameter values
    true_preval <- abc_parameters$true_preval[i]
    true_crate <- abc_parameters$true_crate[i]
    true_ptran <- abc_parameters$true_ptran[i]
    
    abc_preval <- abc_parameters$abc_preval[i]
    abc_crate <- abc_parameters$abc_crate[i]
    abc_ptran <- abc_parameters$abc_ptran[i]
    
    # Run simulations with true parameters
    cat("  Running simulations with true parameters...\n")
    model_true <- ModelSIRCONN(
      name = paste0("True_", i),
      n = n_population,
      prevalence = true_preval,
      contact_rate = true_crate,
      transmission_rate = true_ptran,
      recovery_rate = recovery_rate
    )
    
    # Create saver for true parameters
    saver_true <- make_saver("total_hist", "reproductive")
    
    # Run the simulations with true parameters
    run_multiple(model_true, ndays = model_ndays, nsims = n_simulations, 
                 saver = saver_true, nthread = 2)
    
    # Retrieve the results for true parameters
    true_results <- run_multiple_get_results(model_true)
    
    # Run simulations with ABC parameters
    cat("  Running simulations with ABC parameters...\n")
    model_abc <- ModelSIRCONN(
      name = paste0("ABC_", i),
      n = n_population,
      prevalence = abc_preval,
      contact_rate = abc_crate,
      transmission_rate = abc_ptran,
      recovery_rate = recovery_rate
    )
    
    # Create saver for ABC parameters
    saver_abc <- make_saver("total_hist", "reproductive")
    
    # Run the simulations with ABC parameters
    run_multiple(model_abc, ndays = model_ndays, nsims = n_simulations, 
                 saver = saver_abc, nthread = 2)
    
    # Retrieve the results for ABC parameters
    abc_results <- run_multiple_get_results(model_abc)
    
    # Process results for this parameter set
    true_data <- process_simulation_results(true_results, "true", i)
    abc_data <- process_simulation_results(abc_results, "abc", i)
    
    # Store results
    true_results_list[[i]] <- true_data
    abc_results_list[[i]] <- abc_data
  }
  
  # Combine all results
  all_true_results <- do.call(rbind, true_results_list)
  all_abc_results <- do.call(rbind, abc_results_list)
  
  # Combine true and ABC results
  all_results <- rbind(all_true_results, all_abc_results)
  
  return(all_results)
}

# Function to process simulation results
process_simulation_results <- function(results, param_type, param_set_id) {
  # Create a list to store results for each simulation
  processed_data <- list()
  
  # First check the structure of the total_hist data
  if (is.list(results$total_hist) && length(results$total_hist) > 0) {
    for (sim_idx in 1:length(results$total_hist)) {
      # Get the data for this simulation
      sim_result <- results$total_hist[[sim_idx]]
      
      # Print the structure of the first result for debugging
      if (sim_idx == 1 && param_set_id == 1) {
        cat("Structure of total_hist data:\n")
        print(str(sim_result))
      }
      
      # Extract infected counts - this will depend on the exact structure
      # For data frames with 'state' column:
      if (is.data.frame(sim_result) && "state" %in% colnames(sim_result)) {
        infected_data <- sim_result[sim_result$state == "Infected", ]
        
        processed_data[[sim_idx]] <- data.frame(
          param_set = param_set_id,
          param_type = param_type,
          sim_id = sim_idx,
          day = infected_data$step,
          infected = infected_data$counts
        )
      } 
      # For lists with 'Infected' element:
      else if (is.list(sim_result) && !is.null(sim_result$Infected)) {
        infected_counts <- sim_result$Infected
        
        processed_data[[sim_idx]] <- data.frame(
          param_set = param_set_id,
          param_type = param_type,
          sim_id = sim_idx,
          day = 0:(length(infected_counts)-1),
          infected = infected_counts
        )
      }
      # For arrays or matrices with named dimensions:
      else if ((is.array(sim_result) || is.matrix(sim_result)) && 
               !is.null(dimnames(sim_result)) && 
               !is.null(dimnames(sim_result)[[1]]) && 
               "Infected" %in% dimnames(sim_result)[[1]]) {
        infected_row <- which(dimnames(sim_result)[[1]] == "Infected")
        infected_counts <- sim_result[infected_row, ]
        
        processed_data[[sim_idx]] <- data.frame(
          param_set = param_set_id,
          param_type = param_type,
          sim_id = sim_idx,
          day = 0:(length(infected_counts)-1),
          infected = infected_counts
        )
      }
      # For a simple vector (unlikely but possible):
      else if (is.vector(sim_result) && !is.list(sim_result)) {
        processed_data[[sim_idx]] <- data.frame(
          param_set = param_set_id,
          param_type = param_type,
          sim_id = sim_idx,
          day = 0:(length(sim_result)-1),
          infected = sim_result
        )
      }
      # Default case - try to handle unknown structure:
      else {
        cat("WARNING: Unknown data structure for simulation", sim_idx, 
            "in parameter set", param_set_id, "\n")
        
        # Try to extract using dplyr if available:
        if (requireNamespace("dplyr", quietly = TRUE)) {
          if (is.data.frame(sim_result)) {
            infected_data <- dplyr::filter(sim_result, state == "Infected")
            if (nrow(infected_data) > 0) {
              processed_data[[sim_idx]] <- data.frame(
                param_set = param_set_id,
                param_type = param_type,
                sim_id = sim_idx,
                day = infected_data$step,
                infected = infected_data$counts
              )
            }
          }
        }
      }
    }
    
    # Combine all simulations for this parameter set
    if (length(processed_data) > 0) {
      combined_data <- do.call(rbind, processed_data)
      return(combined_data)
    } else {
      cat("ERROR: Could not process data for parameter set", param_set_id, 
          "with type", param_type, "\n")
      return(NULL)
    }
  } else {
    cat("ERROR: total_hist data not available or empty for parameter set", 
        param_set_id, "with type", param_type, "\n")
    return(NULL)
  }
}

# --------------------------
# Run simulations and process results
# --------------------------
cat("Running simulations for", nrow(abc_parameters), "parameter sets, with", 
    n_simulations, "replications each...\n")

# Run all simulations
all_results <- run_and_compare()

# Save the full results
saveRDS(all_results, "simulation_comparison_results.rds")

# --------------------------
# Analyze results
# --------------------------
# Calculate summary statistics by parameter set, day, and type
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming all_results is already loaded in the environment
# If not, load it from the saved RDS file:
# all_results <- readRDS("simulation_comparison_results.rds")

# Check what parameter sets are available
unique_params <- unique(all_results$param_set)
cat("Available parameter sets:", paste(unique_params, collapse = ", "), "\n")

# --------------------------
# Calculate summary statistics
# --------------------------

# Calculate overall statistics across all simulations for each day and parameter type
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming all_results is already loaded in the environment
# If not, load it from the saved RDS file:
# all_results <- readRDS("simulation_comparison_results.rds")

# Check what parameter sets are available
unique_params <- unique(all_results$param_set)
cat("Available parameter sets:", paste(unique_params, collapse = ", "), "\n")

# --------------------------
# Calculate summary statistics
# --------------------------

# Calculate overall statistics across all simulations for each day and parameter type
daily_stats <- all_results %>%
  group_by(day, param_type) %>%
  summarize(
    mean_infected = mean(infected, na.rm = TRUE),
    median_infected = median(infected, na.rm = TRUE),
    # Use more robust approach for confidence intervals
    lower_ci = if(length(infected) > 1) max(min(infected), mean(infected) - 1.96 * sd(infected)) else mean(infected),
    upper_ci = if(length(infected) > 1) min(max(infected), mean(infected) + 1.96 * sd(infected)) else mean(infected),
    sd_infected = if(length(infected) > 1) sd(infected) else 0,
    n_sims = n()
  ) %>%
  ungroup()

# Calculate statistics for each parameter set separately
param_set_stats <- all_results %>%
  group_by(param_set, day, param_type) %>%
  summarize(
    mean_infected = mean(infected, na.rm = TRUE),
    median_infected = median(infected, na.rm = TRUE),
    # Use more robust approach for confidence intervals
    lower_ci = if(length(infected) > 1) max(min(infected), mean(infected) - 1.96 * sd(infected)) else mean(infected),
    upper_ci = if(length(infected) > 1) min(max(infected), mean(infected) + 1.96 * sd(infected)) else mean(infected),
    sd_infected = if(length(infected) > 1) sd(infected) else 0,
    n_sims = n()
  ) %>%
  ungroup()

# Save the summary statistics
write.csv(daily_stats, "overall_daily_stats.csv", row.names = FALSE)
write.csv(param_set_stats, "parameter_set_daily_stats.csv", row.names = FALSE)

# --------------------------
# Create visualizations
# --------------------------

# 1. Overall plot with confidence intervals
overall_plot <- ggplot(daily_stats, aes(x = day, y = mean_infected, color = param_type)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = param_type), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("true" = "blue", "abc" = "red"),
                     labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
  scale_fill_manual(values = c("true" = "blue", "abc" = "red"),
                    labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
  labs(title = "Average Infected Counts Over Time with 95% CI",
       x = "Day", y = "Mean Number of Infected",
       color = "Parameter Type", fill = "Parameter Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the overall plot
ggsave("overall_epidemic_curves.png", plot = overall_plot, width = 10, height = 6)

# 2. Individual plots for each parameter set
for (ps in unique_params) {
  # Filter data for this parameter set
  ps_data <- param_set_stats %>% filter(param_set == ps)
  
  # Create plot
  ps_plot <- ggplot(ps_data, aes(x = day, y = mean_infected, color = param_type)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = param_type), alpha = 0.2, color = NA) +
    scale_color_manual(values = c("true" = "blue", "abc" = "red"),
                       labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
    scale_fill_manual(values = c("true" = "blue", "abc" = "red"),
                      labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
    labs(title = paste("Parameter Set", ps, "- Average Infected Counts Over Time with 95% CI"),
         x = "Day", y = "Mean Number of Infected",
         color = "Parameter Type", fill = "Parameter Type") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Save parameter set specific plot
  ggsave(paste0("parameter_set_", ps, "_epidemic_curves.png"), plot = ps_plot, width = 10, height = 6)
}

# 3. Create error analysis if both true and abc are present
if (all(c("true", "abc") %in% unique(all_results$param_type))) {
  # Calculate errors
  error_data <- param_set_stats %>%
    select(param_set, day, param_type, mean_infected) %>%
    pivot_wider(names_from = param_type, values_from = mean_infected) %>%
    filter(!is.na(true) & !is.na(abc)) %>%
    mutate(
      absolute_error = abc - true,
      relative_error = ifelse(true > 0, (abc - true) / true * 100, NA)
    )
  
  # Create absolute error plot
  abs_error_plot <- ggplot(error_data, aes(x = day, y = absolute_error, group = param_set)) +
    geom_line(alpha = 0.3, color = "gray") +
    geom_smooth(aes(group = 1), color = "red", method = "loess", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Absolute Error: ABC vs True Parameters",
         x = "Day", y = "Error (ABC - True)") +
    theme_minimal()
  
  # Create relative error plot
  rel_error_plot <- ggplot(error_data, aes(x = day, y = relative_error, group = param_set)) +
    geom_line(alpha = 0.3, color = "gray") +
    geom_smooth(aes(group = 1), color = "red", method = "loess", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Relative Error: ABC vs True Parameters",
         x = "Day", y = "Relative Error (%)") +
    theme_minimal() +
    ylim(-100, 100)  # Limit y-axis for better visualization
  
  # Save error plots
  ggsave("absolute_error.png", plot = abs_error_plot, width = 10, height = 6)
  ggsave("relative_error.png", plot = rel_error_plot, width = 10, height = 6)
}

# 4. Violin plot of infected distribution by day and parameter type
# Only include a subset of days to make the plot clearer
if (length(unique(all_results$day)) > 10) {
  days_to_plot <- sort(unique(all_results$day))[seq(1, length(unique(all_results$day)), by = 10)]
} else {
  days_to_plot <- sort(unique(all_results$day))
}

violin_data <- all_results %>%
  filter(day %in% days_to_plot)

# Check if there's enough variation to create a violin plot
if (length(unique(violin_data$infected)) > 1) {
  violin_plot <- ggplot(violin_data, aes(x = factor(day), y = infected, fill = param_type)) +
    geom_violin(position = position_dodge(width = 0.7), alpha = 0.6) +
    scale_fill_manual(values = c("true" = "blue", "abc" = "red"),
                      labels = c("true" = "True Parameters", "abc" = "ABC Parameters")) +
    labs(title = "Distribution of Infected Counts on Selected Days",
         x = "Day", y = "Number of Infected",
         fill = "Parameter Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save violin plot
  ggsave("infected_distribution_violin.png", plot = violin_plot, width = 12, height = 7)
} else {
  cat("Not enough variation in infected counts to create a meaningful violin plot.\n")
}

cat("Analysis complete. Results and plots saved to current directory.\n")