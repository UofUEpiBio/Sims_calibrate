# Install required packages if not already installed
# install.packages(c("reticulate", "epiworldR", "data.table", "tidyverse", "gridExtra", "cowplot"))

# Load required libraries
library(reticulate)
library(epiworldR)
library(data.table)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(cowplot)

# --------------------------
# Set up Python environment for RNN
# --------------------------
# Initialize Python environment
torch <- import("torch")
nn <- torch$nn
joblib <- import("joblib")
np <- import("numpy")
os <- import("os")

# Define paths (adjust to your actual paths)
model_path <- normalizePath("~/Desktop/Sims_calibrate/model4_bilstm_new.pt")
scaler_additional_path <- normalizePath("~/Desktop/Sims_calibrate/scaler_additional.pkl")
scaler_targets_path <- normalizePath("~/Desktop/Sims_calibrate/scaler_targets.pkl")

# Set up Python environment with model
py_run_string('
import torch
import torch.nn as nn
import joblib
import numpy as np

# Define Model Architecture (BiLSTM)
class BiLSTMModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_layers, additional_dim, output_dim, dropout):
        super(BiLSTMModel, self).__init__()
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        
        # Bidirectional LSTM
        self.bilstm = nn.LSTM(input_dim, hidden_dim, num_layers, batch_first=True, 
                             dropout=dropout, bidirectional=True)
        
        # Fully connected layers
        self.fc1 = nn.Linear(2 * hidden_dim + additional_dim, 64)
        self.fc2 = nn.Linear(64, output_dim)

        self.sigmoid = nn.Sigmoid()
        self.softplus = nn.Softplus()

    def forward(self, x, additional_inputs):
        _, (h_n, _) = self.bilstm(x)
        h_n = torch.cat((h_n[-2], h_n[-1]), dim=1)
        combined = torch.cat((h_n, additional_inputs), dim=1)
        x = self.fc1(combined)
        out = self.fc2(x)
        out = torch.stack([
            self.sigmoid(out[:, 0]),
            self.softplus(out[:, 1]),
            self.softplus(out[:, 2])
        ], dim=1)
        return out
')

# Load model and scalers
py_run_string(sprintf('
# Model parameters
input_dim = 1
hidden_dim = 128
num_layers = 2
additional_dim = 2
output_dim = 3
dropout = 0.3

# Set device
device = torch.device("cpu")

# Initialize model
model = BiLSTMModel(input_dim, hidden_dim, num_layers, additional_dim, output_dim, dropout)

# Load model weights
model.load_state_dict(torch.load("%s", map_location=device))
model.to(device)
model.eval()

# Load scalers
scaler_additional = joblib.load("%s")
scaler_targets = joblib.load("%s")
', model_path, scaler_additional_path, scaler_targets_path))

# Create RNN prediction function
predict_with_bilstm <- function(time_series, n, recov) {
  # Check inputs
  if (!is.numeric(time_series) || !is.numeric(n) || !is.numeric(recov)) {
    stop("Inputs must be numeric")
  }
  
  # Reshape time series to match expected input (batch_size, seq_len, features)
  time_series <- array(time_series, dim = c(1, length(time_series), 1))
  
  # Pass all data to Python
  py$time_series <- time_series
  py$n_value <- n
  py$recov_value <- recov
  
  # Run the prediction in one clean Python block
  py_run_string("
# Scale additional inputs
new_add = np.array([[n_value, recov_value]])
new_add_scaled = scaler_additional.transform(new_add)

# Convert inputs to tensors
new_seq_tensor = torch.tensor(time_series, dtype=torch.float32).to(device)
new_add_tensor = torch.tensor(new_add_scaled, dtype=torch.float32).to(device)

# Run model
with torch.no_grad():
    pred = model(new_seq_tensor, new_add_tensor).cpu().numpy()

# Inverse-transform predictions
pred_unscaled = scaler_targets.inverse_transform(pred)
")
  
  # Extract predictions
  predictions <- py$pred_unscaled[1, ]
  names(predictions) <- c("ptran", "crate", "R0")
  
  return(predictions)
}

# --------------------------
# Global Simulation Settings
# --------------------------
model_ndays <- 50   # simulation duration (days)
model_seed  <- 122  # seed for reproducibility
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
# ABC (LFMCMC) Calibration Functions - MODIFIED to include R0
# --------------------------

# This function simulates using the proposed parameters and returns all daily infected counts
simulation_fun <- function(params, lfmcmc_obj) {
  # Extract parameters for simulation
  # We have: prevalence, R0, recovery_rate, transmission_prob, contact_rate
  # Need to pass to simulate_epidemic_calib: prevalence, contact_rate, recovery_rate, transmission_prob
  
  sim_params <- c(
    params[1],  # prevalence
    params[5],  # contact_rate (calculated from R0, recovery, transmission)
    params[3],  # recovery_rate
    params[4]   # transmission_prob
  )
  
  # Run simulation with the parameters
  infected_counts <- simulate_epidemic_calib(sim_params, ndays = model_ndays)
  return(as.numeric(infected_counts))
}

# This returns the observed data (all daily infected counts)
summary_fun <- function(data, lfmcmc_obj) {
  return(as.numeric(data))
}

# Generate new parameter proposals - MODIFIED to include R0
proposal_fun <- function(old_params, lfmcmc_obj) {
  # Proposals with appropriate step sizes for each parameter
  # old_params contains: prevalence, R0, recovery_rate, transmission_prob, contact_rate
  
  new_preval <- plogis(qlogis(old_params[1]) + rnorm(1, sd = 0.1))
  new_R0     <- old_params[2] * exp(rnorm(1, sd = 0.1))  # Log-normal proposal for R0
  new_recov  <- old_params[3] * exp(rnorm(1, sd = 0.1))  # Log-normal proposal
  new_ptran  <- plogis(qlogis(old_params[4]) + rnorm(1, sd = 0.1))
  
  # Calculate contact rate based on R0, recovery rate, and transmission rate
  # R0 = (contact_rate * transmission_rate) / recovery_rate
  # Therefore: contact_rate = (R0 * recovery_rate) / transmission_rate
  new_crate <- (new_R0 * new_recov) / new_ptran
  
  return(c(new_preval, new_R0, new_recov, new_ptran, new_crate))
}

# Kernel function using sum of squared differences across all days
kernel_fun <- function(simulated_stat, observed_stat, epsilon, lfmcmc_obj) {
  # Calculate the sum of squared differences across all days
  diff <- sum((simulated_stat - observed_stat)^2)
  return(exp(-diff / (2 * epsilon^2)))
}

# --------------------------
# Function: Simulation–Calibration Study with RNN Comparison - MODIFIED
# --------------------------
simulate_and_calibrate <- function(true_params, sim_id) {
  cat("Running simulation", sim_id, "of", N_SIMS, "\n")
  
  # Calculate true R0 from true parameters
  # R0 = (contact_rate * transmission_rate) / recovery_rate
  true_R0 <- (true_params[3] * true_params[5]) / true_params[4]
  
  # Step 1: Simulate the observed epidemic using the true parameters
  observed_infected <- simulate_epidemic_observed(true_params, ndays = model_ndays, 
                                                  seed = model_seed + sim_id)
  
  # Step 2: ABC calibration with R0
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
  
  # Initial parameters for ABC: 
  # prevalence, R0, recovery rate, transmission probability, contact rate
  init_params <- c(
    true_params[2],              # prevalence
    true_R0,                     # R0
    true_params[4],              # recovery rate
    true_params[5],              # transmission probability
    true_params[3]               # contact rate (derived from the above)
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
    abc_preval <- calibrated_params_raw[1]
    abc_R0 <- calibrated_params_raw[2]
    abc_recov <- calibrated_params_raw[3]
    abc_ptran <- calibrated_params_raw[4]
    abc_crate <- calibrated_params_raw[5]  # This is derived from R0, recov, ptran
    
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
  
  # Step 4: RNN prediction using observed time series
  time_series_for_rnn <- observed_infected
  n_rnn <- global_n
  recov_rnn <- true_params[4]  # Use true recovery rate as input to RNN
  
  # Get RNN predictions
  tryCatch({
    rnn_predictions <- predict_with_bilstm(time_series_for_rnn, n_rnn, recov_rnn)
    
    # Extract RNN-predicted parameters
    rnn_ptran <- rnn_predictions["ptran"]
    rnn_crate <- rnn_predictions["crate"]  # This is the predicted contact rate directly from RNN
    rnn_R0 <- rnn_predictions["R0"]
    
    # Calculate TRUE contact rate based on R0, recovery rate, and transmission rate
    # R0 = (contact_rate * transmission_rate) / recovery_rate
    # Therefore: contact_rate = (R0 * recovery_rate) / transmission_rate
    rnn_true_crate <- (rnn_R0 * recov_rnn) / rnn_ptran
    
    # Calculate prevalence from RNN using initial values
    initial_infected <- observed_infected[1]
    rnn_preval <- initial_infected / n_rnn
    
    # Create RNN parameters vector with the derived true contact rate
    rnn_params <- c(rnn_preval, rnn_true_crate, recov_rnn, rnn_ptran)
    
    # Step 5: Simulate predictions using RNN-predicted parameters with calculated contact rate
    rnn_predicted_infected <- simulate_epidemic_calib(rnn_params, ndays = model_ndays, 
                                                      seed = model_seed + sim_id + 300)
  }, error = function(e) {
    cat("Error in RNN prediction for simulation", sim_id, ":", e$message, "\n")
    rnn_preval <- NA
    rnn_crate <- NA  # Direct from RNN
    rnn_ptran <- NA
    rnn_R0 <- NA
    rnn_true_crate <- NA  # Derived from R0
    rnn_predicted_infected <- rep(NA, model_ndays + 1)
  })
  
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
    # RNN predicted parameters
    rnn_preval = rnn_preval,
    rnn_crate = rnn_crate,  # Direct from RNN
    rnn_recov = recov_rnn,
    rnn_ptran = rnn_ptran,
    rnn_R0 = rnn_R0,
    rnn_true_crate = rnn_true_crate,  # Derived from R0
    # Observed and predicted counts for each day
    observed_infected = observed_infected,
    abc_predicted_infected = abc_predicted_infected,
    rnn_predicted_infected = rnn_predicted_infected,
    # Bias calculation for both methods
    abc_bias = abc_predicted_infected - observed_infected,
    rnn_bias = rnn_predicted_infected - observed_infected,
    # Relative bias (to account for scale differences)
    abc_rel_bias = ifelse(observed_infected > 0, 
                          (abc_predicted_infected - observed_infected) / observed_infected, 
                          NA),
    rnn_rel_bias = ifelse(observed_infected > 0, 
                          (rnn_predicted_infected - observed_infected) / observed_infected, 
                          NA)
  )
  
  # Parameter comparison for both methods - MODIFIED to include R0 and derived contact rate
  param_comparison <- data.frame(
    sim_id = sim_id,
    parameter = c("prevalence", "contact_rate", "recovery_rate", "transmission_prob", "R0", "rnn_true_contact_rate"),
    true_value = c(true_params[2:5], true_R0, true_params[3]),
    abc_calibrated_value = c(calibrated_params, abc_R0, abc_crate),
    rnn_predicted_value = c(rnn_preval, rnn_crate, recov_rnn, rnn_ptran, rnn_R0, rnn_true_crate),
    abc_error_pct = c(
      (calibrated_params[1] - true_params[2]) / true_params[2] * 100,
      (calibrated_params[2] - true_params[3]) / true_params[3] * 100,
      (calibrated_params[3] - true_params[4]) / true_params[4] * 100,
      (calibrated_params[4] - true_params[5]) / true_params[5] * 100,
      (abc_R0 - true_R0) / true_R0 * 100,
      (abc_crate - true_params[3]) / true_params[3] * 100
    ),
    rnn_error_pct = c(
      (rnn_preval - true_params[2]) / true_params[2] * 100,
      (rnn_crate - true_params[3]) / true_params[3] * 100,
      (recov_rnn - true_params[4]) / true_params[4] * 100,
      (rnn_ptran - true_params[5]) / true_params[5] * 100,
      (rnn_R0 - true_R0) / true_R0 * 100,
      (rnn_true_crate - true_params[3]) / true_params[3] * 100
    )
  )
  
  return(list(
    daily_results = daily_results,
    param_comparison = param_comparison
  ))
}

# --------------------------
# Run Simulations in Parallel
# --------------------------
# Determine how many cores to use
n_cores <- min(detectCores() - 1, N_SIMS)
if (n_cores < 1) n_cores <- 1
n_cores=18  # Adjust as needed
cat("Running", N_SIMS, "simulations using", n_cores, "cores...\n")
set.seed(model_seed)  # Set global seed for reproducibility

# Run simulations in parallel if available
if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
  results_list <- mclapply(1:N_SIMS, function(i) {
    simulate_and_calibrate(as.numeric(theta_use[i]), i)
  }, mc.cores = n_cores)
} else {
  # Sequential processing
  results_list <- list()
  for (i in 1:N_SIMS) {
    results_list[[i]] <- simulate_and_calibrate(as.numeric(theta_use[i]), i)
  }
}

# --------------------------
# Combine and Analyze Results
# --------------------------
# Extract and combine all daily results
daily_results <- bind_rows(lapply(results_list, function(x) x$daily_results))

# Extract and combine all parameter comparisons
param_comparisons <- bind_rows(lapply(results_list, function(x) x$param_comparison))
write.csv(daily_stats, "abc_rnn_comparison_daily_stats.csv", row.names = FALSE)
write.csv(param_comparisons, "abc_rnn_param_comparisons.csv", row.names = FALSE)


# Save the raw daily_results and param_comparisons dataframes
saveRDS(daily_results, "daily_results.rds")
saveRDS(param_comparisons, "param_comparisons.rds")

# Additionally save as CSV for easier access if needed
write.csv(daily_results, "daily_results.csv", row.names = FALSE)
write.csv(param_comparisons, "param_comparisons.csv", row.names = FALSE)

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
  rnn_predicted_values <- day_data$rnn_predicted_infected
  
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
  
  # Calculate confidence intervals for RNN predicted values
  rnn_pred_ci <- calc_CI(rnn_predicted_values)
  coverage_stats$rnn_pred_mean[coverage_stats$day == day_num] <- rnn_pred_ci["mean"]
  coverage_stats$rnn_pred_lower_ci[coverage_stats$day == day_num] <- rnn_pred_ci["lower_ci"]
  coverage_stats$rnn_pred_upper_ci[coverage_stats$day == day_num] <- rnn_pred_ci["upper_ci"]
  
  # Calculate coverage for ABC
  abc_coverage <- mean(abc_predicted_values >= obs_ci["lower_ci"] & 
                         abc_predicted_values <= obs_ci["upper_ci"], na.rm = TRUE) * 100
  coverage_stats$abc_coverage[coverage_stats$day == day_num] <- abc_coverage
  
  # Calculate coverage for RNN
  rnn_coverage <- mean(rnn_predicted_values >= obs_ci["lower_ci"] & 
                         rnn_predicted_values <= obs_ci["upper_ci"], na.rm = TRUE) * 100
  coverage_stats$rnn_coverage[coverage_stats$day == day_num] <- rnn_coverage
}

# Calculate daily statistics across simulations
daily_stats <- daily_results %>%
  group_by(day) %>%
  summarize(
    mean_observed = mean(observed_infected),
    abc_mean_predicted = mean(abc_predicted_infected),
    rnn_mean_predicted = mean(rnn_predicted_infected, na.rm = TRUE),
    abc_mean_bias = mean(abc_bias),
    rnn_mean_bias = mean(rnn_bias, na.rm = TRUE),
    abc_mean_abs_bias = mean(abs(abc_bias)),
    rnn_mean_abs_bias = mean(abs(rnn_bias), na.rm = TRUE),
    n_sims = n()
  )

# Join with coverage statistics
daily_stats <- left_join(daily_stats, coverage_stats, by = "day")
write.csv(daily_stats, "daily_stats.csv", row.names = FALSE)
saveRDS(daily_stats, "daily_stats.rds")

# --------------------------
# Enhanced Visualization - Compare ABC and RNN - MODIFIED to include R0
# --------------------------

# 1. Combined mean bias comparison plot
bias_comparison_plot <- ggplot(daily_stats, aes(x = day)) +
  geom_line(aes(y = abc_mean_bias, color = "ABC"), size = 1) +
  geom_line(aes(y = rnn_mean_bias, color = "RNN"), size = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("ABC" = "blue", "RNN" = "red")) +
  labs(title = "Mean Bias Comparison: ABC vs RNN",
       subtitle = paste(N_SIMS, "simulations"),
       x = "Day", y = "Mean Bias (Predicted - Observed)") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 2. Epidemic curves with 3 methods
epidemic_curve_comparison <- ggplot(daily_stats, aes(x = day)) +
  geom_line(aes(y = obs_mean, color = "Observed"), size = 1) +
  geom_line(aes(y = abc_pred_mean, color = "ABC"), size = 1, linetype = "dashed") +
  geom_line(aes(y = rnn_pred_mean, color = "RNN"), size = 1, linetype = "dotted") +
  geom_ribbon(aes(ymin = abc_pred_lower_ci, ymax = abc_pred_upper_ci), fill = "blue", alpha = 0.1) +
  geom_ribbon(aes(ymin = rnn_pred_lower_ci, ymax = rnn_pred_upper_ci), fill = "red", alpha = 0.1) +
  scale_color_manual(values = c("Observed" = "black", "ABC" = "blue", "RNN" = "red")) +
  labs(title = "Epidemic Curves: Observed vs ABC vs RNN",
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
    geom_point(aes(y = abc_calibrated_value, color = "ABC"), alpha = 0.6) +
    geom_point(aes(y = rnn_predicted_value, color = "RNN"), alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
    scale_color_manual(values = c("ABC" = "blue", "RNN" = "red")) +
    labs(title = str_to_title(gsub("_", " ", param)),
         x = "True Value", 
         y = "Estimated Value") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  param_plots[[param]] <- p
}

# 4. New plot for comparing derived contact rate
derived_contact_rate_plot <- ggplot(
  param_comparisons[param_comparisons$parameter == "rnn_true_contact_rate", ], 
  aes(x = true_value)
) +
  geom_point(aes(y = abc_calibrated_value, color = "ABC"), alpha = 0.6) +
  geom_point(aes(y = rnn_predicted_value, color = "RNN Derived"), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("ABC" = "blue", "RNN Derived" = "red")) +
  labs(title = "Derived Contact Rate Comparison",
       subtitle = "Contact rate derived from R0, recovery rate, and transmission rate",
       x = "True Contact Rate", 
       y = "Derived Contact Rate") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 5. Comprehensive comparison dashboard
dashboard <- cowplot::plot_grid(
  bias_comparison_plot,
  epidemic_curve_comparison,
  param_plots$prevalence,
  param_plots$R0,
  param_plots$transmission_prob,
  derived_contact_rate_plot,
  ncol = 2,
  nrow = 3,
  labels = "AUTO"
)

# Calculate aggregate statistics for parameter estimation
param_stats <- param_comparisons %>%
  group_by(parameter) %>%
  summarize(
    mean_true = mean(true_value, na.rm = TRUE),
    mean_abc = mean(abc_calibrated_value, na.rm = TRUE),
    mean_rnn = mean(rnn_predicted_value, na.rm = TRUE),
    mean_abc_error_pct = mean(abc_error_pct, na.rm = TRUE),
    mean_rnn_error_pct = mean(rnn_error_pct, na.rm = TRUE),
    rmse_abc = sqrt(mean((abc_calibrated_value - true_value)^2, na.rm = TRUE)),
    rmse_rnn = sqrt(mean((rnn_predicted_value - true_value)^2, na.rm = TRUE))
  )

# Create a comparison of ABC and RNN for R0 specifically
r0_comparison <- ggplot(
  param_comparisons[param_comparisons$parameter == "R0", ],
  aes(x = sim_id)
) +
  geom_point(aes(y = true_value, color = "True R0"), size = 3) +
  geom_point(aes(y = abc_calibrated_value, color = "ABC R0"), size = 2) +
  geom_point(aes(y = rnn_predicted_value, color = "RNN R0"), size = 2) +
  scale_color_manual(values = c("True R0" = "black", "ABC R0" = "blue", "RNN R0" = "red")) +
  labs(title = "R0 Estimation Comparison Across Simulations",
       x = "Simulation ID", 
       y = "R0 Value") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Display main plots
print(bias_comparison_plot)
print(epidemic_curve_comparison)
print(r0_comparison)
print(derived_contact_rate_plot)
print(dashboard)

# Print parameter statistics table
print(param_stats)

# --------------------------
# Save Results and Plots
# --------------------------
# Save data
write.csv(daily_stats, "abc_rnn_comparison_daily_stats.csv", row.names = FALSE)
write.csv(param_comparisons, "abc_rnn_param_comparisons.csv", row.names = FALSE)
write.csv(param_stats, "abc_rnn_parameter_statistics.csv", row.names = FALSE)

# Save plots
ggsave("bias_comparison.png", plot = bias_comparison_plot, width = 10, height = 6)
ggsave("epidemic_curve_comparison.png", plot = epidemic_curve_comparison, width = 10, height = 6)
ggsave("r0_comparison.png", plot = r0_comparison, width = 12, height = 6)
ggsave("derived_contact_rate.png", plot = derived_contact_rate_plot, width = 10, height = 6)
ggsave("comprehensive_comparison.png", plot = dashboard, width = 14, height = 16)

# --------------------------
# Additional Analysis - R0 and Contact Rate Relationship
# --------------------------

# Create a dataframe to examine the relationship between R0 and contact rate
r0_crate_relationship <- daily_results %>%
  filter(day == 0) %>%  # Just use one row per simulation
  select(sim_id, true_R0, true_crate, 
         calib_R0, calib_crate, 
         rnn_R0, rnn_true_crate) %>%
  distinct()

# Plot the relationship
r0_crate_plot <- ggplot(r0_crate_relationship, aes(x = true_R0)) +
  geom_point(aes(y = true_crate, color = "True"), size = 3) +
  geom_point(aes(y = calib_crate, color = "ABC"), size = 2) +
  geom_point(aes(y = rnn_true_crate, color = "RNN Derived"), size = 2) +
  scale_color_manual(values = c("True" = "black", "ABC" = "blue", "RNN Derived" = "red")) +
  labs(title = "Relationship: R0 and Contact Rate",
       x = "R0", 
       y = "Contact Rate") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(r0_crate_plot)
ggsave("r0_contact_rate_relationship.png", plot = r0_crate_plot, width = 10, height = 6)

cat("\nSimulation and analysis complete with ABC and RNN comparison including R0 and derived contact rates!\n")
