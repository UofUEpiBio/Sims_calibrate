# ───────────────────────────────────
# Complete code to compare epidemic curves for ONE parameter set
# ───────────────────────────────────

# 1) libraries
library(tidyverse)
library(glue)
source("~/Desktop/Sims_calibrate/01a-bilstm.R")

library(epiworldR)
# ──────────────────────────────────────────────────────────────────
# 2) Load ABC predictions and run LSTM predictions
# ──────────────────────────────────────────────────────────────────

filename <- paste0("~/Desktop/Sims_calibrate/abc_parameters.csv")
ndays    <- 60
nsims    <- 1
nthreads <- 10
abc_pred <- read_csv(filename)

# Run LSTM predictions for all parameter sets
lstm_results <- vector("list", nrow(abc_pred))
# Your for-loop code here
for (i in seq_len(nrow(abc_pred))) {
  sim_id <- abc_pred$sim_id[i]
  cat(sprintf("Processing sim_id %d (%d/%d)\n", sim_id, i, nrow(abc_pred)))
  
  tryCatch({
    # True SIR run
    m0 <- ModelSIRCONN(
      name              = paste0("true_sim_", sim_id),
      n                 = abc_pred$true_n[i],
      prevalence        = abc_pred$true_preval[i],
      contact_rate      = abc_pred$true_crate[i],
      transmission_rate = abc_pred$true_ptran[i],
      recovery_rate     = abc_pred$true_recov[i]
    )
    run(m0, ndays = ndays)
    incidence <- plot_incidence(m0, plot = FALSE)
    ts        <- incidence[, 1]
    
    # LSTM prediction
    joblib        <- import("joblib")
    library(reticulate)
    lstm_out      <- predict_with_bilstm(ts, abc_pred$true_n[i], abc_pred$true_recov[i])
    ptran_lstm    <- lstm_out[1]
    crate_lstm    <- lstm_out[2]
    R0_lstm       <- lstm_out[3]
    prevalence_lstm    <- abc_pred$true_preval[i]
    contact_rate_lstm  <- R0_lstm * abc_pred$true_recov[i] / ptran_lstm
    transmission_rate_lstm <- ptran_lstm
    
    lstm_results[[i]] <- tibble(
      sim_id               = sim_id,
      prevalence_lstm      = prevalence_lstm,
      contact_rate_lstm    = contact_rate_lstm,
      transmission_rate_lstm = transmission_rate_lstm,
      R0_lstm              = R0_lstm
    )
    cat(sprintf("  ✓ Success for sim_id %d\n", sim_id))
    
  }, error = function(e) {
    cat(sprintf("  ✗ Error for sim_id %d: %s\n", sim_id, e$message))
    lstm_results[[i]] <- tibble(
      sim_id               = sim_id,
      prevalence_lstm      = NA,
      contact_rate_lstm    = NA,
      transmission_rate_lstm = NA,
      R0_lstm              = NA
    )
  })
}

lstm_preds <- bind_rows(lstm_results) %>%
  filter(!is.na(prevalence_lstm))

cat(sprintf("\nSummary:\n"))
cat(sprintf("Total simulations: %d\n", nrow(abc_pred)))
cat(sprintf("Successful LSTM predictions: %d\n", nrow(lstm_preds)))
cat(sprintf("Success rate: %.2f%%\n", 100 * nrow(lstm_preds) / nrow(abc_pred)))

# ──────────────────────────────────────────────────────────────────
# 3) Stack true / abc / lstm parameters
# ──────────────────────────────────────────────────────────────────

params_all <- abc_pred %>%
  select(sim_id,
         true_n,
         prevalence_true        = true_preval,
         contact_rate_true      = true_crate,
         transmission_rate_true = true_ptran,
         recovery_rate_true     = true_recov,
         R0_true                = true_R0,
         R0_abc                 = abc_R0,
         prevalence_abc         = abc_preval,
         contact_rate_abc       = abc_crate,
         transmission_rate_abc  = abc_ptran,
         recovery_rate_abc      = abc_recov) %>%
  inner_join(lstm_preds, by = "sim_id") %>%
  mutate(
    recovery_rate_lstm = recovery_rate_true
  ) %>%
  pivot_longer(
    cols = -c(sim_id, true_n),
    names_to  = c(".value", "param_type"),
    names_pattern = "(.*)_(true|abc|lstm)"
  )

write.csv(params_all, "params_all.csv", row.names = FALSE)

# ──────────────────────────────────────────────────────────────────
# 4) Pick one parameter set for comparison
# ──────────────────────────────────────────────────────────────────
# Identify the sim_ids you want to keep (those with abc values not both over 30)
set.seed(120)
#params_all_filtered%>% filter(sim_id==17)
chosen_sim_id <- sample.int(1000, size = 1)

single_param_set <- params_all %>%
  filter(sim_id == chosen_sim_id)

cat(sprintf("Comparing epidemic curves for sim_id: %d\n", chosen_sim_id))
print(single_param_set)


###combined plot with CI 
n_reps <- 100
single_param_set$run_id <- seq_len(nrow(single_param_set))
single_param_set[2,7] <- single_param_set[1,7]

# Container to collect results
all_avg_series <- tibble()

# Loop over parameter sets
for (i in seq_len(nrow(single_param_set))) {
  row <- single_param_set[i, ]
  
  cat(sprintf("Running %s parameters (rep %d/%d)...\n", row$param_type, row$run_id, n_reps))
  
  # 🔒 Fix random seed for reproducibility
  set.seed(123 + row$sim_id + i)  # ensures same seed for same param row
  
  # Build the model
  calibrated_model <- ModelSIRCONN(
    name              = paste0("single_", row$param_type, "_", row$sim_id),
    n                 = row$true_n,
    prevalence        = row$prevalence,
    contact_rate      = row$contact_rate,
    transmission_rate = row$transmission_rate,
    recovery_rate     = row$recovery_rate
  )
  
  # Run simulations and get results
  saver <- make_saver("total_hist")
  run_multiple(calibrated_model, ndays = ndays, nsims = n_reps, saver = saver, nthreads = 18)
  ans <- run_multiple_get_results(calibrated_model,nthreads=18)
  #incidence <- plot_incidence(ans, plot = FALSE)
  # Extract total hist
  df <- ans$total_hist
  
  # Calculate mean and quantile-based confidence intervals across simulations
  avg_states <- df %>%
    group_by(date, state) %>%
    summarize(
      mean_count = mean(counts),
      # 95% confidence intervals using quantiles
      ci_lower = quantile(counts, 0.025),
      ci_upper = quantile(counts, 0.975),
      .groups = "drop"
    ) %>%
    mutate(
      param_type = row$param_type,
      sim_id     = row$sim_id
    )
  
  # Collect
  all_avg_series <- bind_rows(all_avg_series, avg_states)
}

# === Generate Plots with Confidence Intervals ===
library(ggplot2)
library(patchwork)

# Helper function to create plots with CI
create_sir_plot <- function(data, state_name, title) {
  plot_data <- data %>% filter(state == state_name)
  
  p <- ggplot(plot_data, aes(x = date, color = param_type)) +
    # Confidence interval ribbon only for non-"true" methods
    geom_ribbon(
      data = plot_data %>% filter(param_type != "true"),
      aes(ymin = ci_lower, ymax = ci_upper, fill = param_type), 
      alpha = 0.2, color = NA
    ) +
    # Mean lines for all methods
    geom_line(aes(y = mean_count), size = 1) +
    labs(
      title = title,
      x = "Day", 
      y = paste0("Average ", state_name, " (95% CI)"), 
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12, hjust = 0.5)
    )
  
  return(p)
}

# Create individual plots
p_infected <- create_sir_plot(all_avg_series, "Infected", "Infected Over Time by Method")
p_susceptible <- create_sir_plot(all_avg_series, "Susceptible", "Susceptible Over Time by Method") 
p_recovered <- create_sir_plot(all_avg_series, "Recovered", "Recovered Over Time by Method")

ggsave("p_infected.png", p_infected, width = 12, height = 10, dpi = 300)
ggsave("p_susceptible.png", p_susceptible, width = 12, height = 10, dpi = 300)
ggsave("p_recovered.png", p_recovered, width = 12, height = 10, dpi = 300)
# Update the plotting function to increase y-axis breaks
library(ggplot2)
library(patchwork)
library(dplyr)


# -- Plotting Function with CI and axis breaks --
create_sir_plot_detailed <- function(data, state_name, subtitle) {
  plot_data <- data %>% filter(state == state_name)
  
  y_max <- max(plot_data$ci_upper, na.rm = TRUE)
  y_breaks <- pretty(c(0, y_max), n = 6)
  
  ggplot(plot_data, aes(x = date, color = param_type)) +
    geom_ribbon(
      data = plot_data %>% filter(param_type != "true"),
      aes(ymin = ci_lower, ymax = ci_upper, fill = param_type),
      alpha = 0.2, color = NA
    ) +
    geom_line(aes(y = mean_count), size = 1) +
    scale_y_continuous(breaks = y_breaks) +
    labs(
      title = subtitle,
      x = "Day",
      y = paste(state_name, "(mean ± 95% CI)"),
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12)
    )
}

# -- Create 3 SIR plots --
p1 <- create_sir_plot_detailed(all_avg_series, "Infected", "Infected Over Time")
p2 <- create_sir_plot_detailed(all_avg_series, "Susceptible", "Susceptible Over Time")
p3 <- create_sir_plot_detailed(all_avg_series, "Recovered", "Recovered Over Time")

# -- Combine Horizontally with Shared Legend and Caption --
final_plot <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") +  # <-- Collect legends into one
  plot_annotation(
    title = paste("Epidemic Curves for Randomly Chosen sim_id:", chosen_sim_id),
    caption = "Legend:\nBlue Line = true (no CI)\nRed Line + Light Red = abc (95% CI)\nGreen Line + Light Cyan = lstm (95% CI)",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 11, hjust = 0, margin = margin(t = 10))
    )
  )

# -- Save and Show --
ggsave("epidemic_curves_annotated.png", final_plot, width = 12, height = 7, dpi = 300)
print(final_plot)

create_sir_plot <- function(data, state_name, title) {
  # Filter data to only include days 0-40 and the specified state
  plot_data <- data %>% 
    filter(state == state_name, date <= 40)
  
  p <- ggplot(plot_data, aes(x = date, color = param_type)) +
    # Confidence interval ribbon only for non-"true" methods
    geom_ribbon(
      data = plot_data %>% filter(param_type != "true"),
      aes(ymin = ci_lower, ymax = ci_upper, fill = param_type), 
      alpha = 0.2, color = NA
    ) +
    # Mean lines for all methods
    geom_line(aes(y = mean_count), size = 1) +
    # Set x-axis limits to ensure consistent range across plots
    scale_x_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
    labs(
      title = title,
      x = "Day", 
      y = paste0("Average ", state_name, " (95% CI)"), 
      color = "Method",
      fill = "Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12, hjust = 0.5)
    )
  
  return(p)
}

# Create the plots with the modified function
p_infected <- create_sir_plot(all_avg_series, "Infected", "Infected Over Time by Method")
p_susceptible <- create_sir_plot(all_avg_series, "Susceptible", "Susceptible Over Time by Method") 
p_recovered <- create_sir_plot(all_avg_series, "Recovered", "Recovered Over Time by Method")

# Optional: Create combined plots
library(patchwork)

# Horizontal layout (plots arranged in rows)
combined_horizontal <- p_susceptible / p_infected / p_recovered + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Add title
combined_horizontal <- combined_horizontal + 
  plot_annotation(
    title = "SIR Model Comparison (Days 0-40)",
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
  )

# Display plots
print(combined_horizontal)

# Save plot
ggsave("sir_comparison_day40_horizontal.png", combined_horizontal, width = 10, height = 10, dpi = 900)


####Running over all 1000 simulations 100 time for each 
library(tidyvers)
library(patchwork)
library(ggplot2)
library(dplyr)
library(slurmR)

# Parameters
n_runs <- 100         # Runs per parameter set per method
ndays <- 60          # Days to simulate
n_param_sets <- max(params_all$sim_id)  # Should be 1000
cat("Starting bias and coverage analysis for", n_param_sets, "parameter sets...\n")
cat("Each parameter set will be run", n_runs, "times for", ndays, "days\n")
set.seed(123)

# Function to run simulation for a single parameter set with specific method
run_parameter_set_method <- function(sim_idx, method_name) {
  
  # Load required packages on worker nodes
  library(epiworldR)
  library(dplyr)
  
  # Get parameters for current method - REMOVED dplyr::
  method_params <- params_all %>% 
    filter(sim_id == sim_idx, param_type == method_name)
  
  if (nrow(method_params) == 0) {
    cat("Warning: No parameters found for sim_id", sim_idx, "method", method_name, "\n")
    return(tibble())
  }
  
  # Create model with current parameters
  model <- ModelSIRCONN(
    name = paste0("sim_", sim_idx, "_", method_name),
    n = method_params$true_n,
    prevalence = method_params$prevalence,
    contact_rate = method_params$contact_rate,
    transmission_rate = method_params$transmission_rate,
    recovery_rate = method_params$recovery_rate
  )
  
  # Run multiple simulations
  saver <- epiworldR::make_saver("total_hist")
  epiworldR::run_multiple(model, ndays = ndays, nsims = n_runs, saver = saver, nthreads = 18)  # Reduced threads
  ans <- epiworldR::run_multiple_get_results(model)
  
  # Extract results and focus on Infected state
  raw_hist <- ans$total_hist
  
  # Calculate statistics across the runs for each day - REMOVED dplyr::
  infected_stats <- raw_hist %>%
    filter(state == "Infected") %>%
    group_by(date) %>%
    summarize(
      mean_infected = mean(counts),
      std_infected = sd(counts),
      q025 = quantile(counts, 0.025),  # Lower CI
      q975 = quantile(counts, 0.975),  # Upper CI
      .groups = "drop"
    ) %>%
    mutate(
      sim_id = sim_idx,
      param_type = method_name
    )
  
  # Force garbage collection to free memory
  gc()
  
  return(infected_stats)
}

# Create wrapper function for each method to avoid variable scoping issues
run_true_method <- function(sim_idx) {
  return(run_parameter_set_method(sim_idx, "true"))
}

run_abc_method <- function(sim_idx) {
  return(run_parameter_set_method(sim_idx, "abc"))
}

run_lstm_method <- function(sim_idx) {
  return(run_parameter_set_method(sim_idx, "lstm"))
}

# Run separate SLURM jobs for each method individually
cat("Submitting SLURM jobs for method: true\n")
ans <- Slurm_lapply(
  X = 1:n_param_sets,
  FUN = run_true_method,
  job_name = "bias_coverage_true",
  njobs = 100,
  overwrite = TRUE,
  plan = "submit",
  sbatch_opt = list(
    partition = "vegayon-shared-np",
    account = "vegayon-np",
    time = "02:00:00",
    `mem-per-cpu` = "8G",
    `cpus-per-task` = 1
  ),
  export = c(
    "run_parameter_set_method",
    "run_true_method",
    "params_all",
    "n_runs",
    "ndays"
  )
)

cat("Submitting SLURM jobs for method: abc\n")
ans1 <- Slurm_lapply(
  X = 1:n_param_sets,
  FUN = run_abc_method,
  job_name = "bias_coverage_abc",
  njobs = 100,
  overwrite = TRUE,
  plan = "submit",
  sbatch_opt = list(
    partition = "vegayon-shared-np",
    account = "vegayon-np",
    time = "02:00:00",
    `mem-per-cpu` = "8G",
    `cpus-per-task` = 1
  ),
  export = c(
    "run_parameter_set_method",
    "run_abc_method",
    "params_all",
    "n_runs",
    "ndays"
  )
)

cat("Submitting SLURM jobs for method: lstm\n")
ans2 <- Slurm_lapply(
  X = 1:n_param_sets,
  FUN = run_lstm_method,
  job_name = "bias_coverage_lstm",
  njobs = 100,
  overwrite = TRUE,
  plan = "submit",
  sbatch_opt = list(
    partition = "vegayon-shared-np",
    account = "vegayon-np",
    time = "02:00:00",
    `mem-per-cpu` = "8G",
    `cpus-per-task` = 1
  ),
  export = c(
    "run_parameter_set_method",
    "run_lstm_method",
    "params_all",
    "n_runs",
    "ndays"
  )
)

# Collect results for each method separately
cat("Collecting results from SLURM jobs for method: true\n")
results_true <- Slurm_collect(ans)
method_results_true <- bind_rows(results_true)
# Remove failed results
results_true_clean <- results_true[!sapply(results_true, function(x) inherits(x$res, "error"))]

# Then bind
method_results_true <- dplyr::bind_rows(results_true_clean)

saveRDS(method_results_true, "results_true2.rds")
cat("Saved results for method: true to results_true.rds\n")

cat("Collecting results from SLURM jobs for method: abc\n")
results_abc <- Slurm_collect(ans1)
results_abc_clean <- results_abc[!sapply(results_abc, function(x) inherits(x$res, "error"))]

method_results_abc <- bind_rows(results_abc_clean)
saveRDS(method_results_abc, "results_abc2.rds")
cat("Saved results for method: abc to results_abc.rds\n")

cat("Collecting results from SLURM jobs for method: lstm\n")
results_lstm <- Slurm_collect(ans2)
results_lstm_clean <- results_lstm[!sapply(results_lstm, function(x) inherits(x$res, "error"))]

method_results_lstm <- bind_rows(results_lstm_clean)
saveRDS(method_results_lstm, "results_lstm2.rds")
cat("Saved results for method: lstm to results_lstm.rds\n")

# Combine all method results
cat("Combining all method results...\n")
method_results_lstm=read_rds("results_lstm2.rds")
method_results_abc=read_rds("results_abc2.rds")
method_results_true=read_rds("results_true2.rds")
all_simulation_results <- bind_rows(method_results_true, method_results_abc, method_results_lstm)

# Save combined results
saveRDS(all_simulation_results, "all_simulation_results2.rds")
cat("Saved combined results to all_simulation_results.rds\n")

# Summary statistics
cat("Analysis complete! Total results:", nrow(all_simulation_results), "rows\n")
cat("Unique parameter sets:", length(unique(all_simulation_results$sim_id)), "\n")
cat("Methods processed:", paste(unique(all_simulation_results$param_type), collapse = ", "), "\n")
# ========== BIAS AND COVERAGE CALCULATION ==========
# STEP 1: Filter sim_ids present in all three methods


wide_results <- all_simulation_results %>%
  select(sim_id, date, param_type, mean_infected, q025, q975) %>%
  pivot_wider(
    names_from = param_type,
    values_from = c(mean_infected, q025, q975),
    names_sep = "_"
  )

# Calculate bias and coverage
bias_coverage_results <- wide_results %>%
  mutate(
    # Bias calculations
    bias_abc = mean_infected_abc - mean_infected_true,
    bias_lstm = mean_infected_lstm - mean_infected_true,
    abs_bias_abc = abs(bias_abc),
    abs_bias_lstm = abs(bias_lstm),
    
    # Coverage: Is true value within ABC/LSTM confidence intervals?
    coverage_abc = (mean_infected_true >= q025_abc) & (mean_infected_true <= q975_abc),
    coverage_lstm = (mean_infected_true >= q025_lstm) & (mean_infected_true <= q975_lstm)
  )

# ========== SUMMARY STATISTICS ==========

# Overall statistics across all parameter sets and days
overall_stats <- bias_coverage_results %>%
  summarize(
    # ABC Statistics
    mean_bias_abc = mean(bias_abc, na.rm = TRUE),
    median_bias_abc = median(bias_abc, na.rm = TRUE),
    rmse_abc = sqrt(mean(bias_abc^2, na.rm = TRUE)),
    mean_abs_bias_abc = mean(abs_bias_abc, na.rm = TRUE),
    
    # LSTM Statistics  
    mean_bias_lstm = mean(bias_lstm, na.rm = TRUE),
    median_bias_lstm = median(bias_lstm, na.rm = TRUE),
    rmse_lstm = sqrt(mean(bias_lstm^2, na.rm = TRUE)),
    mean_abs_bias_lstm = mean(abs_bias_lstm, na.rm = TRUE),
    
    # Coverage Statistics (will be calculated separately using quantiles)
    n_observations = n()
  )

# Calculate overall coverage statistics using the corrected method
overall_coverage_stats <- all_simulation_results %>%
  select(sim_id, date, param_type, mean_infected) %>%
  pivot_wider(names_from = param_type, values_from = mean_infected) %>%
  summarize(
    # Calculate overall 95% CI from ABC and LSTM predictions
    abc_q025 = quantile(abc, 0.025, na.rm = TRUE),
    abc_q975 = quantile(abc, 0.975, na.rm = TRUE),
    lstm_q025 = quantile(lstm, 0.025, na.rm = TRUE),
    lstm_q975 = quantile(lstm, 0.975, na.rm = TRUE),
    
    # Overall coverage
    coverage_abc = mean((true >= abc_q025) & (true <= abc_q975), na.rm = TRUE) * 100,
    coverage_lstm = mean((true >= lstm_q025) & (true <= lstm_q975), na.rm = TRUE) * 100
  )

# Combine overall stats
overall_stats <- bind_cols(overall_stats, overall_coverage_stats %>% select(coverage_abc, coverage_lstm))
print(overall_stats)

# Calculate coverage properly using quantiles across parameter sets
coverage_by_day <- all_simulation_results %>%
  select(sim_id, date, param_type, mean_infected) %>%
  pivot_wider(names_from = param_type, values_from = mean_infected) %>%
  group_by(date) %>%
  summarize(
    # Calculate 95% CI from ABC and LSTM predictions across all parameter sets
    abc_q025 = quantile(abc, 0.025, na.rm = TRUE),
    abc_q975 = quantile(abc, 0.975, na.rm = TRUE),
    lstm_q025 = quantile(lstm, 0.025, na.rm = TRUE),
    lstm_q975 = quantile(lstm, 0.975, na.rm = TRUE),
    
    # Coverage: Check if true values fall within these CIs
    coverage_abc = mean((true >= abc_q025) & (true <= abc_q975), na.rm = TRUE) * 100,
    coverage_lstm = mean((true >= lstm_q025) & (true <= lstm_q975), na.rm = TRUE) * 100,
    
    .groups = "drop"
  )

# Daily statistics (averaged across all 1000 parameter sets)
daily_stats <- bias_coverage_results %>%
  group_by(date) %>%
  summarize(
    # Mean bias by day
    mean_bias_abc = mean(bias_abc, na.rm = TRUE),
    mean_bias_lstm = mean(bias_lstm, na.rm = TRUE),
    
    # Absolute bias by day
    mean_abs_bias_abc = mean(abs_bias_abc, na.rm = TRUE),
    mean_abs_bias_lstm = mean(abs_bias_lstm, na.rm = TRUE),
    
    # Standard errors for plotting
    se_bias_abc = sd(bias_abc, na.rm = TRUE) / sqrt(n()),
    se_bias_lstm = sd(bias_lstm, na.rm = TRUE) / sqrt(n()),
    
    .groups = "drop"
  ) %>%
  # Add the corrected coverage calculations
  left_join(coverage_by_day, by = "date")

# ========== VISUALIZATIONS ==========

# Plot 1: Mean Bias over Time (Days 0-60)
bias_long <- daily_stats %>%
  select(date, mean_bias_abc, mean_bias_lstm, se_bias_abc, se_bias_lstm) %>%
  pivot_longer(
    cols = c(mean_bias_abc, mean_bias_lstm),
    names_to = "method", 
    values_to = "bias"
  ) %>%
  mutate(
    method = case_when(
      method == "mean_bias_abc" ~ "ABC",
      method == "mean_bias_lstm" ~ "LSTM"
    ),
    se = case_when(
      method == "ABC" ~ se_bias_abc,
      method == "LSTM" ~ se_bias_lstm
    )
  ) %>%
  select(-se_bias_abc, -se_bias_lstm)

p1_bias <- ggplot(bias_long, aes(x = date, y = bias, color = method)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_ribbon(aes(ymin = bias - 1.96*se, ymax = bias + 1.96*se, fill = method), 
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "Mean Bias Over Time (1000 Parameter Sets)",
    subtitle = "Bias = Predicted - True, with 95% confidence bands",
    x = "Day (0-60)", 
    y = "Mean Bias (Infected Count)",
    color = "Method", 
    fill = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB")) +
  scale_fill_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB"))

# Plot 2: Coverage over Time
coverage_long <- daily_stats %>%
  select(date, coverage_abc, coverage_lstm) %>%
  pivot_longer(
    cols = c(coverage_abc, coverage_lstm),
    names_to = "method", 
    values_to = "coverage"
  ) %>%
  mutate(method = case_when(
    method == "coverage_abc" ~ "ABC",
    method == "coverage_lstm" ~ "LSTM"
  ))

p2_coverage <- ggplot(coverage_long, aes(x = date, y = coverage, color = method)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray50", alpha = 0.7) +
  labs(
    title = "Coverage Over Time (95% Confidence Intervals)",
    subtitle = "Percentage of true values within predicted confidence intervals",
    x = "Day (0-60)", 
    y = "Coverage (%)",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB")) +
  ylim(0, 100)

# Plot 3: Absolute Bias over Time
abs_bias_long <- daily_stats %>%
  select(date, mean_abs_bias_abc, mean_abs_bias_lstm) %>%
  pivot_longer(
    cols = c(mean_abs_bias_abc, mean_abs_bias_lstm),
    names_to = "method", 
    values_to = "abs_bias"
  ) %>%
  mutate(method = case_when(
    method == "mean_abs_bias_abc" ~ "ABC",
    method == "mean_abs_bias_lstm" ~ "LSTM"
  ))

p3_abs_bias <- ggplot(abs_bias_long, aes(x = date, y = abs_bias, color = method)) +
  geom_line(size = 1.2, alpha = 0.8) +
  labs(
    title = "Mean Absolute Bias Over Time",
    x = "Day (0-60)", 
    y = "Mean Absolute Bias",
    color = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB"))

# Plot 4: Summary Statistics
summary_table <- data.frame(
  Method = c("ABC", "LSTM"),
  Mean_Bias = c(overall_stats$mean_bias_abc, overall_stats$mean_bias_lstm),
  RMSE = c(overall_stats$rmse_abc, overall_stats$rmse_lstm),
  Mean_Abs_Bias = c(overall_stats$mean_abs_bias_abc, overall_stats$mean_abs_bias_lstm),
  Coverage = c(overall_stats$coverage_abc, overall_stats$coverage_lstm)
)

# Summary metrics plot
summary_long <- summary_table %>%
  select(-Coverage) %>%
  pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value")

p4_summary <- ggplot(summary_long, aes(x = Method, y = Value, fill = Method)) +
  geom_col(alpha = 0.7) +
  facet_wrap(~Metric, scales = "free_y") +
  labs(
    title = "Summary Statistics Comparison",
    x = "Method", y = "Value"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB"))

# Coverage summary plot
p5_coverage_summary <- ggplot(summary_table, aes(x = Method, y = Coverage, fill = Method)) +
  geom_col(alpha = 0.7) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray50") +
  labs(
    title = "Overall Coverage Comparison",
    subtitle = "Target: 95%",
    x = "Method", y = "Coverage (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB")) +
  ylim(0, 100)

# Combine all plots
final_plot <- (p1_bias | p2_coverage)/(p4_summary| p5_coverage_summary)  +
  plot_annotation(
    title = "Bias and Coverage Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

final_plot <- ( p4_summary| p5_coverage_summary)  +
  plot_annotation(
    title = "Bias and Coverage Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Display the final plot
print(final_plot)

ggsave("final_plot.png", plot = final_plot, width = 10, height = 6, dpi = 300)
# Print detailed summary
cat("\n=== DETAILED SUMMARY TABLE ===\n")
print(summary_table)


###
# Calculate bias and coverage with PERCENTAGE BIAS
# Calculate bias as (Pred - True) / Pred
bias_coverage_results <- wide_results %>%
  mutate(
    # Percentage bias based on predicted values
    pct_bias_abc = ((mean_infected_abc - mean_infected_true) / mean_infected_abc) * 100,
    pct_bias_lstm = ((mean_infected_lstm - mean_infected_true) / mean_infected_lstm) * 100,
    
    # Coverage
    coverage_abc = (mean_infected_true >= q025_abc) & (mean_infected_true <= q975_abc),
    coverage_lstm = (mean_infected_true >= q025_lstm) & (mean_infected_true <= q975_lstm)
  )

# Daily statistics: mean and SE
daily_stats_pct <- bias_coverage_results %>%
  group_by(date) %>%
  summarize(
    mean_pct_bias_abc = mean(pct_bias_abc, na.rm = TRUE),
    mean_pct_bias_lstm = mean(pct_bias_lstm, na.rm = TRUE),
    se_pct_bias_abc = sd(pct_bias_abc, na.rm = TRUE) / sqrt(n()),
    se_pct_bias_lstm = sd(pct_bias_lstm, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Prepare long format for plotting
pct_bias_long <- daily_stats_pct %>%
  pivot_longer(
    cols = c(mean_pct_bias_abc, mean_pct_bias_lstm),
    names_to = "method", 
    values_to = "pct_bias"
  ) %>%
  mutate(
    method = case_when(
      method == "mean_pct_bias_abc" ~ "ABC",
      method == "mean_pct_bias_lstm" ~ "LSTM"
    ),
    se = case_when(
      method == "ABC" ~ se_pct_bias_abc,
      method == "LSTM" ~ se_pct_bias_lstm
    )
  )

# Plot
p_pct_bias <- ggplot(pct_bias_long, aes(x = date, y = pct_bias, color = method)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_ribbon(aes(ymin = pct_bias - 1.96 * se, ymax = pct_bias + 1.96 * se, fill = method), 
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "Percentage Bias Over Time",
    subtitle = "",
    x = "Day", 
    y = "Mean Percentage Bias (%)",
    color = "Method", 
    fill = "Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB")) +
  scale_fill_manual(values = c("ABC" = "#E74C3C", "LSTM" = "#3498DB"))

p_pct_bias

final_plot <- (p_pct_bias | p2_coverage)/(p4_summary| p5_coverage_summary)  +
  plot_annotation(
    title = "Bias and Coverage Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )
print(final_plot)
