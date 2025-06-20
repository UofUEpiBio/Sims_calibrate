# ───────────────────────────────────
# Complete code to compare epidemic curves for ONE parameter set
# ───────────────────────────────────

# 1) libraries
library(tidyverse)
library(glue)
source("~/Desktop/Sims_calibrate/loading_lstm_model_003.R")

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
    ts        <- incidence[, 1]        # drop day0
             # exactly 60 days
    
    # LSTM prediction
    joblib        <- import("joblib")
    library(reticulate)
    lstm_out      <- predict_with_bilstm(ts, abc_pred$true_n[i], abc_pred$true_recov[i])
    ptran_lstm    <- lstm_out[1]
    crate_lstm    <- lstm_out[2]
    R0_lstm       <- lstm_out[3]
    prevalence_lstm    <- abc_pred$true_preval[i]
    contact_rate_lstm  <- R0_lstm*abc_pred$true_recov[i]/ptran_lstm
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

#params_all_filtered%>% filter(sim_id==17)
chosen_sim_id <- 4

single_param_set <- params_all %>%
  filter(sim_id == chosen_sim_id)

cat(sprintf("Comparing epidemic curves for sim_id: %d\n", chosen_sim_id))
print(single_param_set)

# ──────────────────────────────────────────────────────────────────
# 5) Run MULTIPLE simulations for the single parameter set with CI
# ──────────────────────────────────────────────────────────────────
# 
# n_reps <- 10
# cat(sprintf("\nRunning %d replications for single parameter set...\n", n_reps))
# 
# single_sims <- single_param_set %>%
#   mutate(run_id = row_number()) %>%
#   pmap_dfr(function(run_id, sim_id, true_n, param_type,
#                     prevalence, contact_rate,
#                     transmission_rate, recovery_rate, R0) {
#     cat(sprintf("Running %s parameters (rep %d/%d)...\n", param_type, run_id, n_reps))
#     calibrated_model <- ModelSIRCONN(
#       name              = paste0("single_", param_type, "_", sim_id),
#       n                 = true_n,
#       prevalence        = prevalence,
#       contact_rate      = contact_rate,
#       transmission_rate = transmission_rate,
#       recovery_rate     = recovery_rate
#     )
#     saver <- make_saver("total_hist", "reproductive")
#     run_multiple(calibrated_model, ndays = ndays, nsims = n_reps, saver = saver, nthreads = 4)
#     incidence_data <- plot_incidence(calibrated_model, plot = FALSE)
#     
#     map_dfr(seq_len(min(n_reps, ncol(incidence_data))), function(rep_num) {
#       ts <- incidence_data[-1, rep_num]
#       tibble(
#         time             = seq_along(ts),
#         I                = ts,
#         param_type       = param_type,
#         sim_id           = sim_id,
#         rep              = rep_num,
#         prevalence       = prevalence,
#         contact_rate     = contact_rate,
#         transmission_rate= transmission_rate,
#         recovery_rate    = recovery_rate,
#         R0               = R0
#       )
#     })
#   })
# 
# single_param_set$run_id <- seq_len(nrow(single_param_set))
# 
# single_sims <- tibble()
# 
# for (i in seq_len(nrow(single_param_set))) {
#   row <- single_param_set[i, ]
#   
#   cat(sprintf("Running %s parameters (rep %d/%d)...\n", row$param_type, row$run_id, n_reps))
#   
#   calibrated_model <- ModelSIRCONN(
#     name              = paste0("single_", row$param_type, "_", row$sim_id),
#     n                 = row$true_n,
#     prevalence        = row$prevalence,
#     contact_rate      = row$contact_rate,
#     transmission_rate = row$transmission_rate,
#     recovery_rate     = row$recovery_rate
#   )
#   
#   saver <- make_saver("total_hist")
#   run_multiple(calibrated_model, ndays = ndays, nsims = 10, saver = saver, nthreads = 4)
#   get_hist_total(calibrated_model)
#   raw_hist <- run_multiple_get_results(calibrated_model)
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   raw_hist <- get_hist_total(calibrated_model)
#   df=raw_hist$total_hist
#   # Filter Infected counts only
#   infected_ts <- df %>%
#     filter(state == "Infected") %>%
#     select(sim_num, date, counts)
#   avg_infected <- df %>%
#     filter(state == "Infected") %>%
#     group_by(date) %>%
#     summarize(mean_infected = mean(counts), .groups = "drop")
#   
#   # Plot
#   ggplot(avg_infected, aes(x = date, y = mean_infected)) +
#     geom_line(size = 1) +
#     labs(
#       title = "Average Infected Over Time",
#       x = "Day",
#       y = "Average Number of Infected Individuals"
#     ) +
#     theme_minimal()
#   
#   incidence_data <- plot_incidence(calibrated_model, plot = FALSE,all = TRUE)
#   
#   for (rep_num in seq_len(min(n_reps, ncol(incidence_data)))) {
#     ts <- incidence_data[-1, rep_num]
#     sim_data <- tibble(
#       time              = seq_along(ts),
#       I                 = ts,
#       param_type        = row$param_type,
#       sim_id            = row$sim_id,
#       rep               = rep_num,
#       prevalence        = row$prevalence,
#       contact_rate      = row$contact_rate,
#       transmission_rate = row$transmission_rate,
#       recovery_rate     = row$recovery_rate,
#       R0                = row$R0
#     )
#     single_sims <- bind_rows(single_sims, sim_data)
#   }
# }
# 
# # ──────────────────────────────────────────────────────────────────
# # 6) Calculate summary statistics with confidence intervals
# # ──────────────────────────────────────────────────────────────────
# 
# summary_stats <- single_sims %>%
#   group_by(param_type, time) %>%
#   summarise(
#     median_I = median(I),
#     lower95  = quantile(I, 0.025),
#     upper95  = quantile(I, 0.975),
#     .groups = "drop"
#   )
# 
# # ──────────────────────────────────────────────────────────────────
# # 7) Create comparison plot with a single grey 95% ribbon
# # ──────────────────────────────────────────────────────────────────
# 
# comparison_plot <- ggplot(summary_stats, aes(x = time, color = param_type)) +
#   
#   # one light-grey 95% CI under everything
#   geom_ribbon(
#     aes(ymin = lower95, ymax = upper95),
#     fill  = "grey80",
#     alpha = 0.3
#   ) +
#   
#   # median trajectories
#   geom_line(aes(y = median_I), size = 1.2) +
#   
#   # labels
#   labs(
#     x        = "Day",
#     y        = "Active Infected (I)",
#     title    = sprintf("Epidemic Curves (95%% CI) – Sim ID %d", chosen_sim_id),
#     subtitle = sprintf("True vs ABC vs LSTM (%d reps each)", n_reps),
#     color    = "Parameter Source"
#   ) +
#   
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom") +
#   
#   scale_color_manual(
#     values = c(
#       "true" = "#2E8B57",
#       "abc"  = "#FF6347",
#       "lstm" = "#4169E1"
#     ),
#     labels = c("True", "ABC", "LSTM")
#   )
# 
# print(comparison_plot)
# 
# # ──────────────────────────────────────────────────────────────────
# # 8) (Optional) Spaghetti plot
# # ──────────────────────────────────────────────────────────────────
# 
# spaghetti_plot <- ggplot(single_sims, aes(x = time, y = I, color = param_type)) +
#   geom_line(aes(group = interaction(param_type, rep)), alpha = 0.3, size = 0.5) +
#   geom_line(data = summary_stats, aes(y = median_I), size = 1.5) +
#   labs(
#     x     = "Day",
#     y     = "Active Infected (I)",
#     title = sprintf("Individual Simulation Trajectories - Sim ID %d", chosen_sim_id),
#     subtitle = sprintf("All %d reps with median overlay", n_reps),
#     color = "Parameter Source"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "bottom") +
#   scale_color_manual(
#     values = c("true" = "#2E8B57", "abc" = "#FF6347", "lstm" = "#4169E1"),
#     labels = c("True", "ABC", "LSTM")
#   )
# 
# print(spaghetti_plot)
# 
# # ──────────────────────────────────────────────────────────────────
# # 9) Save figures
# # ──────────────────────────────────────────────────────────────────
# 
# ggsave(
#   sprintf("~/Desktop/Sims_calibrate/figures/single_comparison_CI_sim_%d.png", chosen_sim_id),
#   comparison_plot, width = 12, height = 8, dpi = 300
# )
# 
# ggsave(
#   sprintf("~/Desktop/Sims_calibrate/figures/single_comparison_spaghetti_sim_%d.png", chosen_sim_id),
#   spaghetti_plot, width = 12, height = 8, dpi = 300
# )
# 

n_reps=100
single_param_set$run_id <- seq_len(nrow(single_param_set))
single_param_set[2,7]=single_param_set[1,7]
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
  
  # Run simulations
  saver <- make_saver("total_hist")
  run_multiple(calibrated_model, ndays = ndays, nsims = n_reps, saver = saver, nthreads = 8)
  
  # Extract results
  raw_hist <- get_hist_total(calibrated_model)
  df <- raw_hist
  
  # Average across simulations
  avg_states <- df %>%
    group_by(date, state) %>%
    summarize(mean_count = mean(counts), .groups = "drop") %>%
    mutate(
      param_type = row$param_type,
      sim_id     = row$sim_id
    )
  
  # Collect
  all_avg_series <- bind_rows(all_avg_series, avg_states)
}

# === Generate Plots ===

# Infected
p_infected <- all_avg_series %>%
  filter(state == "Infected") %>%
  ggplot(aes(x = date, y = mean_count, color = param_type)) +
  geom_line(size = 1) +
  labs(
    title = "Infected Over Time by Method",
    x = "Day", y = "Average Infected", color = "Method"
  ) +
  theme_minimal()

# Susceptible
p_susceptible <- all_avg_series %>%
  filter(state == "Susceptible") %>%
  ggplot(aes(x = date, y = mean_count, color = param_type)) +
  geom_line(size = 1) +
  labs(
    title = "Susceptible Over Time by Method",
    x = "Day", y = "Average Susceptible", color = "Method"
  ) +
  theme_minimal()

# Recovered
p_recovered <- all_avg_series %>%
  filter(state == "Recovered") %>%
  ggplot(aes(x = date, y = mean_count, color = param_type)) +
  geom_line(size = 1) +
  labs(
    title = "Recovered Over Time by Method",
    x = "Day", y = "Average Recovered", color = "Method"
  ) +
  theme_minimal()
library(patchwork)
# Combine plots
combined_plot <- (p_infected / p_susceptible / p_recovered) +
  plot_annotation(title = "SIR Dynamics Averaged Over Simulations by Method")

# Show the final plot
print(combined_plot)

# Bias and Coverage Analysis for Infected Counts: ABC vs LSTM
# Using params_all data structure (1000 parameter sets × 3 methods)

library(tidyverse)
library(patchwork)
library(ggplot2)

# Assuming you have params_all loaded
# params_all should have 3000 rows (1000 sim_id × 3 param_type)

# Parameters
n_runs <- 100         # Runs per parameter set per method
ndays <- 60          # Days to simulate
n_param_sets <- max(params_all$sim_id)  # Should be 1000

cat("Starting bias and coverage analysis for", n_param_sets, "parameter sets...\n")
cat("Each parameter set will be run", n_runs, "times for", ndays, "days\n")

# Container to collect all results
all_simulation_results <- tibble()

# Loop through all parameter sets
for (sim_idx in 1:n_param_sets) {
  
  if (sim_idx %% 100 == 0) {
    cat("Processing parameter set", sim_idx, "of", n_param_sets, "...\n")
  }
  
  # Get current parameter set (all 3 methods)
  current_sim_params <- params_all %>% 
    filter(sim_id == sim_idx)
  
  # Run simulations for each method (true, abc, lstm)
  for (method in c("true", "abc", "lstm")) {
    
    # Get parameters for current method
    method_params <- current_sim_params %>% 
      filter(param_type == method)
    
    if (nrow(method_params) == 0) {
      cat("Warning: No parameters found for sim_id", sim_idx, "method", method, "\n")
      next
    }
    
    # Create model with current parameters
    model <- ModelSIRCONN(
      name = paste0("sim_", sim_idx, "_", method),
      n = method_params$true_n,
      prevalence = method_params$prevalence,
      contact_rate = method_params$contact_rate,
      transmission_rate = method_params$transmission_rate,
      recovery_rate = method_params$recovery_rate
    )
    
    # Set seed for reproducibility
    set.seed(123 + sim_idx * 1000 + match(method, c("true", "abc", "lstm")))
    
    # Run multiple simulations
    saver <- make_saver("total_hist")
    run_multiple(model, ndays = ndays, nsims = n_runs, saver = saver, nthreads = 4)
    
    # Extract results and focus on Infected state
    raw_hist <- get_hist_total(model)
    
    # Calculate statistics across the 100 runs for each day
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
        param_type = method
      )
    
    # Store results
    all_simulation_results <- bind_rows(all_simulation_results, infected_stats)
  }
}

cat("Simulations complete! Calculating bias and coverage...\n")

# ========== BIAS AND COVERAGE CALCULATION ==========

# Reshape data to have true, abc, lstm in separate columns
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
    
    # Coverage Statistics
    coverage_abc = mean(coverage_abc, na.rm = TRUE) * 100,  # Convert to percentage
    coverage_lstm = mean(coverage_lstm, na.rm = TRUE) * 100
  )

print("=== OVERALL BIAS AND COVERAGE STATISTICS ===")
print(overall_stats)

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
    
    # Coverage by day
    coverage_abc = mean(coverage_abc, na.rm = TRUE) * 100,
    coverage_lstm = mean(coverage_lstm, na.rm = TRUE) * 100,
    
    # Standard errors for plotting
    se_bias_abc = sd(bias_abc, na.rm = TRUE) / sqrt(n()),
    se_bias_lstm = sd(bias_lstm, na.rm = TRUE) / sqrt(n()),
    
    .groups = "drop"
  )

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
final_plot <- (p1_bias | p2_coverage) / (p3_abs_bias | p5_coverage_summary) / p4_summary +
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

cat("\n=== KEY FINDINGS ===\n")
cat("Sample Size: 1000 parameter sets × 100 runs × 60 days =", 
    format(1000 * 100 * 60, big.mark = ","), "data points per method\n\n")

cat("BIAS RESULTS:\n")
cat("ABC  - Mean Bias:", sprintf("%.3f", overall_stats$mean_bias_abc), 
    "| RMSE:", sprintf("%.3f", overall_stats$rmse_abc), "\n")
cat("LSTM - Mean Bias:", sprintf("%.3f", overall_stats$mean_bias_lstm), 
    "| RMSE:", sprintf("%.3f", overall_stats$rmse_lstm), "\n\n")

cat("COVERAGE RESULTS (Target: 95%):\n")
cat("ABC  Coverage:", sprintf("%.1f%%", overall_stats$coverage_abc), "\n")
cat("LSTM Coverage:", sprintf("%.1f%%", overall_stats$coverage_lstm), "\n\n")

if (overall_stats$rmse_abc < overall_stats$rmse_lstm) {
  cat("🏆 ABC shows lower RMSE (better accuracy)\n")
} else {
  cat("🏆 LSTM shows lower RMSE (better accuracy)\n")
}

if (abs(overall_stats$coverage_abc - 95) < abs(overall_stats$coverage_lstm - 95)) {
  cat("🎯 ABC shows coverage closer to target 95%\n")
} else {
  cat("🎯 LSTM shows coverage closer to target 95%\n")
}

# Save results
write_csv(bias_coverage_results, "bias_coverage_detailed_results.csv")
write_csv(daily_stats, "daily_bias_coverage_summary.csv")
write_csv(summary_table, "method_comparison_summary.csv")

cat("\n✅ Analysis complete! Results saved to CSV files.\n")
############################################################################



write_csv(all_simulation_results, "all_simulation_results.csv")



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
final_plot <- (p1_bias | p2_coverage) / (p3_abs_bias | p5_coverage_summary) / p4_summary +
  plot_annotation(
    title = "Bias and Coverage Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

# Display the final plot
print(final_plot)

# Print detailed summary
cat("\n=== DETAILED SUMMARY TABLE ===\n")
print(summary_table)

cat("\n=== KEY FINDINGS ===\n")
cat("Sample Size: 1000 parameter sets × 100 runs × 60 days =", 
    format(1000 * 100 * 60, big.mark = ","), "data points per method\n\n")

cat("BIAS RESULTS:\n")
cat("ABC  - Mean Bias:", sprintf("%.3f", overall_stats$mean_bias_abc), 
    "| RMSE:", sprintf("%.3f", overall_stats$rmse_abc), "\n")
cat("LSTM - Mean Bias:", sprintf("%.3f", overall_stats$mean_bias_lstm), 
    "| RMSE:", sprintf("%.3f", overall_stats$rmse_lstm), "\n\n")

cat("COVERAGE RESULTS (Target: 95%):\n")
cat("ABC  Coverage:", sprintf("%.1f%%", overall_stats$coverage_abc), "\n")
cat("LSTM Coverage:", sprintf("%.1f%%", overall_stats$coverage_lstm), "\n\n")

if (overall_stats$rmse_abc < overall_stats$rmse_lstm) {
  cat("🏆 ABC shows lower RMSE (better accuracy)\n")
} else {
  cat("🏆 LSTM shows lower RMSE (better accuracy)\n")
}

if (abs(overall_stats$coverage_abc - 95) < abs(overall_stats$coverage_lstm - 95)) {
  cat("🎯 ABC shows coverage closer to target 95%\n")
} else {
  cat("🎯 LSTM shows coverage closer to target 95%\n")
}

# Save results
write_csv(bias_coverage_results, "bias_coverage_detailed_results.csv")
write_csv(daily_stats, "daily_bias_coverage_summary.csv")
write_csv(summary_table, "method_comparison_summary.csv")

cat("\n✅ Analysis complete! Results saved to CSV files.\n")
ggsave("bias_coverage_analysis_plot.png", final_plot, 
       width = 16, height = 12, dpi = 300, bg = "white")

ggsave("bias_over_time.png", p1_bias, 
       width = 12, height = 6, dpi = 300, bg = "white")

ggsave("coverage_over_time.png", p2_coverage, 
       width = 12, height = 6, dpi = 300, bg = "white")

ggsave("absolute_bias_over_time.png", p3_abs_bias, 
       width = 12, height = 6, dpi = 300, bg = "white")

ggsave("summary_statistics.png", p4_summary, 
       width = 10, height = 6, dpi = 300, bg = "white")

ggsave("coverage_summary.png", p5_coverage_summary, 
       width = 8, height = 6, dpi = 300, bg = "white")
ggsave("final_plot_correct.png", final_plot, 
       width = 10, height = 6, dpi = 300, bg = "white")

