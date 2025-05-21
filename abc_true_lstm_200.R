# ── SCRIPT: compare_true_abc_lstm_fixed.R ───────────────────────────────────

library(tidyverse)
library(epiworldR)
library(glue)
source("~/Desktop/Sims_calibrate/loading_lstm_model.R")
# 1) load your ABC predictions
abc_pred <- readRDS("~/Desktop/Sims_calibrate/predicted_parameters/abc_predicted_parameters_200_fixedrecn.rds")

# 2) constants
n        <- 5000
recov    <- 0.1
ndays    <- 60
nsims    <- 1
nthreads <- 10

# 3) Instead of using rowwise(), let's use a for loop for better control
lstm_results <- vector("list", nrow(abc_pred))

for (i in 1:nrow(abc_pred)) {
  sim_id <- abc_pred$sim_id[i]
  
  cat(sprintf("Processing sim_id %d (%d/%d)\n", sim_id, i, nrow(abc_pred)))
  
  tryCatch({
    # a) run a single-rep true SIR
    m0 <- ModelSIRCONN(
      name              = paste0("true_sim_", sim_id),
      n                 = n,
      prevalence        = abc_pred$true_preval[i],
      contact_rate      = abc_pred$true_crate[i],
      transmission_rate = abc_pred$true_ptran[i],
      recovery_rate     = recov
    )
    run(m0,ndays = 60)
    incidence <- plot_incidence(m0,plot=FALSE)
    all_counts=incidence[,1]
    
    ts <- all_counts[-1]           # drop day0
    ts <- ts[1:ndays]              # force exactly 60 values
    
    if (length(ts) != ndays) {
      stop(glue("After dropping day0 and subsetting, length(ts) = {length(ts)}, not {ndays}!"))
    }
    joblib <- import("joblib")
    library(reticulate)
    
  
    # c) predict parameters
    lstm_out <- predict_with_bilstm(ts, n, recov)
    
    # Extract values from LSTM output (ptran, crate, R0)
    ptran_lstm <- lstm_out[1]
    crate_lstm <- lstm_out[2]  
    R0_lstm <- lstm_out[3]
    
    # Use fixed values and calculate contact rate from R0
    # Formula: R0 = (contact_rate * transmission_rate) / recovery_rate
    # So: contact_rate = (R0 * recovery_rate) / transmission_rate
    prevalence_lstm <- abc_pred$true_preval[i]  # Use true prevalence
    contact_rate_lstm <- (R0_lstm * recov) / ptran_lstm  # Calculate from R0
    transmission_rate_lstm <- ptran_lstm  # Use directly from LSTM
    
    # Store results
    lstm_results[[i]] <- data.frame(
      sim_id = sim_id,
      prevalence_lstm = prevalence_lstm,
      contact_rate_lstm = contact_rate_lstm,
      transmission_rate_lstm = transmission_rate_lstm,
      R0_lstm=R0_lstm,
      stringsAsFactors = FALSE
    )
    
    cat(sprintf("  ✓ Success for sim_id %d\n", sim_id))
    
  }, error = function(e) {
    cat(sprintf("  ✗ Error for sim_id %d: %s\n", sim_id, e$message))
    # Store NA results for failed predictions
    lstm_results[[i]] <- data.frame(
      sim_id = sim_id,
      prevalence_lstm = NA,
      contact_rate_lstm = NA,
      transmission_rate_lstm = NA,
      stringsAsFactors = FALSE
    )
  })
}

# Combine results
lstm_preds <- bind_rows(lstm_results)

# Remove failed predictions
lstm_preds <- lstm_preds %>%
  filter(!is.na(prevalence_lstm))

cat(sprintf("\nSummary:\n"))
cat(sprintf("Total simulations: %d\n", nrow(abc_pred)))
cat(sprintf("Successful LSTM predictions: %d\n", nrow(lstm_preds)))
cat(sprintf("Success rate: %.2f%%\n", 100 * nrow(lstm_preds) / nrow(abc_pred)))

# 4) stack true / abc / lstm
params_all <- abc_pred %>%
  select(sim_id,
         prevalence_true        = true_preval,
         contact_rate_true      = true_crate,
         transmission_rate_true = true_ptran,
         recovery_rate_true     = true_recov,
         R0_true=true_R0,
         R0_abc=abc_R0,
         prevalence_abc         = abc_preval,
         contact_rate_abc       = abc_crate,
         transmission_rate_abc  = abc_ptran,
         recovery_rate_abc      = abc_recov) %>%
  inner_join(lstm_preds, by = "sim_id") %>% 
  # Add recovery rate for LSTM (should be 0.1)
  mutate(recovery_rate_lstm = 0.1)%>%mutate(recovery_rate_abc=0.1) %>%
  pivot_longer(
    cols = -sim_id,
    names_to     = c(".value","param_type"),
    names_pattern= "(.*)_(true|abc|lstm)"
  )

# 5) run 20 sims × 60 days for each param_type
cat("\nRunning forward simulations...\n")
# Modify your pmap_dfr function to include R0 as a parameter
# 5) run a single simulation for each param_type and extract incidence
cat("\nRunning simulations...\n")
# 5) run a single simulation for each param_type and extract incidence
cat("\nRunning simulations...\n")
all_sims <- params_all %>%
  mutate(run_id = row_number()) %>%
  pmap_dfr(function(run_id, sim_id, param_type,
                    prevalence, contact_rate,
                    transmission_rate, recovery_rate, R0) {  # Add R0 parameter here
    
    # Create model with the given parameters
    calibrated_model <- ModelSIRCONN(
      name              = paste0("run_", run_id),
      n                 = n,
      prevalence        = prevalence,
      contact_rate      = contact_rate,
      transmission_rate = transmission_rate,
      recovery_rate     = recovery_rate
    )
    
    # Run a single simulation
    run(calibrated_model, ndays = 60)
    
    # Extract incidence data
    incidence_calib <- plot_incidence(calibrated_model, plot = FALSE)
    infected_calibrated <- incidence_calib[, 1]
    infected_calibrated <- infected_calibrated[-1]  # remove day 0
    
    # Create a data frame with time and infected count
    tibble(
      time = 1:length(infected_calibrated),
      I = infected_calibrated,
      rep = 1,  # Since we're only doing one rep per parameter set
      param_type = param_type,
      sim_id = sim_id
    )
  })
# 6) summarize
summary_df <- all_sims %>%
  group_by(param_type, time) %>%
  summarise(
    median_I = median(I),
    lower95  = quantile(I, 0.025),
    upper95  = quantile(I, 0.975),
    .groups  = "drop"
  )
View(summary_df)
# 7) plot all three
library(ggplot2)
ggplot(summary_df,
       aes(x = time, y = median_I,
           color = param_type, fill = param_type)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95),
              alpha = 0.2, color = NA) +
  labs(
    x      = "Day",
    y      = "Active Infected (I)",
    title  = "True vs. ABC vs. LSTM‐Estimated SIRCONN Trajectories",
    color  = "Parameter set",
    fill   = "Parameter set"
  ) +
  theme_minimal(base_size = 14)

# Save the plot
ggsave("comparison_plot.png", width = 10, height = 8, dpi = 300)
cat("\nPlot saved as 'comparison_plot.png'\n")

# Save the results
results_summary <- list(
  abc_predictions = abc_pred,
  lstm_predictions = lstm_preds,
  combined_parameters = params_all,
  simulation_results = all_sims,
  summary_statistics = summary_df
)

saveRDS(results_summary, "comparison_results.rds")
cat("Results saved as 'comparison_results.rds'\n")

# ── end script ───────────────────────────────────────────────────────────────

#biases:
# Calculate bias between ABC/LSTM parameters and true parameters
cat("\nCalculating parameter biases...\n")

# First, ensure we have R0 for all parameter sets if not already calculated
params_all <- params_all %>%
  mutate(
    # Calculate R0 as contact_rate * transmission_rate / recovery_rate if not already there
    R0 = contact_rate * transmission_rate / recovery_rate
  )

# Reshape the data to wide format to make comparisons easier
params_wide <- params_all %>%
  select(sim_id, param_type, prevalence, contact_rate, transmission_rate, recovery_rate, R0) %>%
  pivot_wider(
    id_cols = sim_id,
    names_from = param_type,
    values_from = c(prevalence, contact_rate, transmission_rate, recovery_rate, R0)
  )

# Calculate absolute and relative biases for each parameter
param_biases <- params_wide %>%
  mutate(
    # Absolute biases
    abs_bias_prevalence_abc = prevalence_abc - prevalence_true,
    abs_bias_contact_rate_abc = contact_rate_abc - contact_rate_true,
    abs_bias_transmission_rate_abc = transmission_rate_abc - transmission_rate_true,
    abs_bias_recovery_rate_abc = recovery_rate_abc - recovery_rate_true,
    abs_bias_R0_abc = R0_abc - R0_true,
    
    # Relative biases (as proportions, not percentages)
    rel_bias_prevalence_abc = (prevalence_abc - prevalence_true) / prevalence_true,
    rel_bias_contact_rate_abc = (contact_rate_abc - contact_rate_true) / contact_rate_true,
    rel_bias_transmission_rate_abc = (transmission_rate_abc - transmission_rate_true) / transmission_rate_true,
    rel_bias_recovery_rate_abc = (recovery_rate_abc - recovery_rate_true) / recovery_rate_true,
    rel_bias_R0_abc = (R0_abc - R0_true) / R0_true,
    
    # Also calculate for LSTM if needed
    abs_bias_prevalence_lstm = prevalence_lstm - prevalence_true,
    abs_bias_contact_rate_lstm = contact_rate_lstm - contact_rate_true,
    abs_bias_transmission_rate_lstm = transmission_rate_lstm - transmission_rate_true,
    abs_bias_recovery_rate_lstm = recovery_rate_lstm - recovery_rate_true,
    abs_bias_R0_lstm = R0_lstm - R0_true,
    
    rel_bias_prevalence_lstm = (prevalence_lstm - prevalence_true) / prevalence_true,
    rel_bias_contact_rate_lstm = (contact_rate_lstm - contact_rate_true) / contact_rate_true,
    rel_bias_transmission_rate_lstm = (transmission_rate_lstm - transmission_rate_true) / transmission_rate_true,
    rel_bias_recovery_rate_lstm = (recovery_rate_lstm - recovery_rate_true) / recovery_rate_true,
    rel_bias_R0_lstm = (R0_lstm - R0_true) / R0_true
  )

# Summarize the biases
bias_summary <- param_biases %>%
  summarise(
    # ABC Absolute bias statistics
    mean_abs_bias_prevalence_abc = mean(abs_bias_prevalence_abc),
    mean_abs_bias_contact_rate_abc = mean(abs_bias_contact_rate_abc),
    mean_abs_bias_transmission_rate_abc = mean(abs_bias_transmission_rate_abc),
    mean_abs_bias_recovery_rate_abc = mean(abs_bias_recovery_rate_abc),
    mean_abs_bias_R0_abc = mean(abs_bias_R0_abc),
    
    median_abs_bias_prevalence_abc = median(abs_bias_prevalence_abc),
    median_abs_bias_contact_rate_abc = median(abs_bias_contact_rate_abc),
    median_abs_bias_transmission_rate_abc = median(abs_bias_transmission_rate_abc),
    median_abs_bias_recovery_rate_abc = median(abs_bias_recovery_rate_abc),
    median_abs_bias_R0_abc = median(abs_bias_R0_abc),
    
    # ABC Relative bias statistics (as proportions)
    mean_rel_bias_prevalence_abc = mean(rel_bias_prevalence_abc),
    mean_rel_bias_contact_rate_abc = mean(rel_bias_contact_rate_abc),
    mean_rel_bias_transmission_rate_abc = mean(rel_bias_transmission_rate_abc),
    mean_rel_bias_recovery_rate_abc = mean(rel_bias_recovery_rate_abc),
    mean_rel_bias_R0_abc = mean(rel_bias_R0_abc),
    
    median_rel_bias_prevalence_abc = median(rel_bias_prevalence_abc),
    median_rel_bias_contact_rate_abc = median(rel_bias_contact_rate_abc),
    median_rel_bias_transmission_rate_abc = median(rel_bias_transmission_rate_abc),
    median_rel_bias_recovery_rate_abc = median(rel_bias_recovery_rate_abc),
    median_rel_bias_R0_abc = median(rel_bias_R0_abc),
    
    # LSTM statistics - if needed
    mean_abs_bias_prevalence_lstm = mean(abs_bias_prevalence_lstm),
    mean_abs_bias_contact_rate_lstm = mean(abs_bias_contact_rate_lstm),
    mean_abs_bias_transmission_rate_lstm = mean(abs_bias_transmission_rate_lstm),
    mean_abs_bias_recovery_rate_lstm = mean(abs_bias_recovery_rate_lstm),
    mean_abs_bias_R0_lstm = mean(abs_bias_R0_lstm),
    
    median_abs_bias_prevalence_lstm = median(abs_bias_prevalence_lstm),
    median_abs_bias_contact_rate_lstm = median(abs_bias_contact_rate_lstm),
    median_abs_bias_transmission_rate_lstm = median(abs_bias_transmission_rate_lstm),
    median_abs_bias_recovery_rate_lstm = median(abs_bias_recovery_rate_lstm),
    median_abs_bias_R0_lstm = median(abs_bias_R0_lstm),
    
    mean_rel_bias_prevalence_lstm = mean(rel_bias_prevalence_lstm),
    mean_rel_bias_contact_rate_lstm = mean(rel_bias_contact_rate_lstm),
    mean_rel_bias_transmission_rate_lstm = mean(rel_bias_transmission_rate_lstm),
    mean_rel_bias_recovery_rate_lstm = mean(rel_bias_recovery_rate_lstm),
    mean_rel_bias_R0_lstm = mean(rel_bias_R0_lstm),
    
    median_rel_bias_prevalence_lstm = median(rel_bias_prevalence_lstm),
    median_rel_bias_contact_rate_lstm = median(rel_bias_contact_rate_lstm),
    median_rel_bias_transmission_rate_lstm = median(rel_bias_transmission_rate_lstm),
    median_rel_bias_recovery_rate_lstm = median(rel_bias_recovery_rate_lstm),
    median_rel_bias_R0_lstm = median(rel_bias_R0_lstm)
  )

# Print the summary
print(bias_summary)

# Create a more readable table for the biases
bias_table <- data.frame(
  Parameter = c("Prevalence", "Contact Rate", "Transmission Rate", "Recovery Rate", "R0"),
  
  ABC_Mean_Abs_Bias = c(
    bias_summary$mean_abs_bias_prevalence_abc,
    bias_summary$mean_abs_bias_contact_rate_abc,
    bias_summary$mean_abs_bias_transmission_rate_abc,
    bias_summary$mean_abs_bias_recovery_rate_abc,
    bias_summary$mean_abs_bias_R0_abc
  ),
  
  LSTM_Mean_Abs_Bias = c(
    bias_summary$mean_abs_bias_prevalence_lstm,
    bias_summary$mean_abs_bias_contact_rate_lstm,
    bias_summary$mean_abs_bias_transmission_rate_lstm,
    bias_summary$mean_abs_bias_recovery_rate_lstm,
    bias_summary$mean_abs_bias_R0_lstm
  )
)

# Print the formatted table
print(bias_table)
