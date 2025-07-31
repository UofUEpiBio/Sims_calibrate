# ───────────────────────────────────

library(tidyverse)
library(epiworldR)
library(glue)
source("~/Desktop/Sims_calibrate/01a-bilstm.R")
# 1) load your ABC predictions

# 2) constants

filename <- paste0("~/Desktop/Sims_calibrate/abc_parameters.csv")
ndays=60
nsims    <- 1
nthreads <- 10
abc_pred <- read_csv(filename)
# 3) Instead of using rowwise(), let's use a for loop for better control
lstm_results <- vector("list", nrow(abc_pred))

for (i in 1:nrow(abc_pred)) {
  sim_id <- abc_pred$sim_id[i]
  
  cat(sprintf("Processing sim_id %d (%d/%d)\n", sim_id, i, nrow(abc_pred)))
  
  tryCatch({
    # a) run a single-rep true SIR
    m0 <- ModelSIRCONN(
      name              = paste0("true_sim_", sim_id),
      n                 = (abc_pred$true_n[i]),
      prevalence        = abc_pred$true_preval[i],
      contact_rate      = abc_pred$true_crate[i],
      transmission_rate = abc_pred$true_ptran[i],
      recovery_rate     = abc_pred$true_recov[i]
    )
    run(m0,ndays = 60)
    incidence <- plot_incidence(m0,plot=FALSE)
    all_counts=incidence[,1]
    
    ts <- all_counts[-1]           # drop day0
    ts <- ts[0:ndays]              # force exactly 60 values
    
    if (length(ts) != ndays) {
      stop(glue("After dropping day0 and subsetting, length(ts) = {length(ts)}, not {ndays}!"))
    }
    joblib <- import("joblib")
    library(reticulate)
    
  
    # c) predict parameters
    lstm_out <- predict_with_bilstm(ts, abc_pred$true_n[i], abc_pred$true_recov[i])
    
    # Extract values from LSTM output (ptran, crate, R0)
    ptran_lstm <- lstm_out[1]
    crate_lstm <- lstm_out[2]  
    R0_lstm <- lstm_out[3]
    
    # Use fixed values and calculate contact rate from R0
    # Formula: R0 = (contact_rate * transmission_rate) / recovery_rate
    # So: contact_rate = (R0 * recovery_rate) / transmission_rate
    prevalence_lstm <- abc_pred$true_preval[i]  # Use true prevalence
    contact_rate_lstm <- (R0_lstm * abc_pred$true_recov[i]) / ptran_lstm  # Calculate from R0
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
    recovery_rate_lstm = recovery_rate_true                # dummy tag to keep true_n
  ) %>%
  pivot_longer(
    cols = -c(sim_id, true_n),  # exclude true_n from pivot
    names_to = c(".value", "param_type"),
    names_pattern = "(.*)_(true|abc|lstm)"
  ) 

write.csv(params_all, "params_all.csv", row.names = FALSE)
# 5) run 20 sims × 60 days for each param_type
cat("\nRunning forward simulations...\n")
# Modify your pmap_dfr function to include R0 as a parameter
# 5) run a single simulation for each param_type and extract incidence
cat("\nRunning simulations...\n")
# 5) run a single simulation for each param_type and extract incidence
cat("\nRunning simulations...\n")
all_sims <- params_all %>%
  mutate(run_id = row_number()) %>%
  pmap_dfr(function(run_id, sim_id, true_n,param_type,
                    prevalence, contact_rate,
                    transmission_rate, recovery_rate, R0) {  # Add R0 parameter here
    
    # Create model with the given parameters
    calibrated_model <- ModelSIRCONN(
      name              = paste0("run_", run_id),
      n                 = true_n,
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
ggsave("~/Desktop/Sims_calibrate/figures/comparison_plot.png", width = 10, height = 8, dpi = 300)
cat("\nPlot saved as 'comparison_plot.png'\n")

# Save the results
results_summary <- list(
  abc_predictions = abc_pred,
  lstm_predictions = lstm_preds,
  combined_parameters = params_all,
  simulation_results = all_sims,
  summary_statistics = summary_df
)

saveRDS(results_summary, "~/Desktop/Sims_calibrate/figures/comparison_results.rds")
cat("Results saved as 'comparison_results.rds'\n")

# ── end script ───────────────────────────────────────────────────────────────

#biases:
# Calculate bias between ABC/LSTM parameters and true parameters
cat("\nCalculating parameter biases...\n")




