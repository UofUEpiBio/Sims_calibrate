# ── SCRIPT: compare_true_abc_lstm_fixed.R ───────────────────────────────────

library(tidyverse)
library(epiworldR)
library(glue)

# 1) load your ABC predictions
abc_pred <- readRDS("~/Desktop/Sims_calibrate/abc_predicted_parameters_200_fixedrecn.rds")

# 2) constants
n        <- 5000
recov    <- 0.1
ndays    <- 60
nsims    <- 40
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
    sv <- make_saver("total_hist")
    run_multiple(m0, ndays = ndays, nsims = 1, saver = sv, nthreads = 1)
    df0 <- run_multiple_get_results(m0)$total_hist
    
    # b) pull out the counts vector, drop the first (day 0), then take days 1:60
    all_counts <- df0 %>%
      filter(state == "Infected") %>%
      arrange(date) %>%
      pull(counts)
    
    ts <- all_counts[-1]           # drop day0
    ts <- ts[1:ndays]              # force exactly 60 values
    
    if (length(ts) != ndays) {
      stop(glue("After dropping day0 and subsetting, length(ts) = {length(ts)}, not {ndays}!"))
    }
    
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
         prevalence_abc         = abc_preval,
         contact_rate_abc       = abc_crate,
         transmission_rate_abc  = abc_ptran,
         recovery_rate_abc      = abc_recov) %>%
  inner_join(lstm_preds, by = "sim_id") %>%
  # Add recovery rate for LSTM (should be 0.1)
  mutate(recovery_rate_lstm = 0.1) %>%
  pivot_longer(
    cols = -sim_id,
    names_to     = c(".value","param_type"),
    names_pattern= "(.*)_(true|abc|lstm)"
  )

# 5) run 20 sims × 60 days for each param_type
cat("\nRunning forward simulations...\n")
all_sims <- params_all %>%
  mutate(run_id = row_number()) %>%
  pmap_dfr(function(run_id, sim_id, param_type,
                    prevalence, contact_rate,
                    transmission_rate, recovery_rate) {
    m <- ModelSIRCONN(
      name              = paste0("run_", run_id),
      n                 = n,
      prevalence        = prevalence,
      contact_rate      = contact_rate,
      transmission_rate = transmission_rate,
      recovery_rate     = recovery_rate
    )
    sv <- make_saver("total_hist")
    run_multiple(m,
                 ndays    = ndays,
                 nsims    = nsims,
                 saver    = sv,
                 nthreads = nthreads)
    res <- run_multiple_get_results(m)$total_hist
    
    res %>%
      filter(state == "Infected") %>%
      rename(time = date, rep = sim_num, I = counts) %>%
      select(time, rep, I) %>%
      mutate(param_type = param_type,
             sim_id     = sim_id)
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




# ── SCRIPT: compare_true_abc_lstm_fixed_incidence.R ───────────────────────

library(tidyverse)
library(epiworldR)
library(glue)

# 1) load your ABC predictions
abc_pred <- readRDS("~/Desktop/Sims_calibrate/abc_predicted_parameters_200_fixedrecn.rds")

# 2) constants
n        <- 5000
recov    <- 0.1
ndays    <- 60
nsims    <- 20
nthreads <- 10

# 3) Instead of using rowwise(), let's use a for loop for better control
lstm_results <- vector("list", nrow(abc_pred))

for (i in 1:nrow(abc_pred)) {
  sim_id <- abc_pred$sim_id[i]
  
  cat(sprintf("Processing sim_id %d (%d/%d)\n", sim_id, i, nrow(abc_pred)))
  
  tryCatch({
    # a) run a single-rep true SIR to get incidence
    model_sir <- ModelSIRCONN(
      name              = paste0("true_sim_", sim_id),
      n                 = n,
      prevalence        = abc_pred$true_preval[i],
      contact_rate      = abc_pred$true_crate[i],
      transmission_rate = abc_pred$true_ptran[i],
      recovery_rate     = recov
    )
    
    run(model_sir, ndays = ndays)
    incidence <- plot_incidence(model_sir, plot = FALSE)
    infected_original <- incidence[, 1]
    
    # b) drop the first day (day 0) and take days 1:60
    ts <- infected_original[-1]    # drop day0
    ts <- ts[1:ndays]              # force exactly 60 values
    
    if (length(ts) != ndays) {
      stop(glue("After dropping day0 and subsetting, length(ts) = {length(ts)}, not {ndays}!"))
    }
    
    # c) predict parameters using LSTM with incidence data
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
         prevalence_abc         = abc_preval,
         contact_rate_abc       = abc_crate,
         transmission_rate_abc  = abc_ptran,
         recovery_rate_abc      = abc_recov) %>%
  inner_join(lstm_preds, by = "sim_id") %>%
  # Add recovery rate for LSTM (should be 0.1)
  mutate(recovery_rate_lstm = 0.1) %>%
  pivot_longer(
    cols = -sim_id,
    names_to     = c(".value","param_type"),
    names_pattern= "(.*)_(true|abc|lstm)"
  )

# 5) run 20 sims × 60 days for each param_type
cat("\nRunning forward simulations...\n")
all_sims <- params_all %>%
  mutate(run_id = row_number()) %>%
  pmap_dfr(function(run_id, sim_id, param_type,
                    prevalence, contact_rate,
                    transmission_rate, recovery_rate) {
    m <- ModelSIRCONN(
      name              = paste0("run_", run_id),
      n                 = n,
      prevalence        = prevalence,
      contact_rate      = contact_rate,
      transmission_rate = transmission_rate,
      recovery_rate     = recovery_rate
    )
    sv <- make_saver("total_hist")
    run_multiple(m,
                 ndays    = ndays,
                 nsims    = nsims,
                 saver    = sv,
                 nthreads = nthreads)
    res <- run_multiple_get_results(m)$total_hist
    
    res %>%
      filter(state == "Infected") %>%
      rename(time = date, rep = sim_num, I = counts) %>%
      select(time, rep, I) %>%
      mutate(param_type = param_type,
             sim_id     = sim_id)
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
    title  = "True vs. ABC vs. LSTM‐Estimated SIRCONN Trajectories (Using Incidence)",
    color  = "Parameter set",
    fill   = "Parameter set"
  ) +
  theme_minimal(base_size = 14)

# Save the plot
ggsave("comparison_plot_incidence.png", width = 10, height = 8, dpi = 300)
cat("\nPlot saved as 'comparison_plot_incidence.png'\n")

# Save the results
results_summary <- list(
  abc_predictions = abc_pred,
  lstm_predictions = lstm_preds,
  combined_parameters = params_all,
  simulation_results = all_sims,
  summary_statistics = summary_df
)

saveRDS(results_summary, "comparison_results_incidence.rds")
cat("Results saved as 'comparison_results_incidence.rds'\n")

# ── end script ───────────────────────────────────────────────────────────────