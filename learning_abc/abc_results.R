library(epiworldR)
library(tidyverse)
library(ggplot2)
summary(abc_predicted_parameters_200)

abc_predicted_parameters_200 <- readRDS("~/Desktop/Sims_calibrate/abc_predicted_parameters_1000_fixedrecn.rds")
# 1) Stack true vs abc into one tibble (as before)
params_all <- abc_predicted_parameters_200 %>%
  select(sim_id,
         prevalence_true    = true_preval,
         contact_rate_true  = true_crate,
         transmission_rate_true = true_ptran,
         recovery_rate_true     = true_recov,
         prevalence_abc     = abc_preval,
         contact_rate_abc   = abc_crate,
         transmission_rate_abc  = abc_ptran,
         recovery_rate_abc      = abc_recov) %>%
  pivot_longer(
    cols = -sim_id,
    names_to = c(".value","param_type"),
    names_pattern = "(.*)_(true|abc)"
  )

# 2) For each parameter set, run 20 sims × 60 days, saving total_hist
all_sims <- params_all %>%
  mutate(run_id = row_number()) %>%
  pmap_dfr(function(run_id, sim_id, param_type,
                    prevalence, contact_rate,
                    transmission_rate, recovery_rate) {
    
    # Build model
    m <- ModelSIRCONN(
      name              = "COVID-19",
      n                 = 5000,
      prevalence        = prevalence,
      contact_rate      = contact_rate,
      transmission_rate = transmission_rate,
      recovery_rate     = recovery_rate
    )
    
    # Saver for full state counts
    saver <- make_saver("total_hist")  # only "total_hist" is supported :contentReference[oaicite:1]{index=1}
    
    # Run
    run_multiple(m,
                 ndays    = 60,
                 nsims    = 20,
                 saver    = saver,
                 nthreads = 10)
    
    # Pull results
    res <- run_multiple_get_results(m)  # returns a named list with your saver output :contentReference[oaicite:2]{index=2}
    
    # Extract *only* the Infected counts
    res$total_hist %>%
      filter(state == "Infected") %>%           # keep just the Infected rows
      rename(
        time = date,                            # day → time
        rep  = sim_num,                         # replicate → rep
        I    = counts                           # counts → I
      ) %>%
      select(time, rep, I) %>%                  # we just need these three
      mutate(
        param_type = param_type,
        sim_id     = sim_id
      )
  })

# 3) Summarize: median + 95% CI of I over sims
summary_df <- all_sims %>%
  group_by(param_type, time) %>%
  summarize(
    median_I = median(I),
    lower95  = quantile(I, 0.025),
    upper95  = quantile(I, 0.975),
    .groups  = "drop"
  )
saveRDS(summary_df, file = "summary_df.rds")

# Later, to read it back in:
summary_df_loaded <- readRDS("summary_df.rds")
summary_df=summary_df_loaded
# 4) Plot true vs abc curves with ribbons
ggplot(summary_df,
       aes(x = time, y = median_I,
           color = param_type, fill = param_type)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95),
              alpha = 0.2, color = NA) +
  labs(
    x     = "Day",
    y     = "Active Infected (I)",
    title = "True vs. ABC‐Estimated SIRCONN Trajectories",
    color = "Parameter set",
    fill  = "Parameter set"
  ) +
  theme_minimal(base_size = 14)


####################
#──────────────────────────────────────────────────────────────────────────────
# Complete end-to-end R script:
#   • Reads your 200×11 tibble of true vs ABC parameters
#   • Simulates a single “true” trajectory for each sim_id
#   • Calibrates via BiLSTM to get a third (“calibrated”) parameter set
#   • Runs 20 stochastic sims for each of the three methods
#   • Summarizes median ± 95% CI
#   • Plots all three on one figure
#──────────────────────────────────────────────────────────────────────────────

# 0) load libraries -----------------------------------------------------------
library(epiworldR)
library(tidyverse)
library(reticulate)
library(ggplot2)

# 0a) source your BiLSTM‐predictor (adjust path if needed)
#    This makes predict_with_bilstm() available in R.
#source_python("predict_with_bilstm.py")

#──────────────────────────────────────────────────────────────────────────────
# 1) helper to run one replicate and return the I(t) vector (drop day 0)
#──────────────────────────────────────────────────────────────────────────────
simulate_single <- function(preval, crate, ptran, recov, n = 5000, ndays = 60) {
  m <- ModelSIRCONN(
    name              = "COVID-19",
    n                 = n,
    prevalence        = preval,
    contact_rate      = crate,
    transmission_rate = ptran,
    recovery_rate     = recov
  )
  run_multiple(m,
               ndays    = ndays,
               nsims    = 1,
               saver    = make_saver("total_hist"),
               nthreads = 10)
  df <- run_multiple_get_results(m)$total_hist
  infected <- df %>%
    filter(state == "Infected") %>%
    arrange(date) %>%
    pull(counts)
  infected[-1]
}

#──────────────────────────────────────────────────────────────────────────────
# 2) read in your 200×11 tibble
#──────────────────────────────────────────────────────────────────────────────
params_raw <- readRDS("~/Desktop/Sims_calibrate/abc_predicted_parameters_200_fixedrecn.rds")

#──────────────────────────────────────────────────────────────────────────────
# 3) pivot true vs abc into long form: sim_id × type ∈ {true,abc}
#──────────────────────────────────────────────────────────────────────────────
params_long <- params_raw %>%
  select(sim_id, matches("^(true|abc)_(preval|crate|ptran|recov)$")) %>%
  pivot_longer(
    cols      = -sim_id,
    names_to  = c("type","param"),
    names_sep = "_"
  ) %>%
  pivot_wider(
    names_from  = param,
    values_from = value
  )
# now: 400 rows (200 sims × 2 types) with columns:
#   sim_id, type ∈ {true, abc}, preval, crate, ptran, recov

#──────────────────────────────────────────────────────────────────────────────
# 4) simulate one‐rep for each “true” to get infected_true curves
#──────────────────────────────────────────────────────────────────────────────
true_curves <- params_long %>%
  filter(type == "true") %>%
  mutate(
    infected = pmap(
      list(preval, crate, ptran, recov),
      simulate_single
    )
  )

#──────────────────────────────────────────────────────────────────────────────
# 5) calibrate each curve via BiLSTM → “calibrated” params
#──────────────────────────────────────────────────────────────────────────────
calibrated_params <- true_curves %>%
  mutate(
    out       = map(infected, ~ predict_with_bilstm(.x, n = 5000, recov = recov)),
    ptran_cal = map_dbl(out, 1),
    crate_cal = map_dbl(out, 2),
    r0_cal    = map_dbl(out, 3),
    # recompute contact_rate so that R0 = (crate_cal × ptran_cal) / recov
    crate_cal = (r0_cal * recov) / ptran_cal,
    type      = "calibrated"
  ) %>%
  transmute(
    sim_id,
    type,
    preval = preval,
    crate  = crate_cal,
    ptran  = ptran_cal,
    recov
  )

#──────────────────────────────────────────────────────────────────────────────
# 6) bind together the three method‐sets: true, abc, calibrated
#──────────────────────────────────────────────────────────────────────────────
params_all <- bind_rows(
  params_long   %>% transmute(sim_id, type, preval, crate, ptran, recov),
  calibrated_params
)

#──────────────────────────────────────────────────────────────────────────────
# 7) run 20 sims for each parameter set, extract infected counts
#──────────────────────────────────────────────────────────────────────────────
all_sims <- params_all %>%
  arrange(sim_id, type) %>%
  mutate(uid = row_number()) %>%
  pmap_dfr(function(sim_id, type, preval, crate, ptran, recov, uid) {
    m <- ModelSIRCONN(
      name              = paste0("sim", sim_id, "-", type),
      n                 = 5000,
      prevalence        = preval,
      contact_rate      = crate,
      transmission_rate = ptran,
      recovery_rate     = recov
    )
    run_multiple(m,
                 ndays    = 60,
                 nsims    = 20,
                 saver    = make_saver("total_hist"),
                 nthreads = 4)
    res <- run_multiple_get_results(m)$total_hist
    res %>%
      filter(state == "Infected") %>%
      rename(time = date, rep = sim_num, I = counts) %>%
      select(time, rep, I) %>%
      mutate(sim_id = sim_id, param_type = type)
  })

#──────────────────────────────────────────────────────────────────────────────
# 8) summarize median + 95% CI over all sims for each method
#──────────────────────────────────────────────────────────────────────────────
summary_df <- all_sims %>%
  group_by(param_type, time) %>%
  summarize(
    median_I = median(I),
    lower95  = quantile(I, 0.025),
    upper95  = quantile(I, 0.975),
    .groups  = "drop"
  )

#──────────────────────────────────────────────────────────────────────────────
# 9) plot true vs ABC vs calibrated
#──────────────────────────────────────────────────────────────────────────────
ggplot(summary_df, aes(x = time, y = median_I,
                       color = param_type, fill = param_type)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.2) +
  labs(
    x     = "Day",
    y     = "Active Infected (I)",
    title = "True vs ABC vs BiLSTM‐Calibrated SIRCONN Trajectories",
    color = "Method",
    fill  = "Method"
  ) +
  theme_minimal(base_size = 14)

#──────────────────────────────────────────────────────────────────────────────
# End of script
#──────────────────────────────────────────────────────────────────────────────
