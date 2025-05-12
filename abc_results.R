library(epiworldR)
library(tidyverse)
library(ggplot2)


abc_predicted_parameters_200 <- readRDS("~/Desktop/Sims_calibrate/abc_predicted_parameters_200_fixedrecn.rds")
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
      n                 = 1000,
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
