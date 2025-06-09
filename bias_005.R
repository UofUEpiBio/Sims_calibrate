library(tidyverse)
df=params_all
# Step 1: Pivot to wide format
df_wide <- df %>%
  pivot_wider(
    id_cols = sim_id,
    names_from = param_type,
    values_from = c(contact_rate, transmission_rate)
  ) %>%
  rename(
    true_crate  = contact_rate_true,
    abc_crate   = contact_rate_abc,
    lstm_crate  = contact_rate_lstm,
    true_ptran  = transmission_rate_true,
    abc_ptran   = transmission_rate_abc,
    lstm_ptran  = transmission_rate_lstm
  )

# Step 2: Compute biases
bias_df <- df_wide %>%
  mutate(
    bias_crate_abc  = abc_crate  - true_crate,
    bias_crate_lstm = lstm_crate - true_crate,
    bias_ptran_abc  = abc_ptran  - true_ptran,
    bias_ptran_lstm = lstm_ptran - true_ptran
  ) %>%
  select(sim_id, starts_with("bias_")) %>%
  arrange(sim_id)

# Step 3: Convert to long format
bias_long <- bias_df %>%
  pivot_longer(
    cols = -sim_id,
    names_to = c("param", "method"),
    names_pattern = "bias_(.*)_(abc|lstm)"
  )

# Step 4: Line plot — trends over sim_id
ggplot(bias_long, aes(x = sim_id, y = value, color = method, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 1.2) +
  facet_wrap(~param, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(crate = "Contact Rate Bias",
                                      ptran = "Transmission Rate Bias"))) +
  scale_color_manual(values = c("abc" = "darkred", "lstm" = "orange")) +
  labs(
    title = "Bias of ABC and LSTM Over Simulations",
    x = "Simulation ID",
    y = "Bias (Estimate − True)",
    color = "Method"
  ) +
  theme_minimal(base_size = 14)


bias_long_filtered <- bias_long %>%
  filter(sim_id <= 200)
ggplot(bias_long_filtered, aes(x = sim_id, y = value, color = method, group = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 1.2) +
  facet_wrap(~param, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(crate = "Contact Rate Bias",
                                      ptran = "Transmission Rate Bias"))) +
  scale_color_manual(values = c("abc" = "darkred", "lstm" = "orange")) +
  labs(
    title = "Bias of ABC and LSTM Estimates (Simulations 1–200)",
    x = "Simulation ID",
    y = "Bias (Estimate − True)",
    color = "Method"
  ) +
  theme_minimal(base_size = 14)

