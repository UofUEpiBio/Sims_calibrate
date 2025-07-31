library(tidyverse)
library(ggplot2)
library(patchwork)
source("~/Desktop/Sims_calibrate/abc_true_lstm_004.R")
# Assuming you have params_all from your main script
df <- read_csv("params_all.csv")

# Step 1: Pivot to wide format and calculate R0
df_wide <- df %>%
  pivot_wider(
    id_cols = sim_id,
    names_from = param_type,
    values_from = c(contact_rate, transmission_rate, recovery_rate)
  ) %>%
  rename(
    true_crate  = contact_rate_true,
    abc_crate   = contact_rate_abc,
    lstm_crate  = contact_rate_lstm,
    true_ptran  = transmission_rate_true,
    abc_ptran   = transmission_rate_abc,
    lstm_ptran  = transmission_rate_lstm,
    true_recov  = recovery_rate_true,
    abc_recov   = recovery_rate_abc,
    lstm_recov  = recovery_rate_lstm
  ) %>%
  mutate(
    true_R0  = (true_crate * true_ptran) / true_recov,
    abc_R0   = (abc_crate * abc_ptran) / abc_recov,
    lstm_R0  = (lstm_crate * lstm_ptran) / lstm_recov
  )

# Step 2: Compute biases
bias_df <- df_wide %>%
  mutate(
    bias_crate_abc  = abc_crate  - true_crate,
    bias_crate_lstm = lstm_crate - true_crate,
    bias_ptran_abc  = abc_ptran  - true_ptran,
    bias_ptran_lstm = lstm_ptran - true_ptran,
    bias_R0_abc     = abc_R0     - true_R0,
    bias_R0_lstm    = lstm_R0    - true_R0
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

# Step 4: Calculate MAE
mae_stats <- bias_long %>%
  group_by(param, method) %>%
  summarise(
    mae = mean(abs(value), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mae_label = sprintf("MAE: %.4f", mae)
  )

# Labels
param_labels <- c(
  "crate" = "Contact Rate",
  "ptran" = "Transmission Rate", 
  "R0" = "Basic Reproduction Number (R₀)"
)


# Bias summary table
bias_summary_table <- bias_long %>%
  group_by(param, method) %>%
  summarise(
    MAE = mean(abs(value), na.rm = TRUE),
    RMSE = sqrt(mean(value^2, na.rm = TRUE)),
    Mean_Bias = mean(value, na.rm = TRUE),
    Median_Bias = median(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(param, method)

write_csv(bias_summary_table, "~/Desktop/Sims_calibrate/figures/bias_summary_table.csv")
saveRDS(bias_summary_table, "~/Desktop/Sims_calibrate/figures/bias_summary_table.rds")

# Box plots


# 1. Compute mean bias per param × method
mean_bias_df <- bias_long %>%
  group_by(param, method) %>%
  summarise(
    mean_bias = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

param_labels <- c(
  "crate" = "Contact Rate",
  "ptran" = "Transmission Rate", 
  "R0" = "Basic Reproduction Number (R₀)"
)

# ✅ Updated method colors
method_colors <- c("abc" = "#d62728", "lstm" = "#2ca02c")  # red and green

# 2. Helper to build one plot at a time
make_param_plot <- function(param_name, param_label) {
  df_sub  <- bias_long %>% filter(param == param_name)
  mean_df <- mean_bias_df %>% filter(param == param_name)
  
  # determine where to put the text (5% below the min bias)
  y_min  <- min(df_sub$value, na.rm = TRUE)
  y_max  <- max(df_sub$value, na.rm = TRUE)
  y_text <- y_min - 0.05 * (y_max - y_min)
  
  ggplot(df_sub, aes(x = method, y = value, fill = method)) +
    geom_hline(yintercept = 0,
               linetype = "dashed", color = "gray50") +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_violin(alpha = 0.3) +
    
    # mean‐bias annotation
    geom_text(
      data = mean_df,
      aes(x = method, y = y_text, 
          label = sprintf("Mean: %.4f", mean_bias)),
      inherit.aes = FALSE,
      vjust = 1, size = 4, fontface = "bold"
    ) +
    
    scale_fill_manual(
      values = method_colors,
      labels = c("abc" = "ABC", "lstm" = "LSTM"),
      name = "Method"
    ) +
    scale_x_discrete(labels = c("abc" = "ABC", "lstm" = "LSTM")) +
    
    labs(
      title = param_label,
      x     = "Method",
      y     = "Bias"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

# 3. Build each plot
p_crate <- make_param_plot("crate", "Contact Rate")
p_ptran <- make_param_plot("ptran", "Transmission Rate")
make_param_plot <- function(param_name, param_label) {
  df_sub <- bias_long %>%
    filter(param == param_name) %>%
    filter(!(param == "R0" & value > 15))  # remove R0 values > 15
  
  mean_df <- mean_bias_df %>%
    filter(param == param_name)
  
  # determine where to put the text (5% below the min bias, still safe for y-limits)
  y_min  <- min(df_sub$value, na.rm = TRUE)
  y_max  <- max(df_sub$value, na.rm = TRUE)
  y_text <- max(-10, y_min - 0.05 * (y_max - y_min))  # keep above -10
  
  ggplot(df_sub, aes(x = method, y = value, fill = method)) +
    geom_hline(yintercept = 0,
               linetype = "dashed", color = "gray50") +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_violin(alpha = 0.3) +
    
    # mean-bias annotation
    geom_text(
      data = mean_df,
      aes(x = method, y = y_text, 
          label = sprintf("Mean: %.4f", mean_bias)),
      inherit.aes = FALSE,
      vjust = 1, size = 4, fontface = "bold"
    ) +
    
    scale_fill_manual(
      values = method_colors,
      labels = c("abc" = "ABC", "lstm" = "LSTM"),
      name = "Method"
    ) +
    scale_x_discrete(labels = c("abc" = "ABC", "lstm" = "LSTM")) +
    
    labs(
      title = param_label,
      x     = "Method",
      y     = "Bias"
    ) +
    
    coord_cartesian(ylim = c(-10, 10)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
}

p_R0 <- make_param_plot("R0", "Basic Reproductive number")

library(patchwork)

# Combine the three plots side by side (or use / for stacking)
combined_plot <- p_crate + p_ptran + p_R0 +
  plot_layout(ncol = 3) +
  plot_annotation(title = "Comparison of Parameter Recovery for ABC vs LSTM")

# Display the plot
print(combined_plot)

library(ggplot2)
library(dplyr)
library(patchwork)


# 4. Save or display
ggsave("bias_box_crate.png", p_crate, width = 4, height = 6, dpi = 300)
ggsave("bias_box_ptran.png", p_ptran, width = 4, height = 6, dpi = 300)

# to view in the console:
p_crate
p_ptran

make_param_plot <- function(param_name, param_label, log_y = FALSE) {
  
  # 1) Subset & truncate R₀ if requested
  df_sub <- bias_long %>% 
    filter(param == param_name)
  if (param_name == "R0") {
    df_sub <- df_sub %>% filter(abs(value) <= 15)
  }
  
  # 2) Compute MAE and Mean Bias per method
  stats_df <- df_sub %>%
    group_by(method) %>%
    summarise(
      mae        = mean(abs(value), na.rm = TRUE),
      mean_bias  = mean(value, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      label = sprintf("MAE: %.4f\nMean: %.4f", mean_bias)
    )
  
  # 3) Figure out where to place text (5% below the min)
  y_min   <- min(df_sub$value, na.rm = TRUE)
  y_max   <- max(df_sub$value, na.rm = TRUE)
  y_text  <- y_min - 0.05 * (y_max - y_min)
  
  # 4) Build the plot
  p <- ggplot(df_sub, aes(x = method, y = value, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    geom_violin(alpha = 0.3) +
    
    # annotate MAE and Mean Bias underneath
    geom_text(
      data = stats_df,
      aes(x = method, y = y_text, label = label),
      inherit.aes = FALSE,
      vjust = 1, size = 4, fontface = "bold",
      show.legend = FALSE
    ) +
    
    scale_fill_manual(
      values = method_colors,
      labels = c("abc" = "ABC", "lstm" = "LSTM"),
      name   = "Method"
    ) +
    scale_x_discrete(labels = c("abc" = "ABC", "lstm" = "LSTM")) +
    
    labs(
      title    = param_label,
      subtitle = if (param_name == "R0") "" else NULL,
      x        = "Method",
      y        = if (log_y) "Bias (log₁₀)" else "Bias"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "gray40"),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  # 5) optionally log‐scale Y
  if (log_y) {
    p <- p  # no log transformation applied in this code, but placeholder kept
  }
  
  return(p)
}

p_R0    <- make_param_plot("R0",    "Basic Reproduction Number (R₀) Bias", log_y = FALSE)

# 4. Save or display
ggsave("bias_box_R0.png", p_R0, width = 4, height = 6, dpi = 300)
