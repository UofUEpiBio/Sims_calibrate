# Enhanced color palettes and styling
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Define consistent color palettes
# Option 1: Professional palette (recommended)
colors_professional <- c("ABC" = "#2E86AB", "LSTM" = "#A23B72")
colors_fill_professional <- c("ABC" = "#2E86AB", "LSTM" = "#A23B72")

# Option 2: High contrast palette
colors_high_contrast <- c("ABC" = "#1f77b4", "LSTM" = "#ff7f0e")
colors_fill_high_contrast <- c("ABC" = "#1f77b4", "LSTM" = "#ff7f0e")

# Option 3: Colorblind-friendly palette
colors_colorblind <- c("ABC" = "#0173B2", "LSTM" = "#DE8F05")
colors_fill_colorblind <- c("ABC" = "#0173B2", "LSTM" = "#DE8F05")

# Choose your preferred palette (change this to switch between options)
method_colors <- colors_professional
method_fills <- colors_fill_professional

# Enhanced theme for consistency
enhanced_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.3),
    strip.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(15, 15, 15, 15)
  )

# Plot 1: Mean Bias over Time (Enhanced)
p1_bias <- ggplot(bias_long, aes(x = date, y = bias, color = method, fill = method)) +
  geom_ribbon(aes(ymin = bias - 1.96*se, ymax = bias + 1.96*se), 
              alpha = 0.15, color = NA) +
  geom_line(size = 1.3, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.6, size = 0.8) +
  labs(
    title = "Mean Bias Over Time",
    subtitle = "Bias = Predicted - True | 1000 Parameter Sets with 95% Confidence Bands",
    x = "Day (0-60)", 
    y = "Mean Bias (Infected Count)",
    color = "Method", 
    fill = "Method"
  ) +
  enhanced_theme +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_fills) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    fill = guide_legend(override.aes = list(alpha = 0.3))
  )

# Plot 2: Coverage over Time (Enhanced)
p2_coverage <- ggplot(coverage_long, aes(x = date, y = coverage, color = method)) +
  geom_line(size = 1.3, alpha = 0.9) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray30", alpha = 0.8, size = 0.8) +
  geom_point(size = 1.5, alpha = 0.7) +
  annotate("text", x = max(coverage_long$date) * 0.8, y = 97, 
           label = "Target: 95%", color = "gray30", size = 3.5, fontface = "italic") +
  labs(
    title = "Coverage Over Time",
    subtitle = "Percentage of True Values Within 95% Confidence Intervals",
    x = "Day (0-60)", 
    y = "Coverage (%)",
    color = "Method"
  ) +
  enhanced_theme +
  scale_color_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# Plot 3: Absolute Bias over Time (Enhanced)
p3_abs_bias <- ggplot(abs_bias_long, aes(x = date, y = abs_bias, color = method)) +
  geom_line(size = 1.3, alpha = 0.9) +
  geom_point(size = 1.2, alpha = 0.6) +
  labs(
    title = "Mean Absolute Bias Over Time",
    subtitle = "Absolute Deviation from True Values",
    x = "Day (0-60)", 
    y = "Mean Absolute Bias (Infected Count)",
    color = "Method"
  ) +
  enhanced_theme +
  scale_color_manual(values = method_colors) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# Plot 4: Summary Statistics (Enhanced)
# Create a more detailed summary plot
summary_long_enhanced <- summary_table %>%
  select(-Coverage) %>%
  pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = case_when(
      Metric == "Mean_Bias" ~ "Mean Bias",
      Metric == "RMSE" ~ "RMSE",
      Metric == "Mean_Abs_Bias" ~ "Mean Absolute Bias",
      TRUE ~ Metric
    )
  )

p4_summary <- ggplot(summary_long_enhanced, aes(x = Method, y = Value, fill = Method)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_text(aes(label = round(Value, 3)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  labs(
    title = "Summary Statistics Comparison",
    subtitle = "Lower Values Indicate Better Performance",
    x = "Method", 
    y = "Value"
  ) +
  enhanced_theme +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    panel.spacing = unit(1, "lines")
  ) +
  scale_fill_manual(values = method_fills)

# Improved Plot 5: Coverage Summary (Much Cleaner)
p5_coverage_summary_improved <- ggplot(summary_table, aes(x = Method, y = Coverage, fill = Method)) +
  geom_col(alpha = 0.8, width = 0.6, color = "white", size = 1) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "black", size = 1, alpha = 0.7) +
  geom_text(aes(label = paste0(round(Coverage, 1), "%")), 
            vjust = -0.7, size = 5, fontface = "bold", color = "black") +
  annotate("text", x = 0.6, y = 96.5, 
           label = "Target: 95%", color = "black", size = 4, fontface = "bold") +
  labs(
    title = "Overall Coverage Comparison",
    subtitle = "Closer to 95% is Better",
    x = "", 
    y = "Coverage (%)"
  ) +
  enhanced_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = method_fills) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(85, 105))

# Alternative: Horizontal bar chart for better readability
p5_coverage_horizontal <- ggplot(summary_table, aes(x = Coverage, y = Method, fill = Method)) +
  geom_col(alpha = 0.8, width = 0.6, color = "white", size = 1) +
  geom_vline(xintercept = 95, linetype = "dashed", color = "black", size = 1, alpha = 0.7) +
  geom_text(aes(label = paste0(round(Coverage, 1), "%")), 
            hjust = -0.1, size = 5, fontface = "bold", color = "black") +
  annotate("text", x = 95, y = 2.3, 
           label = "Target: 95%", color = "black", size = 4, fontface = "bold") +
  labs(
    title = "Overall Coverage Comparison",
    subtitle = "Closer to 95% is Better",
    x = "Coverage (%)", 
    y = ""
  ) +
  enhanced_theme +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = method_fills) +
  scale_x_continuous(limits = c(85, 105), breaks = seq(85, 105, 5))

# Alternative: Simple dot plot
p5_coverage_dots <- ggplot(summary_table, aes(x = Method, y = Coverage, color = Method)) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray40", size = 1, alpha = 0.8) +
  geom_point(size = 8, alpha = 0.8) +
  geom_text(aes(label = paste0(round(Coverage, 1), "%")), 
            color = "white", size = 4, fontface = "bold") +
  annotate("text", x = 0.6, y = 96.5, 
           label = "Target: 95%", color = "gray40", size = 4, fontface = "italic") +
  labs(
    title = "Overall Coverage Comparison",
    subtitle = "Percentage of True Values Within Confidence Intervals",
    x = "", 
    y = "Coverage (%)"
  ) +
  enhanced_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = method_colors) +
  scale_y_continuous(limits = c(90, 100), breaks = seq(90, 100, 2))

# Updated comprehensive plot with improved coverage chart
final_plot_comprehensive_fixed <- (pct_bias_long | p2_coverage) / 
  (p4_summary | p5_coverage_summary_improved) +
  plot_annotation(
    title = "Comprehensive Performance Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    caption = "Top row: Time series analysis | Bottom row: Summary statistics",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray40", margin = margin(b = 12)),
      plot.caption = element_text(size = 12, hjust = 1, color = "gray50", margin = margin(t = 15)),
      plot.margin = margin(25, 25, 25, 25)
    )
  )

# Alternative with horizontal coverage plot
final_plot_horizontal_coverage <- (p_pct_bias | p2_coverage) / 
  (p4_summary | p5_coverage_horizontal) +
  plot_annotation(
    title = "Comprehensive Performance Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    caption = "Top row: Time series analysis | Bottom row: Summary statistics",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray40", margin = margin(b = 12)),
      plot.caption = element_text(size = 12, hjust = 1, color = "gray50", margin = margin(t = 15)),
      plot.margin = margin(25, 25, 25, 25)
    )
  )

# Alternative with dot plot
final_plot_dots_coverage <- (p1_bias | p2_coverage) / 
  (p4_summary | p5_coverage_dots) +
  plot_annotation(
    title = "Comprehensive Performance Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    caption = "Top row: Time series analysis | Bottom row: Summary statistics",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray40", margin = margin(b = 12)),
      plot.caption = element_text(size = 12, hjust = 1, color = "gray50", margin = margin(t = 15)),
      plot.margin = margin(25, 25, 25, 25)
    )
  )

# Show the individual improved coverage plots
print("Improved vertical bar chart:")
print(p5_coverage_summary_improved)

print("Horizontal bar chart alternative:")
print(p5_coverage_horizontal)

print("Dot plot alternative:")
print(p5_coverage_dots)

# Show the complete fixed plots
print("Fixed comprehensive plot:")
print(final_plot_comprehensive_fixed)
# Plot 5: Coverage Summary (Enhanced)
p5_coverage_summary <- ggplot(summary_table, aes(x = Method, y = Coverage, fill = Method)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "gray30", size = 0.8) +
  geom_text(aes(label = paste0(round(Coverage, 1), "%")), 
            vjust = -0.5, size = 4, fontface = "bold", color = "black") +
  annotate("text", x = 1.5, y = 97, 
           label = "Target: 95%", color = "gray30", size = 4, fontface = "italic") +
  labs(
    title = "Overall Coverage Comparison",
    subtitle = "Closer to 95% is Better",
    x = "Method", 
    y = "Coverage (%)"
  ) +
  enhanced_theme +
  theme(legend.position = "none") +
  scale_fill_manual(values = method_fills) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))

# Create the final combined plots with enhanced styling
final_plot_v1 <- (p1_bias | p2_coverage) +
  plot_annotation(
    title = "Bias and Coverage Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    caption = "Shaded areas represent 95% confidence intervals",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40", margin = margin(b = 10)),
      plot.caption = element_text(size = 11, hjust = 1, color = "gray50", margin = margin(t = 10)),
      plot.margin = margin(20, 20, 20, 20)
    )
  )

final_plot_v2 <- (p4_summary | p5_coverage_summary) +
  plot_annotation(
    title = "Summary Performance Metrics: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    caption = "Values shown on bars for precise comparison",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40", margin = margin(b = 10)),
      plot.caption = element_text(size = 11, hjust = 1, color = "gray50", margin = margin(t = 10)),
      plot.margin = margin(20, 20, 20, 20)
    )
  )

# Alternative: Comprehensive 2x2 layout
final_plot_comprehensive <- (p1_bias | p2_coverage) / (p4_summary | p5_coverage_summary) +
  plot_annotation(
    title = "Comprehensive Performance Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    caption = "Top row: Time series analysis | Bottom row: Summary statistics",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray40", margin = margin(b = 12)),
      plot.caption = element_text(size = 12, hjust = 1, color = "gray50", margin = margin(t = 15)),
      plot.margin = margin(25, 25, 25, 25)
    )
  )

# Additional plot: Absolute Bias if needed
final_plot_with_abs_bias <- (p1_bias | p2_coverage) / (p3_abs_bias | p5_coverage_summary) +
  plot_annotation(
    title = "Detailed Bias and Coverage Analysis: ABC vs LSTM Methods",
    subtitle = "1000 Parameter Sets × 100 Runs Each × 60 Days | Infected Counts",
    caption = "Comprehensive view including absolute bias trends",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray40", margin = margin(b = 12)),
      plot.caption = element_text(size = 12, hjust = 1, color = "gray50", margin = margin(t = 15)),
      plot.margin = margin(25, 25, 25, 25)
    )
  )


# Print the plots
print(final_plot_v1)
print(final_plot_v2)
print(final_plot_comprehensive)

# Save high-resolution versions
ggsave("bias_coverage_analysis_enhanced.png", final_plot_v1, 
       width = 14, height = 7, dpi = 300, bg = "white")
ggsave("summary_analysis_enhanced.png", final_plot_v2, 
       width = 14, height = 7, dpi = 300, bg = "white")
ggsave("comprehensive_analysis_enhanced.png", final_plot_comprehensive, 
       width = 14, height = 12, dpi = 300, bg = "white")

# Color palette reference for easy switching
cat("Available color palettes:\n")
cat("Professional: ABC =", method_colors["ABC"], "| LSTM =", method_colors["LSTM"], "\n")
cat("High Contrast: ABC = #1f77b4 | LSTM = #ff7f0e\n")
cat("Colorblind-friendly: ABC = #0173B2 | LSTM = #DE8F05\n")