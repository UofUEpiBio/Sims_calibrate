# Collect results
results_list <- list()
for (i in 1:N_SIMS) {
  result_file <- paste0("result_", i, ".rds")
  if (file.exists(result_file)) {
    results_list[[i]] <- readRDS(result_file)
  }
}

# Extract and combine all daily results
daily_results <- bind_rows(lapply(results_list, function(x) x$daily_results))

# Extract and combine all parameter comparisons
param_comparisons <- bind_rows(lapply(results_list, function(x) x$param_comparison))

