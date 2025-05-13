
# Load required libraries
library(epiworldR)
library(data.table)
library(parallel)

# Source the main simulation function
source('simulation_functions.R')

# Get the simulation IDs for this job
sim_ids <- c(3)

# Run simulations for this chunk
results_list <- list()
for (i in seq_along(sim_ids)) {
  results_list[[i]] <- simulate_and_calibrate_slurm(sim_ids[i])
}

# Save results
saveRDS(results_list, paste0('chunk_', 3, '_results.rds'))

