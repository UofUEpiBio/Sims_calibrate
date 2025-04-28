#!/bin/bash
#SBATCH --job-name=epidemic_sim
#SBATCH --output=sim_%A_%a.out
#SBATCH --error=sim_%A_%a.err
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --account=vegayon-np
#SBATCH --partition=notchpeak
#SBATCH --qos=vegayon-np

# Load R module
module load R/4.2.2

# Run the R script 
Rscript run_single_sim.R


