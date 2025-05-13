#!/bin/bash
#SBATCH --job-name=epidemic_sim_9
#SBATCH --partition=notchpeak-freecycle
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=slurm_output/job_9.out
#SBATCH --error=slurm_output/job_9.err

# Load R module if needed (uncomment and modify as needed)
# module load R

# Run the R script
cd $SLURM_SUBMIT_DIR
Rscript slurm_scripts/job_9.R

