#!/bin/bash
#SBATCH --job-name=epidemic_sim_6
#SBATCH --partition=notchpeak-freecycle
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=slurm_output/job_6.out
#SBATCH --error=slurm_output/job_6.err

# Load R module if needed (uncomment and modify as needed)
# module load R

# Run the R script
cd $SLURM_SUBMIT_DIR
Rscript slurm_scripts/job_6.R

