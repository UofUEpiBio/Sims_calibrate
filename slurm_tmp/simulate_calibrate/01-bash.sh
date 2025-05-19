#!/bin/sh
#SBATCH --job-name=simulate_calibrate
#SBATCH --output=slurm_tmp/simulate_calibrate/02-output-%A-%a.out
#SBATCH --array=1-10
#SBATCH --job-name=simulate_calibrate
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  slurm_tmp/simulate_calibrate/00-rscript.r
