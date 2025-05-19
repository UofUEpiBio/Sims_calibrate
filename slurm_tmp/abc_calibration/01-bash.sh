#!/bin/sh
#SBATCH --job-name=abc_calibration
#SBATCH --output=slurm_tmp/abc_calibration/02-output-%A-%a.out
#SBATCH --array=1-10
#SBATCH --job-name=abc_calibration
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  slurm_tmp/abc_calibration/00-rscript.r
