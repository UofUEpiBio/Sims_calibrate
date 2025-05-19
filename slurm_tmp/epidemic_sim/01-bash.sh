#!/bin/sh
#SBATCH --job-name=epidemic_sim
#SBATCH --output=slurm_tmp/epidemic_sim/02-output-%A-%a.out
#SBATCH --array=1-100
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=epidemic_sim
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  slurm_tmp/epidemic_sim/00-rscript.r
