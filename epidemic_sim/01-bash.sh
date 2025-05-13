#!/bin/sh
#SBATCH --job-name=epidemic_sim
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/epidemic_sim/02-output-%A-%a.out
#SBATCH --array=1-100
#SBATCH --job-name=epidemic_sim
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  /uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/epidemic_sim/00-rscript.r
