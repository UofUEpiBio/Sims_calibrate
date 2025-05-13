#!/bin/sh
#SBATCH --job-name=slurmr-job-323ecd6e5a6039
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/slurmr-job-323ecd6e5a6039/02-output-%A-%a.out
#SBATCH --array=1-100
#SBATCH --job-name=slurmr-job-323ecd6e5a6039
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  /uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/slurmr-job-323ecd6e5a6039/00-rscript.r
