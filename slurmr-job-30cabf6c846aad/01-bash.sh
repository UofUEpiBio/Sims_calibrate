#!/bin/sh
#SBATCH --job-name=slurmr-job-30cabf6c846aad
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/slurmr-job-30cabf6c846aad/02-output-%A-%a.out
#SBATCH --array=1-100
#SBATCH --time=12:00:00
#SBATCH --mem_per_cpu=4G
#SBATCH --job_name=epi_sim
#SBATCH --job-name=slurmr-job-30cabf6c846aad
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  /uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/slurmr-job-30cabf6c846aad/00-rscript.r
