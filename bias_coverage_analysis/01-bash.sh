#!/bin/sh
#SBATCH --job-name=bias_coverage_analysis
#SBATCH --output=/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/bias_coverage_analysis/02-output-%A-%a.out
#SBATCH --array=1-10
#SBATCH --partition=vegayon-shared-np
#SBATCH --account=vegayon-np
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bias_coverage_analysis
#SBATCH --ntasks=1
/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/bin/Rscript  /uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/bias_coverage_analysis/00-rscript.r
