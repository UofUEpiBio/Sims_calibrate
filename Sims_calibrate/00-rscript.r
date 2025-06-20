.libPaths(c("/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4", "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0", "/uufs/chpc.utah.edu/sys/installdir/r8/R/4.4.0/lib64/R/library"))
message("[slurmR info] Loading variables and functions... ", appendLF = FALSE)
Slurm_env <- function (x = "SLURM_ARRAY_TASK_ID") 
{
    y <- Sys.getenv(x)
    if ((x == "SLURM_ARRAY_TASK_ID") && y == "") {
        return(1)
    }
    y
}
ARRAY_ID  <- as.integer(Slurm_env("SLURM_ARRAY_TASK_ID"))

# The -snames- function creates the write names for I/O of files as a 
# function of the ARRAY_ID
snames    <- function (type, array_id = NULL, tmp_path = NULL, job_name = NULL) 
{
    if (length(array_id) && length(array_id) > 1) 
        return(sapply(array_id, snames, type = type, tmp_path = tmp_path, 
            job_name = job_name))
    type <- switch(type, r = "00-rscript.r", sh = "01-bash.sh", 
        out = "02-output-%A-%a.out", rds = if (missing(array_id)) "03-answer-%03i.rds" else sprintf("03-answer-%03i.rds", 
            array_id), job = "job.rds", stop("Invalid type, the only valid types are `r`, `sh`, `out`, and `rds`.", 
            call. = FALSE))
    sprintf("%s/%s/%s", tmp_path, job_name, type)
}
TMP_PATH  <- "/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate"
JOB_NAME  <- "Sims_calibrate"

# The -tcq- function is a wrapper of tryCatch that on error tries to recover
# the message and saves the outcome so that slurmR can return OK.
tcq <- function (...) 
{
    ans <- tryCatch(..., error = function(e) e)
    if (inherits(ans, "error")) {
        ARRAY_ID. <- get("ARRAY_ID", envir = .GlobalEnv)
        msg <- paste0("[slurmR info] An error has ocurred while evualting the expression:\n[slurmR info]   ", 
            paste(deparse(match.call()[[2]]), collapse = "\n[slurmR info]   "), 
            "\n[slurmR info] in ", "ARRAY_ID # ", ARRAY_ID., 
            "\n[slurmR info] The error will be saved and quit R.\n")
        message(msg, immediate. = TRUE, call. = FALSE)
        ans <- list(res = ans, array_id = ARRAY_ID., job_name = get("JOB_NAME", 
            envir = .GlobalEnv), slurmr_msg = structure(msg, 
            class = "slurm_info"))
        saveRDS(list(ans), snames("rds", tmp_path = get("TMP_PATH", 
            envir = .GlobalEnv), job_name = get("JOB_NAME", envir = .GlobalEnv), 
            array_id = ARRAY_ID.))
        message("[slurmR info] job-status: failed.\n")
        q(save = "no")
    }
    invisible(ans)
}
message("done loading variables and functions.")
message("[slurmR info] Loading packages ... ")
tcq({
  library(tidyverse, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(ggplot2, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(tibble, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(tidyr, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(readr, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(purrr, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(dplyr, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(stringr, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(forcats, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(lubridate, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(glue, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(reticulate, lib.loc = "/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4")
  library(epiworldR, lib.loc = "/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4")
  library(gridExtra, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(patchwork, lib.loc = "/uufs/chpc.utah.edu/common/home/u1418987/R/x86_64-pc-linux-gnu-library/4.4")
  library(data.table, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
  library(cowplot, lib.loc = "/uufs/chpc.utah.edu/sys/installdir/r8/RLibs/4.4.0")
})
message("[slurmR info] done loading packages.")
tcq({
  INDICES <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/INDICES.rds")
})
tcq({
  X <- readRDS(sprintf("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/X_%04d.rds", ARRAY_ID))
})
tcq({
  FUN <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/FUN.rds")
})
tcq({
  mc.cores <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/mc.cores.rds")
})
tcq({
  simulate_and_calibrate <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/simulate_and_calibrate.rds")
})
tcq({
  simulate_epidemic_observed <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/simulate_epidemic_observed.rds")
})
tcq({
  simulate_epidemic_calib <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/simulate_epidemic_calib.rds")
})
tcq({
  simulation_fun <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/simulation_fun.rds")
})
tcq({
  summary_fun <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/summary_fun.rds")
})
tcq({
  proposal_fun <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/proposal_fun.rds")
})
tcq({
  kernel_fun <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/kernel_fun.rds")
})
tcq({
  theta_use <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/theta_use.rds")
})
tcq({
  model_ndays <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/model_ndays.rds")
})
tcq({
  model_seed <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/model_seed.rds")
})
tcq({
  global_n <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/global_n.rds")
})
tcq({
  N_SIMS <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/N_SIMS.rds")
})
tcq({
  seeds <- readRDS("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/seeds.rds")
})
set.seed(seeds[ARRAY_ID], kind = NULL, normal.kind = NULL)
tcq({
  ans <- parallel::mclapply(
    X                = X,
    FUN              = FUN,
    mc.cores         = mc.cores
)
})
saveRDS(ans, sprintf("/uufs/chpc.utah.edu/common/home/u1418987/Desktop/Sims_calibrate/Sims_calibrate/03-answer-%03i.rds", ARRAY_ID), compress = TRUE)
message("[slurmR info] job-status: OK.\n")
