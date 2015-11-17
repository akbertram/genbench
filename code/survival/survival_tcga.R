# Survival analysis
# andre.verissimo@tecnico.ulisboa.pt
# November 2015

### set up session
rm(list=ls())
# reproducibility
set.seed(8008)
stopifnot(file.exists(file.path("..", "..", "data"))) # data path is relative

## packages
library(glmnet)
library(survival)

## global vars
VERBOSE <- TRUE # print progress?
INPUT <- "../../data/survival/tcga_ov.csv"

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))
# holder for results
BENCHMARK <- "survival"
RESULTS <- results(benchmark_name = BENCHMARK)
TIMES   <- timings(benchmark_name = BENCHMARK)

# parameters for survival analysis
params = list();
params$nlambda <- 2
params$lam_max <- 1e-7
params$lam_rat <- 0.001
params$lambda  <- seq(params$lam_max, params$lam_rat * params$lam_max, length.out=params$nlambda)
params$alpha   <- seq(0,1,0.1)
params$alpha   <- c(0.1, 0.5)
params$nalpha  <- length(params$alpha)

source("survival_functions.R")

### run functions and time them
## load data and compute matrix
TIMES <- addRecord(TIMES, record_name = "load",
                   record = system.time(gcFirst = T,
                                        surv_data <- do.load(INPUT)
                   )
)
TIMES <- addRecord(TIMES, record_name = "calc",
                   record = system.time(gcFirst = T,
                                        results <- do.calc(surv_data)
                   )
)

## output results for comparison
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()
