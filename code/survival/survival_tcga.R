# Survival analysis
# andre.verissimo@tecnico.ulisboa.pt
# November 2015

##### set up session #####
rm(list=ls())
# reproducibility
set.seed(8008)
stopifnot(file.exists(file.path("..", "..", "data"))) # data path is relative

##### loading packages #####
library(glmnet)
library(survival)

##### Set global vars #####
VERBOSE <- TRUE # print progress?
INPUT <- "../../data/survival/tcga_ov.csv"

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))
# holder for results
BENCHMARK <- "survival"
RESULTS <- genbench_results(benchmark_name = BENCHMARK)
TIMES   <- genbench_timings(benchmark_name = BENCHMARK)

# parameters for survival analysis
params = list();
params$nlambda <- 10000 # number of diferent lambdas to be tested
#
# commented out as they are not being used
#params$lam_max <- 1e-7
#params$lam_rat <- 0.001
#params$lambda  <- seq(params$lam_max, params$lam_rat * params$lam_max, length.out=params$nlambda)
#
params$alpha   <- seq(0,1,0.1)
params$nalpha  <- length(params$alpha)

##### Blocks for timing #####

do.load <- function(INPUT){
  print('Loading all necessary data...')
  ### import data
  # read data into data.frames
  temp_dat <- read.csv(INPUT)
  dat <- data.frame( temp_dat[,c(2:(dim(temp_dat)[2]))], row.names = temp_dat[,1] )
  #
  var_len = dim(dat)[2]
  var_arr = 3:var_len
  # get the variables from data
  #  ydata is a data.frame keeps the status of the patient and time of last follow-up
  ydata     <- cbind( time=dat$time, status=dat$status)
  #  xdata keeps the gene expression for each patient
  xdata     <- dat[,var_arr]
  #
  surv_data <- list(ydata = ydata, xdata=xdata)
  return(surv_data)
}

do.calc <- function(surv_data){
  print('Starting calculation...')
  #
  alpha_vec  <- array(0,params$nalpha)
  #
  xdata <- as.matrix(surv_data$xdata)
  ydata <- surv_data$ydata
  #
  # create results structure, as it was not possible to determine before the first loop
  my_results = list()
  # for each alpha and lambda, determine the glmnet
  for (mm in 1:length(params$alpha)) {
    # set the alpha value
    alpha_v = params$alpha[mm]
    # get local results
    temp_results = glmnet( xdata, ydata, family='cox', alpha=alpha_v, nlambda=params$nlambda, standardize=FALSE )
    # save results
    item = list()
    item$lambda = temp_results$lambda
    item$beta   = temp_results$beta
    my_results[[mm]] <- item
  }
  alpha_vec  <- params$alpha
  return( list( my_results=my_results, alphas=alpha_vec))
}

############################################################################
################### TIMING AND REPORTING ###################################
############################################################################

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
