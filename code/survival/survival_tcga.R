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

#### functions

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
  # parameters for survival analysis
  params.nlambda <- 1000
  params.lam_max <- 1e-7
  params.lam_rat <- 0.001
  params.lambda  <- seq(params.lam_max, params.lam_rat * params.lam_max, length.out=params.nlambda)
  params.alpha   <- seq(0,1,0.1)
  params.nalpha  <- length(params.alpha)
  results        <- NULL
  #
  lambda_vec <- array(0,params.nlambda)
  alpha_vec  <- array(0,params.nalpha)
  #
  xdata <- as.matrix(surv_data$xdata)
  ydata <- surv_data$ydata
  #
  # create results structure, as it was not possible to determine before the first loop
  results <- array(0,c(params.nlambda, length(params.alpha), dim(xdata)[2]))
  #
  # for each alpha and lambda, determine the glmnet
  for (mm in 1:length(params.alpha)) {
    # set the alpha value
    alpha_v = params.alpha[mm]
    # get local results
    temp_results = glmnet( xdata, ydata, family='cox', alpha=alpha_v, nlambda=params.nlambda, standardize=FALSE )
    # save results
    for (nn in 1:length(temp_results$lambda)) {
      results[nn,mm,] <- temp_results$beta[,nn]
    }
  }
  lambda_vec <- temp_results$lambda
  alpha_vec  <- params.alpha
  
  return( list( results=results, lambdas=lambda_vec, alphas=alpha_vec))
}

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
