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
  #
  alpha_vec  <- array(0,params$nalpha)
  #
  xdata <- as.matrix(surv_data$xdata)
  ydata <- surv_data$ydata
  #
  # create results structure, as it was not possible to determine before the first loop
  results = list()
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
    results[[mm]] <- item
  }
  alpha_vec  <- params$alpha
  return( list( results=results, alphas=alpha_vec))
}
