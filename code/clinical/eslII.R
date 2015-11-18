# Supervised learning on clinical datasets
# reproducing examples in ESLII
# ieuan.clay@gmail.com
# May 2015

### set up session
rm(list=ls())

# reproducibility
set.seed(8008)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

## packages
library(ncvreg) # source datasets from http://cran.r-project.org/web/packages/ncvreg/ncvreg.pdf
library(boot)
library(lars)
library(lasso2)
library(mda)
library(leaps)
#library(survival)

## global vars
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "eslII"

# holder for results
RESULTS <- genbench_results(benchmark_name = BENCHMARK)
TIMES <- genbench_timings(benchmark_name = BENCHMARK)

### functions
## Utility

## Timings blocks
do.load <- function(){

  # run locally to avoid having extra datasets loaded
  e <- new.env()
  ### heart dataset
  # Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer.
  # Rousseauw, J., et al. (1983). Coronary risk factor screening in three rural communities. South African Medical Journal, 64, 430-436.
  data(heart, envir = e)

  ### prostate dataset
  # Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer.
  # Stamey, T., et al. (1989). Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy treated patients. Journal of Urology, 16: 1076-1083.
  data(prostate, envir = e)

  ### lung dataset
  # http://CRAN.R-project.org/package=survival
  # Kalbfleisch D and Prentice RL (1980), The Statistical Analysis of Failure Time Data. Wiley, New York.
  data(Lung, envir = e)

  alldata <- list(heart=e$heart, lung=e$Lung, prostate=e$prostate)

  return(alldata)

}

do.varselect <- function(data, plot_results=FALSE){

  ### variable selection using coordinate descent
  ### on prostate and heart datasets from ncvreg (see do.load)
  # expects input from do.load
  # see: http://myweb.uiowa.edu/pbreheny/publications/Breheny2011.pdf

  # capture
  results <- list()

  ## prostate
  # cross validation and model fitting
  X <- as.matrix(data$prostate[,1:8])
  y <- data$prostate$lpsa
  cvfit <- cv.ncvreg(X, y, penalty="lasso", seed = 8008, nfolds=10 )

  if(plot_results){
    plot(cvfit)
    summary(cvfit)
  }
  results <- append(results,
                    list(data.frame(
                      dat="prostate",
                      var=rownames(cvfit$fit$beta),
                      coeff=cvfit$fit$beta[,as.character(round(cvfit$lambda.min, digits = 4))]
                    ))
  )

  ## heart
  # cross validated model fitting
  X2 <- as.matrix(data$heart[,sapply(data$heart, is.numeric)])
  y2 <- data$heart$transplant
  cvfit <- cv.ncvreg(X2, y2, penalty="lasso", seed = 8008, nfolds=100 )

  if(plot_results){
    plot(cvfit)
    summary(cvfit)
  }
  results <- append(results,
                    list(data.frame(
                      dat="heart",
                      var=rownames(cvfit$fit$beta),
                      coeff=cvfit$fit$beta[,as.character(round(cvfit$lambda.min, digits = 4))]
                    ))
  )

  return(do.call("rbind",results))

}

do.prostate <- function(data, plot_results=FALSE){

  ### some modelling on prostate dataset from ncvreg (see do.load)
  # expects input from do.load

  # see http://www-stat.stanford.edu/ElemStatLearn
  # 3.2.1 Example: Prostate Cancer
  ## code is adapted from
  # http://cran.r-project.org/web/packages/ElemStatLearn/ElemStatLearn.pdf

  # placeholder
  results <- list()

  # examine data
  if(plot_results){
    cor( data$prostate[,1:8] )
    pairs( data$prostate[,1:9], col="violet" )
  }

  # set test/train
  traintest <- rep(TRUE, nrow(data$prostate))
  traintest[sample(1:length(traintest), size = 30, replace = FALSE)] <- FALSE

  train <- data$prostate[traintest,1:9]
  test <- data$prostate[!traintest,1:9]

  # The book (page 56) uses only train subset, so we the same:
  prostate.leaps <- regsubsets( lpsa ~ . , data=train, nbest=70,
                                really.big=TRUE )
  prostate.leaps.sum <- summary( prostate.leaps )
  prostate.models <- prostate.leaps.sum$which
  prostate.models.size <- as.numeric(attr(prostate.models, "dimnames")[[1]])
  if(plot_results){hist( prostate.models.size )}
  prostate.models.rss <- prostate.leaps.sum$rss
  prostate.models.best.rss <-
    tapply( prostate.models.rss, prostate.models.size, min )

  # add results for the only intercept model
  prostate.dummy <- lm( lpsa ~ 1, data=train )
  prostate.models.best.rss <- c(
    sum(resid(prostate.dummy)^2),
    prostate.models.best.rss)

  if (plot_results){
    plot( 0:8, prostate.models.best.rss,
          type="b", xlab="subset size", ylab="Residual Sum Square",
          col="red2" )
    points( prostate.models.size, prostate.models.rss, pch=17, col="brown",cex=0.7 )
  }
  # capture some results
  results <- append(results,
                    list(data.frame(
                      dat="rss",
                      var=names(prostate.models.best.rss),
                      coeff=prostate.models.best.rss
                    ))
  )

  ## Calculations for the lasso:
  prostate.lasso <- l1ce( lpsa ~ ., data=train, trace=TRUE, sweep.out=~1,
                          bound=seq(0,1,by=0.1) )
  prostate.lasso.coef <- sapply(prostate.lasso, function(x) x$coef)
  colnames(prostate.lasso.coef) <- seq( 0,1,by=0.1 )
  if(plot_results){
  matplot( seq(0,1,by=0.1), t(prostate.lasso.coef[-1,]), type="b",
           xlab="shrinkage factor", ylab="coefficients",
           xlim=c(0, 1.2), col="blue", pch=17 )
  }
  results <- append(results,
                    list(data.frame(
                      dat="lasso",
                      var=rownames(prostate.lasso.coef),
                      coeff=prostate.lasso.coef[,"1"]
                    ))
  )

  ## lasso with lars:
  prostate.lasso.lars <- lars( as.matrix(train[,1:8]), train[,9],
                               type="lasso", trace=TRUE )
  prostate.lasso.larscv <- cv.lars( as.matrix(train[,1:8]), train[,9], plot.it=plot_results,
           type="lasso", trace=TRUE, K=10 )
  results <- append(results,
                    list(data.frame(
                      dat="lars",
                      var=colnames(prostate.lasso.lars$beta),
                      coeff=prostate.lasso.lars$lambda
                    ))
  )

  ## CV (cross-validation) using package boot:
  prostate.glm <- glm( lpsa ~ ., data=train )
  # repeat this some times to make clear that cross-validation is
  # a random procedure
  prostate.glmcv <- cv.glm( train, prostate.glm, K=10 )

  results <- append(results,
                    list(data.frame(
                      dat="glm",
                      var=names(prostate.glm$coefficients),
                      coeff=prostate.glm$coefficients
                    ))
  )

  return(do.call("rbind",results))
}

# do.lung(data, plot_results=FALSE){
#   ### not yet implemented
#   ### survey data using Lung dataset from ncvreg
#   # expects output from do.load
#
# }

# do.survival <- function(data, plot_results=FALSE){
#   ### survival analysis on Lung dataset from ncvreg (see do.load)
#   # expects input from do.load
#
#   results <- list()
#
#   ## kaplan-meier
#   survfit()
#   results <- append(results,
#                     list(data.frame(
#                       dat="kaplan-meier",
#                       var=,
#                       coeff=
#                     ))
#   )
#
#   ## cox stats
#   results <- append(results,
#                     list(data.frame(
#                       dat="cox",
#                       var=,
#                       coeff=
#                     ))
#   )
#
#   ## return
#   return(do.call("rbind",results))
# }

### reporting
# load data
data <- do.load()
# score on sliding window
TIMES <- addRecord(TIMES, record_name = "varselect",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="varselect",
                                                             record=do.varselect(data)
                                        )
                   ))
TIMES <- addRecord(TIMES, record_name = "prostate",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="prostate",
                                                             record=do.prostate(data)
                                        )
                   ))

# not yet written
# TIMES <- addRecord(TIMES, record_name = "lung",
#                    record = system.time(gcFirst = T,
#                                         RESULTS <- addRecord(RESULTS, record_name="lung",
#                                                              record=do.lung(data)
#                                         )
#                    ))

# not yet written
# TIMES <- addRecord(TIMES, record_name = "survival",
#                    record = system.time(gcFirst = T,
#                                         RESULTS <- addRecord(RESULTS, record_name="survival",
#                                                              record=do.survival(data)
#                                         )
#                    ))

## output results for comparison
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()
