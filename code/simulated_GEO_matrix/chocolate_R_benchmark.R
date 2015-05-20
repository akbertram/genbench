# this is a plain-R version of vanilla_R_benchmark.R, without extra data management packages
# only works with NGENES and NPATIENTS >= 1000

## set up session
rm(list=ls())
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative
set.seed(8008)
# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

# packages
library(biclust)
library(s4vd)
library(irlba)

# data collection
RESULTS <- list()
TIMES <- list()
BENCHMARK <- "chocolate_geo"

# needs info about path and what size of data to run on
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  PATH <- args[1]
  NGENES <- args[2]
  NPATIENTS <- args[3]
} else {
  PATH <- "../../data/simulated_GEO_matrix/"
  NGENES <- "500"
  NPATIENTS <- "500"
}
GEO <-      file.path(PATH, paste('GEO-', NGENES, '-', NPATIENTS, '.rds', sep=""))
GO <-       file.path(PATH, paste('GO-', NGENES, '-', NPATIENTS, '.rds', sep=""))
GENES <-    file.path(PATH, paste('GeneMetaData-', NGENES, '-', NPATIENTS, '.rds', sep=""))
PATIENTS <- file.oath(PATH, paste('PatientMetaData-', NGENES, '-', NPATIENTS, '.rds', sep=""))


# plain-R q&d replacement for acast(A, list(names(A)[1], names(A)[2]))
df2mxc <- function(df) {
  d1 <- factor(df[,1])
  d2 <- factor(df[,2])
  m <- matrix(data=NA, nrow=length(levels(d1)), 
    ncol=length(levels(d2)), dimnames=list(levels(d1), levels(d2)))
  m[cbind(d1, d2)] <- df[,3]
  m
}

# plain-R q&d replacement for sparseMatrix(go[,1], go[,2], x=go[,3])
df2mxs <- function(df) {
  d1 <- df[,1]
  d2 <- df[,2]
  m <- matrix(data=NA, nrow=max(d1), 
    ncol=max(d2))
  m[cbind(d1, d2)] <- df[,3]
  m
}


regression <- function() {
  ptm = proc.time()

  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  genes    <- readRDS(GENES)
  patients <- readRDS(PATIENTS)

  # subset
  sub_gmd = genes[genes$func < 250,]
  colnames(sub_gmd)[1] = "geneid"

  response = patients[,"drug.response"]

  # join
  A = merge(geo, sub_gmd)[,c("patientid", "geneid", "expression.value")]
  
  # matrix cast
  A <- df2mxc(A)

  ### Data management ops end ###
  TIMES$regression.data_management <<- proc.time() - ptm
  ptm = proc.time()

  # run regression
  RESULTS$ regression <<- lm.fit(x=A, y=response) # todo: report data.frame
  TIMES$regression.analytics <<- proc.time() - ptm
}

covariance <- function() {
  ptm <- proc.time()

  ### Data Management ops start ###

  geo      <- readRDS(GEO)
  genes    <- readRDS(GENES)
  patients <- readRDS(PATIENTS)

  sub_pmd <- patients[patients$disease==5,]

  # convert to data tables
  colnames(sub_pmd)[1] = "patientid"

  # join
  A <- merge(geo, sub_pmd)[,c("patientid", "geneid", "expression.value")]
  
  # convert to matrix
  A <- df2mxc(A)

  midtm <- proc.time() - ptm
  ptm <- proc.time()  

  # calculate covariance
  covar <- cov(A)
  TIMES$covariance.analytics <<- proc.time() - ptm
  ptm <- proc.time()

  covar <- which(covar>0.01*(max(covar)), arr.ind=T)
  res <- merge(covar, genes, by.x='row', by.y='id')
  RESULTS$covariance <<- merge(res, genes, by.x='col', by.y='id')  
 
  ### Data management ops end ###
  TIMES$covariance.data_management <<- (proc.time() - ptm) + midtm
}

biclustering<-function() {
  ptm = proc.time()

  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  patients <- readRDS(PATIENTS)

  sub_pmd <- patients[patients$gender == 1 & patients$age <= 40, ]
  colnames(sub_pmd)[1] <- "patientid"
  A <- merge(geo, sub_pmd)[,c("patientid", "geneid", "expression.value")]
  A <- df2mxc(A)

  ### Data management ops end ###

  TIMES$biclust.data_management <<- proc.time() - ptm
  ptm <- proc.time()
  
  # run biclustering
  RESULTS$biclust <<- biclust(A, method=BCssvd, K=1) # todo: report a data.frame
  TIMES$biclust.analytics <<- proc.time() - ptm
} 

svd_irlba <- function() {
  ptm <- proc.time()

  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  genes    <- readRDS(GENES)

  sub_gmd <- genes[genes$func < 250,]

  # convert to data tables
  colnames(sub_gmd)[1] = "geneid"
  # join
  A <- merge(geo, sub_gmd)[,c("patientid", "geneid", "expression.value")]

  # store as matrix
  A <- df2mxc(A)

  ### Data management ops end ###
  TIMES$svd.data_management <<- proc.time() - ptm
  ptm <- proc.time()

  # run svd
  RESULTS$svd <- irlba(A, nu=50, nv=50, sigma="ls")
  TIMES$svd.analytics <<- proc.time() - ptm
}

stats <- function() {
  ptm <- proc.time()

  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  go       <- readRDS(GO)

  # update code to start all ids at 1
  geo[,1] <- geo[,1]+1
  geo[,2] <- geo[,2]+1

  # select subset of patients, but breaks if we do. why?? too few
  #geo <- geo[geo$patientid < 0.0025 * max(geo$patientid),]
  A <- df2mxc(geo)

  go[,1] <- go[,1] + 1
  go[,2] <- go[,2] + 1
  go <- df2mxs(go)

  ### Data management ops end ###
  TIMES$stats.data_management <<- proc.time() - ptm
  ptm <- proc.time()

  for   (ii in 1:dim(go)[2]) {
    for (jj in 1:dim(A) [1]) {
      set1 <- A[jj, go[,ii] == 1]
      set2 <- A[jj, go[,ii] == 0]
      print(str(set1))
      print(str(set2))

      w < wilcox.test(set1, set2, alternative="less")
    }
  }
  # todo: capture results
  
  TIMES$stats.analytics <<- proc.time() - ptm
}

### reporting of timings
TIMES$regression <- system.time(regression(),   gcFirst=T)
TIMES$svd <- system.time(svd_irlba(),    gcFirst=T)
TIMES$covariance <- system.time(covariance(),   gcFirst=T)
TIMES$biclust <- system.time(biclustering(), gcFirst=T)
#TIMES$stats <- system.time(stats(),        gcFirst=T) 

## output results for comparison
# check output directories exist
check_generated()
# write results to file
report_results(RESULTS = RESULTS, BENCHMARK = BENCHMARK)

# timings
report_timings(TIMES = TIMES, BENCHMARK = BENCHMARK)

# final clean up
rm(list=ls())
gc()