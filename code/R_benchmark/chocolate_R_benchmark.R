# this is a plain-R version of vanilla_R_benchmark.R, without extra data management packages
# only works with NGENES and NPATIENTS >= 1000

# needs info about path and what size of data to run on
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  PATH <- args[1]
  NGENES <- args[2]
  NPATIENTS <- args[3]
} else {
  PATH <- "."
  NGENES <- "500"
  NPATIENTS <- "500"
}
GEO <-      paste(PATH, '/GEO-', NGENES, '-', NPATIENTS, '.rds', sep="")
GO <-       paste(PATH, '/GO-', NGENES, '-', NPATIENTS, '.rds', sep="")
GENES <-    paste(PATH, '/GeneMetaData-', NGENES, '-', NPATIENTS, '.rds', sep="")
PATIENTS <- paste(PATH, '/PatientMetaData-', NGENES, '-', NPATIENTS, '.rds', sep="")


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
  cat(sprintf('Regression data management: %f\n', (proc.time() - ptm)['elapsed']))
  ptm = proc.time()

  # run regression
  lm.fit(x=A, y=response)
  cat(sprintf('Regression analytics: %f\n', (proc.time() - ptm)['elapsed']))
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

  midtm <- (proc.time() - ptm)['elapsed']
  ptm <- proc.time()  

  # calculate covariance
  covar <- cov(A)
  cat(sprintf('Covariance analytics: %f\n', (proc.time() - ptm)['elapsed']))
  ptm <- proc.time()

  covar <- which(covar>0.01*(max(covar)), arr.ind=T)
  res <- merge(covar, genes, by.x='row', by.y='id')
  res <- merge(res, genes, by.x='col', by.y='id')  
 
  ### Data management ops end ###
  cat(sprintf('Covariance data management: %f\n', (proc.time() - ptm)['elapsed'] + midtm))
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

  cat(sprintf('Biclust data management: %f\n', (proc.time() - ptm)['elapsed']))
  ptm <- proc.time()
  
  # run biclustering
  library(biclust)
  library(s4vd)
  biclust(A, method=BCssvd, K=1)
  cat(sprintf('Biclust analytics: %f\n', (proc.time() - ptm)['elapsed']))
} 

svd_irlba <- function() {
  library(irlba)
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
  cat(sprintf('SVD data management: %f\n', (proc.time() - ptm)['elapsed']))
  ptm <- proc.time()

  # run svd
  irlba(A, nu=50, nv=50, sigma="ls")
  cat(sprintf('SVD analytics: %f\n', (proc.time() - ptm)['elapsed']))
}

stats <- function() {
  ptm <- proc.time()

  ### Data Management ops start ###
  geo      <- readRDS(GEO)
  go       <- readRDS(GO)

  # update code to start all ids at 1
  geo[,1] <- geo[,1]+1
  geo[,2] <- geo[,2]+1

  # select subset of patients, but breaks if we do. why??
  #geo <- geo[geo$patientid < 0.0025 * max(geo$patientid),]
  A <- df2mxc(geo)

  go[,1] <- go[,1] + 1
  go[,2] <- go[,2] + 1
  go <- df2mxs(go)

  ### Data management ops end ###
  cat(sprintf('Stats data management: %f\n', (proc.time() - ptm)['elapsed']))
  ptm <- proc.time()

  for   (ii in 1:dim(go)[2]) {
    for (jj in 1:dim(A) [1]) {
      wilcox.test(A[jj, go[,ii] == 1], A[jj, go[,ii] == 0], alternative="less")
    }
  }
  cat(sprintf('Stats analytics: %f\n', (proc.time() - ptm)['elapsed']))
}

print(paste('Regression: ', system.time(regression(), gcFirst=T)['elapsed'], sep=''));
print(paste('SVD: ', system.time(svd_irlba(), gcFirst=T)['elapsed'], sep=''));
print(paste('Covariance: ', system.time(covariance(), gcFirst=T)['elapsed'], sep=''));
print(paste('Biclustering: ', system.time(biclustering(), gcFirst=T)['elapsed'], sep=''));
print(paste('Stats: ', system.time(stats(), gcFirst=T)['elapsed'], sep='')); 
