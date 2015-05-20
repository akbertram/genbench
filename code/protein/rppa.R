# RPPA classification
# ieuan.clay@gmail.com
# April 2015

### set up session
rm(list=ls())

# reproducibility
set.seed(8008)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

## packages
library(stats)

## global vars
# https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/
INPUT <- "http://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/RPPA_input.csv"
OUTPUT <- "http://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/RPPA_output.csv"
VERBOSE <- TRUE # print progress?
BENCHMARK <- "rppa" # name of benchmark
DOWNLOAD <- FALSE

# holder for results
RESULTS <- list()
TIMES <- list()


### functions
do.download <- function(csv){
  ## load data
  if (VERBOSE){print("Loading Data")}
  
  # samples x features matrix, including some sample metadata
  try(download.file(csv, destfile = file.path("..","..","data","protein","rppa.csv"), method="internal"))
  
  return(TRUE)
}
do.load <- function(){
  #try(rppa <- read.csv(csv, header=T, stringsAsFactors=F))
  rppa <- read.csv(file.path("..","..","data","protein","rppa.csv"), 
                   header=T, stringsAsFactors=F, 
                   sep=",", row.names=NULL,
                   blank.lines.skip = TRUE
                   )
  
  # drop non-numeric columns
  rows <- rppa$TCGA_ID
  rppa <- rppa[, sapply(rppa, is.numeric)]
  row.names(rppa) <- rows
  
  # check it
  # str(rppa)
  # dim(rppa) # 3467  131
  
  if (VERBOSE){print("Loading Data: complete")}
  return(rppa)
}

### classification
## original description of method here: http://www.cell.com/cms/attachment/2019543870/2039643570/mmc1.pdf
# Unsupervised Clustering: 
#  - unsupervised clustering
#  - Pearson correlation was used as the distance metric
#  - Ward was used as the linkage algorithm
#  - We identified eight robust clusters
#  - The input data matrix for RPPA clustering is available in Synapse at syn1759392 and the subtype assignments are at syn1756922.

do.dist <- function(input_data){
  ## compute distance matrix
  if (VERBOSE){print("Calculating Distance Matrix")}
  # transpose input data to get distances between samples, not features
  # convert pearson correlation to distance (i.e. bound 0-1, 0 is close)
  dist_mat <- as.dist((1-cor(t(input_data), method="pearson"))/2)
  
  if (VERBOSE){print("Calculating Distance Matrix: complete")}
  return(dist_mat)
}

## unsupervised clustering
do.within.ss <- function(d = NULL, clustering){
  # cut from 'fpc' package function cluster.stats()
  #   a generalisation of the within clusters sum of squares 
  #   (k-means objective function), which is obtained if d is a Euclidean 
  #   distance matrix. For general distance measures, this is half the sum of the 
  #   within cluster squared dissimilarities divided by the cluster size.
  
  # variables
  cluster.size <- numeric(0)
  dmat <- as.matrix(d)
  within.cluster.ss <- 0
  di <- list()
  
  # iterate thought distance matrix
  for (i in 1:max(clustering)) {
    cluster.size[i] <- sum(clustering == i)
    di <- as.dist(dmat[clustering == i, clustering == i])
    if (i <= max(clustering)) {
      within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
      
    }
    
  }
  
  return(within.cluster.ss)
  
}

do.elbow <- function(df){
  ## return optimal cluster size according to elbow method
  # this is a simple method, not nessecarily recommended,
  # which looks to find a compromise between the minimal number 
  # of clusters and within group variance
  # see: http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
  # see: http://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set
  # EXPECTS: names(df) == c("cluster", "within_ss")
  
  # check ordering of cluster sizes, small to large
  df <- df[ order(df$cluster, decreasing = FALSE) ,]
  
  # calculate delta to next cluster size
  df$delta <- c(
                df[1:(nrow(df) -1),"within_ss"] - df[2:nrow(df),"within_ss"],
                0
  )
  
  # fit curve, to remove small local variability and predict k=1 to k=2
  df$k <- df$cluster-1
  df$smooth <- predict(lm(delta ~ poly((cluster), 4, raw=TRUE), data=df), data.frame(cluster=df$k))
  
  # return first cluster number (k) where
  # the difference between the change in within_ss for (k) and (k+1)
  # is less than 10% of the starting (i.e. k=1 to k=2)
  return(min(which(df$smooth < df[1,"smooth"]/10)))
  
}
# hierarchical
do.hc <- function(dist_mat){
  if (VERBOSE){print("Calculating Hierarchical clustering")}
  ptm <- proc.time()
  
  require(stats)
  
  # hierarchical clustering, WARD as linkage
  res <- hclust(d = dist_mat, method="ward")
  
  TIMES$hc.clust <<- proc.time() - ptm
  
  # determine clustering statistics (within cluster SS), for a range of 'cuts'
  ptm <- proc.time() # reset clock
  if (VERBOSE){print("Calculating Hierarchical clustering: cutting tree")}
  cuts <- lapply(2:25, FUN = function(i){ # note: 1 cluster => 'Inf' error
        if (VERBOSE){print(paste("    >", i, "clusters"))}
        return(data.frame(
            cluster = i,
            within_ss = do.within.ss(dist_mat, cutree(res, k=i))
          ))
    }  
  )
  
  TIMES$hc.cut <<- proc.time() - ptm
  # determine optimal cut and return labels
  # use within cluster SS to be consistent with kmeans
  best_cut <- do.elbow(do.call("rbind", cuts)) # 8 according to paper
  res <- cutree(res, best_cut)
  res <- data.frame(id=attr(dist_mat, which = "Labels"), cluster=res)
  if (VERBOSE){print("Calculating Hierarchical clustering: complete")}
  return(res)
}

# kmeans
do.km <- function(dist_mat){
  ptm <- proc.time()
  if (VERBOSE){print("Calculating K-means clustering")}
  require(stats)
  
  # kmeans clustering for a range of cluster numbers
  res <- lapply(2:25, FUN = function(i){
      if (VERBOSE){print(paste("    >", i, "clusters"))}
      kmeans(dist_mat, algorithm="Hartigan-Wong", centers=i)
    }
  )
  TIMES$km.clust <<- proc.time() - ptm
  
  # determine clustering statistics for a range of 'cuts'
  ptm <- proc.time() # reset
  cuts <- lapply(res, function(x){data.frame(
                                              cluster = max(x$cluster),
                                              within_ss = sum(x$withinss)
                                              )})
  TIMES$km.cut <<- proc.time() - ptm
  
  # determine optimal cut and return labels for this
  best_cut <- do.elbow(do.call("rbind", cuts)) # 8 according to paper
  res <- res[[best_cut]]
  res <- data.frame(id=attr(dist_mat, which = "Labels"), cluster=res$cluster)
  if (VERBOSE){print("Calculating K-means clustering: complete")}
  return(res)
}

# random forrest
do.rf <- function(dist_mat){
  # TODO: package randomForest not yet implemented in renjin
  
  # return 2 column dataframe of tumour ID and cluster id
}

# SVM
do.svm <- function(dist_mat){
  # TODO: package e1071 not yet implemented in renjin
  
  # return 2 column dataframe of tumour ID and cluster id
}

# bayesian
do.bayes <- function(dist_mat){
  # TODO: package e1071 not yet implemented in renjin
  
  # return 2 column dataframe of tumour ID and cluster id
}

### reporting
## load data and compute matrix
if(DOWNLOAD){
  TIMES$download <- system.time(gcFirst = T,
                            do.download(csv=INPUT))
}
TIMES$load <- system.time(gcFirst = T,
  rppa <- do.load()
)
TIMES$dist <- system.time(gcFirst = T,
  dist_mat <- do.dist(input_data=rppa)
)
## clustering
# hierarchical
TIMES$hc <- system.time(gcFirst = T,
  RESULTS$hc <- do.hc(dist_mat=dist_mat)
)
# kmeans
TIMES$km <- system.time(gcFirst = T,
  RESULTS$km <- do.km(dist_mat=dist_mat)
)

# random forrest

# bayes

### compare results to each other and to published
# collect published results
# RESULTS$pub <- read.csv(OUTPUT, header = TRUE, stringsAsFactors=FALSE)
# names(RESULTS$pub) <- c("id", "cluster")
# # replace dots with dashes in names to make comparible to names given in input
# RESULTS$pub$id <-
#   unlist(lapply(strsplit(
#     # split on dots
#     RESULTS$pub$id, split = "\\.", perl=T), 
#     # rejoin with dashes
#     function(x) paste(x,collapse="-"))
#   )

# compare stability of number of clusters or membership to expected results
# TBD


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
