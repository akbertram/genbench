# RPPA classification
# ieuan.clay@gmail.com
# April 2015

### set up session
## packages
library(stats)

## global vars
# https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/
INPUT <- "https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/RPPA_input.csv"
OUTPUT <- "https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/RPPA_output.csv"

# holder for results
RESULTS <- list()
TIMES <- list()

# reproducibility
set.seed(8008)

### functions
do.load <- function(csv){
  ## load data
  
  # samples x features matrix, including some sample metadata
  rppa <- read.csv(csv, header=T, stringsAsFactors=F)
  
  # drop non-numeric columns
  rows <- rppa$TCGA_ID
  rppa <- rppa[, sapply(rppa, is.numeric)]
  row.names(rppa) <- rows
  
  # check it
  # str(rppa)
  # dim(rppa) # 3467  131
  
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
  
  # transpose input data to get distances between samples, not features
  # convert pearson correlation to distance (i.e. bound 0-1, 0 is close)
  dist_mat <- as.dist((1-cor(t(input_data), method="pearson"))/2)
  
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
# hierarchical
do.hc <- function(dist_mat){
  ptm <- proc.time()
  
  require(stats)
  
  # hierarchical clustering, WARD as linkage
  res <- hclust(d = dist_mat, method="ward.D2")
  
  TIMES$hc.clust <<- proc.time() - ptm
  
  # determine optimal clustering using within cluster SS, for a range of 'cuts'
  ptm <- proc.time() # reset clock
  cuts <- lapply(2:25, FUN = function(i){ # note: 1 cluster => 'Inf' error
    
      cluster.stats(dist_mat, cutree(res, k=i))
    }  
  )
  
  TIMES$hc.cut <<- proc.time() - ptm
  # determine optimal cut and return labels
  # TODO: use within cluster SS to be consistent with kmeans?
  best_cut <- 8 # according to paper
  res <- cutree(res, best_cut)
  res <- data.frame(id=attr(dist_mat, which = "Labels"), cluster=res)

  return(res)
}

# kmeans
do.km <- function(dist_mat){
  ptm <- proc.time()
  
  require(stats)
  
  # kmeans clustering for a range of cluster numbers
  res <- lapply(2:25, FUN = function(i){
      kmeans(dist_mat, algorithm="Hartigan-Wong", centers=i)
    }
  )
  TIMES$km.clust <<- proc.time() - ptm
  
  # determine optimal clustering using cluster.stats, for a range of 'cuts'
  ptm <- proc.time() # reset
  cuts <- lapply(res, function(x){sum(x$withinss)})
  TIMES$km.cut <<- proc.time() - ptm
  
  # TODO: determine optimal cut and return labels for this
  best_cut <- 8 # according to paper
  res <- res[best_cut]
  res <- data.frame(id=attr(dist_mat, which = "Labels"), cluster=res$cluster)
  
  return(res)
}

# random forrest
do.rf <- function(dist_mat){
  # TBD
}

# bayesian
do.bayes <- function(dist_mat){
 # TBD  
}

### reporting
## load data and compute matrix
TIMES$load <- system.time(gcFirst = T,
  rppa <- do.load(csv=INPUT)
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
RESULTS$pub <- read.csv(OUTPUT, header = TRUE, stringsAsFactors=FALSE)
names(RESULTS$pub) <- c("id", "cluster")
# replace dots with dashes in names to make comparible to names given in input
RESULTS$pub$id <-
  unlist(lapply(strsplit(
    # split on dots
    RESULTS$pub$id, split = "\\.", perl=T), 
    # rejoin with dashes
    function(x) paste(x,collapse="-"))
  )

# compare stability of number of clusters or membership to expected results
# TBD

# output results for comparison
# todo: redirect to given file or similar
write.table(file="", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
  do.call("rbind", 
          lapply(names(RESULTS), function(x){
            cbind(RESULTS[[x]], data.frame(res=x))
              }
            )
          )
  )

# timings
# todo: redirect to file or other
write.table(file="", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE,
  format(do.call("rbind", TIMES), digits=5)
)
