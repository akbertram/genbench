# microarray differential expression analysis with limma
# ieuan.clay@gmail.com
# April 2015

### code based on examples found the limma user guide:
# http://www.bioconductor.org/packages/release/bioc/html/limma.html

### set up session
rm(list=ls())
# reproducibility
set.seed(8008)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

## (bioconductor) packages
library(Biobase)
library(affy) # reading and normalising microarray data
library(hgu133plus2cdf) # platform annotations
library(limma) # differential expression

## global vars
VERBOSE <- TRUE # print progress?
DATA_DIR <- file.path("..", "..", "data", "microarray")
DOWNLOAD <- FALSE
INPUT <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45417/suppl/GSE45417_RAW.tar"
BENCHMARK <- "limma"

# holder for results
RESULTS <- genbench_results(benchmark_name = BENCHMARK)
TIMES <- genbench_timings(benchmark_name = BENCHMARK)

#### functions

do.download <- function(INPUT, DATA_DIR){

  ## download CEL files from [INPUT] to [DATA_DIR]
  # download and unpack data

  # make sure DATA_DIR exists
  if(!file.exists(DATA_DIR)){   dir.create(DATA_DIR, recursive = TRUE)    }

  # download files from repo to local data directory and unpack
  download.file(INPUT,destfile = file.path(DATA_DIR, basename(INPUT), method="internal"))

  untar(file.path(DATA_DIR, basename(INPUT), exdir = DATA_DIR))

  return(TRUE)

}

do.load <- function(DATA_DIR){

  ## load downloaded data to ExpressionSet instance and return
  # collect names of downloaded files
  cel.files <- dir(DATA_DIR, pattern = ".CEL", full.names = TRUE)

  # construct "phenodata", i.e. metadata data.frame instance
  # for annotating CEL files with experimental groups, etc
  # by parsing sample names
  pd <- data.frame(Name=basename(cel.files), FileName=cel.files, stringsAsFactors = FALSE)
  pd$Group <-
    do.call("rbind", strsplit(pd$Name, split = "[_\\.]", perl = TRUE))[,2]
  pd$Treatment <-
    sapply(
      strsplit(
        do.call("rbind", strsplit(pd$Name, split = "[_\\.]", perl = TRUE))[,3]
        ,""), function(x) paste(x[1:3],collapse = ''))
  pd$Replicate <-
    sapply(
      strsplit(
        do.call("rbind", strsplit(pd$Name, split = "[_\\.]", perl = TRUE))[,3]
        ,""), function(x) x[4])


  # convert "phenodata" (i.e. sample metadata)
  # into AnnotatedDataFrame (class expected by import functions)
  pd <- new("AnnotatedDataFrame", data=pd)

  # read affymetrix files, and attach phenodata
  cel.files <- read.affybatch(phenoData = pd, filenames = pData(pd)$FileName)

  return(cel.files)
}

do.qc <- function(cel.files){
  ## run basic qc measures on cel file info
  ## [cel.files] must be a affybatch instance

  # check RNA degredation and filter any low quality samples
  cel.qc <- affy::AffyRNAdeg(cel.files, log.it = TRUE)
  # drop samples with abnormal RNAdegredation slope
  cel.files <- cel.files[,(cel.qc$slope >= 3 & cel.qc$slope <= 4.5)]

  return(cel.files)
}

do.norm <- function(cel.files){
  ## normalise and scale affybatch instance ([cel.files])
  ## filter out non-detected probes
  ## return ExpressionSet instance, ready for linear modeling

  # extract expression values using
  # robust multi-array average (RMA) method
  # note: rma returns expression values in log2 scale,
  # if using other expression meaasures, convert to log2
  # (limma expects data to be log2 transformed)
  cel.rma <- rma(cel.files, normalize = TRUE, background = TRUE)


  # the below is not nessecary, but just to demonstrate the other common
  # expression method
  cel.mas <- mas5(cel.files, normalize = TRUE, sc = 150)
  exprs(cel.mas) <- log2(exprs(cel.mas))

  # filter probes
  # remove probes which do not have expression
  # greater than 50 in 25% of the samples
  filterfun <- function(eset, threshold=150, fraction=0.25, test_only=FALSE){
    # returns a logical vector, per probeset
    # TRUE : probeset is expressed at or above [threshold]
    # in ([fraction] * 100)% of the samples

    filter.vec <- apply(exprs(eset), MARGIN = 1, # row-wise (probesets)
                        function(x) (sum(x >= log2(threshold)) / length(x)) >= fraction )
    if(test_only){
      # just test filter settings, return nothing
      cat(sprintf(
        "%i probes, (%0.2f percent), would be retained.\n",
        sum(filter.vec), (sum(filter.vec)/dim(exprs(eset))[1])*100
        ))
      return(sum(filter.vec))
    } else{
    return(
        filter.vec
      )
    }
  }

  # for simplicity continue only with MAS5 expression values
  filterfun(cel.mas, threshold = 50, test_only = TRUE)
  cel.filtered <- cel.mas[ filterfun(cel.mas, threshold = 50)  ,]

  return(cel.filtered)
}

do.limma <- function(cel.filtered){
  ### differential expression using limma (linear models for microarrays)

  ## use factorial design matrix to extract comparisons
  # create factors for:
  # group ("DOX" : depletion of ZXDC)
  # treatment ("PMA" : induction of differentiation)
  # re-order factor such that it reflects the comparisons we want
  # i.e. the "control" always comes first!
  group <- factor(pData(cel.filtered)$Group,
                  levels=unique(rev(sort(pData(cel.filtered)$Group)))
                  )
  treatment <- factor(pData(cel.filtered)$Treatment,
                      levels=unique(rev(sort(pData(cel.filtered)$Treatment)))
                      )

  ## construct design matrix
  # key questions:
  # 1. which genes respond to stimulation in wild-type cells,
  # 2. which genes respond to stimulation in depleted cells, and
  # 3. which genes respond differently in depleted compared to wild-type cells.

  design <- model.matrix(~group*treatment)
  # This creates a design matrix which defines four coefficients with the following interpretations:
  # Coefficient Comparison Interpretation
  # Intercept               |   VEH.VEH                             |   Baseline level of unstimulated WT
  # groupDOX                |   DOX.VEH-VEH.VEH                     |   Difference between unstimulated
  # treatmentPMA            |   VEH.PMA-VEH.VEH                     |   Stimulation effect for non-depleted cells
  # groupDOX:treatmentPMA   |   (DOX.PMA-DOX.VEH)-(VEH.PMA-VEH.VEH) |   Interaction

  # fit model: note that question (2) is not present in the above design,
  # so we need to construct a specific contrast matrix
  fit <- lmFit(cel.filtered, design)
  cont.matrix <- cbind(treatment_on_WT=c(0,0,1,0),treatment_on_depleted=c(0,0,1,1),Diff=c(0,0,0,1))
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)

  ## output results of linear model, BH correction
  # topTable(fit2, adjust="BH")

  # output results (top 10 differentially expressed genes,
  # use number = Inf to return all) for each contrast
  results <-
    do.call("rbind",
          lapply(
            dimnames(cont.matrix)[[2]],
            function(x){
              tt <-topTable(fit2, adjust="BH", coef=x, number = 10)
              tt$contrast <- x
              tt$ID <- rownames(tt)
              return(tt)
            }
            ))

  return(results[,c("ID", "contrast", "adj.P.Val", "logFC")]) # return 4 column dataframe
}

do.geneset.examples <- function(  row_dim = 1e4, col_dim = 20, set_size=40){
  # run gene set testing examples as provided
  # in LIMMA manual
  # http://www.bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf

  # results collection
  results <- list()

  # CAMERA
#   y <- matrix(rnorm(1000*6),1000,6)
#   design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
#
#   # First set of 20 genes are genuinely differentially expressed
#   index1 <- 1:20
#   y[index1,4:6] <- y[index1,4:6]+1
#
#   # Second set of 20 genes are not DE
#   index2 <- 21:40
#
#   camera(y, index1, design)
#   camera(y, index2, design)
#
#   camera(y, list(set1=index1,set2=index2), design)

  ## romer
  # construct random matrix
  set.seed(888) # ensure reproducibility
  y <- matrix(rnorm(row_dim*col_dim),row_dim,col_dim)
  # arbitrary design
  design <- cbind(Intercept=1,Group=round(rnorm(col_dim, 0, 1) > 0))
  # set1 of 40 genes that are genuinely differential between the groups
  iset1 <- 1:set_size
  y[iset1,design[,"Group"]==1] <- y[iset1,design[,"Group"]==1]+rnorm(40, mean=3, sd=2)
  iset <- list(
    iset1,
    # second set where only half are differential
    iset1 + round(length(iset1)/2)
    )
  for(N in 1:round(row_dim/20)){
    # append additional (hopefully) non significant sets by shifting along the rows
    tmp <- iset1+(length(iset1) + N)
    # replace any values outside of dimensions
    tmp <- unlist(lapply(tmp, function(x) ifelse(x > row_dim, x %% row_dim, x)))
    iset <- append(iset, list(tmp))
  }
  names(iset) <- paste("iset", 1:length(iset), sep='')
  # run simulation
  r <- romer(
              index=iset,
              y=y,design=design,contrast=2,nrot=99
              )
  # reformat results
  r <- data.frame(topRomer(r))
  r$id <- row.names(r)
  r$sim <- "romer"
  r <- r[,c("id", "sim", "Up", "Mixed")]
  results <- append(results, list(r))

  ## roast
  # re-use datasets from romer
  r <- mroast(y=y, index=iset,design,contrast=2)
  r <- r[ order(r$PValue.Mixed, decreasing = FALSE), ]
  r <- subset(r, Direction=="Up")
  r <- r[1:10,] # just take top 10
  r$id <- row.names(r)
  r$sim <- "roast"
  r <- r[,c("id", "sim", "PValue", "PValue.Mixed")]
  names(r) <- c("id", "sim", "Up", "Mixed")

  results <- append(results, list(r))

  # combine and reformat results
  results <- do.call("rbind", results)
  return(results)
}

do.geneset.real <- function(){
  # competitive gene set enrichment testing
  # controls for type I error by accounting for for inter-gene correlation
  # paper: http://nar.oxfordjournals.org/content/early/2012/05/24/nar.gks461.full

  # gene sets
  #http://bioinf.wehi.edu.au/software/MSigDB/

  # TODO

}


### run and time code
## get data
if(DOWNLOAD){
  TIMES <- addRecord(TIMES, record_name = "download",
                     record = system.time(gcFirst = T,
                                          do.download(INPUT, DATA_DIR)
                     )
  )
}

## load data and compute matrix
TIMES <- addRecord(TIMES, record_name = "load",
                   record = system.time(gcFirst = T,
                          cel.files <- do.load(DATA_DIR)
))
## run QC
TIMES <- addRecord(TIMES, record_name = "qc",
                   record = system.time(gcFirst = T,
                          cel.files <- do.qc(cel.files)
))
## normalise and scale
TIMES <- addRecord(TIMES, record_name = "norm",
                   record = system.time(gcFirst = T,
                          eset <- do.norm(cel.files)
))
## linear modelling
TIMES <- addRecord(TIMES, record_name = "limma",
                   record = system.time(gcFirst = T,
                          RESULTS <- addRecord(RESULTS, record_name = "limma",
                                               record = do.limma(eset))
))

if(as.integer(R.version$major) >= 3){
  ## simulated genesets
  TIMES <- addRecord(TIMES, record_name = "genesets.examples",
                     record = system.time(gcFirst = T,
                            RESULTS <- addRecord(RESULTS, record_name = "genesets.examples",
                                                 record = do.geneset.examples())

  ))
}

### reporting
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()
