# microarray differential expression analysis with limma
# ieuan.clay@gmail.com
# April 2015

### code based on examples found the limma user guide:
# http://www.bioconductor.org/packages/release/bioc/html/limma.html

### set up session
rm(list=ls())

## (bioconductor) packages
library(affy) # reading and normalising microarray data
library(hgu133plus2cdf) # platform annotations
library(limma) # differential expression

## global vars
VERBOSE <- TRUE # print progress?
DATA_DIR <- "../../data/microarray"
INPUT <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45417/suppl/GSE45417_RAW.tar"

# holder for results
RESULTS <- list()
TIMES <- list()

# reproducibility
set.seed(8008)
stopifnot(all(rev(strsplit(getwd(), "/")[[1]])[1:3] == c("microarray","code","genbase"))) # must be run from directory containing rppa.R

#### functions

do.load <- function(INPUT, DATA_DIR){
  ## download CEL files from [INPUT] to [DATA_DIR]
  ## load downloaded data to ExpressionSet instance and return
  
  # make sure DATA_DIR exists
  if(!file.exists(DATA_DIR)){   dir.create(DATA_DIR, recursive = TRUE)    }
  
  # download files from repo to local data directory and unpack
  ptm <- proc.time()
  download.file(INPUT,destfile = file.path(DATA_DIR, basename(INPUT), method="internal"))
  TIMES$load.download <<- proc.time() - ptm
  ptm <- proc.time()
  untar(file.path(DATA_DIR, basename(INPUT), exdir = DATA_DIR))
  TIMES$load.untar <<- proc.time() - ptm
  ptm <- proc.time()
    
  # collect names of downloaded files
  cel.files <- dir(DATA_DIR, pattern = ".CEL", full.names = TRUE)
  TIMES$load.read <<- proc.time() - ptm
  
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
  
  ptm <- proc.time()
  
  # extract expression values using
  # robust multi-array average (RMA) method
  # note: rma returns expression values in log2 scale,
  # if using other expression meaasures, convert to log2 
  # (limma expects data to be log2 transformed)
  cel.rma <- rma(cel.files, normalize = TRUE, background = TRUE)
  TIMES$norm.rma <<- proc.time() - ptm
  ptm <- proc.time()
  
  
  # the below is not nessecary, but just to demonstrate the other common
  # expression method
  cel.mas <- mas5(cel.files, normalize = TRUE, sc = 150)
  exprs(cel.mas) <- log2(exprs(cel.mas))
  TIMES$norm.mas5 <<- proc.time() - ptm
  ptm <- proc.time()
  
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
  TIMES$norm.filter <<- proc.time() - ptm
  
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
              return(tt)
            } 
            ))
  
  return(results)
}


### run and time code
## load data and compute matrix
TIMES$load <- system.time(gcFirst = T,
                          cel.files <- do.load(INPUT, DATA_DIR)
)
TIMES$qc <- system.time(gcFirst = T,
                          cel.files <- do.qc(cel.files)
)
TIMES$norm <- system.time(gcFirst = T,
                          eset <- do.norm(cel.files)
)
TIMES$limma <- system.time(gcFirst = T,
                          RESULTS$limma <- do.limma(eset)
)

### retporting
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
