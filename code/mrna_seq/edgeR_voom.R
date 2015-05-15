# mRNAseq with edgeR and limma
# ieuan.clay@gmail.com
# April 2015

### code based on examples found in:
# http://www.bioconductor.org/packages/release/bioc/html/edgeR.html
# http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# http://bioinf.wehi.edu.au/RNAseqCaseStudy/

### set up session
rm(list=ls())
## (bioconductor) packages
library(limma)
library(edgeR)

## global vars
VERBOSE <- TRUE # print progress?
INPUT <- "http://bowtie-bio.sourceforge.net/recount/countTables/gilad_count_table.txt"
PDATA <- "http://bowtie-bio.sourceforge.net/recount/phenotypeTables/gilad_phenodata.txt"

# holder for results
RESULTS <- list()
TIMES <- list()
BENCHMARK <- "edgeR_voom"

# reproducibility
set.seed(8008)
# stopifnot(all(rev(strsplit(getwd(), "[\\\\/]", perl=TRUE)[[1]])[1:3] == c("mrna_seq","code","genbase"))) # must be run from directory containing rppa.R
stopifnot(file.exists("../../data")) # data path is relative

#### functions

do.load <- function(INPUT){
  ### import data
  # read data into data.frames
  counts <- read.delim(INPUT)
  row.names(counts) <- counts$gene
  counts <- counts[,2:ncol(counts)]
  pd <- read.delim(PDATA, sep = " ", stringsAsFactors=FALSE)
  
  # convert data to DGEList instance
  dge <- DGEList(counts=counts, group = pd$gender)
  
  return(dge)
}

do.norm <- function(dge){
  
  ## filter non-detected features
  ## carry out normalisation
  ## return "EList" instance (appropriate for limma)
  
  # calculate normalisation factors
  dge <- calcNormFactors(dge)
  
  ### filter non-detected features
  filterfun <- function(dgelist, threshold=1, fraction=0.25, test_only=FALSE){ 
    # returns a logical vector, per row
    # TRUE : row (feature) is expressed at or above [threshold]
    # in ([fraction] * 100)% of the samples
    
    filter.vec <- apply(dgelist$counts, MARGIN = 1, # row-wise (probesets)
                        function(x) (sum(x >= threshold) / length(x)) >= fraction )
    if(test_only){
      # just test filter settings, return nothing
      cat(sprintf(
        "Threshold : %i, Fraction ; %0.2f \n\t\t= %i features, (%0.2f percent), would be retained.\n", 
        threshold, fraction, sum(filter.vec), (sum(filter.vec)/dim(dgelist$counts)[1])*100
      ))
      return(sum(filter.vec))
    } else{
      return( 
        filter.vec
      )  
    }
  }
  
  # for simplicity continue only with MAS5 expression values
  filterfun(dge, threshold = 1, fraction = 0.1, test_only = TRUE) # i.e. keep if any signal seen in any sample
  dge <- dge[filterfun(dge, threshold = 1, fraction = 0.1),]
  
  ### normalise (TMM / skedacity compensation)
  design <- model.matrix(~0+dge$samples$group) # levels = F, M, i.e. female is reference
  colnames(design) <- levels(dge$samples$group)
  
  ptm <- proc.time()
  vm <- voom(dge,design,plot=FALSE)
  TIMES$norm.voom <<- proc.time() - ptm

  return(vm)
}

do.limma <- function(vm){
  ### fit linear model, within simple 2 group contrast
  fit <- lmFit(vm,vm$design)
  cont.matrix <- makeContrasts(MvF=M-F, levels=vm$design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  ### export results
  return(topTable(fit2, adjust="BH"))
}


### run functions and time them
## load data and compute matrix
TIMES$load <- system.time(gcFirst = T,
                          dge <- do.load(INPUT)
)
TIMES$norm <- system.time(gcFirst = T,
                          vm <- do.norm(dge)
)
TIMES$limma <- system.time(gcFirst = T,
                          RESULTS$limma <- do.limma(vm)
)


## output results for comparison
# check output directories exist
if(!file.exists(file.path("..", "..", "generated", "results", basename(getwd())))){
  
  dir.create(file.path("..", "..", "generated", "results", basename(getwd())), recursive = TRUE)
  
}
if(!file.exists(file.path("..", "..", "generated", "timings", basename(getwd())))){
  
  dir.create(file.path("..", "..", "generated", "timings", basename(getwd())), recursive = TRUE)
  
}
# write results to file
write.table(file=file.path("..", "..", "generated", "results", 
                           basename(getwd()), paste(BENCHMARK, "results", "tsv", sep = '.')), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
            do.call("rbind", 
                    lapply(names(RESULTS), function(x){
                      cbind(RESULTS[[x]], data.frame(res=x))
                    }
                    )
            )
)

# timings
write.table(file=file.path("..", "..", "generated", "timings", 
                           basename(getwd()), paste(BENCHMARK, format(Sys.time(), "%Y%m%d%H%M%S"), "tsv", sep = '.')), 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE,
            format(do.call("rbind", TIMES), digits=5)
)

# final clean up
rm(list=ls())
gc()
