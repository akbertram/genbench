# mass spectrometry
# ieuan.clay@gmail.com
# Thanks to [Yann Abraham](https://github.com/yannabraham) for suggestions and contributions
# May 2015

### code based on examples found the Bioconductor proteomics workflow:
# http://www.bioconductor.org/help/workflows/proteomics/
# http://arxiv.org/pdf/1305.6559.pdf
# library("RforProteomics")
# This is the ’RforProteomics’ version 1.0.3.
# Run ’RforProteomics()’ in R or visit
# ’http://lgatto.github.com/RforProteomics/’ to get started.
# source("http://bioconductor.org/biocLite.R")
# biocLite("RforProteomics")
# https://github.com/lgatto/RforProteomics/blob/master/vignettes/vigsrc/RforProteomics.R

### set up session
rm(list=ls())
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative

# reproducibility
set.seed(8008)

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

## (bioconductor) packages
# library(rpx) # requires bioc > v3.0
library(MSnbase)

## global vars
VERBOSE <- TRUE # print progress?
DATA_DIR <- file.path("..", "..", "data", "protein")
DOWNLOAD <- FALSE
INPUT <- "TODO"
BENCHMARK <- "mass_spec"

# holder for results
RESULTS <- results(benchmark_name = BENCHMARK)
TIMES <- timings(benchmark_name = BENCHMARK)


#### functions

do.download <- function(INPUT, DATA_DIR){
  
  ## download files from [INPUT] to [DATA_DIR]
  # download and unpack data
  
   # PRD000032
  
  
}

do.load <- function(){
  
  MSnSet
  
  40 x 1000 is typical size
  
  qualitative, number of spectra / molecular weight of protein
  
}


### reporting
## load data
if(DOWNLOAD){
  TIMES <- addRecord(TIMES, record_name = "download",
                     record = system.time(gcFirst = T,
                                          do.download(csv=INPUT)
                     ))
}
TIMES <- addRecord(TIMES, record_name = "load",
                   record = system.time(gcFirst = T,
                                        rppa <- do.load()
                   ))


## output results for comparison
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()
