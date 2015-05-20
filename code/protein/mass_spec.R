# mass spectrometry
# ieuan.clay@gmail.com
# Thanks to [Yann Abraham](https://github.com/yannabraham) for suggestions and contributions
# May 2015

### code based on examples found the Bioconductor proteomics workflow:
# http://www.bioconductor.org/help/workflows/proteomics/

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

# holder for results
RESULTS <- list()
TIMES <- list()
BENCHMARK <- "mass_spec"

#### functions

do.download <- function(INPUT, DATA_DIR){
  
  ## download files from [INPUT] to [DATA_DIR]
  # download and unpack data
  
   # PRD000032
  
  
}

do.load <- function(){
  
  MSnSet
  
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
