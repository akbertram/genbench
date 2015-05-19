# mass spectrometry
# ieuan.clay@gmail.com
# Thanks to [Yann Abraham](https://github.com/yannabraham) for suggestions and contributions
# May 2015

### code based on examples found the Bioconductor proteomics workflow:
# http://www.bioconductor.org/help/workflows/proteomics/

### set up session
rm(list=ls())

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

# reproducibility
set.seed(8008)
stopifnot(file.exists("../../data")) # data path is relative

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