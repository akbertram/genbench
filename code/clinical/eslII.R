# Supervised learning on clinical datasets
# reproducing examples in ESLII
# ieuan.clay@gmail.com
# May 2015

### set up session
rm(list=ls())

# reproducibility
set.seed(8008)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

## packages
library(ncvreg) # source datasets from http://cran.r-project.org/web/packages/ncvreg/ncvreg.pdf

## global vars
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "eslII"

# holder for results
RESULTS <- results(benchmark_name = BENCHMARK)
TIMES <- timings(benchmark_name = BENCHMARK)

### functions
## Utility

## Timings blocks
do.load <- function(){
  ### heart dataset
  # Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer. 
  # Rousseauw, J., et al. (1983). Coronary risk factor screening in three rural communities. South African Medical Journal, 64, 430-436.
  data(heart)
  
  ### prostate dataset
  # Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer. 
  # Stamey, T., et al. (1989). Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy treated patients. Journal of Urology, 16: 1076-1083.
  data(prostate)
  
  ### lung dataset
  # http://CRAN.R-project.org/package=survival
  # Kalbfleisch D and Prentice RL (1980), The Statistical Analysis of Failure Time Data. Wiley, New York.
  data(Lung)
  
  return(list(heart=heart, lung=Lung, prostate=prostate))
  
}

### reporting
# score on sliding window
TIMES <- addRecord(TIMES, record_name = "fam.window",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="fam.window",
                                                             record=do.ibd.window(fam.scores = fam.scores))
                   ))

## output results for comparison
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()
