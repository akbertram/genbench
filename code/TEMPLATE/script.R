# TEMPLATE
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
# CRAN
# Bioc

## global vars
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "TEMPLATE"

# holder for results
RESULTS <- results(benchmark_name = BENCHMARK)
TIMES <- timings(benchmark_name = BENCHMARK)

### functions
## Utility

## Blocks for timing
do.TEMPLATE <-function(plot_results=FALSE){
  ### this funciton is a template
  # it does some stuff
  # it expects some stuff
  
  df <- data.frame(x=rnorm(10000), y=rnorm(10000))
  df$d <- sqrt(df$x**2 + df$y**2)
  df$col <- factor(round(df$d))
  df <- df[ order(df$d, decreasing = TRUE), ]
  
  if(plot_results){
    
    plot(x=df$x, y=df$y, 
         type="p", 
         main="evil eye template", xlab="X", ylab="Y",
         col=df$col, pch=20, 
         )
    
  }
  
  return(head(df))
    
}

### reporting
# TEMPLATE
TIMES <- addRecord(TIMES, record_name = "TEMPLATE",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="TEMPLATE",
                                                             record=do.TEMPLATE()
                   ))
)

## output results for comparison
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()
