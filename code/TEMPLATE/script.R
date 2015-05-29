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
TEMPLATE

## global vars
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- TEMPLATE

# holder for results
RESULTS <- results(benchmark_name = BENCHMARK)
TIMES <- timings(benchmark_name = BENCHMARK)

### functions
## Utility

## Timing blocks


### reporting
# TEMPLATE
TIMES <- addRecord(TIMES, record_name = TEMPLATE,
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name=TEMPLATE,
                                                             record=do.TEMPATE()
                   ))

## output results for comparison
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()
