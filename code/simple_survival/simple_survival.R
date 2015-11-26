# Code by Parham Solaimani
# Test-case workflow for TCGA Browser Shiny app which is developed in Mitch Levesque lab by Phil Cheng
# Analysis code provided by Phil Cheng.

##### set up session #####
rm(list=ls())
set.seed(1000)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative
.libPaths(file.path("~","R","libs"))
source(file.path("..", "..","benchmark_utilities.R"))

print(.libPaths())

##### loading packages #####
library(survival)

cat("\nLoaded all libraries.....\n") #DEBUG

##### Set global vars #####
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "simple_survival"
DATA_DIR <- file.path("..","..", "data","simple_survival")
files = list.files(path = DATA_DIR, pattern = "txt$")
RESULTS <- genbench_results(benchmark_name = BENCHMARK)
TIMES <- genbench_timings(benchmark_name = BENCHMARK)

##### Functions #####


##### Blocks for timing #####
do.load <- function(DATA_DIR){
  return(readRDS(file.path(DATA_DIR, "pat.gene.rda")))
}

do.analyse <- function(DATA){
  # Performs calculations and plotting
  pat.gene <- DATA
  m.surv <- Surv(pat.gene$pfs_days, pat.gene$pfs)
  sdf <- survdiff(m.surv ~ pat.gene$gene2)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

  survplot <- survfit(Surv(pfs_days, pfs) ~ gene2, data = pat.gene)
  half <- summary(survplot)$table[,"median"]
  print(half)

}

############################################################################
################### TIMING AND REPORTING ###################################
############################################################################

cat("\nTIMES   record=do.load(DATA_DIR, 'RNAseq').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "smpl_surv_load",
                   record = system.time(gcFirst = T,
                                        DATA <- do.load(DATA_DIR)
                                        )
)

cat("\nTIMES   record=do.analyse(DATA_rna, 'RNAseq').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "smpl_surv_analyse",
                   record = system.time(gcFirst = T,
                                        do.analyse(DATA)
                                        )
)

### reporting
# write results to file
cat("\nreportRecords(RESULTS).....\n") #DEBUG
reportRecords(RESULTS)

# timings
cat("\nreportRecords(TIMES).....\n") #DEBUG
reportRecords(TIMES)

# final clean up
cat("\nrm(list=ls()).....\n") #DEBUG
rm(list=ls())
gc()
