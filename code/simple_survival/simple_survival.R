# Code by Parham Solaimani
# Test-case workflow for TCGA Browser Shiny app which is developed in Mitch Levesque lab by Phil Cheng
# Analysis code provided by Phil Cheng.

##### set up session #####
rm(list=ls())
set.seed(1000)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative
.libPaths(file.path("~","R","libs"))
source(file.path("..", "..","benchmark_utilities.R"))

##### loading packages #####
library(data.table)
library(survival)
library(DT)

cat("\nLoaded all libraries.....\n") #DEBUG

##### Set global vars #####
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "TCGAbrowser"
DATA_DIR <- file.path("..","..", "data","simple_survival")
files = list.files(path = DATA_DIR, pattern = "txt$")
RESULTS <- genbench_results(benchmark_name = BENCHMARK)
TIMES <- genbench_timings(benchmark_name = BENCHMARK)

##### Functions #####


##### Blocks for timing #####
do.load <- function(DATA_DIR){
  # Load the required Data
  d1 <- fread(paste0(DATA_DIR,"/TCGA_PRAD_RNAseq.txt"), sep="\t", header=T)
  setkey(d1, Gene)
  pat <- fread(paste0(DATA_DIR,"/TCGA_PRAD_patient.txt"), sep="\t", header=T) #converts all spaces and NA and unknown into NA for R, easier for sorting
  setkey(pat, bcr_patient_barcode, name)


  g1 <- "BAZ2A"
  time <- "pfs"
  gleason <- c("2+4", "3+3", "3+4", "3+5", "4+3", "4+4", "4+5", "5+3", "5+4", "5+5")
  high <- (ncol(d1) - 2 ) * 0.75
  low <- (ncol(d1) - 2) * 0.25
  setkey(pat, gleason)
  pat.d1 <- d1[,c("Gene", pat[gleason, name]), with=F]
  setkey(pat.d1, Gene)
  #selecting row for gene and only retrieving values
  pat.d1.gene <- melt(pat.d1[g1, setdiff(colnames(pat.d1), "Gene"), with=F], id.vars = NULL, measure.vars = colnames(pat.d1)[-1], variable.name = "name", value.name="g1")
  setkey(pat.d1.gene, g1)
  pat.d1.gene[, name := factor(name, levels=name)]
  pat.d1.gene[, ':=' (high = g1 > g1[eval(high)], low = g1 < g1[eval(low)])]
  pat.d1.gene[, gene2 := high*2 + low]
  pat.d1.gene[, gene3 := mgsub2(list(c("0", "middle"), c("1", "low"), c("2", "high")), pat.d1.gene$gene2)]
  phenosgene <- merge(pat, pat.d1.gene, by= "name")
  pat.gene <- phenosgene[gene2 !=0]
  return (pat.gene)
}

do.analyse <- function(DATA){
  # Performs calculations and plotting
  pat.gene <- DATA
  m.surv <- Surv(pat.gene$pfs_days, pat.gene$pfs)
  sdf <- survdiff(m.surv ~ pat.gene$gene2)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

  survplot <- survfit(Surv(pfs_days, pfs) ~ gene2, data = pat.gene)
  half <- summary(survplot)$table[,"median"]

  #tiff("Survival_all_met.tif", height=450, width=800)
  plot(survplot, col=1:2, xlab="Survival Time (Days)", ylab="Survival", lwd=3)
  legend("topright", c(sprintf("%s low, n=%s", g1, sdf$n[1]), sprintf("Median survival %s days", half[1]) ,
                       sprintf("%s high, n=%s", g1, sdf$n[2]), sprintf("Median survival %s days", half[2])), col=c(1,1,2,2), lty=c(1,0,1,0))
  legend("bottomleft", paste("p = ", round(p.val,4)))
  lines(c(0,half[1]), c(0.5, 0.5), lwd=1, lty =2, col=1)
  lines(c(half[1],half[1]), c(0, 0.5), lwd=1, lty =2, col=1)
  lines(c(0,half[2]), c(0.5, 0.5), lwd=1, lty=2, col=2)
  lines(c(half[2],half[2]), c(0, 0.5), lwd=1, lty=2, col=2)
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

cat("\nTIMES   record=do.preprocess(DATA_rna, 'RNAseq').....\n") #DEBUG
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
