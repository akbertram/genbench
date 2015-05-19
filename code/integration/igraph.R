# microarray differential expression analysis with limma
# ieuan.clay@gmail.com
# April 2015

### code based on examples found the limma user guide:
# http://www.bioconductor.org/packages/release/bioc/html/limma.html

### set up session
rm(list=ls())

## (bioconductor) packages
library(Biobase)
library(affy) # reading and normalising microarray data
library(hgu133plus2cdf) # platform annotations
library(limma) # differential expression

## global vars
VERBOSE <- TRUE # print progress?
DATA_DIR <- file.path("..", "..", "data", "microarray")
DOWNLOAD <- FALSE
INPUT <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45417/suppl/GSE45417_RAW.tar"

# holder for results
RESULTS <- list()
TIMES <- list()
BENCHMARK <- "limma"

# reproducibility
set.seed(8008)
#stopifnot(all(rev(strsplit(getwd(), "[\\\\/]", perl=TRUE)[[1]])[1:3] == c("microarray","code","genbase"))) # must be run from directory containing rppa.R
stopifnot(file.exists("../../data")) # data path is relative

#### functions

do.download <- function(INPUT, DATA_DIR){
  
  ## download CEL files from [INPUT] to [DATA_DIR]
  # download and unpack data
