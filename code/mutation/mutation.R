# RPPA classification
# ieuan.clay@gmail.com
# April 2015

### set up session
rm(list=ls())
## packages
library(stats)

## global vars
# https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/
INPUT <- "http://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/RPPA_input.csv"
OUTPUT <- "http://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/RPPA_output.csv"
VERBOSE <- TRUE # print progress?

# holder for results
RESULTS <- list()
TIMES <- list()

# reproducibility
set.seed(8008)
stopifnot(all(rev(strsplit(getwd(), "/")[[1]])[1:3] == c("protein","code","genbase"))) # must be run from directory containing rppa.R

### functions

## mutation data
# prepare for downloading
if(!file.exists("../../data/mutation_pop")){
  # create directory for holding data
  dir.create(recursive = TRUE, "../../data/mutation_pop")
}
download.file(destfile ="../../data/mutation_pop/maf.zip", "http://tcga-data.nci.nih.gov/docs/publications/laml_2012/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.12.0.tar.gz")
unzip("../../data/mutation_pop/maf.zip", exdir="../../data/mutation_pop")



readLines("../../data/mutation.maf", 3) # 2 line header

head(read.delim("../../data/mutation.maf", skip=2, blank.lines.skip = TRUE)$ref_context)

## sample data


download.file(destfile ="../../data/mutation.meta.tar", 
              "http://tcga-data.nci.nih.gov/docs/publications/hnsc_2014/HNSC.Clinical.tar")
untar("../../data/mutation.meta.tar", exdir="../../data/")
readLines("../../data/mutation.meta", 3) # 2 line header

dim(read.delim("../../data/mutation.meta", skip=2, blank.lines.skip = TRUE))



## familial data
# prepare for downloading
if(!file.exists("../../data/mutation_fam")){
  # create directory for holding data
  dir.create(recursive = TRUE, "../../data/mutation_fam")
}
