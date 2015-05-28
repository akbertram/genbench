# adapted from: https://github.com/bedatadriven/activityinfo-R/blob/development/inst/ci-build.R
# ieuan clay
# may 2015

### run all benchmarks in genbench
# needs info about path and what size of data to run on
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  NRUNS <- args[1]

} else {
  cat("Using default number of runs\n")
  NRUNS <- 1
}

## funcion definitions
install.dependencies <- function(cran=c(), bioc=c()){
  # Install required packages
  # Set a CRAN mirror to use
  options(repos=structure(c(CRAN="http://cran.rstudio.com")))
  
  # Install packages to Jenkins' root folder
  .libPaths(new="~/R/libs")
  
  # Install CRAN packages
  for(pkg in cran) {
    if(!(pkg %in% installed.packages())) {
      tryCatch(install.packages(pkg, dependencies = TRUE), 
                error=function(e,pkg=pkg){
                      cat(sprintf("CRAN installation of %s failed with the following errors", pkg))
                      e
                  })
    }
  }
  
  # Install bioconductor packages
  source("http://bioconductor.org/biocLite.R")
  for(pkg in bioc) {
    if(!(pkg %in% installed.packages())) {
      tryCatch(biocLite(pkg, suppressUpdates = TRUE, suppressAutoUpdate = TRUE, ask = FALSE),
                error=function(e, pkg=pkg){
                  cat(sprintf("Bioc installation of %s failed with the following errors", pkg))
                  e
                })
    }
  }
}

## install required packages if any
# find any dependencies
if (FALSE){
  # find and report
  for (SCRIPT in rev(dir(file.path(getwd(), "code"), 
                         full.names = TRUE, recursive = TRUE, pattern = "\\.R$", 
                         ignore.case = TRUE))){
    # run benchmark script
    cat("Packages found in ", SCRIPT,"\n")
    lines <- readLines(SCRIPT)
    print(lines[grepl("library", lines)])
    print(lines[grepl("require", lines)])
  }
  
}
cran <- c(
  # clustering
  'stats','biclust', 's4vd', 'irlba',
  # table processing
  "plyr", "reshape","sqldf",
  # utils
  "utils", "R.utils", "XML",
  # graph models
  "igraph",
  # plotting
  "ggplot2", "dplyr"
)
bioc <- c('Biobase', 'affy', 'hgu133plus2cdf', 'limma', 'edgeR')
install.dependencies(bioc=bioc, cran=cran)

# find and run all benchmark scripts
for (SCRIPT in rev(dir(file.path(getwd(), "code"), 
                  full.names = TRUE, recursive = TRUE, pattern = "\\.R$", 
                  ignore.case = TRUE))){
  # run benchmark script
  cat(timestamp(quiet = TRUE), "Running benchmark at ", SCRIPT,"\n")
  for (x in 1:NRUNS){
      cat(sprintf("\t>>>Run %i\n", x))
      # all scripts assume working dir is same as script
      # each script run in a fresh local environment
      try({source(SCRIPT, chdir = TRUE, local=new.env())})
  }
}

## plot results so far
source("examine_benchmarks.R")
