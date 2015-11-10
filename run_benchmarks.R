# adapted from: https://github.com/bedatadriven/activityinfo-R/blob/development/inst/ci-build.R
# ieuan clay
# may 2015

### run all benchmarks in genbench
# set loca libs for installing any packages
.libPaths(file.path("~","R","libs"))
# needs info about path and what size of data to run on
args <- commandArgs(trailingOnly = TRUE)
# reset timings?
if ("--reset" %in% args){
  # reset (remove) all timings files (tsv)
  lapply(dir(file.path("generated/timings/"), pattern = "tsv$",
                   full.names = TRUE, recursive = TRUE, include.dirs = FALSE
                   ),
        function(FILE){
          try(file.remove(FILE))
        }
  )
  # rename last results summary files (pdf)
  lapply(dir(file.path("generated/timings/"), pattern = "^current.*.pdf$",
                   full.names = TRUE, recursive = TRUE, include.dirs = FALSE),
        function(FILE){
          # rename file to a date stamped version
#           newname <- strtrim(FILE, nchar(FILE) -4)
#           newname <- sprintf("%s.%s.pdf", newname, format(Sys.time(), "%Y%m%d"))
          newname <- strsplit(basename(FILE), "current")[[1]]
          newname <- file.path(dirname(FILE), paste(c(format(Sys.time(), "%Y%m%d"), newname), sep = "", collapse = ""))
          
          try(file.rename(FILE, newname))
        }
  )
  
  # remove reset flag from args
  args <- args[args != "--reset"]
}

if (length(args) > 0) {
  NRUNS <- as.integer(args[1])
  cat(sprintf("Using %i runs per benchmark\n", NRUNS))

} else {
  NRUNS <- 1
  cat(sprintf("Using default number of runs (%i)\n", NRUNS))
}

## funcion definitions
install.dependencies <- function(cran=c(), bioc=c(), github=list()){
  # Install required packages
  # Set a CRAN mirror to use
  options(repos=structure(c(CRAN="http://cran.rstudio.com")))
  # set libPath to local user dir
  .libPaths(file.path("~","R","libs"))
  # Install CRAN packages
  for(pkg in cran) {
    if(!(pkg %in% installed.packages())) {
      tryCatch(install.packages(pkg, lib = file.path("~","R","libs"), dependencies = TRUE), 
                error=function(e,pkg=pkg){
                      cat(sprintf("CRAN installation of %s failed with the following errors", pkg))
                      e
                  })
    }
  }
  
  # Install bioconductor packages
  source("https://bioconductor.org/biocLite.R")
  for(pkg in bioc) {
    if(!(pkg %in% installed.packages())) {
      tryCatch(biocLite(pkg, suppressUpdates = TRUE, lib = file.path("~","R","libs"), #lib.loc = file.path("~","R","libs"), 
                        suppressAutoUpdate = TRUE, ask = FALSE),
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
  # levensque
  'data.table', 'magrittr', 'd3heatmap', 'RColorBrewer', 'DT',
  # clustering
  'stats','biclust', 's4vd', 'irlba',
  # ML
  "ncvreg", "boot", "lars", "lasso2", "mda", "leaps", "e1071", "MASS",
  # survival
  "survival", # not yet implemented in clinical/esrII.R
  # table processing
  "plyr", "reshape","sqldf",
  # utils
  "utils", "R.utils", "XML", #"devtools",
  # graph models
  "igraph",
  # plotting
  "ggplot2", "dplyr", "ggvis", "googleVis", 
  # db stuff and reporting
  "RJDBC", "jsonlite", "RPostgreSQL", "rjson" #, "RJSONIO"
)

bioc <- c('Biobase', 'XVector', 'GenomicRanges', 'affy', 'hgu133plus2cdf', 'limma', 'edgeR', 'gage', 'STRINGdb', 'dplyr', 'Rcpp')
install.dependencies(bioc=bioc, cran=cran)


# find and run all benchmark scripts
for (SCRIPT in rev(dir(file.path(getwd(), "code"), 
                  full.names = TRUE, recursive = TRUE, pattern = "\\.R$", 
                  ignore.case = TRUE))){
  
  if(NRUNS >= 1){
    # run benchmark script
    cat(timestamp(quiet = TRUE), "Running benchmark at ", SCRIPT,"\n")
    for (x in 1:NRUNS){
        cat(sprintf("\t>>>Run %i\n", x))
        # all scripts assume working dir is same as script
        # each script run in a fresh local environment
        try({source(SCRIPT, chdir = TRUE, local=new.env())})
    }
  } else {
    # dry run, install packages only
    cat("Not running benchmark at ", SCRIPT,"\n")
  }
}

