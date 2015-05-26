# adapted from: https://github.com/bedatadriven/activityinfo-R/blob/development/inst/ci-build.R
# ieuan clay
# may 2015

### run all benchmarks in genbench

# Install required packages
install.dependencies <- function(cran=c(), bioc=c() {
  
  # Set a CRAN mirror to use
  options(repos=structure(c(CRAN="http://cran.rstudio.com")))
  
  # Install packages to Jenkins' root folder
  .libPaths(new="~/R/libs")
  
  # Install CRAN packages
  for(pkg in cran) {
    if(!(pkg %in% installed.packages())) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  
  # Install bioconductor packages
  source("http://bioconductor.org/biocLite.R")
  for(pkg in bioc) {
    if(!(pkg %in% installed.packages())) {
      biocLite(pkg, suppressUpdates = TRUE, suppressAutoUpdate = TRUE, ask = FALSE)
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
  "igraph"
)
bioc <- c('Biobase', 'affy', 'hgu133plus2cdf', 'limma', 'edgeR')
install.dependencies(bioc=bioc, cran=cran)

# find and run all benchmark scripts
for (SCRIPT in rev(dir(file.path(getwd(), "code"), 
                  full.names = TRUE, recursive = TRUE, pattern = "\\.R$", 
                  ignore.case = TRUE))){
  # run benchmark script
  cat(timestamp(quiet = TRUE), "Running benchmark at ", SCRIPT,"\n")
  try({source(SCRIPT, chdir = TRUE)}) # all scripts assume working dir is same as script
  
}
# multiruns
# for (x in 1:10){cat(sprintf("run %i\n", x)); source("code//mrna_seq//edgeR_voom.R", chdir = TRUE)}