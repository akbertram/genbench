# adapted from: https://github.com/bedatadriven/activityinfo-R/blob/development/inst/ci-build.R
# ieuan clay
# may 2015

### run all benchmarks in genbench

# Install required packages
install.dependencies <- function() {
  
  # Set a CRAN mirror to use
  options(repos=structure(c(CRAN="http://cran.rstudio.com")))
  
  # Install packages to Jenkins' root folder
  .libPaths(new="~/R/libs")
  
  # Install CRAN packages
  for(pkg in c('stats', 'utils', 'biclust', 's4vd', 'irlba')) {
    if(!(pkg %in% installed.packages())) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  
  # Install bioconductor packages
  source("http://bioconductor.org/biocLite.R")
  for(pkg in c('affy', 'hgu133plus2cdf', 'limma', 'edgeR')) {
    if(!(pkg %in% installed.packages())) {
      biocLite(pkg, suppressUpdates = TRUE, suppressAutoUpdate = TRUE, ask = FALSE)
    }
  }
}

# install required packages if any
install.dependencies()

# find and run all benchmark scripts
for (SCRIPT in dir(file.path(getwd(), "code"), 
                  full.names = TRUE, recursive = TRUE, pattern = "\\.R$", 
                  ignore.case = TRUE)){
  # run benchmark script
  cat(timestamp(quiet = TRUE), "Running benchmark at ", SCRIPT,"\n")
  source(SCRIPT, chdir = TRUE) # all scripts assume working dir is same as script
}