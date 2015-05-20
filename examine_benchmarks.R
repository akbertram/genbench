# ieuan clay
# ieuan.clay@gmail.com
# may 2015

### examine results of benchmark runs in genbench

### set up environment
rm(list=ls())
library(ggplot2)
library(dplyr)

### functions
collect_reports <- function(){
  # read in all reported benchmarks
  reports <- lapply(
    dir(path = file.path("generated", "timings"), pattern = ".tsv$",
        full.names = TRUE, recursive = TRUE),
    function(path){
      # read data and parse file names
      tmp <- read.delim(path)
      tmp$path <- path
      tmp$id <- strtrim(basename(path), nchar(basename(path)) - 4) # drop file extension
      tmp$block <- basename(dirname(path))
      tmp$benchmark <- strsplit(tmp$id, '\\.')[[1]][[1]]
      tmp$section <- rownames(tmp)
      tmp$time <- strsplit(tmp$id, '\\.')[[1]][[2]]
      
      return(tmp)
    }
    )
  
  return(do.call("rbind", reports))
    
  
}

### plot

### summarise

### compare to expected results