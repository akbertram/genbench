# ieuan clay
# ieuan.clay@gmail.com
# may 2015

### examine results of benchmark runs in genbench

### set up environment
rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape)

### functions
collect_reports <- function(){
  # read in all reported benchmarks
  require(reshape)
  reports <- lapply(
    dir(path = file.path("generated", "timings"), pattern = ".tsv$",
        full.names = TRUE, recursive = TRUE),
    function(path){
      # read data and parse file names
      tmp <- read.delim(path)
      tmp$benchmark <- rownames(tmp)
      # melt to tall and skinny for plotting
      tmp <- melt(tmp, id.vars = c("benchmark"), na.rm = TRUE)
      suppressWarnings(tmp$value <- as.numeric(tmp$value))
      tmp <- tmp[!is.na(tmp$value),] # drop empty timings
      
      # parse path to group benchmarks
      tmp$path <- path
      tmp$id <- strtrim(basename(path), nchar(basename(path)) - 4) # drop file extension
      tmp$block <- basename(dirname(path))
      tmp$benchmark_group <- strsplit(tmp$id, '\\.')[[1]][[1]]
      
      tmp$time <- strsplit(tmp$id, '\\.')[[1]][[2]]
       
      return(tmp)
    }
    )
  
  return(do.call("rbind", reports))
    
  
}

############
### MAIN ###
############

### collect data
res <- collect_reports()

### plot
g <- ggplot(data = res)
g +
  geom_jitter(aes(y=value, x=benchmark, colour=variable), 
              position = position_jitter(width = .5)) +
  facet_grid(. ~ block + benchmark_group, scales = "free_x") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))

### summarise
# todo
### compare to expected results
# todo

# add expected values to above plots for visualisation