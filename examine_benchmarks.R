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
  # TODO: update to use genbench classes for collection and handling
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
  reports <- do.call("rbind", reports)
  # reorder variable for what was timed so that
  # figures look nicer
  levels(reports$variable) <- rev(sort(levels(reports$variable)))
  
  
  return(reports)
    
  
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
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  ylab("Time (seconds)") +
  xlab("Benchmark section")

# save results
ggsave(filename = "generated/timings/current.all.pdf", width = 20, height = 10)

### summarise
## total time per benchmark run (i.e. per ID)
g <- ggplot(data = 
              res %>%
                group_by(block, benchmark_group, id, variable) %>%
                summarise(
                  total_time = sum(value), # sum all parts of the benchmark
                  time_stamp = as.numeric(unique(time))
                          )
)
g +
  geom_jitter(aes(y=total_time, x=variable, colour=time_stamp), 
              position = position_jitter(width = .5)) +
  facet_grid(. ~ block + benchmark_group, scales = "free_x") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))

            
# save results
ggsave(filename = "generated/timings/current.summaryperrun.pdf", width = 20, height = 10)

## average time per benchmark file (i.e. per group)
g <- ggplot(data = 
              res %>%
              group_by(block, benchmark_group, id, variable) %>% # total time per run
              summarise(
                value = sum(value)
                ) %>%
              group_by(block, benchmark_group, variable) %>% # now summarise across runs
              summarise(
                max_time = max(value),
                mean_time = mean(value),
                variance = var(value),
                se = sd(value)/sqrt(length(value))
              ),
            aes(x=variable)
)
g +
  geom_errorbar(aes(ymin=mean_time - se, ymax=mean_time+se, colour=variable)) +
  geom_point(aes(y=mean_time, colour=variable)) +
  facet_grid(. ~ block + benchmark_group, scales = "free_x") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))


## save results
ggsave(filename = "generated/timings/current.summaryperbenchmark.pdf", width = 20, height = 10)

### compare to expected results
# todo

# add expected values to above plots for visualisation


### clean up
rm(list=ls())
gc()