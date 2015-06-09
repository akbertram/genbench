# ieuan clay
# ieuan.clay@gmail.com
# may 2015

### examine results of benchmark runs in genbench
USE_DB <- TRUE

### set up environment
rm(list=ls())
library(ggplot2)
library(dplyr)
library(reshape)

### functions
collect_reports <- function(USE_DB=TRUE, usr="foo", pwd="bar"){
  if(USE_DB){
    res <- collect_reports.mysql(usr, pwd, )
  } else {
    
    # use locally cached files instead
    res <- collect_reports.local()
    
  }
  
  return(res)
}
collect_reports.mysql <- function(usr="foo", pwd="bar"){
  
  source(file.path("db", "sql_utilities.R"), chdir = TRUE)
  
  conn <- getConnection(usr, pwd, jdbcdriver = "db//mysql-connector-java-5.1.35.jar") 
  
  res <- dbGetQuery(conn, 
                    "
                      SELECT 
                        t.block, t.variable, t.value,
                        CONCAT(m.benchmark, '.', m.timestamp) AS id,
                        m.benchmark_group, m.benchmark, m.timestamp,
                        m.sys_name, m.sys_release,
                        m.lang, m.lang_major, m.lang_minor
                      FROM
                        timings t
                        JOIN meta m ON t.meta_id=m.meta_id
                      ORDER BY
                        m.benchmark_group, m.benchmark, m.timestamp, t.variable
                      ;
                    "
                    )
  return(res)
  
}

collect_reports.local <- function(){
  # read in all reported benchmarks
  # TODO: update to use genbench classes for collection and handling
  require(reshape)
  require(rjson)
  reports <- lapply(
    dir(path = file.path("generated", "timings"), pattern = ".tsv$",
        full.names = TRUE, recursive = TRUE),
    function(path){
      # read metadata header if it exists
      if(length(grep(readLines(path, 1), pattern = "^\\{"))){
        meta <- fromJSON(readLines(path, 1))
      } else {
        meta <- list("platform","arch","os","system","status","major",
                     "minor","year","month","day","svn rev","language",
                     "version.string","nickname","sysname","release", "version"
                     )
        names(meta) <- meta
        meta <- lapply(meta, function(x) "not set") # set all values to not set
      }
      
      # read data and parse file names
      tmp <- read.delim(path, comment.char = "{")
      tmp$block <- rownames(tmp)
      # melt to tall and skinny for plotting
      tmp <- melt(tmp, id.vars = c("block"), na.rm = TRUE)
      suppressWarnings(tmp$value <- as.numeric(tmp$value))
      tmp <- tmp[!is.na(tmp$value),] # drop empty timings
      
      # parse path to group benchmarks
      tmp$id <- strtrim(basename(path), nchar(basename(path)) - 4) # drop file extension
      tmp$benchmark_group <- basename(dirname(path))
      tmp$benchmark <- strsplit(tmp$id, '\\.')[[1]][[1]]
      tmp$timestamp <- strsplit(tmp$id, '\\.')[[1]][[2]]
      
      # add meta data to results
      tmp$sys_name <- meta$sysname
      tmp$sys_release <- meta$release
      tmp$lang <- meta$language
      tmp$lang_major <- meta$major
      tmp$lang_minor <- meta$minor
       
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
res <- collect_reports(USE_DB = FALSE)

### plot
g <- ggplot(data = res)
g +
  geom_jitter(aes(y=value, x=block, colour=variable), 
              position = position_jitter(width = .5)) +
  facet_grid(sys_name + sys_release + lang + lang_major + lang_minor ~ benchmark_group + benchmark, scales = "free_x") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  ggtitle("Time per benchmark run, by block") +
  ylab("Time (seconds)") +
  xlab("Benchmark section")

# save results
ggsave(filename = "generated/timings/current.all.pdf", width = 20, height = 10)

### summarise
## total time per benchmark run (i.e. per ID)
g <- ggplot(data = 
              res %>%
                group_by(sys_name, sys_release, lang, lang_major, lang_minor, benchmark_group, benchmark, id, variable) %>%
                summarise(
                  total_time = sum(value), # sum all parts of the benchmark
                  time_stamp = as.numeric(unique(timestamp))
                          )
)
g +
  geom_jitter(aes(y=total_time, x=variable, colour=time_stamp), 
              position = position_jitter(width = .5)) +
  facet_grid(sys_name+ sys_release+ lang+ lang_major+ lang_minor ~ benchmark_group + benchmark, scales = "free_x") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  ggtitle("Total time per benchmark run") +
  ylab("Time (seconds)") +
  xlab("Timing type")

            
# save results
ggsave(filename = "generated/timings/current.summaryperrun.pdf", width = 20, height = 10)

## average time per benchmark file (i.e. per id)
g <- ggplot(data = 
              res %>%
              group_by(sys_name, sys_release, lang, lang_major, lang_minor,benchmark_group, benchmark, id, variable) %>% # total time per run
              summarise(
                value = sum(value)
                ) %>%
              group_by(sys_name, sys_release, lang, lang_major, lang_minor, benchmark_group, benchmark, variable) %>% # now summarise across runs
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
  facet_grid(sys_name+ sys_release+ lang+ lang_major+ lang_minor ~ benchmark_group + benchmark, scales = "free_x") +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  ggtitle("Summarised time per benchmark") +
  ylab("Mean time (seconds)") +
  xlab("Timing type")


## save results
ggsave(filename = "generated/timings/current.summaryperbenchmark.pdf", width = 20, height = 10)

### compare to expected results
# todo

# add expected values to above plots for visualisation


### clean up
rm(list=ls())
gc()