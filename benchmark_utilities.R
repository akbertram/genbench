# ieuan clay
# ieuan.clay@gmail.com
# started: may 2015

### utility functions for running and reporting benchmarks
# future: proper classes

report_timings <- function(TIMES, BENCHMARK="NOT_SET"){
  
  out <- file.path("..", "..", "generated", "timings", 
                   basename(getwd()), paste(BENCHMARK, format(Sys.time(), "%Y%m%d%H%M%S"), "tsv", sep = '.'))
  
  # timings
  write.table(file=out, 
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE,
              format(do.call("rbind", TIMES), digits=5)
  )
  
  cat(sprintf("Timings written to %s", out))
  return(out)
  
}

report_results <- function(RESULTS,BENCHMARK="NOT_SET"){
  
  out <- file.path("..", "..", "generated", "results", 
                   basename(getwd()), paste(BENCHMARK, "results", "tsv", sep = '.'))
  
  # write results to file
  write.table(file=out, 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
              do.call("rbind", 
                      lapply(names(RESULTS), function(x){
                        names(RESULTS[[x]]) <- c(1:ncol(RESULTS[[x]])) # standardise col names
                        cbind(RESULTS[[x]], data.frame(res=x))
                      }
                      )
              )
  )
  
  cat(sprintf("Results written to %s", out))
  return(out)
  
}

check_generated <- function(){
  # check output directories exist
  generated_paths <- list(
    results = file.path("..", "..", "generated", "results", basename(getwd())),
    timings = file.path("..", "..", "generated", "timings", basename(getwd()))
  )
  
  lapply(generated_paths, function(path){
    
    if(!file.exists(results)){
      
      dir.create(results, recursive = TRUE)
      
    }
    
  })
  
  return(all(unlist(lapply(generated_paths, file.exists))))
  
}
