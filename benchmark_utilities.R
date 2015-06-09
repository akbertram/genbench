# ieuan clay
# ieuan.clay@gmail.com
# started: may 2015

### utility functions, classes and methods for running and reporting benchmarks

## class definitions
# generic "genbase" class with common aspects of results and timings
# benchmark_name must be a character string, by default "NOT_SET"
# benchmark_group is by default the dir holding the benchmark script
genbench <- function(benchmark_name="NOT_SET", benchmark_group=basename(getwd())){
  # check args
  if (!is.character(benchmark_name)) stop("benchmark_name must be character string")
  structure(
    # basic data structure is just a (named) list
    list(
      data=list(),
      # attributes
      benchmark=benchmark_name, 
      wd=getwd(),
      benchmark_group=benchmark_group,
      env=c(R.Version(), Sys.info()[c("sysname", "release", "version")])
      ),
    class = "genbench") 
}

# timings S3 class
timings <- function(benchmark_name="NOT_SET", benchmark_group=basename(getwd())){
  instance <- genbench(benchmark_name, benchmark_group)
  # inherits from genbench
  # will search classpath left to right for methods 
  # (i.e. will call child method first choice, 
  # and will only call parent method if no child method exists)
  class(instance) <- append("timings", class(instance)) 
  return(instance)
}

# results S3 class
results <- function(benchmark_name="NOT_SET", benchmark_group=basename(getwd())){
  instance <- genbench(benchmark_name, benchmark_group)
  # inherits from genbench
  # will search classpath left to right for methods 
  # (i.e. will call child method first choice, 
  # and will only call parent method if no child method exists)
  class(instance) <- append("results", class(instance))
  return(instance)
}

## generic methods
# print
print.genbench <- function(obj){
  tmp <- data.frame(cbind(obj[!names(obj) %in% c("data")]))
  names(tmp) <- c("value")
  tmp <- rbind(tmp, data.frame(value=length(getRecords(obj)), row.names = c("number_of_records")))
  print(tmp)
}
# get/set group
benchmarkGroup <- function(obj, group=NULL){
  UseMethod("benchmarkGroup", obj)  
}
benchmarkGroup.default <- function(obj, group=NULL){
  warning(sprintf("Instance class cannot be handled."))
  return(obj) 
}
benchmarkGroup.genbench <- function(obj, group=NULL){
  if(is.null(group)){
    # get
    return(obj$benchmark_group)
  } else {
    # set
    obj$benchmark_group <- group
    return(obj)
  }
}
# get/set benchmark name
benchmarkName <- function(obj, benchmark=NULL){
  UseMethod("benchmarkName", obj)  
}
benchmarkName.default <- function(obj, benchmark=NULL){
  warning(sprintf("Instance class cannot be handled."))
  return(obj) 
}
benchmarkName.genbench <- function(obj, benchmark=NULL){
  if(is.null(benchmark)){
    # get
    return(obj$benchmark)
  } else {
    # set
    obj$benchmark <- benchmark
    return(obj)
  }
}

# add data
addRecord <- function(obj, record_name="not_set", record=NULL){
  UseMethod("addRecord", obj)  
}
addRecord.default <- function(obj, record_name="not_set", record=NULL){
  warning(sprintf("Instance class cannot be handled."))
  return(obj) 
}
addRecord.genbench <- function(obj, record_name="not_set", record=NULL){
  if(is.null(record)){
    warning(sprintf("Data to record must be provided."))
    return(obj)
  }
  # check record name does not already exist
  if(!is.null(names(obj$data)) && record_name %in% names(obj$data)){
    record_name <- paste(record_name, length(obj$data), sep = "_")
  }
  # add new record
  record_names <- append(names(obj$data), record_name)
  obj$data <- append(obj$data, list(record))
  names(obj$data) <- record_names
  
  return(obj)
}

# get data
getRecords <- function(obj){
  UseMethod("getRecords", obj)  
}
getRecords.default <- function(obj){
  warning(sprintf("Instance class cannot be handled."))
  return(obj) 
}
getRecords.genbench <- function(obj){
  # return all records in data slot
  return(obj$data)
}

# get environment info
getEnv <- function(obj){
  UseMethod("getEnv", obj)  
}
getEnv.default <- function(obj){
  warning(sprintf("Instance class cannot be handled."))
  return(obj) 
}
getEnv.genbench <- function(obj){
  # return environment info captured at instatiation
  return(obj$env)
}

# get filename for writing captured results to
getOutputFile <- function(obj){
  UseMethod("getOutputFile", obj)  
}
getOutputFile.default <- function(obj){
  warning(sprintf("Instance class cannot be handled."))
  return(obj)
}
getOutputFile.genbench <- function(obj){
  # return general file name for writing to
  return(
    file.path("..", "..", "generated", 
                     paste(benchmarkGroup(obj),benchmarkName(obj), "tsv", sep = '.')
    )
  )
}
# check directory exists to write output files into
checkOutputFile <- function(obj, create=FALSE){
  UseMethod("checkOutputFile", obj)  
}
checkOutputFile.default <- function(obj, create=FALSE){
  warning(sprintf("Instance class cannot be handled."))
  return(obj)
}
checkOutputFile.genbench <- function(obj, create=FALSE){
  if(file.exists(dirname(getOutputFile(obj)))){
    return(TRUE)
  } else {
    if(create){
        warning(sprintf("%s did not exist, creating...", dirname(getOutputFile(obj))))
        dir.create(dirname(getOutputFile(obj)), recursive = TRUE)
        return(TRUE)
      } else {
        return(FALSE)
      }
    
  }
  
}

getOutputFile.timings <- function(obj){
  
  makeOutputFileStamped <- function(obj){
    file.path("..", "..", "generated", "timings", 
              benchmarkGroup(obj), 
              paste(benchmarkName(obj), 
                    format(Sys.time(), "%Y%m%d%H%M%S"), # datestamped to the second
                    "tsv", sep = '.')
    )
  }
  
  out <-  makeOutputFileStamped(obj)
  
  while(file.exists(out)){
    # make sure file does not already exist
    out <-  makeOutputFileStamped(obj)
  }
  
  return( out )
}
getOutputFile.results <- function(obj){
  return(
    file.path("..", "..", "generated", "results", 
                     benchmarkGroup(obj), 
                     paste(benchmarkName(obj), "results", "tsv", sep = '.')
    )
  )
}

# reporting methods
reportRecords <- function(obj){
  UseMethod("reportRecords", obj)
}

reportRecords.default <- function(obj){
  warning(sprintf("Instance class cannot be handled."))
  return(obj)
}
reportRecords.genbench <- function(obj){
  checkOutputFile(obj, create = TRUE)
  lapply(getRecords(obj), function(rec){
    # append each record to output file
    write.table(x = rec, file = getOutputFile(obj), append = TRUE,
                sep = "\t", row.names = TRUE, col.names=TRUE      
      )    
  })
}

# reporting method for timings class
reportRecords.timings <- function(obj){
  require(rjson)
  checkOutputFile(obj, create = TRUE)
  out <- getOutputFile(obj)
  # add 1 line JSON serialized header containing environment info
  writeLines(toJSON(getEnv(obj)), con = out)
  write.table(file=out, append = TRUE,
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE,
              format(do.call("rbind", getRecords(obj)), digits=5)
              
  )
}

reportRecords.results <- function(obj){
  require(rjson)
  # check output file and collect results
  checkOutputFile(obj, create = TRUE)
  RESULTS <- getRecords(obj)
  ## make sure all Records have the same number of columns for combining and printing
  lapply(RESULTS, function(x){
    if(max(sapply(RESULTS, ncol)) > ncol(x)){ # if dataframe has fewer columns than largest DF
      x[,(ncol(x) +1):max(sapply(RESULTS, ncol))] <- NA} # add more until they have the same number
    return(x)
    } )
  # add 1 line JSON serialized header containing environment info
  #writeLines(toJSON(getEnv(obj)), con = getOutputFile(obj))
  write.table(file=getOutputFile(obj), append = TRUE,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
              do.call("rbind", 
                      lapply(names(RESULTS), function(x){
                        names(RESULTS[[x]]) <- c(1:ncol(RESULTS[[x]])) # standardise col names
                        cbind(RESULTS[[x]], data.frame(res=x))
                      }
                      )
              )
  )
}

