library(parallel)

setClass("GenBench",
         slots = c(
           data = list(),
           benchmark_name = "character",
           wd = "character",
           benchmark_group = "character",
           env = character()
         ),

         prototype = list(
           data = list(),
           benchmark_name = "NOT_SET",
           wd = getwd(),
           benchmark_group = benchmark_group,
           env = c(R.Version(), Sys.info()[c("sysname", "release", "version")], nphyscores=detectCores(logical = FALSE), nlogcores=detectCores(logical = TRUE))
         )
)

setGeneric(name="print",
           def=function(obj)
           {
             standardGeneric("print")
           }
)

setMethod(f="print",
          signature="GenBench",
          definition=function(obj)
          {
            tmp <- data.frame(cbind(obj[!names(obj) %in% c("data")]))
            names(tmp) <- c("value")
            tmp <- rbind(tmp, data.frame(value=length(getRecords(obj)), row.names = c("number_of_records")))
            print(tmp)
          }
)

setClass("Timings",
         contains = "GenBench"
)

setClass("Results",
         contains = "GenBench"
)

setGeneric(name="benchmarkGroup",
           def=function(obj, group = NULL)
           {
             standardGeneric("benchmarkGroup")
           }
)
setMethod(f="benchmarkGroup",
          signature="GenBench",
          definition=function(obj, group = NULL)
          {
            if(is.null(group)){
              # get
              return(obj$benchmark_group)
            } else {
              # set
              obj$benchmark_group <- group
              return(obj)
            }
          }
)

setGeneric(name="benchmarkName",
           def=function(obj, benchmark = NULL)
           {
             standardGeneric("benchmarkName")
           }
)
setMethod(f="benchmarkName",
          signature="GenBench",
          definition=function(obj, benchmark = NULL)
          {
            if(is.null(benchmark)){
              # get
              return(obj$benchmark)
            } else {
              # set
              obj$benchmark <- benchmark
              return(obj)
            }
          }
)

# add data
setGeneric(name="addRecord",
           def=function(obj, record_name = "not_set", record = NULL)
           {
             standardGeneric("addRecord")
           }
)
setMethod(f="addRecord",
          signature="GenBench",
          definition=function(obj, record_name = "not_set", benchmark = NULL)
          {
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
)

# get data
setGeneric(name="getRecords",
           def=function(obj)
           {
             standardGeneric("getRecord")
           }
)
setMethod(f="getRecords",
          signature="GenBench",
          definition=function(obj)
          {
            return(obj$data)
          }
)

# get environment information
setGeneric(name="getEnv",
           def=function(obj)
           {
             standardGeneric("getEnv")
           }
)
setMethod(f="getEnv",
          signature="GenBench",
          definition=function(obj)
          {
            return(obj$env)
          }
)

# get filenames for writing captured results to
setGeneric(name="getOutputFile",
           def=function(obj)
           {
             standardGeneric("getOutputFile")
           }
)
setMethod(f="getOutputFile",
          signature="GenBench",
          definition=function(obj)
          {
            return(
              file.path("..", "..", "generated",
                        paste(benchmarkGroup(obj),benchmarkName(obj), "tsv", sep = '.')
              ))
          }
)

# check directory exists to write output files into
setGeneric(name="checkOutputFile",
           def=function(obj, create = FALSE)
           {
             standardGeneric("getOutputFile")
           }
)
setMethod(f="checkOutputFile",
          signature="GenBench",
          definition=function(obj, create = FALSE)
          {
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
)

setMethod(f="getOutputFile",
          signature="Timings",
          definition=function(obj)
          {
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
)

setMethod(f="getOutputFile",
          signature="Results",
          definition=function(obj)
          {
            return(
              file.path("..", "..", "generated", "results",
                        benchmarkGroup(obj),
                        paste(benchmarkName(obj), "results", "tsv", sep = '.')
              )
            )
          }
)


# reporting methods
setGeneric(name="reportRecords",
           def=function(obj)
           {
             standardGeneric("reportRecords")
           }
)
setMethod(f="reportRecords",
          signature="GenBench",
          definition=function(obj)
          {
            checkOutputFile(obj, create = TRUE)
            lapply(getRecords(obj), function(rec){
              # append each record to output file
              write.table(x = rec, file = getOutputFile(obj), append = TRUE,
                          sep = "\t", row.names = TRUE, col.names=TRUE
              )
            })
          }
)

setMethod(f="reportRecords",
          signature="Timings",
          definition=function(obj)
          {
            # replace RJSONIO with jsonlite
            require(RJSONIO)
            checkOutputFile(obj, create = TRUE)
            out <- getOutputFile(obj)
            # add 1 line JSON serialized header containing environment info
            writeLines(toJSON(getEnv(obj)), con = out)
            write.table(file=out, append = TRUE,
                        quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE,
                        format(do.call("rbind", getRecords(obj)), digits=5)

            )
          }
)

setMethod(f="reportRecords",
          signature="Results",
          definition=function(obj)
          {
            # replace RJSONIO with jsonlite
            require(RJSONIO)
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
)
