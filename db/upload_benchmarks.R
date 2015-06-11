# ieuan clay
# may 2015

### upload benchmarks to MySQL instance
rm(list = ls())
# set libPath to local user dir
.libPaths("~/R/libs")
# packages
library(RJDBC)
library(rjson)
library(reshape)

# Options
CREATE_NEW <- FALSE # default: do not recreate all tables
CONN_INFO <- NA
# check args to see if we should use database connection
args <- commandArgs(trailingOnly = TRUE)
# reset timings?
if ("--create-new" %in% args){
  CREATE_NEW <<- TRUE
}
# check for other required arguments
conn_info <- args[args!="--create-new"]
# split characters
conn_info <- lapply(conn_info, function(x){strsplit(x, split = "=")[[1]]})
# separate into key:value named list for ease of access
names(conn_info) <- lapply(conn_info, function(x) x[[1]])
conn_info <- lapply(conn_info, function(x) x[[2]])

#   if(!"--driver" %in% names(conn_info)){
#     # default: use mysql driver
#     conn_info$`--driver` <- file.path(getwd(), dir("db", pattern = "mysql-connector-java-5.1.35.jar", full.names = TRUE))
#   }

if(all(c("--usr","--pwd","--conn") %in% names(conn_info))){
  # all parameters are provided, so proceed with using database
  # todo: more thorough testing of connection
  
  CONN_INFO <<- conn_info
  USE_DB <<- TRUE
  
  cat("Working database credentials supplied, using database for examine_benchmarks.R")
  
} else {
  
  cat("Working database credentials NOT supplied, using local files for examine_benchmarks.R")
  
}



### functions
source("sql_utilities.R")

### connect to DB
conn <- getConnection(usr=CONN_INFO$`--usr`, 
                      pwd=CONN_INFO$`--pwd`, 
                      conn_string=CONN_INFO$`--conn`)
#dbListTables(conn)

### create tables:
## meta table to hold metadata about each benchmark run
## timings table to hold results
if(CREATE_NEW){
  # empty DB
  lapply(readSQL("drop_tables.sql"), 
         function(statement) dbSendUpdate(conn, statement))
  
  # create new tables to hold data
  lapply(readSQL("create_tables.sql"), 
         function(statement) dbSendUpdate(conn, statement))
}

### collect data from existing reports
reports <- lapply(
  dir(path = file.path("..", "generated", "timings"), pattern = ".tsv$",
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
    tmp <- NA
    try(tmp <- read.delim(path, comment.char = "{"))
    if(is.na(tmp)){return(NA)}
    tmp$block <- rownames(tmp)
    # melt to tall and skinny for plotting
    tmp <- melt(tmp, id.vars = c("block"), na.rm = TRUE)
    suppressWarnings(tmp$value <- as.numeric(tmp$value))
    tmp <- tmp[!is.na(tmp$value),] # drop empty timings
    
    # parse path to group benchmarks
    # add additional meta data
    meta$path <- path
    meta$id <- strtrim(basename(path), nchar(basename(path)) - 4) # drop file extension
    meta$benchmark_group <- basename(dirname(path))
    meta$benchmark <- strsplit(meta$id, '\\.')[[1]][[1]]
    
    meta$time <- strsplit(meta$id, '\\.')[[1]][[2]]
    
    return(list(df=tmp, meta=meta))
  }
)


### push data into meta and timings table
lapply(reports[!is.na(reports)], 
       function(report){
         # check if this report has already been inserted
         if (nrow(dbGetQuery(conn, sprintf("
                                  SELECT * 
                                  FROM meta 
                                  WHERE 
                                    benchmark=\'%s\' AND 
                                    benchmark_group=\'%s\' AND 
                                    timestamp=\'%s\'",
                                  report$meta$benchmark, report$meta$benchmark_group, report$meta$time
                                  ))) >=1){
           
           # report completion
           cat(sprintf("Report already loaded from \'%s\'\n", report$meta$path))
           return(report)
         }
         # insert metadata
         dbSendUpdate(conn,
                      sprintf(
                        "
                          INSERT INTO meta (timestamp, benchmark, benchmark_group, sys_name, sys_release, lang, lang_major, lang_minor)
                          VALUES ( \'%s\', \'%s\', \'%s\', \'%s\', \'%s\', \'%s\', \'%s\', \'%s\' );
                        ", 
                            report$meta$time, report$meta$benchmark, report$meta$benchmark_group,
                            report$meta$sysname, report$meta$release, 
                            report$meta$lang, report$meta$major, report$meta$minor
                        )
                      )
         # retrieve insert meta id
         m_id <- dbGetQuery(conn, "SELECT LAST_INSERT_ID() AS meta_id;")

         # insert timings (row by row for report)
         lapply(1:nrow(report$df), function(i){
              dbSendUpdate(conn,
                  # print(
                          sprintf(
                            "
                              INSERT INTO timings (meta_id, block, variable, value)
                              VALUES (%i, \'%s\', \'%s\', %.3f);
                            ",
                            m_id$meta_id, report$df[i,"block"], 
                            report$df[i,"variable"], report$df[i,"value"]
                            
                            )
               
               )}
         )
         
         # report completion
         cat(sprintf("Report successfully loaded from \'%s\'\n", report$meta$path))
         return(report)
  }
)

### check results
dbGetQuery(conn, "SELECT * FROM timings t JOIN meta m ON t.meta_id=m.meta_id WHERE t.variable=\"elapsed\";")
