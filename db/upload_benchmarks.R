# ieuan clay
# may 2015

### upload benchmarks to MySQL instance
# set libPath to local user dir
.libPaths(file.path("~","R","libs"))
# packages
library(RJDBC)
# Use jsonlite instead of RJSONIO
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
(conn <- getConnection(usr=CONN_INFO$`--usr`,
                      pwd=CONN_INFO$`--pwd`,
                      conn_string=CONN_INFO$`--conn`)
)
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
    # holder for data
    tmp <- NA
    # holder for header lines
    header <- readLines(path)
    # read metadata header if it exists
    if(length(grep(header[[1]], pattern = "^\\{"))){
      # cut header to lines enclosed in curly brackets
      header <- header[1:max(sapply(1:length(header), FUN=function(x) if(length(grep(header[[x]], pattern="\\}$"))) return(x) else 0))]
      meta <- fromJSON(paste(header, collapse = "\n"))
    } else {
      header <- list()
      meta <- list("platform","arch","os","system","status","major",
                   "minor","year","month","day","svn rev","language",
                   "version.string","nickname","sysname","release", "version"
      )
      names(meta) <- meta
      meta <- lapply(meta, function(x) "not set") # set all values to not set
    }

    # read data and parse file names
    try(tmp <- read.delim(path, comment.char = "{", skip=length(header)), silent = TRUE)
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
loaded <- lapply(reports[!is.na(reports)],
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
           return(NA)
         } else {
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
         # insert additional meta data from header (row by row for report)
         report$meta <- melt(
                        data.frame(report$meta, meta_id=m_id),
                        id.vars = "meta_id")
         lapply(1:nrow(report$meta), function(i){
           dbSendUpdate(conn,
                        # print(
                        sprintf(
                          "
                              INSERT INTO extra_meta (meta_id, variable, value)
                              VALUES (%i, \'%s\', \'%s\');
                            ",
                          report$meta[i,"meta_id"],
                          report$meta[i,"variable"],
                          as.character(report$meta[i,"value"])

                        )

           )}
         )

         # report completion
         cat(sprintf("Report successfully loaded from \'%s\'\n", report$meta$path))
         return(report)
         }
       }
)
# report loading
cat(sprintf("%i new reports loaded to database\n", sum(!is.na(loaded))))
cat(sprintf("%i reports already present in database\n", sum(is.na(loaded))))

### check results
cat(sprintf("DB now contains %i benchmarks.\n",
            dbGetQuery(conn, "
                                  SELECT count(distinct meta_id) as cnt
                                  FROM meta;")$cnt)
    )

## clean up and get ou
dbDisconnect(conn)
rm(list = ls())
gc()
