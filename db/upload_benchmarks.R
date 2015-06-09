# upload benchmarks to MySQL instance

library(RJDBC)
library(rjson)
library(reshape)

# Options
CREATE_NEW <- FALSE

### functions
# whitespace stripper
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}
# statement reader (returns list, one element per statement in file)
# read statements, dropping any empty lines (note: NO COMMENTS!)
readSQL <- function(path){
  statements <- lapply(
                      strsplit(split = ";", paste(readLines(path), collapse = " ")),
                      function(x) ifelse(trim(x)== "", NA, paste(x,";")))[[1]]
  statements <- statements[!is.na(statements)]
  return(statements)
}

### connect to mysql instance and check connection
# http://mvnrepository.com/artifact/mysql/mysql-connector-java/5.1.35
getConnection <- function(usr="foo", pwd="bar"){
  drv <- JDBC("com.mysql.jdbc.Driver",
              "mysql-connector-java-5.1.35.jar",
              identifier.quote="`"
              )
  conn <- dbConnect(drv, "jdbc:mysql://173.194.246.104/Rbenchmarks",
                    user=usr, password=pwd)
  return(conn)
}
conn <- getConnection()
dbListTables(conn)

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
    tmp <- read.delim(path, comment.char = "{")
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
lapply(reports, 
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
