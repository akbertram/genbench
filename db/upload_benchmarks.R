# upload benchmarks to MySQL instance
# mysql -uieuan -pbiolion -h173.194.246.104 Rbenchmarks

library(RJDBC)
library(rjson)

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
drv <- JDBC("com.mysql.jdbc.Driver",
            "mysql-connector-java-5.1.35.jar",
            identifier.quote="`"
            )
conn <- dbConnect(drv, "jdbc:mysql://173.194.246.104/Rbenchmarks",
                  user="ieuan", password="biolion")

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
### push data into metatable
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
    
    # add meta data to results
    tmp$sysname <- meta$sysname
    tmp$release <- meta$release
    tmp$lang <- meta$language
    tmp$major <- meta$major
    tmp$minor <- meta$minor
    
    return(tmp)
  }
)

LAST_INSERT_ID()



### push data into timings table
