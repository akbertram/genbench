
require(RJDBC)

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
getConnection <- function(usr="foo", pwd="bar", 
                          driver=NULL,
                          conn_string="jdbc:mysql://173.194.246.104/Rbenchmarks"){
  if(is.null(driver)){
    # the default driver is the local mysql.jar
    # this part is required such that sourcing this script from any dir
    # will still correctly load the driver, without having to specify it each time
    # NB: if the user specifies a custom driver, they must pass the full path!
    tryCatch(driver <- file.path(
                          getwd(), 
                          dir(pattern = "mysql-connector-java.*.jar", 
                                     full.names = TRUE, recursive = TRUE, 
                                     include.dirs = FALSE, no.. = TRUE)[[1]]
                        ),
             error = function(e) cat(sprintf("Driver could not be found: %s\n", e)))
  }
  print(driver)
  drv <- JDBC("com.mysql.jdbc.Driver",
              driver,
              identifier.quote="`"
  )
  conn <- dbConnect(drv, conn_string,
                    user=usr, password=pwd)
  return(conn)
}
