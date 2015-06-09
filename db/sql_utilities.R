
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
                          jdbcdriver="mysql-connector-java-5.1.35.jar",
                          connectionstring="jdbc:mysql://173.194.246.104/Rbenchmarks"){
  drv <- JDBC("com.mysql.jdbc.Driver",
              jdbcdriver,
              identifier.quote="`"
  )
  conn <- dbConnect(drv, connectionstring,
                    user=usr, password=pwd)
  return(conn)
}