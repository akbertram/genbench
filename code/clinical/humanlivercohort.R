# human liver cohort
# synapse account
# user=synapsedata
# pwd=synapse
# ieuan.clay@gmail.com
# June 2015

### set up session
rm(list=ls())

# reproducibility
set.seed(8008)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

## packages
# CRAN
library(e1071)
# Bioc

## global vars
VERBOSE <- FALSE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "hlc"
DATA_DIR <- file.path("..","..","data","integration","hlc")

# holder for results
RESULTS <- results(benchmark_name = BENCHMARK)
TIMES <- timings(benchmark_name = BENCHMARK)

### functions
## Utility

## Blocks for timing
do.load <-function(DATA_DIR, DOWNLOAD=FALSE){
  ### download (optional) data and load into R
  # expects a directory containing downloaded data
  
  if(DOWNLOAD){
    ### load data from synapse repository
    # http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060107
    # https://www.synapse.org/#!Synapse:syn4499
    
    # load synapse package and check supplied data directory
    source('http://depot.sagebase.org/CRAN.R')
    pkgInstall(c("synapseClient"))
    
    library(synapseClient)
    synapseLogin('synapsedata','synapse')
    
    # Obtain a pointer and download the data
    if(!file.exists(DATA_DIR)){
      dir.create(DATA_DIR, recursive = T)
    }
    
    ## load phenotype data
    synGet(id='syn4629', load=F, downloadFile = TRUE, downloadLocation = DATA_DIR)
    
  }
  
  ## load zipped file and remove temp unzipped files
  unzip(file.path(DATA_DIR, "curatedPhenotype.zip"), exdir = DATA_DIR)
  data <- lapply(dir(DATA_DIR, pattern = ".txt$", full.names = TRUE), read.delim, stringsAsFactors=FALSE)
  names(data) <- dir(DATA_DIR, pattern = ".txt$", full.names = FALSE)
  
  file.remove(dir(DATA_DIR, pattern = ".txt$", full.names = TRUE))
  
  # check loaded data
  if(VERBOSE){cat(str(data))}
  
  ## reformat loaded data
  # bind columns together (each row in file is one column)
  tmp <- data.frame(t(rbind(data$phenotype.txt)), stringsAsFactors = FALSE)
  tmp$indvidual_id <- rownames(tmp)
  names(tmp)<-tmp[1,] # first line in file is column names
  tmp <- tmp[2:nrow(tmp),] # drop header
  
  # convert numeric columns from character
  numeric_cols <- c(1, 5:18)
  for (numeric_col in numeric_cols){
    tmp[,numeric_col] <- as.numeric(tmp[,numeric_col])
    tmp[is.na(tmp[,numeric_col]), numeric_col] <- 0 # replace all NAs with 0s, i know this is horrible.
  }
  
  return(tmp)
  
}

do.svm <- function(liverdata){
  ### run some simple predictive modelling on liver cohort clinical data
  # see integration/humanLiverCohort.R for more advanced calculations
  
  results <- list() # placeholder
  
  ## SVM on gender and liver stats using liver enzyme data
  # set test/train, sampling 1/3 rows as test
  traintest <- rep(TRUE, nrow(liverdata))
  traintest[sample(1:length(traintest), size = round(nrow(liverdata)/3), replace = FALSE)] <- FALSE
  
  train <- liverdata[traintest,]
  test <- liverdata[!traintest,]
  livercols <- 9:18
  model <- svm(x = as.matrix(train[,livercols]), # all liver enzme activiy stats 
               y=factor(train$GENDER, levels = c("Male", "Female")), 
               scale = TRUE, type = "C")
  # test model on training set
  pred <- predict(model, as.matrix(train[,livercols]))
  res <- classAgreement(table(pred, factor(train$GENDER, levels = c("Male", "Female"))))
  results <- append(results, 
                    list(data.frame(
                      dat="gendertrain",
                      var=names(res),
                      coeff=unlist(res)
                    ))
  )
  
  # and on the testset
  pred <- predict(model,as.matrix(test[,livercols]))
  res <- classAgreement(table(pred, test$GENDER))
  
  results <- append(results, 
                    list(data.frame(
                      dat="gendertest",
                      var=names(res),
                      coeff=unlist(res)
                    ))
  )
  
  ## liver triglycerides
  model <- svm(x = as.matrix(train[,livercols]), # all liver enzme activiy stats 
               y=train$`Liver_Triglyceride_(mg_per_dL)`, 
               scale = TRUE, type = "eps-regression")
  # test model on training set
  train$pred <- predict(model, as.matrix(train[,livercols]))
  # and on the testset
  test$pred <- predict(model,as.matrix(test[,livercols]))
  
  # return results
  results <- append(results, 
                    list(data.frame(
                      dat=c("triglyceridesTrain", "triglyceridesTest"),
                      var="spearman",
                      coeff=c(
                        cor(x=train$`Liver_Triglyceride_(mg_per_dL)`, y=train$pred, method = "spearman"),
                        cor(x=test$`Liver_Triglyceride_(mg_per_dL)`, y=test$pred, method = "spearman")
                        )
                    ))
  )
  
  return(do.call("rbind",results))
  
  
}

### reporting
# load data
liverdata <- do.load(DATA_DIR=DATA_DIR, DOWNLOAD=DOWNLOAD)

# some simple predictions on basic data
TIMES <- addRecord(TIMES, record_name = "svm",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="svm",
                                                             record=do.svm(liverdata)
                                        )))
                   
## output results for comparison
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()

