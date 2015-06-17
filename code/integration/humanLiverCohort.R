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
library(stats)
library(e1071)
library(MASS) # rlm()
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
    
    # synapse dataset ids
    ids <- c(
      'syn4629', #phenotype
      'syn88644', #expression data
      'syn89614' #genotype data
      
      )
    ## download data
    for (id in ids){
      if(VERBOSE){
        cat(sprintf("Downloading synapse dataset %s", id))
      }
      synGet(id=id, load=F, downloadFile = TRUE, downloadLocation = DATA_DIR)
    }
  }
  
  ## unzip and load zipped files and remove temp unzipped files
  data <- lapply(dir(DATA_DIR, include.dirs = FALSE), function(zipfile){
    unzip(file.path(DATA_DIR, zipfile), exdir = DATA_DIR)
    data <- lapply(dir(DATA_DIR, pattern = ".txt$", full.names = TRUE), read.delim, stringsAsFactors=FALSE)
    names(data) <- dir(DATA_DIR, pattern = ".txt$", full.names = FALSE)
    file.remove(dir(DATA_DIR, pattern = ".txt$", full.names = TRUE))
    
    return(data)
  })
  names(data) <- strtrim(dir(DATA_DIR, include.dirs = FALSE), 11)
  
  ### reformat loaded data
  ## phenotype data
  # bind columns together (each row in file is one column)
  tmp <- data.frame(t(rbind(data$curatedPhen$"phenotype.txt")), stringsAsFactors = FALSE)
  names(tmp)<-tmp[1,] # first line in file is column names
  tmp <- tmp[2:nrow(tmp),] # drop header
  
  # focus on male caucasian patients
  if(VERBOSE){
    table(tmp$inferred_population, tmp$self_reported_ethnicity, tmp$GENDER)
  }
  tmp <- subset(tmp, GENDER=="Male" & inferred_population=="Cauc" & self_reported_ethnicity=="W")
  
  
  # convert numeric columns from character
  numeric_cols <- c(1, 6:7, 9:18)
  for (numeric_col in numeric_cols){
    tmp[,numeric_col] <- as.numeric(tmp[,numeric_col])
  }
  
  # complete cases only
  tmp <- tmp[complete.cases(tmp[,numeric_cols]),]
  tmp$indvidual_id <- rownames(tmp)
  
  # overwrite parent dirty data
  data$curatedPhen <- tmp
  
  ## expression data
  tmp <- data$curatedExpr$"expression.txt"
  # drop cases not in phenotype data
  tmp <- tmp[,names(tmp) %in% c("feature_id", data$curatedPhen$indvidual_id)]
  # merge on gene symbol
  tmp <- merge(tmp, data$curatedExpr$"features.txt"[,c("feature_id", "genesymbol")], by = "feature_id")
  # old skool split, apply combine to get mean expression per feature, per patient
  tmp <- split(tmp, as.factor(tmp$genesymbol))
  tmp <- lapply(tmp, function(df){
    if(nrow(df)==1){
      return(
        # if only one row then just return reformatted dataframe
        data.frame(
          # keep gene symbol as feature id
          feature_id=ifelse(df[1,"genesymbol"]=="", sprintf("feat_%.0f",df[1,"feature_id"]), df[1,"genesymbol"]),
          # for other features calculate 5% trimmed mean
          df[1,!names(df) %in% c("genesymbol", "feature_id")])
        
        )
    } else {
      return(
        
      data.frame(
        # keep gene symbol as feature id
        feature_id=ifelse(df[1,"genesymbol"]=="", sprintf("feat_%.0f",df[1,"feature_id"]), df[1,"genesymbol"]),
        # for other features calculate 5% trimmed mean per subject
        lapply(df[,!names(df) %in% c("genesymbol", "feature_id")], 
               function(x) mean(x, na.rm=TRUE, trim=0.025)))
      
      )
    }
  })
  tmp <- do.call("rbind", tmp)
  
  # remove invariant features and features with missing data
  # be conservative - later feature selection will remove further features
  tmp <- tmp[complete.cases(tmp),]
  tmp <- tmp[, c(
    TRUE, # keep feature ID
    # remove lowest quartile (by variance) of columns
    unlist(apply(tmp[,2:ncol(tmp)], 2, var, na.rm = TRUE)) 
      > quantile(unlist(apply(tmp[,2:ncol(tmp)], 2, var, na.rm = TRUE)), na.rm = TRUE, probs = 0.25, names = FALSE)[[1]]
    )]
  
  # pivot to standard form
  tmp <- list(hdr=as.character(tmp$feature_id), 
              mat=t(tmp[,names(tmp) != "feature_id"])
              )
  colnames(tmp$mat) <- tmp$hdr
  
  # overwrite uncleaned data
  data$curatedExpr <- tmp$mat
  
  ## genotype data
  tmp <- data$curatedGeno$"genotype.txt"
  # drop cases not in phenotype data
  tmp <- tmp[,names(tmp) %in% c("feature_id", data$curatedPhen$indvidual_id)]
  # drop incomplete features (rows)
  tmp <- tmp[complete.cases(tmp),]
  # drop invariate features (only one variant)
  # conservative - later feature selection will remove more
  tmp <- tmp[ 
            # keep rows (SNPs) with more than 1 type of call
            apply(tmp[, !names(tmp) %in% "feature_id"], 1, function(x){length(unique(x))}) > 1 
            , ]
  # recode to factors
  genotypes <- as.numeric(
                factor(unlist(tmp[,!names(tmp) %in% c("feature_id")]))
                )
  for(icol in 2:ncol(tmp)){ # skip feature id
    col <- names(tmp)[icol]
    tmp[,col] <- genotypes[ ((icol-2)*(nrow(tmp)) + 1):((icol-1)*(nrow(tmp))) ]
    
  }
  
  # pivot to standard form
  tmp <- list(hdr=as.character(tmp$feature_id), mat=t(tmp[,names(tmp) != "feature_id"]))
  colnames(tmp$mat) <- tmp$hdr
  
  # overwrite old data
  data$curatedGeno <- tmp$mat
  
  # return cleaned data
  return(data)
  
}

do.svm <- function(liverdata){
  ### run some simple predictive modelling on liver cohort clinical data
  
  results <- list() # placeholder
  
  ## SVM on age and liver stats using liver enzyme data
  # set test/train, sampling 1/3 rows as test
  traintest <- rep(TRUE, nrow(liverdata$curatedPhen))
  traintest[sample(1:length(traintest), size = round(nrow(liverdata$curatedPhen)/3), replace = FALSE)] <- FALSE
  
  train <- liverdata$curatedPhen[traintest,]
  test <- liverdata$curatedPhen[!traintest,]
  livercols <- 9:18
  model <- svm(x = as.matrix(train[,livercols]), # all liver enzme activiy stats 
               y=factor(train$"AGE_(YRS)" > 50), 
               scale = TRUE, type = "C")
  # test model on training set
  pred <- predict(model, as.matrix(train[,livercols]))
  res <- classAgreement(table(pred, factor(train$GENDER, levels = c("Male", "Female"))))
  results <- append(results, 
                    list(data.frame(
                      dat="agetrain",
                      var=names(res),
                      coeff=unlist(res)
                    ))
  )
  
  # and on the testset
  pred <- predict(model,as.matrix(test[,livercols]))
  res <- classAgreement(table(pred, test$GENDER))
  
  results <- append(results, 
                    list(data.frame(
                      dat="agetest",
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
  # TODO: calculate error
  
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

do.nb_expr <- function(liverdata){
  ### using naive bayes (assumes features are independent) to
  ### build classifier
  ### classification based on groups from liver activity
  
  results <- list() # placeholder
  
  ## build categorical classification of patients based on liver enzyme activities
  if(VERBOSE){
    ### explore potential groups of patients
    # feature correlations?
    heatmap(cor(
      liverdata$curatedPhen[rownames(liverdata$curatedExpr) ,strtrim(names(liverdata$curatedPhen), 3)=="CYP"])
    ) # 4 groups? but good, not all very high correlation
    # any obvious groups?
    heatmap(cor(
    t(liverdata$curatedPhen[rownames(liverdata$curatedExpr) ,strtrim(names(liverdata$curatedPhen), 3)=="CYP"])
    )) # ~ 5 groups, maybe only two are distinct
    pairs(prcomp(scale. = TRUE, center = TRUE,
      cor(t(liverdata$curatedPhen[rownames(liverdata$curatedExpr),strtrim(names(liverdata$curatedPhen), 3)=="CYP"]
      )))$x[,c(1,2,3)]) # only 2 distinct groups
    # try plotting groups in PCA plot
    pairs(
      prcomp(scale. = TRUE, center = TRUE,
                 cor(t(liverdata$curatedPhen[rownames(liverdata$curatedExpr),strtrim(names(liverdata$curatedPhen), 3)=="CYP"]
      )))$x[,c(1,2,3)],
      col=c("red", "green3", "blue", "black")[cutree(hclust(dist(cor(t(liverdata$curatedPhen[rownames(liverdata$curatedExpr) ,strtrim(names(liverdata$curatedPhen), 3)=="CYP"])))), k = 4)]
    )
    # and adding to heatmap
    heatmap(cor(
      t(liverdata$curatedPhen[rownames(liverdata$curatedExpr) ,strtrim(names(liverdata$curatedPhen), 3)=="CYP"])
    ),
    RowSideColors = c("red", "green3", "blue", "black")[cutree(hclust(dist(cor(t(liverdata$curatedPhen[rownames(liverdata$curatedExpr) ,strtrim(names(liverdata$curatedPhen), 3)=="CYP"])))), k = 4)]
    )
  
    # combine non-"red" groups into a separate group for classification
    grp <- c("red", "notred", "notred", "notred")[
      cutree(hclust(dist(cor(t(
        liverdata$curatedPhen[rownames(liverdata$curatedExpr), # forces row ordering to be the same in pheno data and expr data
                              strtrim(names(liverdata$curatedPhen), 3)=="CYP"]
        )))), k = 4)]
    # check that against actual activities and other pheno data
    pairs(
      liverdata$curatedPhen[rownames(liverdata$curatedExpr),strtrim(names(liverdata$curatedPhen), 3)=="CYP"],
      col=grp
      ) 
    # no obvious patterns, hmmm, maybe groups are not so great
    # also try to use classification targets as in paper (PLOS Biology: Mapping the Genetic Architecture of Gene ...)
    # top quartile of aldehyde activity
    
  }
  
  ## binary category based on aldehyde oxydase activity
  grps <- list()
  grps$ald <- c("lo", "hi")[
    cut(liverdata$curatedPhen[rownames(liverdata$curatedExpr) ,"aldehyde_oxydase"],
        breaks=c(-1,
                 quantile(probs = 0.75, 
                          liverdata$curatedPhen[rownames(liverdata$curatedExpr) ,"aldehyde_oxydase"]),
                 max(liverdata$curatedPhen[rownames(liverdata$curatedExpr) ,"aldehyde_oxydase"]) +1
        )
    )
    ]
  ## and on overall liver enzyme activity levels
  grps$enz <- c("red", "notred", "notred", "notred")[
    cutree(hclust(dist(cor(t(
      liverdata$curatedPhen[rownames(liverdata$curatedExpr), # forces row ordering to be the same in pheno data and expr data
                            strtrim(names(liverdata$curatedPhen), 3)=="CYP"]
    )))), k = 4)]
  
  ## split data into test and train
  # set test/train, sampling 1/3 rows as test
  traintest <- rep(TRUE, nrow(liverdata$curatedExpr))
  traintest[sample(1:length(traintest), size = round(nrow(liverdata$curatedExpr)/3), replace = FALSE)] <- FALSE
  
  train <- liverdata$curatedExpr[traintest,]
  test <- liverdata$curatedExpr[!traintest,]
  
  ## train model on each set of target classes
  results <- lapply(names(grps), function(grptype){
    grp <- grps[[grptype]]
    # do training
    model <- naiveBayes(x = train, y = factor(grp[traintest]))
    
    # check error on training data
    pred <- predict(model, train)
    trainres <- classAgreement(table(pred, grp[traintest]))
    
    # and on the testset
    pred <- predict(model,test)
    res <- classAgreement(table(pred, grp[!traintest]))
    
    return(      list(data.frame(
                        dat=paste("nbtest", grptype, sep = "-"),
                        var=names(res),
                        coeff=unlist(res)
                      ),
                      data.frame(
                        dat=paste("nbtrain", grptype, sep = "-"),
                        var=names(trainres),
                        coeff=unlist(trainres)
                      )
              )
    )
  })
  
  # return results
  return(do.call("rbind", unlist(results, recursive = F)))
}

do.rlm_expr <- function(liverdata){

  #"Expression trait processing. Expression traits were adjusted for age, sex, and medical center. Residuals were computed using rlm function from R statistical package (M-estimation with Tukey's bisquare weights). In examining the distributions of the mean log ratio measures for each expression trait in the HLC set, we noted a high rate of outliers. As a result, we used robust residuals and nonparametric tests to carry out the association analyses in the HLC. For each expression trait, residual values deviating from the median by more than three robust standard deviations were filtered out as outliers."
  # PLOS Biology: Mapping the Genetic Architecture of Gene ...
  # http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060107
  
  ### using robus linear models to build predictive model for
  ### Liver_Triglyceride levels
  
  results <- list() # placeholder
  
  ## split data into test and train
  # set test/train, sampling 1/3 rows as test
  traintest <- rep(TRUE, nrow(liverdata$curatedExpr))
  traintest[sample(1:length(traintest), size = round(nrow(liverdata$curatedExpr)/3), replace = FALSE)] <- FALSE
  
  # make test and train datasets
  train <- liverdata$curatedExpr[traintest,]
  test <- liverdata$curatedExpr[!traintest,]
  
  ## feature selection
  # simplest: most variable features
  feats <- apply(train, MARGIN = 2, var)
  if(VERBOSE){
    # ~waterfall plot
    # clear split between variable and non variable features
    plot(feats, rank(feats), col=factor(rank(feats) > ncol(train) - 50))
    # proceed with top variable features
  }
  feats <- (rank(feats) > (ncol(train) - 50))
  # exclude any exact linear combinations of variables (singularity causes rlm to fail)
  # see: http://stats.stackexchange.com/questions/70899/what-correlation-makes-a-matrix-singular-and-what-are-implications-of-singularit
  feats[feats] <- !(apply(abs(cor(train[,feats])) >= 0.75, MARGIN = 1, sum) > 1)  ## train model against triglyceride data
  if(VERBOSE){
    # make sure determinant is higher than 0
    det(cor(train[,feats]))
  }
  
  ## train model against triglyceride data
  # do training
  model <- rlm( 
                triglyc ~ ., 
                # adding response column as first column
                data = data.frame(triglyc=liverdata$curatedPhen[rownames(liverdata$curatedExpr)[traintest], "Liver_Triglyceride_(mg_per_dL)"],
                                                  train[,feats]),
#                x = train[,feats], 
#                y = liverdata$curatedPhen[rownames(liverdata$curatedExpr)[traintest], "Liver_Triglyceride_(mg_per_dL)"],
                method= "M", psi = psi.bisquare)
  
  # check error on training data
  results <- append(results, 
                    list(data.frame(
                      dat="rlmtrain-topvar",
                      var="cor",
                      coeff=cor(liverdata$curatedPhen[rownames(liverdata$curatedExpr)[traintest], "Liver_Triglyceride_(mg_per_dL)"], 
                                predict(model, data.frame(train[,feats])))
                    ))
  )
  
  # and on the testset
  results <- append(results, 
                    list(data.frame(
                      dat="rlmtest-topvar",
                      var="cor",
                      coeff=cor(liverdata$curatedPhen[rownames(liverdata$curatedExpr)[!traintest], "Liver_Triglyceride_(mg_per_dL)"], 
                                predict(model, data.frame(test[,feats])))
                    ))
  )


  ### results above by heuristic feature selection are really poor,
  # so try feature selection by sampling

  bestcor <- 0
  set.seed(8008)
  for(i in 1:50){
    if(VERBOSE){
      cat(sprintf("starting iteration %i for rlm()\n", i))
    }
    feats <- apply(train, MARGIN = 2, var)
    feats <- (rank(feats) > (ncol(train) - 1000)) # start from top 1000 by variablility
    subfeats <- 1:length(feats)
    subfeats <- subfeats[feats] # all top features
    subfeats <- names(feats[sample(subfeats, 25, replace = F)]) # sample
    # exclude any exact linear combinations of variables (singularity causes rlm to fail)
    subfeats <- subfeats[!(apply(abs(cor(train[,subfeats])) >= 0.75, MARGIN = 1, sum) > 1)]  ## train model against triglyceride data
    if(VERBOSE){
      # make sure determinant is higher than 0
      det(cor(train[,subfeats]))
    }
    
    ## train model against triglyceride data
    # do training with ERROR HANDLING
    possibleError <- tryCatch(error=function(e) e,
      model <- rlm( maxit=50,
        triglyc ~ ., 
        # adding response column as first column
        data = data.frame(triglyc=liverdata$curatedPhen[rownames(liverdata$curatedExpr)[traintest], "Liver_Triglyceride_(mg_per_dL)"],
                          train[,subfeats]),
        #                x = train[,feats], 
        #                y = liverdata$curatedPhen[rownames(liverdata$curatedExpr)[traintest], "Liver_Triglyceride_(mg_per_dL)"],
        method= "M", psi = psi.bisquare))
    # skip if model failed to converge
    if(inherits(possibleError, "error")){
      if(VERBOSE){
        cat(sprintf("\titeration %i failed to converge\n", i))
      }
      next
    }
    
    # check error on training data
    tmp_res <- cor(liverdata$curatedPhen[rownames(liverdata$curatedExpr)[traintest], "Liver_Triglyceride_(mg_per_dL)"], 
                   predict(model, data.frame(train[,subfeats])))
    # stop now if the result is no better than the previous best
    # no point testing on test set if training set is not better than previous best
    if(tmp_res < bestcor){
      next
    }
    # and on the testset
    tmp_res <- cor(liverdata$curatedPhen[rownames(liverdata$curatedExpr)[!traintest], "Liver_Triglyceride_(mg_per_dL)"], 
                   predict(model, data.frame(test[,subfeats])))
    
    if(tmp_res > bestcor){
      bestcor <- tmp_res
    }
    
  }  
  # append the best result so far
  results <- append(results, 
                    list(data.frame(
                      dat="rlm-sample",
                      var="cor",
                      coeff=bestcor
                    ))
  )
  
  # return results
  return(do.call("rbind",results))
  
}

### reporting
# load data
TIMES <- addRecord(TIMES, record_name = "dataload",
                   record = system.time(gcFirst = T,
                                        liverdata <- do.load(DATA_DIR=DATA_DIR, DOWNLOAD=DOWNLOAD)
                                        ))
# some simple predictions on basic data
TIMES <- addRecord(TIMES, record_name = "svm-pheno",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="svm-pheno",
                                                             record=do.svm(liverdata)
                                        )))

# can expression data predict a class of liver enzyme activity?
TIMES <- addRecord(TIMES, record_name = "naivebayes-expr",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="naivebayes-expr",
                                                             record=do.nb_expr(liverdata)
                                        )))

# can expression data predict liver enzyme activity?
TIMES <- addRecord(TIMES, record_name = "rlm-expr",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="rlm-expr",
                                                             record=do.rlm_expr(liverdata)
                                        )))

## output results for comparison
# write results to file
reportRecords(RESULTS)

# timings
reportRecords(TIMES)

# final clean up
rm(list=ls())
gc()

