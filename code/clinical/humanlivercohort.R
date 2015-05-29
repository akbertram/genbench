# human liver cohort
# synapse account
# user=synapsedata
# pwd=synapse

source('http://depot.sagebase.org/CRAN.R')
pkgInstall(c("synapseClient"))

library(synapseClient)
synapseLogin('synapsedata','synapse')

# Obtain a pointer and download the data
DATA_DIR <- "../../data/integration/hlc"
dir.create(DATA_DIR, recursive = T)

## phenotype
synGet(id='syn4629', load=F, downloadFile = TRUE, downloadLocation = DATA_DIR)
unzip("../../data/integration/hlc//curatedPhenotype.zip", exdir = DATA_DIR)
data <- lapply(dir(DATA_DIR, pattern = ".txt$", full.names = TRUE), read.delim, stringsAsFactors=FALSE)
names(data) <- dir(DATA_DIR, pattern = ".txt$", full.names = FALSE)

file.remove(dir(DATA_DIR, pattern = ".txt$", full.names = TRUE))

str(data)

tmp <- data.frame(t(rbind(data$phenotype.txt)), stringsAsFactors = FALSE)
tmp$indvidual_id <- rownames(tmp)
names(tmp)<-tmp[1,]
tmp <- tmp[2:nrow(tmp),]
numeric_cols <- c(1, 5:18)
for (numeric_col in numeric_cols){
  tmp[,numeric_col] <- as.numeric(tmp[,numeric_col])
  tmp[is.na(tmp[,numeric_col]), numeric_col] <- 0 # replace all NAs with 0s, i know this is horrible.
}


library(e1071)
# set test/train
traintest <- rep(TRUE, nrow(tmp))
traintest[sample(1:length(traintest), size = round(nrow(tmp)/3), replace = FALSE)] <- FALSE

train <- tmp[traintest,]
test <- tmp[!traintest,]
model <- svm(x = as.matrix(train[,numeric_cols]), 
             y=factor(train$GENDER, levels = c("Male", "Female")), 
             scale = TRUE, type = "C")
# test model on training set
pred <- predict(model, as.matrix(train[,numeric_cols]))
table(pred, factor(train$GENDER, levels = c("Male", "Female")))

pred <- predict(model,as.matrix(test[,numeric_cols]))

table(pred, test$GENDER) # turns out liver stats suck at predicting gender :)
