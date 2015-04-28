DATA_DIR <- "../../data/microarray"
INPUT <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45417/suppl/GSE45417_RAW.tar"
dir.create(DATA_DIR, recursive = TRUE)

download.file(INPUT,destfile = file.path(DATA_DIR, basename(INPUT), method="internal")

untar(file.path(DATA_DIR, basename(INPUT), exdir = DATA_DIR)

cel.files <- dir(DATA_DIR, pattern = ".CEL", full.names = TRUE)

pd <- data.frame(Name=basename(cel.files), FileName=cel.files, stringsAsFactors = FALSE)
pd$Group <- 
  do.call("rbind", strsplit(pd$Name, split = "[_\\.]", perl = TRUE))[,2]
pd$Treatment <- 
  sapply(
    strsplit(
      do.call("rbind", strsplit(pd$Name, split = "[_\\.]", perl = TRUE))[,3]
      ,""), function(x) paste(x[1:3],collapse = ''))
pd$Replicate <- 
  sapply(
    strsplit(
      do.call("rbind", strsplit(pd$Name, split = "[_\\.]", perl = TRUE))[,3]
      ,""), function(x) x[4])

# 
# write.table(pData, sep = "\t", col.names = TRUE, row.names=FALSE,
#             file = file.path(DATA_DIR, paste(basename(INPUT), "pData.txt", sep=".")))

library(affy) # reading and normalising microarray data
library(hgu133plus2cdf) # platform annotations
library(limma) # differential expression

# convert "phenodata" (i.e. sample metadata) into AnnotatedDataFrame (class expected by import functions)
pd <- new("AnnotatedDataFrame", data=pd)

# read affymetrix files
cel.files <- read.affybatch(phenoData = pd, filenames = pData(pd)$FileName)

# check RNA degredation and filter any low quality samples
cel.qc <- affy::AffyRNAdeg(cel.files, log.it = TRUE)
# drop samples with abnormal RNAdegredation slope
cel.files <- cel.files[,(cel.qc$slope >= 3 & cel.qc$slope <= 4.5)]

# extract expression values using
# robust multi-array average (RMA) method
# note: rma returns expression values in log2 scale,
# if using other expression meaasures, convert to log2 
# (limma expects data to be log2 transformed)
cel.rma <- rma(cel.files, normalize = TRUE, background = TRUE)
# scale to 150? need to do this before log transformation?

# the below is not nessecary, but just to demonstrate the other common
# expression method
cel.mas <- mas5(cel.files, normalize = TRUE, sc = 150)
exprs(cel.mas) <- log2(exprs(cel.mas))


# filter probes
# remove probes which do not have expression
# greater than 50 in 25% of the samples
filterfun <- function(eset, threshold=150, fraction=0.25, test_only=FALSE){ 
  # returns a logical vector, per probeset
  # TRUE : probeset is expressed at or above [threshold]
  # in ([fraction] * 100)% of the samples
  
  filter.vec <- apply(exprs(eset), MARGIN = 1, # row-wise (probesets)
                      function(x) (sum(x >= log2(threshold)) / length(x)) >= fraction )
  if(test_only){
    # just test filter settings, return nothing
    cat(sprintf(
      "%i probes, (%0.2f percent), would be retained.\n", 
      sum(filter.vec), (sum(filter.vec)/dim(exprs(eset))[1])*100
      ))
    return(sum(filter.vec))
  } else{
  return( 
      filter.vec
    )  
  }
}

# for simplicity continue only with MAS5 expression values
filterfun(cel.mas, threshold = 50, test_only = TRUE)
cel.filtered <- cel.mas[ filterfun(cel.mas, threshold = 50)  ,]

# differential expression
# use factorial design matrix to extract comparisons

# output results