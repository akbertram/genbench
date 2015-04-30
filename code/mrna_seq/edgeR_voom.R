library(limma)
library(edgeR)

# http://www.bioconductor.org/packages/release/bioc/html/edgeR.html
# http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# http://bioinf.wehi.edu.au/RNAseqCaseStudy/

INPUT <- "http://bowtie-bio.sourceforge.net/recount/countTables/gilad_count_table.txt"
PDATA <- "http://bowtie-bio.sourceforge.net/recount/phenotypeTables/gilad_phenodata.txt"

### import data
# read data into data.frames
counts <- read.delim(INPUT)
row.names(counts) <- counts$gene
counts <- counts[,2:ncol(counts)]
pd <- read.delim(PDATA, sep = " ", stringsAsFactors=FALSE)

# convert data to DGEList instance
dge <- DGEList(counts=counts, group = pd$gender)

# calculate normalisation factors
dge <- calcNormFactors(dge)

### filter non-detected features
filterfun <- function(dgelist, threshold=1, fraction=0.25, test_only=FALSE){ 
  # returns a logical vector, per row
  # TRUE : row (feature) is expressed at or above [threshold]
  # in ([fraction] * 100)% of the samples
  
  filter.vec <- apply(dgelist$counts, MARGIN = 1, # row-wise (probesets)
                      function(x) (sum(x >= threshold) / length(x)) >= fraction )
  if(test_only){
    # just test filter settings, return nothing
    cat(sprintf(
      "Threshold : %i, Fraction ; %0.2f \n\t\t= %i features, (%0.2f percent), would be retained.\n", 
      threshold, fraction, sum(filter.vec), (sum(filter.vec)/dim(dgelist$counts)[1])*100
    ))
    return(sum(filter.vec))
  } else{
    return( 
      filter.vec
    )  
  }
}

# for simplicity continue only with MAS5 expression values
filterfun(dge, threshold = 1, fraction = 0.1, test_only = TRUE) # i.e. keep if any signal seen in any sample
dge <- dge[filterfun(dge, threshold = 1, fraction = 0.1),]

### normalise
design <- model.matrix(~0+dge$samples$group) # levels = F, M, i.e. female is reference
colnames(design) <- levels(dge$samples$group)

vm <- voom(dge,design,plot=FALSE)

### fit linear model
fit <- lmFit(vm,design)
cont.matrix <- makeContrasts(MvF=M-F, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

### export results
topTable(fit2, adjust="BH")
