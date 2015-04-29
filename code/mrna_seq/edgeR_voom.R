library(limma)
library(edgeR)

# http://www.bioconductor.org/packages/release/bioc/html/edgeR.html
# 
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


### normalise


### fit linear model


### export results