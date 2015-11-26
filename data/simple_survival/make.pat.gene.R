
# Executes the preliminary data processing step and
# saves the results to an RDA file.

library(data.table)

mgsub2 <- function(myrepl, mystring){
  gsub2 <- function(l, x){
    do.call('gsub', list(x = x, pattern = l[1], replacement = l[2]))
  }
  Reduce(gsub2, myrepl, init = mystring, right = T)
}


# Load the required Data
d1 <- fread("TCGA_PRAD_RNAseq.txt", sep="\t", header=T)
setkey(d1, Gene)
pat <- fread("TCGA_PRAD_patient.txt", sep="\t", header=T) #converts all spaces and NA and unknown into NA for R, easier for sorting
setkey(pat, bcr_patient_barcode, name)


g1 <- "BAZ2A"
time <- "pfs"
gleason <- c("2+4", "3+3", "3+4", "3+5", "4+3", "4+4", "4+5", "5+3", "5+4", "5+5")
high <- (ncol(d1) - 2 ) * 0.75
low <- (ncol(d1) - 2) * 0.25
setkey(pat, gleason)
pat.d1 <- d1[,c("Gene", pat[gleason, name]), with=F]
setkey(pat.d1, Gene)
#selecting row for gene and only retrieving values
pat.d1.gene <- melt(pat.d1[g1, setdiff(colnames(pat.d1), "Gene"), with=F], id.vars = NULL, measure.vars = colnames(pat.d1)[-1], variable.name = "name", value.name="g1")
setkey(pat.d1.gene, g1)
pat.d1.gene[, name := factor(name, levels=name)]
pat.d1.gene[, ':=' (high = g1 > g1[eval(high)], low = g1 < g1[eval(low)])]
pat.d1.gene[, gene2 := high*2 + low]
pat.d1.gene[, gene3 := mgsub2(list(c("0", "middle"), c("1", "low"), c("2", "high")), pat.d1.gene$gene2)]
phenosgene <- merge(pat, pat.d1.gene, by= "name")
pat.gene <- phenosgene[gene2 !=0]
saveRDS(as.data.frame(pat.gene), "pat.gene.rda")
