# Genetics
# ieuan.clay@gmail.com
# April 2015

### set up session
rm(list=ls())
## packages
library(stats)
library(utils)

## global vars
VERBOSE <- TRUE # print progress?

# holder for results
RESULTS <- list()
TIMES <- list()

# reproducibility
set.seed(8008)
stopifnot(all(rev(strsplit(getwd(), "/")[[1]])[1:3] == c("mutation","code","genbase"))) # must be run from directory containing rppa.R

### functions

## mutation data
# prepare for downloading
if(!file.exists("../../data/mutation_pop")){
  # create directory for holding data
  dir.create(recursive = TRUE, "../../data/mutation_pop")
}
## supplementary data described here:
# http://www.nejm.org/doi/suppl/10.1056/NEJMoa1301689/suppl_file/nejmoa1301689_appendix.pdf
#  supp table S6: All somatic mutations with annotation and read counts from DNA and RNA sequencing
download.file(destfile ="../../data/mutation_pop/laml.maf", "http://tcga-data.nci.nih.gov/docs/publications/laml_2012/SupplementalTable06.tsv")
readLines("../../data/mutation_pop/laml.maf", 3) # no header info other than col names

maf <- read.delim("../../data/mutation_pop/laml.maf", blank.lines.skip = TRUE, , stringsAsFactors = FALSE)

## sample data
download.file(destfile ="../../data/mutation_pop/laml.meta.tsv", 
              "http://tcga-data.nci.nih.gov/docs/publications/laml_2012/clinical_patient_laml.tsv")
meta <- read.delim("../../data/mutation_pop/laml.meta.tsv", blank.lines.skip = TRUE, stringsAsFactors = FALSE)

stopifnot(length(intersect(unique(maf$TCGA_id), unique(meta$bcr_patient_barcode))) == 197)

## figure 1A: mutations per sample, split by mutation tier and disease status
# split-apply-combine
fig_1a <- do.call("rbind", 
                  lapply(
                    split(maf, 
                          f = maf$TCGA_id), 
                    function(df){
                      tmp <- data.frame(table(df$tier))
                      names(tmp) <- c("tier", "mut_count")
                      tmp$TCGA_id <- df[1,"TCGA_id"]
                      return(tmp)})
                  )
fig_1a <- merge(x=fig_1a, y=meta, 
                by.x = "TCGA_id", by.y = "bcr_patient_barcode", 
                all.x = TRUE)
# subset data
subset(fig_1a, tier=="tier1")
# plot
plot(y=fig_1a$mut_count, 
     x = fig_1a$acute_myeloid_leukemia_calgb_cytogenetics_risk_category)

## figure 1B: samples per mutated gene
fig_1b <- do.call("rbind", 
                  lapply(
                    split(maf, 
                          f = maf$gene_name), 
                    function(df){
                      tmp <- data.frame(table(df$tier))
                      names(tmp) <- c("tier", "sample_count")
                      tmp$gene_name <- df[1,"gene_name"]
                      return(tmp)})
)
fig_1b <- fig_1b[ order(fig_1b$sample_count, decreasing = TRUE) ,]
fig_1b <- head(subset(fig_1b, tier == "tier1"), n = 100) # top 100 tier 1 genes by count

## familial data
# prepare for downloading
if(!file.exists("../../data/mutation_fam")){
  # create directory for holding data
  dir.create(recursive = TRUE, "../../data/mutation_fam")
}
fam_files <- data.frame(do.call("rbind", unlist(lapply(strsplit(split = "\n",
"Daniel MacArthur | DGM001 | http://s3.amazonaws.com/gnz.genotypes/DGM001_genotypes.zip
Luke Jostins | LXJ001 | http://s3.amazonaws.com/gnz.genotypes/LXJ001_genotypes.zip
Dan Vorhaus | DBV001 | http://s3.amazonaws.com/gnz.genotypes/DBV001_genotypes.zip
Caroline Wright | CFW001 | http://s3.amazonaws.com/gnz.genotypes/CFW001_genotypes.zip
Kate Morley | KIM001 | http://s3.amazonaws.com/gnz.genotypes/KIM001_genotypes.zip
Vincent Plagnol | VXP001 | http://s3.amazonaws.com/gnz.genotypes/VXP001_genotypes.zip
Jeff Barrett | JCB001 | http://s3.amazonaws.com/gnz.genotypes/JCB001_genotypes.zip
Jan Aerts | JXA001 | http://s3.amazonaws.com/gnz.genotypes/JXA001_genotypes.zip
Joe Pickrell | JKP001 | http://s3.amazonaws.com/gnz.genotypes/JKP001_genotypes.zip
Don Conrad | DFC001 | http://s3.amazonaws.com/gnz.genotypes/JKP001_genotypes.zip
Carl Anderson | CAA001 | http://s3.amazonaws.com/gnz.genotypes/CAA001_genotypes.zip
Ilana Fisher | IPF001 | http://s3.amazonaws.com/gnz.genotypes/IPF001_genotypes.zip")[[1]], 
function(x) strsplit(x, " | ", fixed = TRUE)), recursive=FALSE)), stringsAsFactors=FALSE)

names(fam_files) <- c("member","dataset id","link")
fam_files$target <- basename(fam_files$link)

lapply(1..nrow(fam_files), function(x) download.file())


