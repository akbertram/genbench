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
DOWNLOAD <- FALSE # download fresh data?

# holder for results
RESULTS <- list()
TIMES <- list()
BENCHMARK <- "mutation"

# reproducibility
set.seed(8008)
stopifnot(all(rev(strsplit(getwd(), "[\\\\/]", perl=TRUE)[[1]])[1:3] == c("mutation","code","genbase"))) # must be run from directory containing rppa.R

### functions
## general
do.download <- function(n=3){
  # pop data
  DATA_DIR <- file.path("..", "..", "data", "mutation")
  
  # prepare for downloading
  if(!file.exists(file.path(DATA_DIR, "pop"))){
    # create directory for holding data
    dir.create(recursive = TRUE, file.path(DATA_DIR, "pop"))
  }
  ## supplementary data described here:
  # http://www.nejm.org/doi/suppl/10.1056/NEJMoa1301689/suppl_file/nejmoa1301689_appendix.pdf
  #  supp table S6: All somatic mutations with annotation and read counts from DNA and RNA sequencing
  download.file(destfile =file.path(DATA_DIR, "pop","laml.maf"), "http://tcga-data.nci.nih.gov/docs/publications/laml_2012/SupplementalTable06.tsv")
  
  download.file(destfile =file.path(DATA_DIR, "pop","laml.meta.tsv"), 
                "http://tcga-data.nci.nih.gov/docs/publications/laml_2012/clinical_patient_laml.tsv")
  
  
  # 'family' data
  
  ## download data for individuals from genomesunzipped.org
  # n : limit number of individuals loaded
  
  # prepare for downloading
  if(!file.exists(file.path(DATA_DIR, "fam"))){
    # create directory for holding data
    dir.create(recursive = TRUE, file.path(DATA_DIR, "fam"))
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
  fam_files$target <- file.path(DATA_DIR, "fam", basename(fam_files$link))
  
  # only take n files: could expand later if needed
  if(n > nrow(fam_files)){ n <- nrow(fam_files)  }
  fam_files <- fam_files[1:n,]
  
  # download and unzip target files
  lapply(1:nrow(fam_files), function(x){
    download.file(fam_files[x,"link"], 
                  destfile = fam_files[x,"target"])
    unzip(fam_files[x,"target"], exdir = file.path(DATA_DIR, "fam"))
    try(file.remove(fam_files[x,"target"]))
    })
  return(TRUE)
}

## mutation data

do.pop.load <- function(){
  
  readLines(file.path(DATA_DIR, "pop","laml.maf"), 3) # no header info other than col names
  
  maf <- read.delim(file.path(DATA_DIR, "pop","laml.maf"), blank.lines.skip = TRUE, , stringsAsFactors = FALSE)
  
  ## sample data
  meta <- read.delim(file.path(DATA_DIR, "pop","laml.meta.tsv"), blank.lines.skip = TRUE, stringsAsFactors = FALSE)
  
  stopifnot(length(intersect(unique(maf$TCGA_id), unique(meta$bcr_patient_barcode))) == 197)

  return(list(maf=maf, meta=meta))
}

do.pop.fig1a <- function(pop){
  ## figure 1A: mutations per sample, split by mutation tier and disease status
  # split-apply-combine
  proc.time() -> ptm
  fig_1a <- do.call("rbind", 
                    lapply(
                      split(pop$maf, 
                            f = pop$maf$TCGA_id), 
                      function(df){
                        tmp <- data.frame(table(df$tier))
                        names(tmp) <- c("tier", "mut_count")
                        tmp$TCGA_id <- df[1,"TCGA_id"]
                        return(tmp)})
                    )
  TIMES$pop.fig1a.sac <<- proc.time() - ptm
  
  # merge metadata
  proc.time() -> ptm
  fig_1a <- merge(x=fig_1a, y=pop$meta, 
                  by.x = "TCGA_id", by.y = "bcr_patient_barcode", 
                  all.x = TRUE)
  TIMES$pop.fig1a.merge <<- proc.time() - ptm
  # subset data
  fig_1a <- subset(fig_1a, tier=="tier1")
  # plot
  plot(y=fig_1a$mut_count, 
       x = factor(fig_1a$acute_myeloid_leukemia_calgb_cytogenetics_risk_category)
       )
  
  # return
  return(fig_1a[,c("TCGA_id", "tier", "mut_count")])
}

do.pop.fig1b <- function(pop){
  ## figure 1B: samples per mutated gene
  proc.time() -> ptm  
  fig_1b <- do.call("rbind", 
                    lapply(
                      split(pop$maf, 
                            f = pop$maf$gene_name), 
                      function(df){
                        tmp <- data.frame(table(df$tier))
                        names(tmp) <- c("tier", "sample_count")
                        tmp$gene_name <- df[1,"gene_name"]
                        return(tmp)})
  )
  
  # reorder
  proc.time() -> ptm
  fig_1b <- fig_1b[ order(fig_1b$sample_count, decreasing = TRUE) ,]
  fig_1b <- head(subset(fig_1b, tier == "tier1"), n = 100) # top 100 tier 1 genes by count
  
  TIMES$pop.fig1b.order <<- proc.time() - ptm
  
  # return
  return(fig_1b[,c("gene_name", "tier", "sample_count")])
}

# "window" apply
# modified from: http://www.r-bloggers.com/wapply-a-faster-but-less-functional-rollapply-for-vector-setups/
wapply <- function(x, width, by = NULL, FUN = NULL, ...){
  FUN <- match.fun(FUN)
  if (is.null(by)) by <- width
  
  starts <- seq(1, length(x) - width + 1, by = by)
  indices <- lapply(starts, function(x) x:(x + width - 1))
  
  windows <- base:::simplify2array(
    lapply(indices, function(a) FUN(x[a], ...)),
    higher = TRUE)
  return(windows)
}


## familial data
do.fam.load <- function(chromosomes=c(10)){
  
  DATA_DIR <- file.path("..", "..", "data", "mutation")
  
  ## load data for individuals from genomesunzipped.org
  # chromosomes : limit chromosomes loaded
  
  # read in each file, and strip down to chromosome 1 (again, could exapnd later if needed)
  fam <- lapply(lapply(dir(file.path(DATA_DIR, "fam"), full.names = TRUE), 
             function(x) read.delim(x, skip= 14, 
                                    blank.lines.skip=TRUE, stringsAsFactors=FALSE)
             ),
             function(x) subset(x, chromosome %in% chromosomes)
             )
  names(fam) <- dir(file.path(DATA_DIR, "fam"), full.names = FALSE)
  
  # return
  return(fam)
}

do.fam.prepare <- function(fam){
  ## scores for individual each SNP when compared between individuals
  # get SNP allele frequencies
  genotypes <- do.call("rbind",lapply(fam, FUN=function(x) table(x$genotype)))
  # create empty score matrix
  scores <- matrix( nrow=ncol(genotypes), ncol=ncol(genotypes) ,
                    dimnames=list(colnames(genotypes), colnames(genotypes)))
  # function for calculating number of alleles
  # allowing for ordering and homozygosity
  do.matching.alleles <- function(i,j){
    
    ## returns number of matching alleles for two allele pairs
    # expects i and j to be 2 character strings, i.e. i="AA", j="GT"
    
    # returns sum of non-zero row sums, 
    # when each character in i is compared to each character in j,
    # such that:  
    # in case of both pairs being heterozygous:
    # - return 2 if both match, 
    # - 1 if one allele matches, 
    # - 0 otherwise
    
    #   > sapply(c("A", "B"), function(x) x == c("A", "B"))
    #           A     B
    #   [A]  TRUE FALSE  -> TRUE
    #   [B] FALSE  TRUE  -> TRUE
    #                        sum = 2
    #   > sapply(c("B", "A"), function(x) x == c("A", "B"))
    #         B     A
    #   [A] FALSE  TRUE  -> TRUE
    #   [B]  TRUE FALSE  -> TRUE
    #                       sum = 2
    #   
    #   > sapply(c("A", "B"), function(x) x == c("B", "C"))
    #         A     B
    #   [B] FALSE  TRUE  -> TRUE
    #   [C] FALSE FALSE  -> FALSE
    #                       sum = 1
    #   
    # in case of both pairs being homozygous:
    # - return 2 if identical,
    # - 0 otherwise
    
    #   > sapply(c("B", "B"), function(x) x == c("B", "B"))
    #         B    B
    #   [1,] TRUE TRUE  -> TRUE
    #   [2,] TRUE TRUE  -> TRUE
    #                     sum = 2
    #   > sapply(c("B", "B"), function(x) x == c("C", "C"))
    #         B     B
    #   [C] FALSE FALSE  -> FALSE
    #   [C] FALSE FALSE  -> FALSE
    #                     sum = 0
    #   
    # in case of one pair being homozygous, the other heterozygous:
    # - return 1 if a heterozygous allele matches the homozygous pair
    # - 0 otherwise
    
    #   > sapply(c("A", "A"), function(x) x == c("A", "B"))
    #           A     A
    #   [A]  TRUE  TRUE  -> TRUE
    #   [B] FALSE FALSE  -> FALSE
    #                     sum = 1
    #   > sapply(c("B", "B"), function(x) x == c("A", "B"))
    #         B     B
    #   [A] FALSE FALSE  -> FALSE
    #   [B]  TRUE  TRUE  -> TRUE
    #                       sum = 1
    #   > sapply(c("A", "A"), function(x) x == c("C", "B"))
    #         A     A
    #   [C] FALSE FALSE  -> FALSE
    #   [B] FALSE FALSE  -> FALSE
    #                       sum = 0
    #   
    #   > sapply(c("A", "B"), function(x) x == c("B", "B"))
    #         A    B
    #   [B] FALSE TRUE     -> TRUE
    #   [B] FALSE TRUE     -> TRUE
    #                         rowsum = 2 
    #                       ***SHOULD EQUAL 1!!***
    #                           take colsum if second match is homozygous   
    #                           colsum = 1
    
    if(nchar(i) == 2 && nchar(j) == 2){
      
      return(
        
        sum(
          apply(
            sapply(strsplit(i, '')[[1]], # for each character in i
                   function(x) x == strsplit(j, '')[[1]]), # compare to each character in j
            ( length(unique(strsplit(j, '')[[1]])) == 1 ) +1, #colsum if homozygous
               function(x) sum(x) > 0) # sum of matches > 0 ?
          )
          
        )
    }
    else { return(0)}
      
  }
  
  # fill score matrix with precalculated match scores 
  # for all possible allele pair comparisons
  lapply(rownames(scores), function(i){
    lapply(colnames(scores), function(j){
      scores[i,j] <<-              # fill score matrix
        do.matching.alleles(i,j) * # number of matching alleles
        (                          # normalised to "liklihood" of seeing this pair of allele pairs
          (sum(genotypes) - sum(genotypes[,c(i,j)])) / sum(genotypes)
          )
    }
    )
  } )
  round(scores, 3)
  
  # return
  return(scores)
}

do.fam.check <- function(fam){
  ## check that all SNPs match
  checks <- c()
  proc.time() -> ptm
  
  # all sets have the same number of SNPs
  if(length(unique(lapply(fam, nrow))) == 1){
    checks <- append(checks, TRUE)
  } else { checks <- append(checks, FALSE)  }
  TIMES$fam.check.unique <<- proc.time() - ptm
  proc.time() -> ptm
  
  # all SNP rsids the same
  if(
    all(
      sapply(fam, 
             # check that the list of SNPs matches the first list of SNPs
             # i.e. intersection length is equal to unintersected
             function(x) length(intersect(fam[[1]]$X..rsid, x$X..rsid)) == length(fam[[1]]$X..rsid)
             )
        )
    ){
    checks <- append(checks, TRUE)
  } else { checks <- append(checks, FALSE)  }
  
  TIMES$fam.check.intersect <<- proc.time() - ptm
  proc.time() -> ptm
  
  # all SNPs in the same order
  if(
    all(
      sapply(fam, 
             # check that the list of SNPs matches the first list of SNPs
             # i.e. intersection length is equal to unintersected
             function(x) all(order(fam[[1]]$X..rsid) == order(x$X..rsid))
      )
    )
  ){
    checks <- append(checks, TRUE)
  } else { checks <- append(checks, FALSE)  }
  
  TIMES$fam.check.order <<- proc.time() - ptm
  
  # all checks good?
  return(all(checks))
}

do.ibd.vector <- function(fam, scores){
  do.ibd <- function(maf1, maf2, scores){
    
    ### IBD
    # http://en.wikipedia.org/wiki/Identity_by_descent
    # (browning and browning, 2007)[http://www.sciencedirect.com/science/article/pii/S0002929707638828]
    # (rapid calculation of IBD matrices)[http://www.biomedcentral.com/content/pdf/1297-9686-33-5-453.pdf]
        
    # combine to 'vector' of paired alleles
    genotype.vec <- cbind(maf1$genotype, maf2$genotype)
    
    # score each pair
    genotype.score <- apply(genotype.vec, 1, function(x) scores[x[1], x[2]])
    
    # return as annotated df for further analysis
    return(
      
      cbind(maf1[,names(maf1) != "genotype"], data.frame(score=genotype.score))
      
      )
    
  }

  fam.scores <- lapply(combn(1:length(fam),2, simplify = FALSE), # all pairwise combinations
                       # calculate ibd score vector
                       function(x) return(
                         list(df=do.ibd(fam[[x[1]]], fam[[x[2]]], scores=scores),
                              pair=names(fam)[x]
                              )
                         )
                       )
  return(fam.scores)
}

do.ibd.window <- function(fam.scores, window.sizes=seq(5, 50, 5)){

  ## calculate some sliding windows summing the scores
  
  # calculate sliding window over each comparison, 
  # calculating average score per SNP in that window
  
  do.call("rbind",
          unlist(recursive = FALSE, 
                 # for a range of window sizes
                 lapply(window.sizes, function(win){ 
                    # run sliding window, summing for each and dividing by number of markers
                    lapply(fam.scores, 
                           function(x) data.frame(width=win, pair=x$pair, q3=summary(wapply(x$df$score, width=win, by=as.integer(win/2), 
                                                  FUN=sum) / win)["3rd Qu."]))
  })))

}


### reporting
## collect data
if(DOWNLOAD){
  TIMES$download <- system.time(gcFirst = T,
                                do.download()
  )
}
## population data
# load data and compute matrix
TIMES$pop.load <- system.time(gcFirst = T,
                          pop <- do.pop.load()
)
# patient level summaries
TIMES$pop.fig1a <- system.time(gcFirst = T,
                          RESULTS$pop.fig1a <- do.pop.fig1a(pop = pop)
)
# gene level summaries
TIMES$pop.fig1b <- system.time(gcFirst = T,
                          RESULTS$pop.fig1b <- do.pop.fig1b(pop = pop)
)
# todo: something more advanced once packages are implemented in renjin
# maybe something from: http://cran.r-project.org/web/views/Genetics.html
# for example something from package "gap", like gcontrol() 
# http://www.inside-r.org/packages/cran/gap/docs/gcontrol


# clean up
TIMES$pop.clean <- system.time(gcFirst = T, expr = {
                                                    rm(pop);
                                                    gc()}
)

## "family" data - comparing between (non related) individuals
# load data and check SNP lists
TIMES$fam.load <- system.time(gcFirst = T,
                              fam <- do.fam.load()
)
TIMES$fam.check <- system.time(gcFirst = T,
                              stopifnot(do.fam.check(fam = fam))
)

# prepare pre-calculated scoring matrix
TIMES$fam.prepare <- system.time(gcFirst = T,
                              scores <- do.fam.prepare(fam = fam)
)

# create vectors comparing pairwise all individuals
TIMES$fam.score <- system.time(gcFirst = T,
                               fam.scores <- do.ibd.vector(fam = fam, scores = scores)
)
# score on sliding window
TIMES$fam.window <- system.time(gcFirst = T,
                               RESULTS$fam.window <- do.ibd.window(fam.scores = fam.scores)
)

# clean up
TIMES$fam.clean <- system.time(gcFirst = T, expr = {rm(fam, scores, fam.scores);
                                                    gc()}
)

## output results for comparison
# check output directories exist
if(!file.exists(file.path("..", "..", "generated", "results", basename(getwd())))){
  
  dir.create(file.path("..", "..", "generated", "results", basename(getwd())), recursive = TRUE)
  
}
if(!file.exists(file.path("..", "..", "generated", "timings", basename(getwd())))){
  
  dir.create(file.path("..", "..", "generated", "timings", basename(getwd())), recursive = TRUE)
  
}
# write results to file
write.table(file=file.path("..", "..", "generated", "results", 
                           basename(getwd()), paste(BENCHMARK, "results", "tsv", sep = '.')), 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE,
            do.call("rbind", 
                    lapply(names(RESULTS), function(x){
                      names(RESULTS[[x]]) <- c(1:ncol(RESULTS[[x]]))
                      cbind(RESULTS[[x]], data.frame(res=x))
                    }
                    )
            )
)

# timings
write.table(file=file.path("..", "..", "generated", "timings", 
                           basename(getwd()), paste(BENCHMARK, format(Sys.time(), "%Y%m%d%H%M%S"), "tsv", sep = '.')), 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE,
            format(do.call("rbind", TIMES), digits=5)
)

# final clean up
rm(list=ls())
gc()

