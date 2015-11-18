# Code by Parham Solaimani
# Test-case workflow for TCGA Browser Shiny app which is developed in Mitch Levesque lab by Phil Cheng
# Analysis code provided by Phil Cheng.

##### set up session #####
rm(list=ls())
set.seed(1000)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative
.libPaths(file.path("~","R","libs"))
source(file.path("..", "..","benchmark_utilities.R"))

##### loading packages #####
library(ncvreg)
library(boot)
library(lars)
library(lasso2)
library(mda)
library(leaps)
library(data.table)
library(reshape2)
library(ggplot2)
library(magrittr)
library(survival)
library(limma)
library(edgeR)
library(googleVis)
library(gage)
library(plyr)
library(reshape2)
library(STRINGdb)
library(grid)
#library(rCharts)
library(d3heatmap)
library(ggvis)
library(RColorBrewer)
library(DT)
library(jsonlite)

cat("\nLoaded all libraries.....\n") #DEBUG

##### Set global vars #####
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "TCGAbrowser"
DATA_DIR <- file.path("..","..", "data","levesque")
files = list.files(path = DATA_DIR, pattern = "txt$")
RESULTS <- results(benchmark_name = BENCHMARK)
TIMES <- timings(benchmark_name = BENCHMARK)

##### Functions #####
mgsub2 <- function(myrepl, mystring){
  gsub2 <- function(l, x){
    do.call('gsub', list(x = x, pattern = l[1], replacement = l[2]))
  }
  Reduce(gsub2, myrepl, init = mystring, right = T)
}
cellinfo <- function(x){
  if(is.null(x)) return(NULL)
  paste(x$callup)
}
f_dowle2 = function(DT) {
  # or by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j, "Wild-type")
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
##### Blocks for timing #####
do.load <- function(DATA_DIR,TARGET){
  # Load the required Data
  DATA <- list()
  if ( TARGET == "RNAseq" ) {
    d1 <- fread(paste0(DATA_DIR,"/TCGA_PRAD_RNAseq.txt"), sep="\t", header=T)
    setkey(d1, Gene)
    pat <- fread(paste0(DATA_DIR,"/TCGA_PRAD_patient.txt"), sep="\t", header=T) #converts all spaces and NA and unknown into NA for R, easier for sorting
    setkey(pat, bcr_patient_barcode, name)
    m1 <- fread(paste0(DATA_DIR,"/TCGA_PRAD_exome.txt"))
    setkey(m1, bcr_patient_barcode)
    gene.name <- d1$Gene
    string_db <- STRINGdb$new(version="10", species=9606, score_threshold=400, input_directory= ".")
    kg.hsa <- kegg.gsets("hsa")
    DATA <- list( d1, pat, m1, gene.name, string_db, kg.hsa )
    return(DATA)
  }

  if ( TARGET == "Survival") {
    d1 <- fread(paste0(DATA_DIR,"/TCGA_PRAD_RNAseq.txt"), sep="\t", header=T)
    setkey(d1, Gene)
    pat <- fread(paste0(DATA_DIR,"/TCGA_PRAD_patient.txt"), sep="\t", header=T) #converts all spaces and NA and unknown into NA for R, easier for sorting
    setkey(pat, bcr_patient_barcode, name)
    DATA <- list(d1,pat)
  }

  if ( TARGET == "Exome") {
    d1 <- fread(paste0(DATA_DIR,"/TCGA_PRAD_RNAseq.txt"), sep="\t", header=T)
    setkey(d1, Gene)
    pat <- fread(paste0(DATA_DIR,"/TCGA_PRAD_patient.txt"), sep="\t", header=T)
    setkey(pat, bcr_patient_barcode, name)
    m1 <- fread(paste0(DATA_DIR,"/TCGA_PRAD_exome.txt"))
    setkey(m1, Hugo_Symbol)
    DATA <- list(d1,pat,m1)
  }
  return(DATA)
}

do.preprocess <- function(DATA,TARGET){
  # Performs rearragement and formatting of Data

  if ( TARGET == "RNAseq" ) {
    d1 	<- DATA[[1]]
    pat 	<- DATA[[2]]
    m1 	<- DATA[[3]]
    gene.name <- DATA[[4]]
    string_db <- DATA[[5]]
    kg.hsa <- DATA[[6]]
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
    d1.gene <- d1[, pat.gene[, name], with=F]
    d1.gene[, Gene := gene.name]
    setkey(d1.gene, Gene)
    DATA <- list(d1, d1.gene, gene.name, pat.gene, kg.hsa, string_db)
  }

  if ( TARGET == 'Survival') {
    d1 	<- DATA[[1]]
    pat 	<- DATA[[2]]
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
    DATA <- list(pat.gene,g1)
  }

  if ( TARGET == 'Exome') {
    d1 <- DATA[[1]]
    pat <- DATA[[2]]
    m1 <- DATA[[3]]
    gene.name <- d1$Gene
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
    glist <- c("FRG1B", "SPOP", "TP53", "ANKRD36C", "KMT2C", "KMT2D", "KRTAP4-11", "SYNE1", "NBPF10", "ATM", "FOXA1", "LRP1B", "OBSCN", "SPTA1", "USH2A", "AHNAK2", "FAT3", "CHEK2", g1)
    glist <- glist[!duplicated(glist)]
    test <- m1[glist, .(bcr_patient_barcode, Hugo_Symbol, Variant_Classification, Amino, Amino2)]
    test <- unique(test, by=c("bcr_patient_barcode", "Hugo_Symbol", "Amino"))
    pat.genem <- pat.gene[pat.gene$name %in% unique(m1$bcr_patient_barcode), .(name, gene2)]
    setkey(pat.genem, gene2)
    setkey(test, bcr_patient_barcode)
    test1 <- test[.(pat.genem[gene2 == 2, name])]

    #for gene high, exome and copy number
    setkey(test1, Hugo_Symbol, Amino)
    test2 <- test1[!is.na(test1$Hugo_Symbol)] #dcast table with NA them remove NA column
    mut2 <- dcast.data.table(test2[!is.na(bcr_patient_barcode)], bcr_patient_barcode ~ Hugo_Symbol)
    mut2 <- rbindlist(list(mut2, data.table(bcr_patient_barcode = test1[is.na(test1$Hugo_Symbol), bcr_patient_barcode])), fill=T)
    for (j in seq_len(ncol(mut2))[-1])  {
      if (any(is.na(mut2[[j]]))) {
        set(mut2, which(mut2[[j]] != 0),j,1)
        set(mut2, which(is.na(mut2[[j]])),j,0)
      } else {
        set(mut2, which(mut2[[j]] != 0),j,1)
      }
    }

    mut2[, glist[!(glist %in% colnames(mut2)[-1])] := 0]
    mut2 <- mut2[, lapply(.SD, as.numeric), by= bcr_patient_barcode]
    setkeyv(mut2, names(sort(apply(mut2[,colnames(mut2)[-1], with=F], 2, sum), decreasing=T)))

    #graphing exome table
    setkeyv(mut2, names(sort(apply(mut2[,colnames(mut2)[-1], with=F], 2, sum), decreasing=T)))

    mut3 <- melt(mut2)
    mut3$bcr_patient_barcode <- factor(mut3$bcr_patient_barcode, levels=rev(mut2$bcr_patient_barcode))
    mut3$variable <- factor(mut3$variable, levels=names(sort(apply(mut2[,colnames(mut2)[-1], with=F], 2, sum), decreasing=F)))
    mut3_gg1 <- mut3

    #using ggvis
    callup <- test1[, list(Amino=list(Amino)), by=c("bcr_patient_barcode", "Hugo_Symbol")] #aggregates multiple mutations per gene per patient into 1 cell
    setkey(callup, bcr_patient_barcode, Hugo_Symbol)
    callup[, Amino3 := sapply(callup$Amino, function(x) paste(unlist(x), collapse = " "))]
    callup2 <- callup[.(mut3$bcr_patient_barcode, mut3$variable)]
    f_dowle2(callup2) #change all NULL to Wild-type

    mut3[, c("Gene", "Amino") := list(callup2$Hugo_Symbol, callup2$Amino3)]
    mut3[, callup := paste(mut3$bcr_patient_barcode, mut3$Hugo_Symbol, mut3$Amino, sep="\n")]

    mut3  %>%
      ggvis(~bcr_patient_barcode, ~variable, fill=~value, key:= ~callup)  %>%
      layer_rects(width=band(), height=band()) %>%
      add_tooltip(cellinfo, "hover")

    #exome graph for gene low group

    test1 <- test[.(pat.genem[gene2 == 1, name])]
    #for gene high, exome and copy number
    setkey(test1, Hugo_Symbol, Amino)
    test2 <- test1[!is.na(test1$Hugo_Symbol)] #dcast table with NA them remove NA column
    mut2 <- dcast.data.table(test2[!is.na(bcr_patient_barcode)], bcr_patient_barcode ~ Hugo_Symbol)
    mut2 <- rbindlist(list(mut2, data.table(bcr_patient_barcode = test1[is.na(test1$Hugo_Symbol), bcr_patient_barcode])), fill=T)
    for (j in seq_len(ncol(mut2))[-1])  {
      if (any(is.na(mut2[[j]]))) {
        set(mut2, which(mut2[[j]] != 0),j,1)
        set(mut2, which(is.na(mut2[[j]])),j,0)
      } else {
        set(mut2, which(mut2[[j]] != 0),j,1)
      }
    }

    #graphing exome table
    mut2[, glist[!(glist %in% colnames(mut2)[-1])] := 0]
    mut2 <- mut2[, lapply(.SD, as.numeric), by= bcr_patient_barcode]
    setkeyv(mut2, names(sort(apply(mut2[,colnames(mut2)[-1], with=F], 2, sum), decreasing=T)))

    mut3 <- melt(mut2)
    mut3$bcr_patient_barcode <- factor(mut3$bcr_patient_barcode, levels=rev(mut2$bcr_patient_barcode))
    mut3$variable <- factor(mut3$variable, levels=names(sort(apply(mut2[,colnames(mut2)[-1], with=F], 2, sum), decreasing=F)))
    DATA = list(mut3_gg1, g1, mut3)
  }
  return(DATA)
}

do.analyse <- function(DATA,TARGET){
  # Performs calculations and plotting
  if ( TARGET == "RNAseq_DEG" ) {
    d1 <- DATA[[1]]
    d1.gene <- DATA[[2]]
    gene.name <- DATA[[3]]
    pat.gene <- DATA[[4]]
    kg.hsa <- DATA[[5]]
    string_db <- DATA[[6]]
    d3 <- DGEList(counts=d1.gene[,setdiff(colnames(d1.gene), "Gene"), with=F], genes=gene.name, group=pat.gene$gene2)
    isexpr <- rowSums(cpm(d3)>1) >= (ncol(d3)/2) #only keeps genes with at least 1 count-per-million in at least half the samples
    d3 <- d3[isexpr,]
    #d4<- calcNormFactors(d3)
    design <- model.matrix(~pat.gene[,gene2])
    v2 <- voom(d3, design, plot=T)
    #plotMDS(v2, top=100, labels=substring(phenosgene$gene3, 1, 5), col=ifelse(phenosgene$gene3 == "high", "blue", "red"), gene.selection="common")
    fit <- lmFit(v2, design)
    fit2 <- eBayes(fit)
    fit3 <- topTable(fit2, coef=2, n=Inf, lfc=1, p.value=0.05)
    head(fit3)
    DATA <- list(d1, d1.gene, gene.name, pat.gene, kg.hsa, string_db, fit3)
  }

  if ( TARGET == "RNAseq_HM" ) {
    d1.gene <- DATA[[2]]
    pat.gene <- DATA[[4]]
    fit3 <- DATA[[7]]

    deg <- fit3$genes[1:100]
    d1.gene2 <- cpm(d1.gene[, colnames(d1.gene)[-length(colnames(d1.gene))], with=F], log=T, normalized.lib.sizes=T)
    rownames(d1.gene2) <- d1.gene$Gene
    map <- d1.gene2[deg,]
    cat("\n before pdf \n")
    pdf(file="RNAseq_DEG_heatmap.pdf",height = 50, width = 50)
    heatmap.2(map, col=redgreen,  trace = "none", keysize=0.6, labRow=F, labCol=F, scale = "row", ColSideColors=as.character(pat.gene[,gene2]+1))
    dev.off()
    cat("\n after pdf \n")
    DATA <- head(fit3)
  }

  if ( TARGET == "RNAseq_KEGG" ) {
    d1 <- DATA[[1]]
    kg.hsa <- DATA[[5]]
    fit3 <- DATA[[7]]

    d1[fit3$genes, Entrez]
    limma.fc <- fit3$logFC
    names(limma.fc) <- d1[fit3$genes, Entrez]
    kf <- "greater"
    fc.kegg.p <- gage(limma.fc, gsets = kg.hsa$kg.sets, ref = NULL, samp = NULL)
    sel <- fc.kegg.p$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.p$greater[,"p.val"])
    greater <- data.frame(cbind(Pathway = rownames(fc.kegg.p$greater[sel,]),round(fc.kegg.p$greater[sel,1:5],5)))
    sel.1 <- fc.kegg.p$less[,"p.val"] < 0.05 & !is.na(fc.kegg.p$less[,"p.val"])
    less <-data.frame(cbind(Pathway = rownames(fc.kegg.p$less[sel.1,]), round(fc.kegg.p$less[sel.1,1:5],5)))
  }

  if ( TARGET == "RNAseq_STRING" ) {
    string_db <- DATA[[6]]
    fit3 <- DATA[[7]]

    gene_mapped <- string_db$map(fit3, "genes", removeUnmappedRows = T)
    gene_mapped_pval05 <- string_db$add_diff_exp_color(subset(gene_mapped, adj.P.Val <0.01), logFcColStr="logFC")
    gene_mapped_pval05 <- gene_mapped_pval05[order(gene_mapped_pval05$adj.P.Val),]
    splot <- gene_mapped_pval05
    sorder <- "adj.P.Val"
    splot <- splot[with(splot, order(abs(get(sorder)), decreasing=T)),]
    hits <- splot$STRING_id[1:50]
    payload_id <- string_db$post_payload(splot$STRING_id, colors=splot$color)
    string_db$plot_network(hits, payload_id, add_link=F)
  }

  if ( TARGET == "RNAseq_rChart" ) {
    pat.gene <- DATA[[4]]

    group <- "clinical_T"
    gr <- pat.gene[, .N , by=.(gene2, with(pat.gene, get(group)))][order(with)]
    setkey(gr, gene2, with)
    setnames(gr, 2, group)
    gr[gene2 == 1, gene3 := "low"]
    gr[gene2 == 2, gene3 := "high"]
    n1 <- nPlot(N ~ gene3, group = group, data = gr, type = "multiBarChart")
    n1$addParams(dom = "bargraph")
    n1$chart(margin = list(left = 100))
    n1$yAxis(axisLabel  ="Number of patients")
    n1$print("test")
    n1$save("test1.html", cdn=T)
  }

  if ( TARGET == 'Survival') {
    pat.gene <- DATA[[1]]
    g1 <- DATA[[2]]

    m.surv <- Surv(pat.gene$pfs_days, pat.gene$pfs)
    sdf <- survdiff(m.surv ~ pat.gene$gene2)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

    survplot <- survfit(Surv(pfs_days, pfs) ~ gene2, data = pat.gene)
    half <- summary(survplot)$table[,"median"]

    #tiff("Survival_all_met.tif", height=450, width=800)
    plot(survplot, col=1:2, xlab="Survival Time (Days)", ylab="Survival", lwd=3)
    legend("topright", c(sprintf("%s low, n=%s", g1, sdf$n[1]), sprintf("Median survival %s days", half[1]) ,
                         sprintf("%s high, n=%s", g1, sdf$n[2]), sprintf("Median survival %s days", half[2])), col=c(1,1,2,2), lty=c(1,0,1,0))
    legend("bottomleft", paste("p = ", round(p.val,4)))
    lines(c(0,half[1]), c(0.5, 0.5), lwd=1, lty =2, col=1)
    lines(c(half[1],half[1]), c(0, 0.5), lwd=1, lty =2, col=1)
    lines(c(0,half[2]), c(0.5, 0.5), lwd=1, lty=2, col=2)
    lines(c(half[2],half[2]), c(0, 0.5), lwd=1, lty=2, col=2)

  }

  if ( TARGET == 'Exome') {
    mut3_gg1 <- DATA[[1]]
    g1 <- DATA[[2]]
    mut3 <- DATA[[3]]

    gg1 <- ggplot(mut3_gg1,aes(bcr_patient_barcode,variable,fill=as.factor(value))) +
                geom_tile(colour=c("white")) +
                labs(x = "Patient", y="Gene") +
                scale_x_discrete(expand=c(0,0)) +
                scale_y_discrete(expand=c(0,0)) +
                scale_fill_brewer(type="qual",name="Legend", palette=6, labels=c("Wild-type", "1 SNV", "2 SNVs", "3 SNVs") ) +
                theme(  title=element_text(size=28),
                          axis.ticks=element_blank(),
                          axis.text.x=element_blank(),
                          axis.text.y=element_text(size=24),
                          axis.text.x=element_text(size=28),
                          axis.text.y=element_text(size=28),
                          legend.text=element_text(size=18)
                      ) +
                ggtitle(paste(g1, "high"))
    gg1
    gg2 <- ggplot(mut3,aes(bcr_patient_barcode,variable,fill=as.factor(value))) +
                geom_tile(colour=c("white")) +
                labs(x = "Patient", y="Gene") +
                scale_x_discrete(expand=c(0,0)) +
                scale_y_discrete(expand=c(0,0)) +
                scale_fill_brewer(type="qual",name="Legend", palette=6, labels=c("Wild-type", "1 SNV", "2 SNVs", "3 SNVs") ) +
                theme(title=element_text(size=28),
                        axis.ticks=element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_text(size=24),
                        axis.text.x=element_text(size=28),
                        axis.text.y=element_text(size=28),
                        legend.text=element_text(size=18)
                      ) +
                ggtitle(paste(g1, "low"))
    gg2
  }
  return(DATA)
}

############################################################################
################### TIMING AND REPORTING ###################################
############################################################################

cat("\nTIMES   record=do.load(DATA_DIR, 'RNAseq').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "load_RNAseq",
                   record = system.time(gcFirst = T,
                                        DATA_rna <- do.load(DATA_DIR, "RNAseq")
                                        )
                   )

cat("\nTIMES   record=do.preprocess(DATA_rna, 'RNAseq').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "proc_RNA",
                   record = system.time(gcFirst = T,
                                        DATA_rna <- do.preprocess(DATA_rna,"RNAseq")
                                        )
                   )

cat("\nTIMES   record=do.analyse(DATA_rna,'RNAseq_DEG').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "compute_RNA_DEG",
                   record = system.time(gcFirst = T,
                                        DATA_rna <- do.analyse(DATA_rna,"RNAseq_DEG")
                                        )
                   )

cat("\nTIMES   record=do.analyse(DATA_rna,'RNAseq_HM').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "compute_RNA_HM",
                   record = system.time(gcFirst = T,
                                        RESULTS <- addRecord(RESULTS, record_name="compute_RNA_HM",
                                                             record = do.analyse(DATA_rna,"RNAseq_HM")
                                        )
                   ))

cat("\nTIMES   record=do.analyse(DATA_rna,'RNAseq_KEGG').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "ml_analyze_RNAKEGG",
                   record = system.time(gcFirst = T,
                                        do.analyse(DATA_rna,"RNAseq_KEGG")
                                        )
                   )

cat("\nTIMES   record=do.analyse(DATA_rna,'RNAseq_STRING').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "compute_RNA_STRING",
                   record = system.time(gcFirst = T,
                                        do.analyse(DATA_rna,"RNAseq_STRING")
                                        )
                   )

# cat("\nTIMES   record=do.analyse(DATA_rna,'RNAseq_rChart').....\n") #DEBUG
# TIMES <- addRecord(TIMES, record_name = "compute_RNAr_Chart",
#                    record = system.time(gcFirst = T,
#                                         do.analyse(DATA_rna,"RNAseq_rChart")
#                                         )
#                    )

#######
cat("\nTIMES   record=do.load(DATA_DIR, 'Exome').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "load_Exome",
                   record = system.time(gcFirst = T,
                                        DATA_exo <- do.load(DATA_DIR, "Exome")
                                        )
                   )

cat("\nTIMES   record=do.preprocess(DATA_exo, 'Exome').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "proc_Exome",
                   record = system.time(gcFirst = T,
                                        DATA_exo <- do.preprocess(DATA_exo,"Exome")
                                        )
                   )

cat("\nTIMES   record=do.analyse(DATA_exo,'Exome').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "compute_Exome",
                   record = system.time(gcFirst = T,
                                        do.analyse(DATA_exo,"Exome")
                                        )
                   )

#######
cat("\nTIMES   record=do.load(DATA_DIR, 'Survival').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "load_Survival",
                   record = system.time(gcFirst = T,
                                        DATA_sur <- do.load(DATA_DIR, "Survival")
                                        )
                   )

cat("\nTIMES   record=do.preprocess(DATA_sur, 'Survival').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "proc_Survival",
                   record = system.time(gcFirst = T,
                                        DATA_sur <- do.preprocess(DATA_sur,"Survival")
                                        )
                   )

cat("\nTIMES   record=do.analyse(DATA_sur,'Survival').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "compute_Survival",
                   record = system.time(gcFirst = T,
                                        do.analyse(DATA_sur,"Survival")
                                        )
                   )

### reporting
# write results to file
cat("\nreportRecords(RESULTS).....\n") #DEBUG
reportRecords(RESULTS)

# timings
cat("\nreportRecords(TIMES).....\n") #DEBUG
reportRecords(TIMES)

# final clean up
cat("\nrm(list=ls()).....\n") #DEBUG
rm(list=ls())
gc()
