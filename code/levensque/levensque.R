# Code by Parham Solaimani
# Code adapted from ...
# Test-case workflow for TCGA Browser Shiny app which is developed in Mitch Levensque lab by Phil Cheng
# Analysis code provided by Phil Cheng.

### set up session
rm(list=ls())
# reproducibility
set.seed(1000)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

## packages used in TCGA prostate RNAseq
library(ncvreg) # source datasets from http://cran.r-project.org/web/packages/ncvreg/ncvreg.pdf
library(boot)
library(lars)
library(lasso2)
library(mda)
library(leaps)
library(data.table)
library(reshape2)
library(ggplot2)
library(magrittr)
library(shiny)
library(survival)
library(limma)
library(edgeR)
library(googleVis)
library(gage)
library(plyr)
library(reshape2)
library(STRINGdb)
library(ggplot2)
library(grid)
library(shinydashboard)
library(rCharts)
library(d3heatmap)
library(ggvis)
library(RColorBrewer)
library(DT)

## global vars
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "TCGAbrowser"
DATA_DIR <- file.path("..", "..", "data", "levensque")


# holder for results
RESULTS <- results(benchmark_name = BENCHMARK)
TIMES <- timings(benchmark_name = BENCHMARK)

#### functions
do.load <- function(DATA_DIR){
  
  # make sure DATA_DIR exists
  if(!file.exists(DATA_DIR)){   dir.create(DATA_DIR, recursive = TRUE)    }
  
  # Get url of all txt files in the folder
  files = list.files(DATA_DIR, pattern = "txt$")
  
  # read in RNAseq data
  d1 = fread(file = files[grep("RNA",files)], sep="\t", header=T)
  setkey(d1, Gene)
  gene.name <- d1$Gene
  
  # read in patient data
  pat = fread(file = files[grep("patient",files)], sep="\t", header=T)
  setkey(pat, bcr_patient_barcode, name)
  
  # read in exome sequencing data
  m1 = fread(file = files[grep("exome",files)], sep="\t", header=T)
  setkey(m1, bcr_patient_barcode)
  
  # place all data in a list to return
  DATA = list()
  DATA$d1 = d1
  DATA$pat = pat
  DATA$m1 = m1
  DATA$gene.name
  rm(d1,pat,m1,files,gene.name)
  
  return(DATA)
}
  
do.multiplot <- function(DATA){
  
  string_db <- STRINGdb$new(version="10", species=9606, score_threshold=400, input_directory= ".")
  kg.hsa <- kegg.gsets("hsa")
  
  #files = list.files(DATA_DIR, pattern = "txt$")
  
  #d1 <- fread(files[grep("RNA",files)], sep="\t", header=T)
  #setkey(d1, Gene)
  #pat <- fread(files[grep("patient",files)], sep="\t", header=T) #converts all spaces and NA and unknown into NA for R, easier for sorting
  #setkey(pat, bcr_patient_barcode, name)
  #gene.name <- d1$Gene
  #m1 <- fread(files[grep("exome",files)])
  #setkey(m1, bcr_patient_barcode)
  
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
  
}


do.survival <- function(DATA){
  
  #d1 <- fread(files[grep("RNA",files)], sep="\t", header=T)
  #setkey(d1, Gene)
  #gene.name <- d1$Gene
  d1 <- DATA$d1
  pat <- DATA$pat
  m1 <- DATA$m1
  gene.name <- DATA$gene.name
  
  g1 <- "BAZ2A"
  time <- "pfs"
  gleason <- c("2+4", "3+3", "3+4", "3+5", "4+3", "4+4", "4+5", "5+3", "5+4", "5+5")
  high <- (ncol(d1) - 2 ) * 0.75
  low <- (ncol(d1) - 2) * 0.25
  setkey(pat, gleason)
  pat.d1 <- d1[,c("Gene", pat[gleason, name]), with=F]
  setkey(pat.d1, Gene)
  
  #selecting row for gene and only retrieving values
  
  mgsub2 <- function(myrepl, mystring){
    gsub2 <- function(l, x){
      do.call('gsub', list(x = x, pattern = l[1], replacement = l[2]))
    }
    Reduce(gsub2, myrepl, init = mystring, right = T) 
  }
  
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
  
  #survival
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

do.deg_hm_degg <- function(){
  
  #differential expression
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
  
  #heatmap
  deg <- fit3$genes[1:100]
  d1.gene2 <- cpm(d1.gene[, colnames(d1.gene)[-length(colnames(d1.gene))], with=F], log=T, normalized.lib.sizes=T)
  rownames(d1.gene2) <- d1.gene$Gene
  map <- d1.gene2[deg,]
  heatmap.2(map, margins=c(5,5),col=redgreen,  trace = "none", keysize=0.6, labRow=F, labCol=F, scale = "row", ColSideColors=as.character(pat.gene[,gene2]+1))   
  
  #kegg factor
  d1[fit3$genes, Entrez]
  limma.fc <- fit3$logFC
  names(limma.fc) <- d1[fit3$genes, Entrez]
  kf <- "greater"
  fc.kegg.p <- gage(limma.fc, gsets = kg.hsa$kg.sets, ref = NULL, samp = NULL)
  sel <- fc.kegg.p$greater[, "p.val"] < 0.05 & !is.na(fc.kegg.p$greater[,"p.val"])
  greater <- data.frame(cbind(Pathway = rownames(fc.kegg.p$greater[sel,]),round(fc.kegg.p$greater[sel,1:5],5)))
  sel.1 <- fc.kegg.p$less[,"p.val"] < 0.05 & !is.na(fc.kegg.p$less[,"p.val"])
  less <-data.frame(cbind(Pathway = rownames(fc.kegg.p$less[sel.1,]), round(fc.kegg.p$less[sel.1,1:5],5)))
  
  #STRINGdb
  gene_mapped <- string_db$map(fit3, "genes", removeUnmappedRows = T)
  gene_mapped_pval05 <- string_db$add_diff_exp_color(subset(gene_mapped, adj.P.Val <0.01), logFcColStr="logFC")
  gene_mapped_pval05 <- gene_mapped_pval05[order(gene_mapped_pval05$adj.P.Val),]
  splot <- gene_mapped_pval05
  sorder <- "adj.P.Val"
  splot <- splot[with(splot, order(abs(get(sorder)), decreasing=T)),]
  hits <- splot$STRING_id[1:50]
  payload_id <- string_db$post_payload(splot$STRING_id, colors=splot$color)
  string_db$plot_network(hits, payload_id, add_link=F)
  
  #rCharts
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

do.exome <- function(){
  
  
  
  #exome data
  glist <- c("FRG1B", "SPOP", "TP53", "ANKRD36C", "KMT2C", "KMT2D", "KRTAP4-11", "SYNE1", "NBPF10", "ATM", "FOXA1", "LRP1B", "OBSCN", "SPTA1", "USH2A", "AHNAK2", "FAT3", "CHEK2", g1)
  glist <- glist[!duplicated(glist)]
  setkey(m1, Hugo_Symbol)
  test <- m1[glist, .(bcr_patient_barcode, Hugo_Symbol, Variant_Classification, Amino, Amino2)]
  test <- unique(test, by=c("bcr_patient_barcode", "Hugo_Symbol", "Amino"))
  pat.genem <- pat.gene[pat.gene$name %in% unique(m1$bcr_patient_barcode), .(name, gene2)]
  setkey(pat.genem, gene2)
  setkey(test, bcr_patient_barcode)
  test1 <- test[.(pat.genem[gene2 == 2, name])]
  #for gene high, exome and copy number
  setkey(test1, Hugo_Symbol, Amino)
  test2 <- test1[!is.na(test1$Hugo_Symbol)] #dcast table with NA them remove NA column
  #test1 <- rbindlist(list(test1[!(c("BRAF", "NRAS"))], test1[.("BRAF", "V600")], test1[.("NRAS", c("G12", "G13", "Q61"))])) #removes all mutations from BRAF and NRAS that are not typical hotspot
  #setkey(test1, Gene, bcr_patient_code)
  mut2 <- dcast.data.table(test2[!is.na(bcr_patient_barcode)], bcr_patient_barcode ~ Hugo_Symbol)
  #e1 <- pat.genem[gene2 == 2, name][!(pat.genem2[gene == 2, name]) %in% mut2$bcr_patient_code]    #checks for all wildtype
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
  gg1 <- ggplot(mut3,aes(bcr_patient_barcode,variable,fill=as.factor(value)))+geom_tile(colour=c("white"))+labs(x = "Patient", y="Gene")+scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+scale_fill_brewer(type="qual",name="Legend", palette=6, labels=c("Wild-type", "1 SNV", "2 SNVs", "3 SNVs") )+theme(title=element_text(size=28), axis.ticks=element_blank(),  axis.text.x=element_blank(), axis.text.y=element_text(size=24), axis.text.x=element_text(size=28), axis.text.y=element_text(size=28), legend.text=element_text(size=18))+ggtitle(paste(g1, "high"))
  gg1
  
  #using ggvis
  callup <- test1[, list(Amino=list(Amino)), by=c("bcr_patient_barcode", "Hugo_Symbol")] #aggregates multiple mutations per gene per patient into 1 cell
  setkey(callup, bcr_patient_barcode, Hugo_Symbol)
  callup[, Amino3 := sapply(callup$Amino, function(x) paste(unlist(x), collapse = " "))]
  callup2 <- callup[.(mut3$bcr_patient_barcode, mut3$variable)]
  f_dowle2 = function(DT) {
    # or by number (slightly faster than by name) :
    for (j in seq_len(ncol(DT)))
      set(DT,which(is.na(DT[[j]])),j, "Wild-type")
  }
  f_dowle2(callup2) #change all NULL to Wild-type
  
  mut3[, c("Gene", "Amino") := list(callup2$Hugo_Symbol, callup2$Amino3)]
  mut3[, callup := paste(mut3$bcr_patient_barcode, mut3$Hugo_Symbol, mut3$Amino, sep="\n")]
  cellinfo <- function(x){
    if(is.null(x)) return(NULL)
    paste(x$callup)
  }
  
  
  mut3  %>% 
    ggvis(~bcr_patient_barcode, ~variable, fill=~value, key:= ~callup)  %>% 
    layer_rects(width=band(), height=band()) %>% 
    add_tooltip(cellinfo, "hover")
  
  #exome graph for gene low group
  
  test1 <- test[.(pat.genem[gene2 == 1, name])]
  #for gene high, exome and copy number
  setkey(test1, Hugo_Symbol, Amino)
  test2 <- test1[!is.na(test1$Hugo_Symbol)] #dcast table with NA them remove NA column
  #test1 <- rbindlist(list(test1[!(c("BRAF", "NRAS"))], test1[.("BRAF", "V600")], test1[.("NRAS", c("G12", "G13", "Q61"))])) #removes all mutations from BRAF and NRAS that are not typical hotspot
  #setkey(test1, Gene, bcr_patient_code)
  mut2 <- dcast.data.table(test2[!is.na(bcr_patient_barcode)], bcr_patient_barcode ~ Hugo_Symbol)
  #e1 <- pat.genem[gene2 == 2, name][!(pat.genem2[gene == 2, name]) %in% mut2$bcr_patient_code]    #checks for all wildtype
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
  gg2 <- ggplot(mut3,aes(bcr_patient_barcode,variable,fill=as.factor(value)))+geom_tile(colour=c("white"))+labs(x = "Patient", y="Gene")+scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+scale_fill_brewer(type="qual",name="Legend", palette=6, labels=c("Wild-type", "1 SNV", "2 SNVs", "3 SNVs") )+theme(title=element_text(size=28), axis.ticks=element_blank(),  axis.text.x=element_blank(), axis.text.y=element_text(size=24), axis.text.x=element_text(size=28), axis.text.y=element_text(size=28), legend.text=element_text(size=18))+ggtitle(paste(g1, "low"))
  gg2
  
  multiplot(gg1,gg2)
  
}
