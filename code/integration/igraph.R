# co-citation analysis with iGraphs
# ieuan.clay@gmail.com
# May 2015

### a minimal script for extracting/cleaning/and visualising entrez RIF networks
# load iGraph network object from falt table
# do some basic analysis/coloring/plotting
# do some subsetting
#  - random 1%
#  - limit genes by GO terms (or anything other annotation)
#  - limit papers by MeSh terms

### set up session
rm(list=ls())
# reproducibility
set.seed(8008)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative

# load utilities
source(file.path("..", "..","benchmark_utilities.R"))

## packages
library(igraph) # requires > R-15.2
library(XML)
library(R.utils)
library(plyr)
library(reshape)
library(utils)
library(sqldf)

## global vars
VERBOSE <- TRUE # print progress?
DATA_DIR <- file.path("..", "..", "data", "integration")
DOWNLOAD <- FALSE
INPUT <- "ftp://ftp.ncbi.nih.gov/gene/GeneRIF/generifs_basic.gz"

# holder for results
RESULTS <- list()
TIMES <- list()
BENCHMARK <- "igraph"

#### functions

do.download <- function(INPUT, DATA_DIR, DOWNLOAD=TRUE){
  
  ### download files from [INPUT] to [DATA_DIR]
  ## RIFs
  # get RIF file from entrez ftp server
  data_path <- file.path(DATA_DIR, "generifs_basic.gz")
  if(DOWNLOAD){
    download.file(url=INPUT, destfile=data_path)
  }

  return(data_path)
}

do.plot <- function(network, title="", layout=igraph::layout.kamada.kawai,
                    # node properties
                    shape="circle", size=1, colour="black", label=NA,
                    # edge properties
                    weight=1                
                    ){
  
  ### for simplicity, basic plotting of iGraph objects
  # provides default arguments for most parameters
  # otherwise can be passed from:
  # V(..igraph instance..)$..property name.. 
  # and
  # E(..igraph instance..)$..property name.. 
  # functions
  
  plot(
    network,                        # the graph to be plotted
    layout=layout,	                # kamada.kawai|fruchterman.reingold|lgl|auto -  the layout method. see the igraph documentation for details
    main=title,	                    # specifies the title
    vertex.shape=shape,             # nodeshape
    vertex.size=size,               # nodesize  
    vertex.color=colour,            # node colour
    edge.width=weight,              # edge weight
    vertex.label.color='black',		  # the color of the name labels
    vertex.label.font=2,			      # the font of the name labels
    vertex.label=label,		          # specifies the lables of the vertices. in this case the 'name' attribute is used
    vertex.label.cex=1			        # specifies the size of the font of the labels. can also be made to vary
    
  )
  
  
}

do.load.edges <- function(PATH){
  
  # unpack data
  tmpfile <- file.path(basename(PATH), "rif.tmp")
  gunzip(filename=PATH, remove=FALSE,
         destname=tmpfile)
  
  ## load and 'query' downloaded data
  # There is more information on sqldf on the sqldf home page:
  # http://sqldf.googlecode.com
  generifs_basic <- read.delim(tmpfile, header=T)
  names(generifs_basic) <- c("tax_id", "gene_id", "pubmed_ids", "timestamp", "annotation")
  # load statement and execute
  statement <- paste(readLines("get_all_RIF.sql"), collapse="\n")
  #sqldf(drv = "SQLite","select * from generifs_basic limit 5") # 'head'
  edges <- sqldf(drv = "SQLite", statement) # 43491 as of 28/06/13
  
  ## remove tmp file
  file.remove(tmpfile)
  
  return(edges)
  
}

do.load <- function(PATH, plot_loaded_graph=TRUE){
  
  edges <- do.load.edges(PATH)
  
  ### build data into a graph object
  ## make a bipartite graph of pubmed ids and genes
  # graph.data.frame() takes a data frame where the first two columns are nodes IDs (so each row is an edge) 
  # and all the other columns are taken as edge attributes, node attributes are loaded separately - see below
  network <- graph.data.frame(
    d=subset(
      x=edges, 
      subset=pubmed_ids %in% (sample(unique(edges$pubmed_ids), size=(abs(length(unique(edges$pubmed_ids))/100))*5)), 
      select=c("pubmed_ids", "gene_id", "annotation")), # ~5% papers
    directed=F)
  # add node attributes
  # V(network) # all node objects
  V(network)$shape <- c('circle', 'square')[1 + V(network)$name %in% edges$pubmed_ids]
  V(network)$color <- c('black', 'red')[1 + V(network)$name %in% edges$pubmed_ids] # black dots - genes, red boxes - papers
  V(network)$isa <- c('gene', 'paper')[1 + V(network)$name %in% edges$pubmed_ids]
  V(network)$degree <- scale(degree(network), center=F, scale=T)
  # degree(network)
  if(plot_loaded_graph){
    do.plot(network, title="5%er", size = V(network)$degree, shape=V(network)$shape, colour = V(network)$color)
  }

  return(network)
  
}

do.decompose <- function(network, plot_results=TRUE){
  ## pick out the biggest component
  network <- decompose.graph(network) # coverts to list of networks (each component as separate element)
  largest <- function(x){
    # return largest sub-component
    sizes <- sapply(
      X=x, 
      FUN=function(x){return(length(V(x)))}
    )
    max.size <- sapply(
      X=sizes, 
      FUN=function(i){return(i==max(sizes))}
    ) # vector of T/F for (potentially more than one) biggest component(s)
    # just take largest sub-list (biggest component)
    if(sum(max.size==1)){
      return(x[max.size][[1]])
    }
    else {
      return(x[max.size])
    }
  }
  network <- largest(network)
  ## plot it out
  V(network)$degree <- scale(degree(network), center=F, scale=T)
  if(plot_results){
    do.plot(network, title="Largest component", 
            size = V(network)$degree, shape=V(network)$shape, colour = V(network)$color
            )
  }
  
  return(network)
  
}
do.cocitation <- function(network, plot_results = TRUE){
  
  ### use cocitation to project graph into single type space
  # colour and shape according to node type
  genes <- get.data.frame(network, what = "vertices")
  genes <- genes[genes$isa == "gene", "name"]
  
  # do cocitation
  network <- graph.adjacency(
    cocitation(
      graph=network, 
      v=V(network)),
    mode="undirected", weighted=TRUE)
  # add annotations
  V(network)$degree <- scale(degree(network), center=F, scale=T)
  # black dots - genes, red boxes - papers
  V(network)$shape <- c('circle', 'square')[1 + !(V(network)$name %in% genes)]
  V(network)$color <- c('black', 'red')[1 + !(V(network)$name %in% genes)] 
  V(network)$isa <- c('gene', 'paper')[1 + !(V(network)$name %in% genes)] 
  
  # dump "paper" nodes to just examine "gene" nodes
  network <- do.subset(network, node_ids=genes)
  
  # plot it
  if(plot_results){
    do.plot(network, title="Genes only", weight = E(network)$weight,
            size = V(network)$degree, shape=V(network)$shape, colour = V(network)$color,
    )
  }
  
  return(network)
  }
  
  
}

do.subset <- function(network, node_ids){
  
  ## return igraph instance containing only the nodes listed
  # copying all properties
  network <- induced.subgraph(network, vids = genes, impl="copy_and_delete")
  
  return(network)
  
}

do.phospho <- function(DATA_DIR, PATH){
  ### construct network, not with random 1%, but with phosphatase/kinase subset
  ## pull annotations for human PPases and kinases from GO, see later for use
  pk.ids <- read.delim(header=T, file=file.path(DATA_DIR, "pk_pp_9606.txt"))
  # how many of each gene type did we get?
  if(VERBOSE){cat(table(pk.ids[,c("go_term")]))}
  # phosphoprotein phosphatase activity             protein kinase activity 
  # 33                                 195 
  
  ## collect a fresh edge list containing all the information
  network <- do.load.edges(PATH)
  # subset to phospho interaction network
  edges <- subset(edges, gene_id %in% pk.ids$gene_id)
  # create graph
  network<-graph.data.frame(
    d=edges,
    directed=F)
  
  ## decompose and run cocitation
  network <- do.decompose(network)
  V(network)$isa <- c("gene", "paper")[1 + V(network) %in% edges$pubmed_ids]
  network <- do.cocitation(network)
  
  ## plot
  # black dots = kinases, red boxes = PPases
  V(pk.gene)$shape <- c('circle', 'square')[1 + V(pk.gene)$name %in% pk.ids[pk.ids$go_term=="phosphoprotein phosphatase activity","gene_id"]]
  V(pk.gene)$color <- c('black', 'red')[1 + V(pk.gene)$name %in% pk.ids[pk.ids$go_term=="phosphoprotein phosphatase activity","gene_id"]]
  V(pk.gene)$degree <- scale(degree(pk.gene), center=F, scale=T)
  E(pk.gene)$weight.scaled <- scale(E(pk.gene)$weight, center=F, scale=T)
  V(pk.gene)$label <- V(pk.gene)$name # default
  V(pk.gene)[rank(degree(pk.gene), ties.method="max") < length(V(pk.gene)) - 10]$label <- NA # replace everything but the top 10 with NA
  V(pk.gene)$shape.2 <- V(pk.gene)$shape
  V(pk.gene)[!is.na(V(pk.gene)$label)]$shape.2 <- 'none' # is it has a label, then don't plot a shape
  V(pk.gene)$label.pp <- V(pk.gene)$name # default
  V(pk.gene)[V(pk.gene)$shape == 'circle']$label.pp <- NA # replace kinases with NA
  V(pk.gene)$shape.pp <- V(pk.gene)$shape
  V(pk.gene)[!is.na(V(pk.gene)$label.pp)]$shape.pp <- 'none' # is it has a label, then don't plot a shape
  # plot
  pdf("output/example6.PK_RIFome.comp1.genes.pdf")
  plot(
    pk.gene,                          #the graph to be plotted
    layout=layout.kamada.kawai,      # kamada.kawai|fruchterman.reingold|lgl|auto -  the layout method. see the igraph documentation for details
    main='largest component, cocitation, kinases - black dots, PPases - red boxes',                    #specifies the title
    vertex.shape=V(pk.gene)$shape,                    #nodeshape
    vertex.size=3,                           #nodesize  
    # vertex.label.dist=0.5,    	            #puts the name labels slightly off the dots
    vertex.color=V(pk.gene)$color,          #node colour
    # vertex.frame.color='blue', 		          #the color of the border of the dots 
    vertex.label.color=V(pk.gene)$color,		          #the color of the name labels
    vertex.label.font=2,			              #the font of the name labels
    vertex.label=NA,		                    #specifies the lables of the vertices. in this case the 'name' attribute is used
    vertex.label.cex=1,			                #specifies the size of the font of the labels. can also be made to vary
    
    edge.width=E(pk.gene)$weight.scaled      #specifies the thickness of the edges
    
  )
  
  return(network)
}

do.mesh <- function(term="Wnt Signaling Pathway", PATH, plot_results=TRUE){
  
  ### as for do.phospho, 
  ### pulling out papers related to a given MeSh term

  # get all pmids linked to mesh term
  res <- xmlParse(paste(c(
    "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=",
    term, 
    "&field=MeSH%20Terms&rettype=uilist&retmode=text&retmax=100000"), sep='', collapse=''))
  res <- xmlToList(res)
  res <- as.integer(as.vector(c(res$IdList, recursive=TRUE)))
  if(VERBOSE){cat(sprintf('%i PMIDS found for term \"%s\"',length(res), term))}
  
  ## load edges, subset to retrieved terms and load graph
  network <- do.load.edges(PATH)
  edges <- subset(edges, gene_id %in% pk.ids$gene_id)
  network<-graph.data.frame(
    d=edges,
    directed=F)
  
  ## decompose and run cocitation
  network <- do.decompose(network)
  V(network)$isa <- c("gene", "paper")[1 + V(network) %in% edges$pubmed_ids]
  network <- do.cocitation(network)
  
  ## add some community detection calls!
  V(network)$pr <- scale(page.rank(network)$vector, center=F, scale=T) # page rank
  V(network)$walktrap <- walktrap.community(network)$membership # community memberships
  V(network)$eigen <- leading.eigenvector.community(network)$membership
  V(network)$spin <- spinglass.community(network)$membership
  
  # for plotting
  V(network)$shape <- 'circle'
  V(network)$color <- 'black'
  V(network)$degree <- scale(degree(network), center=F, scale=T)
  E(network)$weight.scaled <- scale(E(network)$weight, center=F, scale=T)
  V(network)$label <- V(network)$name # default
  V(network)[rank(degree(network), ties.method="max") < length(V(network)) - 10]$label <- NA # replace everything but the top 10 with NA
  V(network)$shape.2 <- V(network)$shape
  V(network)[!is.na(V(network)$label)]$shape.2 <- 'none' # is it has a label, then don't plot a shape
    
  # NB layout.lgl seems to work well with coloring by community
  plot(
    network,                          #the graph to be plotted
    layout=layout.lgl,      # kamada.kawai|fruchterman.reingold|lgl|auto -  the layout method. see the igraph documentation for details
    main='cocitation, top 10 by degree',                    #specifies the title
    vertex.shape=V(network)$shape.2,                    #nodeshape
    vertex.size=V(network)$degree*4,                           #nodesize  
    # vertex.label.dist=0.5,  		            #puts the name labels slightly off the dots
    vertex.color=V(network)$spin,          #node colour
    # vertex.frame.color='blue', 		          #the color of the border of the dots 
    vertex.label.color=V(network)$color,		          #the color of the name labels
    vertex.label.font=2,			              #the font of the name labels
    vertex.label=V(network)$label,		                    #specifies the lables of the vertices. in this case the 'name' attribute is used
    vertex.label.cex=1,			                #specifies the size of the font of the labels. can also be made to vary
    
    edge.width=E(network)$weight.scaled      #specifies the thickness of the edges

  )
  return(network)
}

### calls

PATH <- do.download(INPUT, DATA_DIR, DOWNLOAD=DOWNLOAD)

network <- do.load(PATH)

network <- do.decompose(network)
network <- do.cocitation(network)

do.phosho(DATA_DIR, PATH)

# look up some interesting mesh terms
# http://www.ncbi.nlm.nih.gov/mesh?term=autism

lapply(c("Wnt Signaling Pathway", "Autistic Disorder", "Melanoma", "Hedgehog Proteins"),
       function(x) do.mesh(term=x, PATH)
       )


## output results for comparison
# check output directories exist
check_generated()
# write results to file
report_results(RESULTS = RESULTS, BENCHMARK = BENCHMARK)

# timings
report_timings(TIMES = TIMES, BENCHMARK = BENCHMARK)

# final clean up
rm(list=ls())
gc()
