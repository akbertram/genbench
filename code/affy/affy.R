# Parham Solaimani parham@bedatadriven.com
# Test-case workflow for analysis of Affymatrix microarry data developed by http://www.arrayanalysis.org

##### set up session #####
rm(list=ls())
set.seed(1000)
stopifnot(file.exists(file.path("..","..", "data"))) # data path is relative
.libPaths(file.path("~","R","libs"))
source(file.path("..", "..","benchmark_utilities.R"))

##### loading packages #####
library("ArrayTools")
library("affy", quietly = TRUE)
library("affycomp", quietly = TRUE)
library("affyPLM", quietly = TRUE)
library("affypdnn", quietly = TRUE)
library("bioDist", quietly = TRUE)
library("simpleaffy", quietly = TRUE)
library("affyQCReport", quietly = TRUE)
library("plier", quietly = TRUE)
if(exists("samplePrep")) library("yaqcaffy", quietly = TRUE)
library("gdata", quietly = TRUE) #trim function
library("gplots", quietly = TRUE) #heatmap.2 functions
print("Libraries have been loaded")


##### Set global vars #####
VERBOSE <- TRUE # print progress?
DOWNLOAD <- FALSE # download fresh data?
BENCHMARK <- "affy"
DATA_DIR <- file.path("..","..", "data","affy")
RESULTS <- genbench_results(benchmark_name = BENCHMARK)
TIMES <- genbench_timings(benchmark_name = BENCHMARK)

affyPARAM = list()
affyPARAM$arrayGroup <- "" #see comments below
affyPARAM$reOrder <- TRUE
affyPARAM$layoutPlot <- TRUE
affyPARAM$controlPlot <- TRUE
affyPARAM$samplePrep <- TRUE
affyPARAM$ratio <- TRUE
affyPARAM$degPlot <- TRUE
affyPARAM$hybrid <- TRUE
affyPARAM$percPres <- TRUE
affyPARAM$posnegDistrib <- TRUE
affyPARAM$bgPlot <- TRUE
affyPARAM$scaleFact <- TRUE
affyPARAM$boxplotRaw <- TRUE
affyPARAM$boxplotNorm <- TRUE
affyPARAM$densityRaw <- TRUE
affyPARAM$densityNorm <- TRUE
affyPARAM$MARaw <- TRUE
affyPARAM$MANorm <- TRUE
affyPARAM$MAOption1 <- "dataset" #see comments below
affyPARAM$spatialImage <- TRUE
affyPARAM$PLMimage <- FALSE
affyPARAM$posnegCOI <- TRUE
affyPARAM$Nuse <- TRUE
affyPARAM$Rle <- TRUE
affyPARAM$correlRaw <- TRUE
affyPARAM$correlNorm <- TRUE
affyPARAM$clusterRaw <- TRUE
affyPARAM$clusterNorm <- TRUE
affyPARAM$clusterOption1 <- "Spearman" #see comments below
affyPARAM$clusterOption2 <- "ward" #see comments below
affyPARAM$PCARaw <- TRUE
affyPARAM$PCANorm <- TRUE
affyPARAM$PMAcalls <- FALSE
affyPARAM$normMeth <- "RMA" #see comments below
affyPARAM$normOption1 <- "dataset" #see comments below
affyPARAM$customCDF <- TRUE
affyPARAM$CDFtype <- "ENSG" #see comments below
affyPARAM$species <- "Hs" #see comments below
affyPARAM$version_nb <- "1.0.0"
affyPARAM$files <- c( "GSM1103973_VEH_VEH1.CEL.gz", "GSM1103974_DOX_VEH1.CEL.gz",
                      "GSM1103975_VEH_PMA1.CEL.gz", "GSM1103976_DOX_PMA1.CEL.gz",
                      "GSM1103977_VEH_VEH2.CEL.gz", "GSM1103978_DOX_VEH2.CEL.gz",
                      "GSM1103979_VEH_PMA2.CEL.gz","GSM1103980_DOX_PMA2.CEL.gz",
                      "GSM1103981_VEH_VEH3.CEL.gz","GSM1103982_DOX_VEH3.CEL.gz",
                      "GSM1103983_VEH_PMA3.CEL.gz","GSM1103984_DOX_PMA3.CEL.gz")

##### Functions #####


##### Blocks for timing #####
do.load <- function(DATA_DIR,affyPARAM){
  affyPARAM <- affyPARAM
  # Load the required Data


  getArrayType <- function(Data) {

    aType <- "PMMM"

    # Test whether the dataset is of aType "PM-only"
    mismatches <- mm(Data[,1])

    if(is.null(mismatches)) {
      # mm does not exist
      aType <- "PMonly"
    } else {
      if(sum(is.na(mismatches))>=length(mismatches)/2){
        # mm is always NA or there are more NA values in the mm probes than
  	  # defined values (assuming these would just be controls)
        aType <- "PMonly"
      } else {
        matches <- pm(Data[,1])
        notNA <- !is.na(mismatches) & !is.na(matches)
        if(sum(mismatches[notNA]!=matches[notNA])==0){
          # MM contains a copy of PM, which indicates a PMonly array
          aType <- "PMonly"
        }
      }
    }

    cat("The arrays are determined to contain", ifelse(aType=="PMMM",
    "perfect match and mismatch probes\n", "perfect match probes only\n"))

    return(aType)
  }


  #######################
  ## addStandardCDFenv ##
  #######################

  addStandardCDFenv <- function(Data, overwrite=FALSE) {
    #if overwrite is FALSE a cdf environment will be kept if already loaded
    #if overwrite is TRUE it will always be overwritten (unless none is found)
    #the first option would be the regular one, used to add cdf environments where
    #   automatic loading failed
    #the second could be used to set back updated cdf files to the standard ones

    #start with r0 just in case this could exist
    rev <- 0

    #initial value
    CDFenv <- 0

    # recall which cdfName was added, in case no other one is found (set back even
    #   if it does not exist)
    presetCDF <- Data@cdfName

    #check whether environment is already correct
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))

    #try the annotation plus cdf as a cdf name
    if ((class(CDFenv)!="environment") | overwrite) {
      CDFenv <- 0 #needed for cases where overwrite is TRUE, but CDFenv already
                  # was an environment
      Data@cdfName<-paste(Data@annotation,".cdf",sep="")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
      #if no succes also try without a dot
      if (class(CDFenv)!="environment") {
        Data@cdfName<-paste(Data@annotation,"cdf",sep="")
        suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
      }
    }

    # don't run the loop if CDF is already known up front, or correct one has been
    # found
    while((class(CDFenv)!="environment") & (rev < 10)) {
      Data@cdfName<-paste(Data@annotation,".r",rev,"cdf",sep="")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
      rev <- rev + 1
    }

    if ((class(CDFenv)!="environment")) {
      Data@cdfName <- presetCDF
      warning("could not automatically retrieve CDF environment for this chip type - object kept as is")
    }

    cat("current cdf environment loaded:",Data@cdfName,"\n")
    return(Data)
  }


  ######################
  ## addUpdatedCDFenv ##
  ######################

  addUpdatedCDFenv <- function(Data, species=NULL, type="ENSG") {
    # note: this function will add an updated cdf environment to the data object
    # and will overwrite a possible already loaded environment, unless no updated
    # cdf environment is found

    # developer's note: it may be of interest to find out whether available
    # species and types can be retrieved automatically from the brainarray website

    if(is.null(species) || (species=="")) stop("The species must be provided")

    types <- c("ENTREZG","REFSEQ","ENSG","ENSE","ENST","VEGAG","VEGAE","VEGAT",
            "TAIRG","TAIRT","UG","MIRBASEF","MIRBASEG")
    if(!tolower(type) %in% tolower(types)) {
      stop("selected type not valid, select from ", paste(types, collapse=" "))
    } else {
      type <- types[match(tolower(type),tolower(types))]
    }

    spp <- c("Ag","At","Bt","Ce","Cf","Dr","Dm","Gg","Hs","MAmu","Mm","Os","Rn",
             "Sc","Sp","Ss")
    names(spp) <- c("Anopheles gambiae","Arabidopsis thaliana","Bos taurus",
             "Caenorhabditis elegans","Canis familiaris", "Danio rerio",
             "Drosophila melanogaster","Gallus gallus","Homo sapiens",
             "Macaca mulatta","Mus musculus", "Oryza sativa","Rattus norvegicus",
             "Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa")
    if(tolower(species) %in% tolower(names(spp)))
         species <- spp[tolower(names(spp))==tolower(species)]
    if(!tolower(species) %in% tolower(spp)) {
      stop("selected species not valid, select from:\n",
            paste(names(spp), collapse="\n"), "\nor abbreviated as ",
            paste(spp,collapse=" "))
    } else {
      species <- spp[match(tolower(species),tolower(spp))]
    }

    #initial value
    CDFenv <- 0

    # recall which cdfName was added, in case no updated one is found (set back
    # even if it does not exist)
    presetCDF <- Data@cdfName

    #try to find updated cdf file of choice
    print(Data@cdfName<-paste(Data@annotation,species,type,sep="_"))
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    #try without a version number
    print(Data@cdfName<-paste(gsub("v[0-9]$","",Data@annotation),species,type,sep="_"))
    suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))

    #if it hasn't loaded, try to download
    if ((class(CDFenv)!="environment")) {
      install.packages(tolower(paste(Data@annotation,species,type,"cdf",sep="")),
        repos="http://brainarray.mbni.med.umich.edu/bioc")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    }

    #if it hasn't loaded, try to download without version number
    if ((class(CDFenv)!="environment")) {
      install.packages(tolower(paste(gsub("v[0-9]$","",Data@annotation),species,type,"cdf",sep="")),
        repos="http://brainarray.mbni.med.umich.edu/bioc")
      suppressWarnings(try(CDFenv <- getCdfInfo(Data),TRUE))
    }

    if ((class(CDFenv)!="environment")) {
      Data@cdfName <- presetCDF
      warning("Could not automatically retrieve CDF environment for this chip type - object kept as is")
    }

    cat("current cdf environment loaded:",Data@cdfName,"\n")
    return(Data)
  }


  ####################
  ## colorsByFactor ##
  ####################

  #create colors for the plots and the legends
  #-------------------------------------------

  colorsByFactor <- function(experimentFactor) {

    #check whether a factor has been provided
    if(class(experimentFactor)!="factor") stop("Parameter 'experimentFactor' must be of class 'factor'")

    if(length(levels(experimentFactor))==1) {
      #if there is only one group (or no groups are provided) take equally spread colors over the rainbow palette
      plotColors <- rainbow(length(experimentFactor),s=.8,v=.7)
  	#set group legend color to white, as there is not a specific group color
  	legendColors <- "white"
    } else {
      #compute the number of colors needed for each class
      tab.tmp <- table(experimentFactor)

      #set the two extreme colors for each class
      colors.light <- rainbow(length(levels(experimentFactor)),s=1-sapply(tab.tmp,min,5)*.1)
      colors.dark <- rainbow(length(levels(experimentFactor)),v=1-sapply(tab.tmp,min,5)*.14)

      #create the colors to plot, and colors for the legend (average one per experimental group)
      plotColors <- NULL
      legendColors <- NULL
      for(l in 1:length(levels(experimentFactor))) {
        colorFun <- colorRampPalette(c(colors.light[l],colors.dark[l]))
        tmpColors <- colorFun(tab.tmp[l])
        plotColors[experimentFactor==levels(experimentFactor)[l]] <- tmpColors
        legendColors[l] <- tmpColors[ceiling(length(tmpColors)/2)]
      }
    }
    return(list(plotColors=plotColors,legendColors=legendColors))

  }


  ###################
  ## deduceSpecies ##
  ###################

  #try to find the species of the chiptype
  #---------------------------------------

  deduceSpecies <- function(descr=NULL) {

    organism <- ""

    if(!is.null(descr) && (descr!="")) {
      try({lib <- paste(descr,".db",sep="")
        eval(parse("",-1,paste("require(",lib,")",sep="")))
        organism <- eval(parse("",-1,paste(descr,"ORGANISM",sep="")))},TRUE)
      if(organism=="") {
        descr <- tolower(descr)
        if(length(grep("^nugomm",descr)) > 0) organism <- "Mus musculus"
        if(length(grep("^nugohs",descr)) > 0) organism <- "Homo Sapiens"
        if(length(grep("^hgu133plus2",descr)) > 0) organism <- "Homo Sapiens"
        if(length(grep("^hugene",descr)) > 0) organism <- "Homo Sapiens"
        if(length(grep("^mogene",descr)) > 0) organism <- "Mus musculus"
        if(length(grep("^ragene",descr)) > 0) organism <- "Rattus norvegicus"
      }
    }

    return(organism)
  }


  ###################
  ## normalizeData ##
  ###################

  #normalize the data set
  #----------------------

  normalizeData <- function(Data, normMeth="", perGroup=FALSE, experimentFactor=NULL,
    customCDF=TRUE, species=NULL, CDFtype=NULL, aType=NULL, WIDTH=1000, HEIGHT=1414) {

    if((normMeth=="") || is.null(normMeth)) stop("normMeth, the requested normalization method, must be provided")
    normMeth <- toupper(normMeth)

  	if(customCDF) {
  		if(is.null(species) || species=="") stop("When customCDF is required, the species must be provided")
  		if(is.null(CDFtype) || CDFtype=="") stop("When customCDF is required, the CDFtype must be provided")
  	}
  	if(perGroup) {
  		if(is.null(experimentFactor)) stop("When normalization per group is requested, experimentFactor must be provided")
  	}
    if((normMeth=="PLIER") && (is.null(aType))) stop("When selecting PLIER normalization, aType must be provided")
    if((normMeth=="GCRMA") && (is.null(aType))) stop("When selecting GCRMA normalization, aType must be provided")

  	#if customCDF option is chosen, apply to copy of Data, in order not to change the original data object
  	Data.copy <- Data
  	if(customCDF){
  		print ("Change CDF before pre-processing")
  		Data.copy <- addUpdatedCDFenv(Data.copy, species, CDFtype)
  	}
  	print ("Pre-processing is running")

    nGroups <- 1
    if(perGroup) {
      nGroups <- max(1,length(levels(experimentFactor)))
      if(nGroups==1) warning("normalization per group requested, but no groups indicated in data set")
    }

    #if per group normalization required, or a method selected that does not return an ExpressionSet object,
    #make a model of class ExpressionSet to paste real values in, use the relatively fast RMA method
    #note that binding of ExpressionSet objects is NOT possible
    if((nGroups>1)) { # || (normMeth=="MAS5")) {
      normData <- rma(Data.copy)
      exprs(normData)[] <- NA
    }

    for(group in 1:nGroups) {
      if(nGroups==1) {
        Data.tmp <- Data.copy
      } else {
        Data.tmp <- Data.copy[,experimentFactor==(levels(experimentFactor)[group])]
      }
      switch(normMeth,
        "MAS5" = {
        #doesn't work
        normData.tmp <- mas5(Data.tmp)
        },
        "GCRMA" = {
        if(customCDF) {
          #probe library needed, first try whether this has been intalled, otherwise do so
          probeLibrary <- tolower(paste(Data@annotation,species,CDFtype,"probe",sep=""))
          loaded <- suppressWarnings(try(eval(parse("",-1,paste("library(",probeLibrary,")", sep=""))),TRUE))
          if(class(loaded)=="try-error") {
            install.packages(probeLibrary, repos="http://brainarray.mbni.med.umich.edu/bioc")
          }
        }
        if(aType == "PMMM") ntype = "fullmodel"
        if(aType == "PMonly") ntype = "affinities" # good results if most of the genes are not expressed
  	  normData.tmp <- gcrma(Data.tmp, type=ntype, fast = FALSE)
        },
        "RMA" = {
        normData.tmp <- rma(Data.tmp)
        },
        "PLIER" = {
        if(aType == "PMMM") ntype = "together"
        if(aType == "PMonly") ntype = "pmonly"
        normData.tmp <- justPlier(Data.tmp, normalize=TRUE, norm.type = ntype)
        }
      )
      if(nGroups==1) {
          normData <- normData.tmp
  		if(normMeth=="MAS5") exprs(normData)<-log2(exprs(normData))
      } else {
        try(
  	  if(normMeth=="MAS5"){
  		exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- log2(exprs(normData.tmp))
  	  }else{
  		exprs(normData)[,match(sampleNames(normData.tmp), sampleNames(normData))] <- exprs(normData.tmp)
  	  },TRUE)
      }
      rm(normData.tmp, Data.tmp)
    }

    #create an 'inter sheet'
    png(file="Cover_2.png", width=WIDTH, height=HEIGHT)
    plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE,
       frame.plot = TRUE, xlim = c(0, 2), ylim = c(0,2))
    text(1,1,"Pre-processing of Raw Data\n\n\n",cex=3)
    text(1,1,paste("\n\nMethod: ",normMeth,"\nAnnotation: ",Data.copy@cdfName),cex=2.5)
    if(perGroup) text(1,1,paste("\n\n\n\n\nNormalization per experimental group"),cex=2.5)
    dev.off()

  	rm(Data.copy)
  	return(normData)
  }

  #####################
  ## computePMAtable ##
  #####################

  #prepare and return a table of PMA calls
  #---------------------------------------

  computePMAtable <- function(Data, customCDF=TRUE, species=NULL, CDFtype=NULL) {

  	if(customCDF) {
  		if(is.null(species)) stop("When customCDF is used, the species must be provided")
  		if(is.null(CDFtype)) stop("When customCDF is used, the CDFtype must be provided")
    }

    Data.copy <- Data
    if(customCDF) {
      Data.copy <- addUpdatedCDFenv(Data.copy, species, CDFtype)
  	}

    PMAtable <- NULL

    #ocmpute the calls using the detection.p.val function from simpleaffy
    try(PMAtable <- detection.p.val(Data.copy)$call,TRUE)

    if(!is.null(PMAtable)) {
      #add column of IDs to table
      PMAtable <- cbind(rownames(PMAtable),PMAtable)

      #remove "_at" from custom probeset IDs to get to real ID
      if(customCDF) {
        control_rows <- grep("affx",tolower(PMAtable[,1]))
        if(length(control_rows) > 0) {
          PMAtable[-control_rows,1] <- substring(PMAtable[-control_rows,1],1,nchar(PMAtable[-control_rows,1])-3)
        } else {
          PMAtable[,1] <- substring(PMAtable[,1],1,nchar(PMAtable[,1])-3)
        }
      }

      #add column names to PMAtable
      colnames(PMAtable)[1] <- ifelse(customCDF,paste(CDFtype,"_ID",sep=""),"Probeset_ID")
    } else {
      warning("PMA table could not be computed for this arraytype")
    }

    return(PMAtable)
  }


  #########################
  ## createNormDataTable ##
  #########################

  #prepare and export the normalized data table
  #--------------------------------------------

  createNormDataTable <- function(normData, customCDF=NULL, species=NULL, CDFtype=NULL) {

    if(is.null(customCDF)) stop("The customCDF parameter must be provided")
  	if(customCDF) {
  		if(is.null(CDFtype)) stop("When customCDF is used, the CDFtype must be provided")
  		if(species=="" || is.null(species)) {
  		  warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
  		  species <- deduceSpecies(rawData@annotation)
  		}
  		if(species=="" || is.null(species)){
  			warning("Could not define species; the CDF will not be changed")
  			customCDF<-FALSE
  		}
  	}

    #add column of IDs and normalized data to normDataTable
    normDataTable<-cbind(rownames(exprs(normData)),exprs(normData))

    #remove "_at" from custom probeset IDs to get to real ID
    if(customCDF) {
      control_rows <- grep("affx",tolower(normDataTable[,1]))
      if(length(control_rows)>0) {
        normDataTable[-control_rows,1]<-substring(normDataTable[-control_rows,1],1,nchar(normDataTable[-control_rows,1])-3)
      } else {
        normDataTable[,1]<-substring(normDataTable[,1],1,nchar(normDataTable[,1])-3)
      }
    }

    #add column names to normDataTable
    colnames(normDataTable)<-c(ifelse(customCDF,paste(CDFtype,"_ID",sep=""),"Probeset_ID"),colnames(exprs(normData)))

    #add gene name and description in case ensembl IDs have been used (otherwise there is no 1 to 1 mapping)
    if(customCDF && CDFtype=="ENSG") {
      #load gene name and description annotations
      library(biomaRt)

      spName <- ""
      if(species=="Ag" || species=="Anopheles gambiae") spName <- "agambiae"
      if(species=="At" || species=="Arabidopsis thaliana") spName <- "athaliana"
      if(species=="Bt" || species=="Bos taurus") spName <- "btaurus"
      if(species=="Ce" || species=="Caenorhabditis elegans") spName <- "celegans"
      if(species=="Cf" || species=="Canis familiaris") spName <- "cfamiliaris"
      if(species=="Dr" || species=="Danio rerio") spName <- "drerio"
      if(species=="Dm" || species=="Drosophila melanogaster") spName <- "dmelanogaster"
      if(species=="Gg" || species=="Gallus gallus") spName <- "ggallus"
      if(species=="Hs" || species=="Homo sapiens") spName <- "hsapiens"
      if(species=="MAmu" || species=="Macaca mulatta") spName <- "mmulatta"
      if(species=="Mm" || species=="Mus musculus") spName <- "mmusculus"
      if(species=="Os" || species=="Oryza sativa") spName <- "osativa"
      if(species=="Rn" || species=="Rattus norvegicus") spName <- "rnorvegicus"
      if(species=="Sc" || species=="Saccharomyces cerevisiae") spName <- "scerevisiae"
      if(species=="Sp" || species=="Schizosaccharomyces pombe") spName <- "spombe"
      if(species=="Ss" || species=="Sus scrofa") spName <- "sscrofa"

      try(ensembl <- useMart("ensembl", dataset = paste(spName,"_gene_ensembl",sep="")))
      if(exists("ensembl")) {
        try(annotationTable<-getBM(attributes=c("ensembl_gene_id","external_gene_name","description"), mart=ensembl, uniqueRows=TRUE),TRUE)
      }
      if(exists("annotationTable")) {
        normDataTable <- as.data.frame(normDataTable,stringsAsFactors=FALSE)
        suppressWarnings(normDataTable<-cbind(normDataTable,annotationTable[match(normDataTable[,1],annotationTable[,1]),2:(dim(annotationTable)[2])]))
        normDataTable[,2:(dim(exprs(normData))[2]+1)] <- apply(normDataTable[,2:(dim(exprs(normData))[2]+1),drop=FALSE],2,as.numeric)
      } else {
        warning("No gene names and annotation could be retrieved from BioMart for this species or no connection could be established, gene information not added to normalized data table")
      }
    }

    return(normDataTable)
  }
  #=============================================================================#
  # ArrayAnalysis - affyAnalysisQC                                              #
  # a tool for quality control and pre-processing of Affymetrix array data      #
  #                                                                             #
  # Copyright 2010-2011 BiGCaT Bioinformatics                                   #
  #                                                                             #
  # Licensed under the Apache License, Version 2.0 (the "License");             #
  # you may not use this file except in compliance with the License.            #
  # You may obtain a copy of the License at                                     #
  #                                                                             #
  # http://www.apache.org/licenses/LICENSE-2.0                                  #
  #                                                                             #
  # Unless required by applicable law or agreed to in writing, software         #
  # distributed under the License is distributed on an "AS IS" BASIS,           #
  # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
  # See the License for the specific language governing permissions and         #
  # limitations under the License.                                              #
  #=============================================================================#

  #####################
  ## coverAndKeyPlot ##
  #####################

  coverAndKeyPlot <- function(description=NULL, refName="", WIDTH=1000,
  	HEIGHT=1414) {

    if(is.null(description)) stop("The description parameters is required")

    #create a 'cover sheet'
    png(file = "Cover_1.png", width=WIDTH, height=HEIGHT)
      plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE,
         frame.plot = TRUE, xlim = c(0,2), ylim = c(0,2))
      text(1,1,"Quality Control & Pre-processing Evaluation\n\n\n\n",pos=3,cex=3)
  	  if(refName != ""){
  		refName <- sub("(_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}_\\d{2})", "",
  			refName)
  		text(1,1,paste("of\n",refName,"\n\n"),pos=3,cex=2.5)
  	  }
      text(1,1,"\n\n\n\n\n\n\n\n\nREPORT",pos=3,cex=5)
    dev.off()

    # create a page indicating the naming and grouping used (it corresponds to the
    # description file)
    linesPerPage <- 35
    noOfPages <- (dim(description)[1] %/% linesPerPage) +
  	((dim(description)[1] %% linesPerPage) > 0)
    for(j in 1:noOfPages) {
    png(file = paste("Description",ifelse(refName!="","_",""),refName,
  	ifelse(noOfPages>1,paste("_",j,sep=""),""),".png",sep=""),
  	width=WIDTH,height=HEIGHT)
        par(oma=c(0,0,0,0))
        plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE,
           frame.plot = TRUE, xlim = c(0,2), ylim = c(0,2))
        text(1,1.90,"Array names and grouping",cex=1.6,pos=3)

        vertical_pos <- seq(1.8, 0.01, length.out= (linesPerPage +1))
        line_coef <- (vertical_pos[1]-vertical_pos[2])/2
        k <- 1
        for (i in c(0,((j-1)*linesPerPage+1):min(j*(linesPerPage),
  		dim(description)[1]))) {
          if(i==0) {
            value <- colnames(description)
          } else {
            value <- description[i,]
          }
          text(0.0,vertical_pos[k],value[1],pos=4,cex=1)
          text(0.9,vertical_pos[k],value[2],pos=4,cex=1)
          text(1.6,vertical_pos[k],value[3],pos=4,cex=1)
          abline(h=vertical_pos[k]-line_coef)
          k <- k+1
        }
      dev.off()
    }
  }


  #################
  ## QCtablePlot ##
  #################

  QCtablePlot <- function(Data,quality=NULL,sprep=NULL,lys=NULL,
    samplePrep=TRUE,ratio=TRUE,hybrid=TRUE,percPres=TRUE,bgPlot=TRUE,
    scaleFact=TRUE, WIDTH=1000, HEIGHT=1414, POINTSIZE=24) {

    rows <- NULL
    indicator <- NULL
    indicCol <- NULL

    if(sum(ratio,hybrid,percPres,bgPlot,scaleFact) > 0 && is.null(quality)) {
      try(quality <- qc(Data),TRUE) #calculate Affymetrix quality data for PMMM (use MAS5)
      if(is.null(quality)) {
        warning("Plots based on the simpleaffy qc function cannot be created for this chip type")
      }
    }

    if(samplePrep){
      if(is.null(sprep)) {
        # find the data
        try(yack <- yaqc(Data),TRUE)
        if(exists("yack")) {
          spnames<-rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3' !
          rownames(yack@morespikes), ignore.case = TRUE),])
          sprep<-t(yack@morespikes[spnames,])
        } else {
          warning("Values based on the yaqc function cannot be computed for this chip type")
        }
      }


      if(is.null(lys)) {
        try({calls<-detection.p.val(Data)$call
        lys<-calls[rownames(calls)[grep("lys.*3",rownames(calls), ignore.case = TRUE)],]
        rm(calls)},TRUE)
        if(is.null(lys)) {
          warning("Values based on the detection.p.val function cannot be computed for this chip type")
        }else{
  		if(length(lys) > length(sampleNames(Data))) { lys<-lys[1,] }
        }
      }

      if(!is.null(sprep)) {
        print ("indicators for sample prep controls"  )

        # test for Lys < Phe < Thr < Dap
        t1 <- (sprep[,1]<sprep[,2] & sprep[,2]<sprep[,3]
              & sprep[,3]<sprep[,4])
        Col1 <- t1
        Col1[t1==TRUE]<-"blue"
        Col1[t1==FALSE]<-"red"
        Text1 <- t1
        Text1[t1==TRUE]<-"T"
        Text1[t1==FALSE]<-"F"

        if(exists("rows")){
          rows<-c(rows,"Sample Prep\nLys<Phe<Thr<Dap if\ncontrols were spiked-in")
        } else{
          rows<-c("Sample Prep\nLys<Phe<Thr<Dap if\ncontrols were spiked-in")
        }

        if(exists("indicator")){
          indicator<-rbind(indicator,Text1)
        } else{
          indicator<-rbind(Text1)
        }
        if(exists("indicCol")){
          indicCol<-rbind(indicCol,Col1)
        } else{
          indicCol<-rbind(Col1)
        }
      }

      if(!is.null(lys)) {
        # test for Lys presence
        t2 <- lys == 'P'
        Col2 <- t2
        Col2[t2==TRUE]<-"blue"
        Col2[t2==FALSE]<-"red"

        if(exists("rows")){
          rows<-c(rows,"Sample Prep\nLys=Present")
        } else{
          rows<-c("Sample Prep\nLys=Present")
        }
        if(exists("indicator")){
  		indicator<-rbind(indicator,lys)
        } else{
  		indicator<-rbind(lys)
        }
        if(exists("indicCol")){
          indicCol<-rbind(indicCol,Col2)
        } else{
          indicCol<-rbind(Col2)
        }
      }
    }

    if(!is.null(quality)) {
      if(ratio){
        print ("indicators for beta-actin & GAPDH 3'/5' ratio")

        ratio35_beta <- quality@qc.probes[,1]/quality@qc.probes[,3]
        ratio35_GAPDH <- quality@qc.probes[,4]/quality@qc.probes[,6]


        if(exists("rows")){
          rows<-c(rows,"3'/5' beta-actin\n(cutoff=3)",
              "3'/5' GAPDH\n(cutoff=1.25)")
        } else{
          rows<-c("3'/5' beta-actin\n(cutoff=3)","3'/5' GAPDH\n(cutoff=1.25)")
        }

        # beta-actin
        t1 <- ratio35_beta <= 3
        Col1 <- t1
        Col1[t1==TRUE]<-"blue"
        Col1[t1==FALSE]<-"red"


        # GAPDH
        t2 <- ratio35_GAPDH <= 1.25
        Col2 <- t2
        Col2[t2==TRUE]<-"blue"
        Col2[t2==FALSE]<-"red"

        if(exists("indicator")){
          indicator<-rbind(indicator,round(ratio35_beta,2),round(ratio35_GAPDH,2))
        } else{
          indicator<-rbind(round(ratio35_beta,2),round(ratio35_GAPDH,2))
        }
        if(exists("indicCol")){
          indicCol<-rbind(indicCol,Col1,Col2)
        } else{
          indicCol<-rbind(Col1,Col2)
        }
      }

      if(hybrid){
        print ("indicators for spike-in hybridization controls"  )

        # test for BioB < BioC < BioD < CreX
        t1 <- (quality@spikes[,1]<quality@spikes[,2]
                     & quality@spikes[,2]<quality@spikes[,3]
                     & quality@spikes[,3]<quality@spikes[,4])
        Col1 <- t1
        Col1[t1==TRUE]<-"blue"
        Col1[t1==FALSE]<-"red"
        Text1 <- t1
        Text1[t1==TRUE]<-"T"
        Text1[t1==FALSE]<-"F"

        # test for BioB presence
        t2 <- quality@bioBCalls == 'P'
        Col2 <- t2
        Col2[t2==TRUE]<-"blue"
        Col2[t2==FALSE]<-"red"

        if(exists("rows")){
          rows<-c(rows,"Hybridization\nBioB<BioC<BioD<CreX",
              "Hybridization\nBioB=Present")
        } else{
          rows<-c("Hybridization\nBioB<BioC<BioD<CreX",
             "Hybridization\nBioB=Present")
        }
        if(exists("indicator")){
  		indicator<-rbind(indicator,Text1,quality@bioBCalls)
        } else{
  		indicator<-rbind(Text1,quality@bioBCalls)
        }
        if(exists("indicCol")){
          indicCol<-rbind(indicCol,Col1,Col2)
        } else{
          indicCol<-rbind(Col1,Col2)
        }
      }

      if(percPres){
        print ("indicators for percent present"  )

  		# Define the reference (it should minimize the number of outliers)
  		ppMin <- min(quality@percent.present)
  		ppMax <- max(quality@percent.present)
  		ppMean<- mean(quality@percent.present)

  		ref<-cbind(c(ppMean-5,ppMin,ppMax-10),c(ppMean+5,ppMin+10,ppMax))
  		t1 <- (quality@percent.present >= (ppMean-5) & quality@percent.present <= (ppMean+5))
  		t2 <- (quality@percent.present >= ppMin & quality@percent.present <= (ppMin+10))
  		t3 <- (quality@percent.present >= (ppMax-10) & quality@percent.present <= ppMax)
  		outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
  		if(outliers[1] == 0){
  			ref<- ref[1,]
  		}else{
  			ref<- ref[outliers==min(outliers),] # less outliers
  			if(length(ref) > 3) { ref<-ref[1,] }
  		}
  		testMin <- ref[1]
  		testMax <- ref[2]
  		t1 <- (quality@percent.present >= testMin & quality@percent.present <= testMax)

        Col <- t1
        Col[t1==TRUE]<-"blue"
        Col[t1==FALSE]<-"red"

        if(exists("rows")){
          rows<-c(rows,"Percent Present\nspread<=10%")
        } else{
          rows<-c("Percent Present\nspread<=10%")
        }
        if(exists("indicator")){
          indicator<-rbind(indicator,paste(round(quality@percent.present,0),"%"))
        } else{
          indicator<-paste(round(quality@percent.present,0),"%")
        }
        if(exists("indicCol")){
          indicCol<-rbind(indicCol,Col)
        } else{
          indicCol<-rbind(Col)
        }

      }

      if(bgPlot){
        print ("indicators for background intensities"  )

  		# Define the reference for the grey rectangle (it should minimize the number of outliers)
  		bgMean<-mean(quality@average.background)
  		bgMin <- min(quality@average.background)
  		bgMax <- max(quality@average.background)
  		ref<-cbind(c(bgMean-10,bgMin,bgMax-20),c(bgMean+10,bgMin+20,bgMax))
  		t1 <- (quality@average.background >= (bgMean-10) & quality@average.background <= (bgMean+10))
  		t2 <- (quality@average.background >= bgMin & quality@average.background <= (bgMin+20))
  		t3 <- (quality@average.background >= (bgMax-20) & quality@average.background <= bgMax)
  		outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
  		if(outliers[1] == 0){
  			ref<- ref[1,]
  		}else{
  			ref<- ref[outliers==min(outliers),] # less outliers
  			if(length(ref) > 3) { ref<-ref[1,] }
  		}
  		testMin <- ref[1]
  		testMax <- ref[2]

  		t1 <- (quality@average.background >= testMin & quality@average.background <= testMax)

        Col <- t1
        Col[t1==TRUE]<-"blue"
        Col[t1==FALSE]<-"red"

        if(exists("rows")){
          rows<-c(rows,"Background\nspread<=20%")
        } else{
          rows<-c("Background\nspread<=20%")
        }
        if(exists("indicator")){
          indicator<-rbind(indicator,round(quality@average.background,0))
        } else{
          indicator<-round(quality@average.background,0)
        }
        if(exists("indicCol")){
          indicCol<-rbind(indicCol,Col)
        } else{
          indicCol<-rbind(Col)
        }
      }

  	if(scaleFact){
  		  print ("indicators for scale factors")
  		# Define the reference for the grey rectangle (it should minimize the number of outliers)
  		sfMin <- min(log2(quality@scale.factors))
  		sfMax <- max(log2(quality@scale.factors))
  		sfMean<- mean(log2(quality@scale.factors))
  		ref<-cbind(c(sfMean-1.5,sfMin,sfMax-3),c(sfMean+1.5,sfMin+3,sfMax))
  		t1 <- (log2(quality@scale.factors) >= (sfMean-1.5) & log2(quality@scale.factors) <= (sfMean+1.5))
  		t2 <- (log2(quality@scale.factors) >= sfMin & log2(quality@scale.factors) <= (sfMin+3))
  		t3 <- (log2(quality@scale.factors) >= (sfMax-3) & log2(quality@scale.factors) <= sfMax)
  		outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
  		if(outliers[1] == 0){
  			ref<- ref[1,]
  		}else{
  			ref<- ref[outliers==min(outliers),] # less outliers
  			if(length(ref) > 3) { ref<-ref[1,] }
  		}
  		testMin <- ref[1]
  		testMax <- ref[2]

  		t1 <- (log2(quality@scale.factors) >= testMin & log2(quality@scale.factors) <= testMax)

  		Col <- t1
  		Col[t1==TRUE]<-"blue"
  		Col[t1==FALSE]<-"red"

      #  if(sfMax - sfMin <= 3){
      #    Col <- rep("blue", length(sampleNames(Data)))
      #  }else{
      #    Col <- rep("red", length(sampleNames(Data)))
      #  }
        if(exists("rows")){
          rows<-c(rows,"Log Scale Factor\nspread<=3")
        } else{
          rows<-c("Log Scale Factor\nspread<=3")
        }
        if(exists("indicator")){
          indicator<-rbind(indicator,round(log2(quality@scale.factors),2))
        } else{
          indicator<-round(log2(quality@scale.factors),2)
        }
        if(exists("indicCol")){
          indicCol<-rbind(indicCol,Col)
        } else{
          indicCol<-rbind(Col)
        }
      }
    }

    #par(parStart)
    if(!is.null(rows)) {#and as such indicator and indicCol
      png(file="QCtable.png",width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
      par(oma=c(2,4.5,0,0))
      plot(c(0,length(rows)), type = 'n', ann = FALSE, axes = FALSE,
        frame.plot = FALSE, ylim = c(0, length(sampleNames(Data))),
        xlim = c(0,length(rows)))
      grid(ny=length(sampleNames(Data))+1, nx = length(rows)+1, equilogs = FALSE)
      for (i in 1:length(rows)){
        for (j in 1:length(sampleNames(Data))) {
          if(length(rows)>1 && length(sampleNames(Data))>0) {
            text(i,length(sampleNames(Data))+1-j,as.character(indicator[i,j]),col=indicCol[i,j],cex=.6)
          } else {
            text(i,length(sampleNames(Data))+1-j,as.character(indicator[ifelse(length(rows)>1,j,i)]),col=indicCol[ifelse(length(rows)>1,j,i)],cex=.6)
          }
        }
      }
      title(main="Summary of raw data quality indicators", cex=0.8)
      mtext("blue = \"within\" / red = \"out of\" recommended cut-off", side=3,
          cex=0.6)
      axis(1, at=1:length(rows),las=2,labels = rows, cex.axis=0.5)
      axis(2,at=1:length(sampleNames(Data)),las=2,labels=rev(substr(sampleNames(Data),1,25)),
          cex.axis=.5)
      dev.off()
    }
  }


  ####################
  ## samplePrepPlot ##
  ####################

  samplePrepPlot <- function(Data,sprep=NULL,lys=NULL,plotColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")

  	if(is.null(lys)) {
  	  try({calls<-detection.p.val(Data)$call
  	  lys<-calls[rownames(calls)[grep("lys.*3",rownames(calls), ignore.case = TRUE)],]
  	  rm(calls)},TRUE)
  	  if(is.null(lys)) {
  		warning("Plots based on the detection.p.val function cannot be created for this chip type")
  	  }else{
  		if(length(lys) > length(sampleNames(Data))) { lys<-lys[1,] }
  	  }
  	}

    if(!is.null(lys) && is.null(sprep)) {
      try(yack <- yaqc(Data),TRUE)
      if(exists("yack")) {
        spnames<-rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3' !
        rownames(yack@morespikes), ignore.case = TRUE),])
        sprep<-t(yack@morespikes[spnames,])
      }else {
        warning("Values based on the yaqc function cannot be computed for this chip type")
      }
    }

    if(!is.null(sprep) && !is.null(lys)) {
      #par(parStart)
      lmin<-min(sprep)
      lmax<-max(sprep)+50

      # create the plot

      png(file = "RawDataSamplePrepControl.png",width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
     if(length(sampleNames(Data)) < MAXARRAY){
      	par(mfrow=c(1,1),oma=c(13,0,0,0))
      }else{
      	par(mfrow=c(1,1),oma=c(0.1,0,0,0))
      }
      plot(c(lmin,lmax), type = 'n', ann = FALSE, axes = FALSE,
             frame.plot = TRUE, ylim = c(lmin, lmax), xlim = c(1,5))

      dummy<-apply(cbind(sprep,plotColors,1:length(sampleNames(Data))),1,
          function(x) {
          par(new=TRUE);
          plot(as.numeric(x[1:(length(x)-2)]),
             type = 'b', ann=FALSE, axes=FALSE, pch='.', cex=4, lwd=2,
             col=x[length(x)-1], lty = as.numeric(x[length(x)]),
             ylim=c(lmin,lmax), xlim=c(1,5))})

      title(main="Spike-in Sample Prep controls intensities and calls",
          ylab="Intensity",cex.lab=0.8)

      label3<-""
      title1<-"Intensity: OK (Lys < Phe < Thr < Dap for all arrays)"
      title2<-"\nLys Present calls: OK (all Lys are called present)"

      # test for Lys < Phe < Thr < Dap
      bad1 <-  colnames(t(sprep))[!(sprep[,1]<sprep[,2]
            & sprep[,2]<sprep[,3]
            & sprep[,3]<sprep[,4])]
      if(length(bad1)>0){
        title1<-paste("Intensity: not OK (some array",ifelse(length(bad1)>1,"s",""),
  		" do",ifelse(length(bad1)>1,"","es")," not follow the rule Lys < Phe < Thr < Dap;",
  		"\nmaybe no Sample Prep controls were spiked on these arrays.)",sep="")
      }

      # test for Lys presence
  	bad2 <- names(lys[lys != 'P'])
      bad2 <- gsub(".present","",bad2)
      if(length(bad2)>0){
        title2<-paste("Lys Present calls: ",length(bad2),
            " Lys not called present",sep="")
        }
       if(length(lys[lys=='P'])==0){
         label3 <-
            "\nApparently no Sample Prep controls were spiked on these arrays."
       }

      title(xlab = c(title1,title2), cex.lab = 0.8)
      legend("topleft", paste("Lys present calls = ",
           length(lys[lys == "P"]),"/", length(lys),label3), bty = "n", cex=0.7)
      if(length(sampleNames(Data)) < (MAXARRAY+20)){
      	cexval <- 0.65
      }else{
      	cexval <- 0.45
      }
      legend("topright", substr(sampleNames(Data),1,20), lwd=2,
  		   col = plotColors, cex =cexval, bty = "n", lt = 1:length(sampleNames(Data)))
      par(cex.axis = 0.8)
  	axis(2)
      axis(1, at=1:4, labels = c("Lys", "Phe", "Thr", "Dap"))

      if(length(bad2)>0){
      text(c(rep(0.87,length(lys[lys!='P']))),sprep[lys!='P',1],
           lys[lys!='P'], pos=4,offset=0.2, cex=0.7, col="red")
      }
      dev.off()
    } else {
      warning("Spike-in sample prep plot is cannot be computed for this chip type")
    }
  }


  ###############
  ## ratioPlot ##
  ###############

  ratioPlot <- function(Data,quality=NULL,experimentFactor=NULL,plotColors=NULL,legendColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
    if(is.null(quality)) try(quality <- qc(Data),TRUE)

    if(!is.null(quality)) {
      #par(parStart)

      plotFun <- function(i,j,cutoff,Cname){
        ratio35 <- quality@qc.probes[,i]/quality@qc.probes[,j]
        ratioM <- quality@qc.probes[,i]/quality@qc.probes[,i+1]

        cMin <- 0
        cMax <- max(max(cbind(ratio35, ratioM)) + 2, 5)
        if(Cname=="beta-actin"){
          symbol=c(17,2) # triangles for beta-actin
        }else{
          symbol=c(19,1) # circles for GAPDH (as in the simpleaffy QC report)
        }

        par(mfrow=c(1,2),oma=c(17,0.5,0,0.5),cex.axis=0.8)
        plot(ratio35, type='n', ann=FALSE, axes=FALSE,
             frame.plot=TRUE, pch='.', cex=10,  ylim=c(cMin,cMax))
        rect(0, 0, length(ratio35)+1, cutoff, col=gray(0.9), border=FALSE)

        for (k in 0:cMax){
        abline(h=k,lty=2,col=gray(0.8))
        }

        par(new=TRUE)
        plot(ratio35, type='h', ann=FALSE, axes=FALSE,
               frame.plot=TRUE, lwd=3, col=plotColors, ylim=c(cMin,cMax))
        par(new=TRUE)
        plot(ratio35,type='p', ann=FALSE, axes=FALSE,
               frame.plot=TRUE, pch=symbol[1], col=plotColors, ylim=c(cMin,cMax))
        par(new=TRUE)
        plot(ratioM,type='p', ann=FALSE, axes=FALSE,
               frame.plot=TRUE, pch=symbol[2], col="black", ylim=c(cMin,cMax))

        title(main= paste("RNA degradation of", Cname),ylab="3'/5' and 3'/M ratios")
        axis(2)
        par(cex.axis=0.65)
  	  if(length(sampleNames(Data))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
  	  	  axis(1,at=1:length(ratio35),las=2,labels=sampleNames(Data))
        }
  	  if(length(levels(experimentFactor))>1){
          legend("topright", levels(experimentFactor),
               col = legendColors, fill = legendColors, bty = "n",cex = 0.55)
        }
  	  legend("topleft", c(paste("3\'/5\' ratio (max=", round(max(ratio35),2),")",
  		sep=""), paste("3'/M ratio (max=", round(max(ratioM),2),")",sep="")),
  		col = c(gray(0.3),"black"),pch=symbol, cex=0.55, bty = "o")
  	  t1 <- ratio35 <= cutoff
  	  if(length(t1[t1==FALSE])>0 && length(sampleNames(Data))>=(MAXARRAY/2)){
  		  if(length(t1[t1==FALSE])<25){
  			textlab<-sampleNames(Data)[t1==FALSE]
  			legend("topleft", c(rep("",4),"Outliers, from left to right",textlab),
  			  text.col = c(rep(1,5), rep("red",length(textlab))), bty = "n",cex=0.38)
  		  }else{
  			mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
  		  }
  	  }
        mtext(paste("Ratios should stand within the grey rectangle [0,",cutoff,"]"),
               side=4, font=1, cex=0.7)

        par(cex.axis=0.8)
        boxplot(cbind(ratio35, ratioM), axes=FALSE, frame.plot=TRUE)
        title(main=paste("Boxplot of", Cname,"ratios"))
        par(cex.lab=0.8)
        if(max(ratio35 ) < cutoff){
           title(xlab=paste(Cname," QC: OK (all 3'/5' ratios < ",cutoff,")",sep=""))
        }  else{
           title(xlab=paste(Cname," QC: not OK (some 3'/5' ratios >",cutoff,
  		 ")\nNote that the threshold of ",cutoff,
  		" \nwas determined for Homo Sapiens.",sep=""), line=4)
        }

        axis(1,at=1:2,labels=c("3'/5' ratio","3'/M ratio"))
        axis(2)
      }

      # Graph creation: two separated graphs
  	png(file = "RawData53ratioPlot_beta-actin.png",width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
      plotFun(1,3,3,"beta-actin")
      dev.off()
      png(file = "RawData53ratioPlot_GAPDH.png",width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
      plotFun(4,6,1.25,"GAPDH")
      dev.off()

    } else {
      warning("3'/5' ratio plot is not computed for this chip type")
    }
  }

  ################
  ## RNAdegPlot ##
  ################

  RNAdegPlot <- function(Data, Data.rnadeg=NULL, plotColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(Data.rnadeg)) Data.rnadeg <- AffyRNAdeg(Data)

    png(file = "RawDataRNAdegradation.png", width=WIDTH, height=HEIGHT, pointsize=POINTSIZE)

    if(length(sampleNames(Data))<MAXARRAY){
  	cexval <- 0.65
  	par(lwd = 2, oma=c(13,0,0,0))
    } else{
  	cexval <- 0.45
  	par(lwd = 2, oma=c(0.1,0,0,0))
    }
    par(mar=c(4,4,4,0), cex.axis=0.6, cex.lab=0.75)
    layout(matrix(c(1,2),1,2,byrow=TRUE), c(2,1), 1, FALSE)
    plotAffyRNAdeg(Data.rnadeg, col = plotColors, lty = 1:length(sampleNames(Data)))
    par(mar=c(4,0,4,0), cex.axis=0.01, cex.lab=0.01)
    plot(1,type = 'n', ann=FALSE, axes=FALSE, frame.plot=FALSE)
    legend("topleft",sampleNames(Data), lwd=2, col = plotColors, cex = cexval, bty = "n", lty = 1)
    dev.off()
  }


  ################
  ## hybridPlot ##
  ################

  hybridPlot <- function(Data,quality=NULL,plotColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(quality)) try(quality <- qc(Data),TRUE)

    if(!is.null(quality)) {

      png(file = "RawDataSpikeinHybridControl.png",width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
     if(length(sampleNames(Data)) < MAXARRAY){
      	par(mfrow=c(1,1),oma=c(13,0,0,0))
      }else{
      	par(mfrow=c(1,1),oma=c(0.1,0,0,0))
      }
      plot(c(min(quality@spikes),max(quality@spikes)),
             type = 'n', ann = FALSE, axes = FALSE, frame.plot = TRUE,
             ylim = c(min(quality@spikes), max(quality@spikes)),
             xlim = c(1,5))

      dummy<-apply(cbind(quality@spikes,plotColors,1:length(sampleNames(Data))),1,
             function(x) {
          par(new=TRUE);
          plot(as.numeric(x[1:(length(x)-2)]),
             type = 'b', ann=FALSE, axes=FALSE, pch='.', cex=4, lwd=2,
             col=x[length(x)-1], lty = as.numeric(x[length(x)]),
             ylim=c(min(quality@spikes),max(quality@spikes)), xlim = c(1,5))})

      title(main = "Spike-in Hybridization controls intensities and calls",
             ylab="Intensity",cex.lab=0.8)
      label3<-""
      title1<-"Intensities: OK (bioB < bioC < bioD < creX for all arrays)"
      title2<-"BioB Present calls: OK (indeed all bioB are called present)"

      # test for bioB<bioC<bioD<creX
      bad1 <- colnames(t(quality@spikes))[!(quality@spikes[,1]<quality@spikes[,2]
                   & quality@spikes[,2]<quality@spikes[,3]
                   & quality@spikes[,3]<quality@spikes[,4])]
      if(length(bad1)>0){
       title1<-paste("Intensity: not OK (some array",ifelse(length(bad1)>1,"s",""),
  		" do",ifelse(length(bad1)>1,"","es")," not follow the rule bioB < bioC < bioD < creX)",sep="")
      }

      # test for bioB presence
      bad2<-gsub(".present","",names(quality@bioBCalls[quality@bioBCalls != 'P']))
      if(length(bad2)>0){
        title2<-paste("BioB present calls: not OK (",length(bad2),"/",
                length(quality@bioBCalls)," bioB not called present)")
        }
       if(length(quality@bioBCalls[quality@bioBCalls =='P'])==0){
         label3 <-
            "\nApparently no Hybridization controls were spiked on these arrays."
        }

        title(xlab = c(title1,title2), cex.lab = 0.8)
        legend("topleft", paste("bioB present calls = ",
                length(quality@bioBCalls[quality@bioBCalls == "P"]),"/",
                length(quality@bioBCalls),label3), bty = "n", cex=0.7)

  	    if(length(sampleNames(Data)) < (MAXARRAY+20)){
  	    	cexval <- 0.65
  	    }else{
  	    	cexval <- 0.45
  	    }

        legend("bottomright", substr(sampleNames(Data),1,20), lwd=2,
        	   col = plotColors, cex = cexval, bty = "n", lt = 1:length(sampleNames(Data)))
  	  par(cex.axis = 0.8)
        axis(2)
        axis(1, at=1:4, labels = names(colnames(quality@spikes)))

        if(length(bad2)>0){
          text(c(rep(0.89,length(quality@bioBCalls[quality@bioBCalls!='P']))),
                quality@spikes[quality@bioBCalls!='P',1],
                quality@bioBCalls[quality@bioBCalls!='P'], pos=4,offset=0.2,
                cex=0.8, col="red")
        }
      dev.off()
    } else {
      warning("Spike-in hybridization plot is not computed for this chip type")
    }
  }


  ##################
  ## percPresPlot ##
  ##################

  percPresPlot <- function(Data,quality=NULL,experimentFactor=NULL,plotColors=NULL,legendColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
    if(is.null(quality)) try(quality <- qc(Data),TRUE)

    if(!is.null(quality)) {
   	# Define the reference for the grey rectangle (it should minimize the number of outliers)
      ppMin <- min(quality@percent.present)
      ppMax <- max(quality@percent.present)
      ppMean<- mean(quality@percent.present)

  	ref<-cbind(c(ppMean-5,ppMin,ppMax-10),c(ppMean+5,ppMin+10,ppMax))
  	reftext<-c("[mean-5% ; mean+5%]", "[min ; min+10%]", "[max-10% ; max]")
  	t1 <- (quality@percent.present >= (ppMean-5) & quality@percent.present <= (ppMean+5))
  	t2 <- (quality@percent.present >= ppMin & quality@percent.present <= (ppMin+10))
  	t3 <- (quality@percent.present >= (ppMax-10) & quality@percent.present <= ppMax)
  	outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
  	if(outliers[1] == 0){
  		ref<- ref[1,]
  		reftext<- reftext[1]
  	}else{
  		ref<- ref[outliers==min(outliers),] # less outliers
  		if(length(ref) > 3) { ref<-ref[1,] }
  		reftext<- reftext[outliers==min(outliers)] # less outliers
  		if(length(reftext) > 1) { reftext<-reftext[1] }
  	}
  	testMin <- ref[1]
  	testMax <- ref[2]

     png(file = "RawDataPercentPresent.png",width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
      par(mfrow=c(1,2),oma=c(17,0.5,0,0.5),cex.axis=0.8)
      plot(quality@percent.present, type='n', ann=FALSE, axes=FALSE,
           frame.plot=TRUE, pch='.', cex=10,  ylim=c(0,100))
      rect(0, testMin, length(quality@scale.factors)+1,
           testMax, col=gray(0.9), border=FALSE)

      for(i in seq(from=0, to=100, by=20)){
          abline(h=i,lty=2,col=gray(0.8))
      }
      abline(h=50,lty=1,col=gray(0.5))

      par(new=TRUE)
      plot(quality@percent.present, type='h', ann=FALSE, axes=FALSE,
             frame.plot=TRUE, pch='.',lwd=3, col=plotColors, ylim=c(0,100))
      par(new=TRUE)
      plot(quality@percent.present,type='p', ann=FALSE, axes=FALSE,
             frame.plot=TRUE, pch='.',cex=10, col=plotColors, ylim=c(0,100))

  	t1 <- (quality@percent.present >= testMin & quality@percent.present <= testMax)
  	if(length(t1[t1==FALSE])>0 && length(sampleNames(Data))>=(MAXARRAY/2)){
  		  if(length(t1[t1==FALSE])<25){
  			textlab<-sampleNames(Data)[t1==FALSE]
  			legend("topleft", c("Outliers, from left to right",textlab),
  			  text.col = c(1, rep("red",length(textlab))), bty = "n",cex=0.38)
  		  }else{
  			mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
  		  }
  	}

  	title(main="Plot of percent present",xlab="",ylab="percentage")
      axis(2)
      par(cex.axis=0.65)
  	if(length(sampleNames(Data))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
  	  	axis(1,at=1:length(quality@percent.present),las=2,labels=sampleNames(Data))
      }

      if(length(levels(experimentFactor))>1){
        legend("topright", levels(experimentFactor),
      	   col = legendColors, fill = legendColors, cex = 0.55, bty = "n")
      }

      legend("bottomright", c(paste("min = ", round(ppMin,2), sep="" ),
          paste("max = ", round(ppMax,2), sep="" ),
          paste("max-min = ", round(ppMax-ppMin,2), sep="")), bty = "o",cex = 0.55)

      mtext(paste(
  		 "Data should be in the grey rectangle representing a spread of 10%"
           #"Data should stand within the grey rectangle",reftext
  		 ),
           side=4, font=1, cex=0.7)
      par(cex.axis=0.8, cex.lab=0.8)
      boxplot(quality@percent.present)
      title(main="Boxplot of percent present")
      if(ppMax-ppMin <= 10){
         title(xlab="Percent present QC: OK (spread <= 10%)")
      } else{
         title(xlab="Percent present QC: not OK (spread > 10%)")
      }

      dev.off()
    } else {
      warning("Percent present plot is not computed for this chip type")
    }
  }


  #################
  ## PNdistrPlot ##
  #################

  PNdistrPlot <- function(Data, WIDTH=1000, HEIGHT=1414, POINTSIZE=24) {
    png(file = "RawDataPosNegDistribution.png",width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(oma=c(17,0,0,0),srt=90)
    borderQC1(Data)
    mtext(paste("All distributions should be similar and","extreme values should not be reached\n"), side=3, cex=0.7)
    dev.off()
  }


  ############
  ## bgPlot ##
  ############

  backgroundPlot <- function(Data,quality=NULL,experimentFactor=NULL,plotColors=NULL,legendColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
    if(is.null(quality)) try(quality <- qc(Data),TRUE)

    if(!is.null(quality)) {

  	# Define the reference for the grey rectangle (it should minimize the number of outliers)
  	bgMean<-mean(quality@average.background)
  	bgMin <- min(quality@average.background)
  	bgMax <- max(quality@average.background)
  	ref<-cbind(c(bgMean-10,bgMin,bgMax-20),c(bgMean+10,bgMin+20,bgMax))
  	reftext<-c("[mean-10 ; mean+10]", "[min ; min+20]", "[max-20 ; max]")
  	t1 <- (quality@average.background >= (bgMean-10) & quality@average.background <= (bgMean+10))
  	t2 <- (quality@average.background >= bgMin & quality@average.background <= (bgMin+20))
  	t3 <- (quality@average.background >= (bgMax-20) & quality@average.background <= bgMax)
  	outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
  	if(outliers[1] == 0){
  		ref<- ref[1,]
  		reftext<-reftext[1]
  	}else{
  		ref<- ref[outliers==min(outliers),] # less outliers
  		if(length(ref) > 3) { ref<-ref[1,] }
  		reftext<- reftext[outliers==min(outliers)] # less outliers
  		if(length(reftext) > 1) { reftext<-reftext[1] }
  	}
  	testMin <- ref[1]
  	testMax <- ref[2]

      ylimit = c(0, max(max(quality@maximum.background),
            testMax)+5)

      png(file = "RawDataBackground.png",width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
      par(mfrow=c(1,2),oma=c(17,0.5,0,0.5),cex.axis=0.8)
      plot(c(quality@minimum.background, quality@maximum.background),
             type = 'n', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', cex=10,
             xlim = c(1,length(quality@maximum.background+1)),
             ylim = ylimit)
      rect(0,testMin,
             length(quality@average.background)+1,
             testMax,
             col=gray(0.8), border=TRUE, density=10, lty=6)
      par(new=TRUE)
      plot(quality@maximum.background,
             type = 'h', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', lwd=3,
             col=plotColors, ylim = ylimit)
      par(new=TRUE)
      plot(quality@minimum.background,
             type = 'h', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch='.', lwd=4,
             col='white', ylim = ylimit)
      par(new=TRUE)
      plot(quality@maximum.background,
             type = 'p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=6, cex=0.5,
             ylim = ylimit)
      par(new=TRUE)
      plot(quality@minimum.background,
             type = 'p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=2, cex=0.5,
             ylim=ylimit)
      par(new=TRUE)
      plot(quality@average.background,
             type='p', ann=FALSE, axes=FALSE, frame.plot=TRUE, pch=3, cex=1,
             ylim=ylimit)
      abline(h=testMin, col=gray(0.8), ylim=ylimit)
      abline(h=testMax, col=gray(0.8), ylim=ylimit)
      title(main="Plot of background intensity",xlab="",ylab="background intensity")
      axis(2)
      par(cex.axis=0.65)
  	if(length(sampleNames(Data))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
  	  	 axis(1,at=1:length(quality@average.background),las=2,labels=sampleNames(Data))
      }
      if(length(levels(experimentFactor))>1){
      legend("bottomright", c(levels(experimentFactor),"max bg","average bg",
      	   "min bg"), col = c(legendColors, 1, 1, 1),
      	   pch = c(rep(15, length(legendColors)),6,3,2), bty = "n", cex=0.55)
      }else{
      legend("bottomright", c("max bg","average bg",
      	   "min bg"), col = c(1, 1, 1), pch = c(6,3,2), bty = "n", cex=0.55)
      }
  	t1 <- (quality@average.background >= testMin & quality@average.background <= testMax)
  	if(length(t1[t1==FALSE])>0 && length(sampleNames(Data))>=(MAXARRAY/2)){
  		  if(length(t1[t1==FALSE])<25){
  			textlab<-sampleNames(Data)[t1==FALSE]
  			legend("bottomleft", c("Outliers, from left to right",textlab),
  				  text.col = c(1, rep("red",length(textlab))), bty = "n",cex=0.38)
  		  }else{
  			mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
  		  }
  	}

      mtext(paste(
  	    "Data should be in the grey rectangle representing a spread of 20"
          #"Data should stand within the dashed grey rectangle",reftext
          ), side=4, font=1, cex=0.7)
      par(cex.axis=0.8, cex.lab=0.8)
      boxplot(quality@average.background)
      title(main="Average background intensity")

      if(bgMax - bgMin <= 20){
         title(xlab="Background QC: OK (spread <= 20)")
      } else{
         title(xlab="Background QC: not OK (spread > 20)")
      }

      legend("bottomleft", c(paste("min = ", round(bgMin,2), sep="" ),
          paste("max = ", round(bgMax,2), sep="" ),
          paste("max-min = ", round(bgMax-bgMin,2), sep="")), cex=0.7)

      dev.off()
    } else {
      warning("Background intensity plot is not computed for this chip type")
    }
  }


  ###################
  ## scaleFactPlot ##
  ###################

  scaleFactPlot <- function(Data,quality=NULL,experimentFactor=NULL,plotColors=NULL,legendColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
    if(is.null(quality)) try(quality <- qc(Data),TRUE)

    if(!is.null(quality)) {

  	# Define the reference for the grey rectangle (it should minimize the number of outliers)
      sfMin <- min(log2(quality@scale.factors))
      sfMax <- max(log2(quality@scale.factors))
  	sfMean<- mean(log2(quality@scale.factors))
  	ref<-cbind(c(sfMean-1.5,sfMin,sfMax-3),c(sfMean+1.5,sfMin+3,sfMax))
  	t1 <- (log2(quality@scale.factors) >= (sfMean-1.5) & log2(quality@scale.factors) <= (sfMean+1.5))
  	t2 <- (log2(quality@scale.factors) >= sfMin & log2(quality@scale.factors) <= (sfMin+3))
  	t3 <- (log2(quality@scale.factors) >= (sfMax-3) & log2(quality@scale.factors) <= sfMax)
  	outliers <- c(length(t1[t1==FALSE]),length(t2[t2==FALSE]), length(t3[t3==FALSE]))
  	if(outliers[1] == 0){
  		ref<- ref[1,]
  	}else{
  		ref<- ref[outliers==min(outliers),] # less outliers
  		if(length(ref) > 3) { ref<-ref[1,] }
  	}
  	testMin <- ref[1]
  	testMax <- ref[2]

      ymin <- min(-3.5, sfMin-1.5, testMin)
      ymax <- max(+3.5, sfMax+1.5, testMax)
      ylimit <-c(ymin, ymax)

      png(file = "RawDataScaleFactors.png",width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)

      par(mfrow=c(1,2), oma=c(17,0.5,0,0.5), cex.axis=0.8)
      plot(log2(quality@scale.factors), type='n', axes=FALSE, frame.plot=TRUE,
             ann=FALSE, pch='.', cex=10, ylim = ylimit)
      rect(0, testMin,
             length(quality@scale.factors)+1,
             testMax, col=gray(0.9), border=FALSE)

      for(i in seq(from=ceiling(ymin), to=floor(ymax), by=1)){
          abline(h=i,lty=2,col=gray(0.8))
      }
      abline(h=0,lty=1,col=gray(0.5))

      par(new=TRUE)
      plot(log2(quality@scale.factors),
             type = 'h', axes=FALSE, frame.plot=TRUE, ann=FALSE,
             pch='.', lwd=3, col=plotColors, ylim = ylimit)
      par(new=TRUE)
      plot(log2(quality@scale.factors),
            type = 'p', axes=FALSE, frame.plot=TRUE, ann=FALSE,pch='.',
            cex=10, col=plotColors, ylim = ylimit)
  	t1 <- (log2(quality@scale.factors) >= testMin & log2(quality@scale.factors) <= testMax)
  	if(length(t1[t1==FALSE])>0 && length(sampleNames(Data))>=(MAXARRAY/2)){
  		  if(length(t1[t1==FALSE])<16){
  			textlab<-sampleNames(Data)[t1==FALSE]
  			legend("bottomleft", c("Outliers, from left to right",textlab),
  			  text.col = c(1, rep("red",length(textlab))), bty = "n",cex=0.38)
  		  }else{
  			mtext("Too many arrays and outliers; \nsee details on the QC table",side=1,cex=0.8,line=3)
  		  }
  	}

      title(main="Plot of Log scale factors", xlab="",ylab="Log2(scale factors)")
      axis(2)
      par(cex.axis=0.65)
  	if(length(sampleNames(Data))<(MAXARRAY/2)){ # array names not reported if more than 20 arrays
  	  	 axis(1,at=1:length(quality@scale.factors),las=2,labels=sampleNames(Data))
      }
      if(length(levels(experimentFactor))>1){
        legend("topright", levels(experimentFactor),
      	   col = legendColors, fill = legendColors, cex = 0.55, bty = "n")
      }
      legend("topleft", c(paste("min = ", round(sfMin,2), sep="" ),
          paste("max = ", round(sfMax,2), sep="" ),
          paste("max - min = ",round((sfMax-sfMin),2),sep="")), bty = "n",
              cex = 0.55)

      mtext(paste(
         "Data should be in the grey rectangle representing 3-fold on a log scale"
         ), side=4,  font=1, cex=0.7)
       par(cex.axis=0.8)
      boxplot(quality@scale.factors)
      title(main="Boxplot of scale factors")
      mtext("(natural scale)\n", side=3, font=1, cex=0.7)
      par(cex=0.8)
      if((sfMax - sfMin) < 3){

         title(xlab="Scale factors QC: OK (spread < 3-fold)")
      } else{
         title(xlab="Scale factors QC: not OK (spread > 3-fold)")
      }
      dev.off()
    } else {
      warning("Scale factor plot is not computed for this chip type")
    }
  }


  ##################
  ## controlPlots ##
  ##################

  controlPlots <- function(Data,plotColors=NULL,experimentFactor=NULL,legendColors=NULL,
  	affxplots=TRUE,boxplots=TRUE, WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

  	if(is.null(plotColors)) stop("The 'plotColors' parameter must be provided")
  	if(is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
      if(is.null(legendColors)) stop("the 'legendColors' parameter is required")

  	affx1 <- unique(probeNames(Data)[grep("AFFX",probeNames(Data),fixed=TRUE)])

  	require("ArrayTools", quietly = TRUE)

  	dataTable <- paste(substr(Data@annotation,1,nchar(Data@annotation)-2),"CONTROL",sep="")
  	suppressWarnings(eval(parse("",-1,paste("data(",dataTable,")",sep="")))) #ArrayTools
  	cntrl <- NULL
  	try(cntrl <- get(dataTable),TRUE)

  	cntrlVars <- NULL
  	if(!is.null(cntrl)) {
  		for(name in names(table(cntrl[,2])[table(cntrl[,2])>0])) {
  		  #often names have a format of xxx->yyy, we only want the last part
  		  name2 <- strsplit(name,"->")[[1]][length(strsplit(name,"->")[[1]])]
  		  #make sure that there are only letters in name2, to be used as variable name
  		  name2 <- gsub("[[:punct:]|[:digit:]]","",name2)
  		  assign(name2,cntrl[grep(name,cntrl[,2]),1])
  		  assign(name2,get(name2)[get(name2) %in% geneNames(Data)])
  		  if(length(get(name2))>0) cntrlVars <- c(cntrlVars,name2)
  		}
  	}

    	if(affxplots) { #plot affx controls

  	  if(length(affx1)>0 | length(grep("affx",cntrlVars))>0) {

  		if(length(grep("affx",cntrlVars))>0) {
  			#affx2 cannot exists as automatically generated variable names are restricted to letters
  			affx2 <- NULL
  			for(a in grep("affx",cntrlVars)) {
  				affx2 <- c(affx2, get(cntrlVars[a]))
  			}
  			if(length(affx1)>0) { #both exist
  				if(sort(affx2)!=sort(affx1)) {
  					#merge both
  					affx2 <- unique(c(affx1,affx2))
  				}
  			}
  		} else { #only length(affx1)>0
  			affx2 <- affx1
  		}
   #		dimnCoef <- 1 # squared layout
  		dimnCoef <- 2
  		dimn<-ceiling(sqrt(length(affx2)/dimnCoef))
  		number_empty <- dimnCoef*dimn*dimn - length(affx2)

  		png(file="RawDataAFFXControlsProfiles.png", width=250*dimn+250, height=dimnCoef*250*dimn, pointsize=POINTSIZE)
  		if(number_empty >= dimn){
  			lineout<-dimn
  			number_empty <- number_empty - dimn
  	    }else{
  			lineout<-0
  		}
  		legcol<-dimnCoef*dimn*dimn-lineout+1
  		layout(cbind(matrix(1:(dimnCoef*dimn*dimn-lineout),ncol=dimn,byrow=TRUE),legcol))
  		par(oma=c(5,3,4,3))

  		dummy<-sapply(
  		  affx2,
  		  function(y) {
  			pm<-log(pm(Data,as.character(y)),2)
  			max_pm<-max(pm)
  			min_pm<-min(pm)
  			plot(0,type='n',ylim=c(min_pm,max_pm),ylab="expression",xlim=c(1,dim(pm)[1]),
  			main=y,xlab=paste("\naverage: ",round(mean(pm))),cex.lab=1.3,cex.main=1.1)
  			assign("col", 0, 1)
  			apply(
  			  pm,
  			  2,
  			  function(z) {
  				assign("col",get("col",1)+1,1)
  				points(z,type='o',col=plotColors[get("col",1)])
  			  }
  			)
  			points(apply(pm,1,mean),type='l',lwd=2,col="black")
  		  }
  		)


  		for(p in 1:number_empty) plot(0,type='n',xaxt='n',yaxt='n',bty='n',xlab="",ylab="",main="")
  		par(mar=c(0,0,0,0))
  		plot(0,type='n',xaxt='n',yaxt='n',bty='n',xlab="",ylab="",main="")
  		legend("topleft",c(sampleNames(Data),"","average"),lwd=c(rep(1,length(sampleNames(Data))),0,2),
  		pch=c(rep(1,length(sampleNames(Data))),-1,-1),col=c(plotColors,"white","black"),cex=1.5,bty="n")

  		mtext("affx control profiles", side = 3, outer = TRUE, font = 2, cex = 2)

  		dev.off()
        } else {
  	    warning("plot of AFFX controls cannot be made, AFFX controls cannot be determined")
  	  }
      }
  	###################
      if(boxplots) {
  		if(length(affx1)>0) {
  			#affx1 will not already be in controlVars as it contains a number
  			cntrlVars <- c("affx1", cntrlVars)
  		}
  		n <- 0
  		for(var in cntrlVars) {
  		    varName <- paste(var,"Data",sep="")
  			assign(varName,NULL)
  			try(assign(varName,log(pm(Data,as.character(get(var))),2)),TRUE)
  			if(!is.null(get(varName))) n<-n+1
  		}
  		#make boxplots of controls
  		if(n > 0) {
  			for(i in 1:n) {
  				if(!is.null(get(varName))) {
  					png(paste("RawData",toupper(cntrlVars[i]),"ControlsBoxplot.png",sep=""), width=WIDTH, height=HEIGHT,pointsize=POINTSIZE)
  					par(oma=c(17,0,0,0))
  					boxplot(get(paste(cntrlVars,"Data",sep="")[i]), main=paste(gsub("1","",
  					cntrlVars[i],fixed=TRUE),"controls"), ylim=c(0,16), col=plotColors, cex=1, axes=FALSE)
  					if(length(sampleNames(Data))<MAXARRAY){
  						cexval <- 0.65
  					}else{
  						cexval <- 0.45
  					}
  				    if(length(levels(experimentFactor))>1){
  						legend("topright", levels(experimentFactor), col = legendColors,
  					         fill = legendColors, cex = 0.7, bg = "white", bty = "o")
  					}
  					axis(1,at=1:length(sampleNames(Data)),las=2,labels=sampleNames(Data), cex.axis=cexval)
  					axis(2, cex.axis=0.7)
  					dev.off()
  				}
  			}
  		} else {
  			warning("no control boxplots can be made, control spots cannot be determined")
  		}
  	}
  }


  ################
  ## boxplotFun ##
  ################

  boxplotFun <- function(Data, experimentFactor=NULL, plotColors=NULL, legendColors=NULL, normMeth="",
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(experimentFactor)) stop("the 'experimentFactor' parameter is required")
    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(legendColors)) stop("the 'legendColors' parameter is required")


    if(class(Data) == "AffyBatch") {
      Type <- "Raw"
      tmain <- "Boxplot of raw intensities"
      tmtext2 <- "Raw log intensity\n\n\n"
    } else {
      if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
      Type <- "Norm"
      tmain <- paste("Boxplot after ", normMeth, sep="")
      tmtext2 <- "Normalized log intensity\n\n\n"
    }

   png(file = paste(Type,"DataBoxplot.png", sep=""),width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
    par(oma=c(17,0,0,0), cex.axis=1)
    suppressWarnings(boxplot(Data, col=plotColors ,main=tmain, axes=FALSE))
    if(length(levels(experimentFactor))>1){
      legend("topright", levels(experimentFactor),
         col=legendColors,fill=legendColors, cex = 0.7, bg = "white", bty = "o")
    }
  	if(length(sampleNames(Data))<MAXARRAY){
  		cexval <- 0.65
  	}else{
  		cexval <- 0.45
  	}
    axis(1,at=1:length(sampleNames(Data)),las=2,labels=sampleNames(Data), cex.axis=cexval)
    axis(2, cex.axis=0.7)
    mtext(tmtext2, side=2, cex=0.8)
    mtext("Distributions should be comparable between arrays\n", side=3, font=1,
        cex=0.7)
    dev.off()
  }


  ################
  ## densityFun ##
  ################

  densityFun <- function(Data, plotColors=NULL, normMeth="",
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41){

    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")

    Type <- ifelse(class(Data) == "AffyBatch","Raw","Norm")

      png(file = paste(Type,"DensityHistogram.png", sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
  	if(length(sampleNames(Data))<MAXARRAY){
  		cexval <- 0.65
  		par(oma=c(12,0,0,0) )
  	}else{
  		cexval <- 0.45
  		par(oma=c(0.1,0,0,0) )
  	}
    if(Type == "Raw"){
      hist(Data, lwd=3, lt = 1:length(sampleNames(Data)), col = plotColors,
        which = "both", main="Density histogram of raw intensities", cex.axis = 0.7, cex.lab=0.8)
    }else{
      if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
      hist(Data, lwd=3, lt = 1:length(sampleNames(Data)), col = plotColors,
        main=paste("Density histogram after ", normMeth,"\n", sep=""), cex.axis = 0.7, cex.lab=0.8)
    }

    legend("topright", substr(sampleNames(Data),1,20), lwd=3, lt = 1:length(sampleNames(Data)),
      col = plotColors, cex = cexval, bty = "n")
    mtext( "Curves should be comparable between arrays\n", side=3, font=1,
      cex=0.7)
    dev.off()
  }


  ##########################
  ## densityFunUnsmoothed ##
  ##########################

  # alternative density function, not called by default by run_affyAnalysisQC.R
  #----------------------------------------------------------------------------

  densityFunUnsmoothed <- function(Data, plotColors=NULL, normMeth="",
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41){

    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")

    #plot(density(exprs(Data)))
    Type <- ifelse(class(Data) == "AffyBatch","Raw","Norm")

    if(Type == "Raw"){
      tmp_data <- log(exprs(Data),2)
      main <- "Density histogram of raw intensities"
    } else { #normalized data
      if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
      tmp_data <- exprs(Data)
      main <- paste("Density histogram after", normMeth)
    }

  	 png(file = paste(Type,"DensityHistogramUnsmoothed.png", sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
  	if(length(sampleNames(Data))<MAXARRAY){
  		cexval <- 0.65
  		par(oma=c(12,0,0,0) )
  	}else{
  		cexval <- 0.45
  		par(oma=c(0.1,0,0,0) )
  	}
    for (i in 1:length(sampleNames(Data))) {
      if (i > 1) par(new=TRUE)
      tmp_hist <- hist(tmp_data[,i],plot=FALSE,breaks=64)
      plot(tmp_hist$mids,tmp_hist$density,type='l',lty=i,col=plotColors[i],xlim=c(min(tmp_data),max(tmp_data)),ylim=c(0,1),
        main = main, xlab='log2 intensity',ylab='density', cex.axis = 0.7, cex.lab=0.8)
      text(tmp_hist$mids[order(tmp_hist$density)[length(tmp_hist$density)]],max(tmp_hist$density),sampleNames(Data)[i],col=i,pos=3,offset=0.2,cex=0.6)
    }

    legend("topright", substr(sampleNames(Data),1,20), lwd=3, lt = 1:length(sampleNames(Data)),
      col = plotColors, cex = cexval, bty = "n")
    mtext( "Curves should be comparable between arrays\n", side=3, font=1,
      cex=0.7)

    dev.off()

  }


  ###########
  ## maFun ##
  ###########

  maFun <- function(Data, experimentFactor=NULL, perGroup=FALSE, normMeth="", aType=NULL,
    WIDTH=1000, HEIGHT=1414, MAXARRAY=41){

    #verify whether raw or normalised data have been provided
    if(class(Data) == "AffyBatch"){
      if(is.null(aType)) stop("When selecting MA plots of raw data, 'aType' must be provided")
      Type <- "Raw"
  	tmain <- paste("MA plots of raw data",ifelse(perGroup,", computed for group",""),sep="")
    } else {
      if(normMeth=="") stop("When selecting MA plots of normalized data, 'normMeth' must be provided")
      Type <- "Norm"
  	tmain <- paste("MA plots after",normMeth,"normalization",ifelse(perGroup,", computed for group ",""),sep="")
    }

    #check whether MA plots have to be created per experimental group or for the whole data set at once
    if(perGroup){
      if(is.null(experimentFactor)) stop("When selecting MA plots per experimental group, 'experimentFactor' must be provided")
    } else {
      experimentFactor <- factor(rep("",length(sampleNames(Data))))
    }

    #plot the images
    for(k in (levels(experimentFactor))) {
      x2 <- Data[,experimentFactor==k]

     #amount of plots per page
  	if(length(sampleNames(x2))<=12) { nPerPage<-6; nCol<- 2; cexmain<-1.9 ; cexval<-1.8} # max 3p of 6 plots
  	if(length(sampleNames(x2))>12) {nPerPage<-15; nCol<- 3; cexmain<-1.3 ; cexval<-1.3}
  	if(length(sampleNames(x2))>45) {nPerPage<-24 ;nCol<- 4; cexmain<-1.4 ; cexval<-1.2}

  	#subselection of data
      if(length(sampleNames(x2)) < 2) {
  	  warning("group ",ifelse(perGroup,k,"")," consists of only one array: no MA plot can be made",
  	    ifelse(perGroup,": consider selecting 'dataset' as option for the MA plots",""))
  	} else {
      nplots <- ceiling(length(sampleNames(x2))/nPerPage)
  	  for(l in 1:nplots){
          png(paste(Type,"DataMAplot", ifelse(k!="","_",""), gsub("[:/|\\?<>*\"]","_",k),
  		  ifelse(nplots>1,"_",""), ifelse(nplots>1,l,""), ".png", sep=""), width=WIDTH, height=HEIGHT)
  		par(mfrow = c((nPerPage/nCol),nCol),oma=c(0,0,6,0),adj=0)
          from <- nPerPage*l - (nPerPage-1)
          to <-   min(nPerPage*l, length(sampleNames(x2)))
          if(Type=="Raw" && aType!="PMMM") {
            MAplot(x2, pairs = FALSE, which= from:to, plot.method = "smoothScatter", lwd=3, type="pm", cex.main=cexmain,cex=cexval)
          } else {
            MAplot(x2, pairs = FALSE, which= from:to, plot.method = "smoothScatter", lwd=3, cex.main=cexmain,cex=cexval)
          }
  		mtext(paste(tmain,ifelse(perGroup,k,""),ifelse(nplots>1,paste(l,"/",nplots),"")),side = 3, outer = TRUE, font = 2, cex = 2)
          dev.off()
        }
      }
    }

  }


  #####################
  ## plotArrayLayout ##
  #####################

  plotArrayLayout <- function(Data,aType=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24) {
    #note that this plot is based on the original Affymetrix cdf file
    #with updated annotations, some originale control probesets could be regular ones, many regular ones will be discarded

    #stop if array type is not provided
    if(is.null(aType)) stop("The aType parameters must be provided")

    #copy first array as a model
    data.tmp <- Data[,1]

    #set every probe to intensity 1
    exprs(data.tmp)[] <- 1

    #determine which are the perfect match rows in exprs
    rn_pm <- sort(as.numeric(rownames(pm(data.tmp))))
    rn_mm <- NULL
    if(aType=="PMMM") {
      #determine which are the mismatch rows in exprs if present
      rn_mm <- sort(as.numeric(rownames(mm(data.tmp))))
      #set intensity of all mismatch probes if present to 2
      exprs(data.tmp)[rn_mm,] <- 2
    }

    #determine control probes
    affx <- grep("AFFX",probeNames(data.tmp),fixed=TRUE)

    require("ArrayTools", quietly = TRUE)

    library <- paste(substr(Data@annotation,1,nchar(Data@annotation)-2),"transcriptcluster.db",sep="")
    suppressWarnings(eval(parse("",-1,paste("require(",library,", quietly = TRUE)",sep=""))))
    #get probe annotations
    all<-NULL
    eval(parse("",-1,paste("try(all2<-",substr(library,1,nchar(library)-3),"ACCNUM,TRUE)",sep="")))
    if(exists("all2")) {
      eval(parse("",-1,paste("try(all<-toTable(all2),TRUE)",sep="")))
    }

    dataTable <- paste(substr(Data@annotation,1,nchar(Data@annotation)-2),"CONTROL",sep="")
    suppressWarnings(eval(parse("",-1,paste("data(",dataTable,")",sep="")))) #ArrayTools
    cntrl <- NULL
    try(cntrl <- get(dataTable),TRUE)

    #controls<-all[all[,2]=="---",1]
    #sum(controls!=sort(cntrl[,1])) #should be 0, so all IDs are equal

    lab_add <- ""
    if(length(affx)>0) {
      #set intensity of control probes to 5 for matches and to 8 for mismatches if present
      pm(data.tmp)[affx,] <- 5
      if(aType=="PMMM") mm(data.tmp)[affx,] <- 8
    }

    #if no cntrl object given: warn and skip plot
    if(!is.null(cntrl)) {
      #set control probes to intensity 2 (affx, intron, exon)
      #these normally are the affx, intron, and exon probes as ag controls are not present in pm
      pm(data.tmp)[probeNames(data.tmp) %in% cntrl[,1],] <- 5
      if(aType=="PMMM") mm(data.tmp)[probeNames(data.tmp) %in% cntrl[,1],] <- 8
    }

    #if given, set probes missing from list of annotated probes to NA
    if(!is.null(all)) {
      missingAnnotation <- !(probeNames(data.tmp) %in% all[,1])
  	if(sum(missingAnnotation)>0) {
  		pm(data.tmp)[missingAnnotation,] <- NA
  		if(aType=="PMMM") mm(data.tmp)[missingAnnotation,] <- NA
  		lab_add <- " - white: other unannotated probe"
  	}
    }

    #set intensity of all probes in exprs, but not in pm (or if present mm) to 15 (all numbers chosen for the colors to work)
    exprs(data.tmp)[setdiff(1:dim(exprs(data.tmp))[1],c(rn_pm,rn_mm)),] <- 15

    #determine colors needed
    colors <- c("black","#222222","blue","lightblue","red")[c(TRUE,aType=="PMMM",TRUE,aType=="PMMM",TRUE)]
    #set up label
    if(aType=="PMMM") {
      xlab <- paste("black/gray: regular match/mismatch probe \n blue/light blue: ",
  	"control match/mismatch probe \n red: unannotated probe (control region)\n",lab_add,sep="")
    } else {
      xlab <- paste("black: regular probe \n blue: control probe \n red: unannotated probe (control region)\n",lab_add,sep="")
    }

    #plot the image
    png("RawDataReferenceArrayLayout.png",width=WIDTH,height=HEIGHT, pointsize=POINTSIZE)
  	par(oma=c(17,2,0,3))
    image(data.tmp,col=colors,xlab=xlab,main="Array reference layout",cex.lab=0.8)
    dev.off()
  }


  ###################
  ## spatialImages ##
  ###################

  spatialImages <- function(Data, Data.pset=NULL, Resid=TRUE, ResSign=TRUE, Raw=TRUE, Weight=TRUE,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24) {

  	if(is.null("Data.pset")&(Resid||ResSign||Weight)) {Data.pset <- fitPLM(Data)}

      splotFun <- function(typeVal){

  	  if(length(sampleNames(Data))<=18) {nPerPage <- 6; nCol<-2; cexval = 1 } # max 3p of 6 plots
  	  if(length(sampleNames(Data))>18) {nPerPage <- 12; nCol<-3; cexval = 0.8 } # max 4p of 12 plots
  	  if(length(sampleNames(Data))>48) {nPerPage <- 20; nCol<-4; cexval = 0.8 }

        nPages <- ceiling(length(sampleNames(Data))/(nPerPage+1))

  	  for(l in 1:nPages){
          from <- nPerPage*l - (nPerPage-1)
          to <-   min(nPerPage*l, length(sampleNames(Data)))
     		png(filename = paste("RawData2DVirtualImage_",typeVal,ifelse(nPages>1,paste("_",l,sep=""),""),".png",sep = ""),
     		               width = WIDTH, height = HEIGHT,pointsize=POINTSIZE)

  		layout(matrix(1:nPerPage,ncol=nCol,byrow=TRUE))
  		par(oma=c(0,0,3,0),mar=c(1,1,3,1))

  		if(typeVal != "raw"){
  	   		for (k in from:to){
  	    		image(Data.pset, which = k, type = typeVal, cex.main = cexval)
  	    	}
  			mtext(paste("2D virtual PLM image for model characteristic:",typeVal,
  			   ifelse(nPages>1,paste(l,"/",nPages),"")),side = 3, outer = TRUE,
  			   font = 2, cex = 1.1)
  	      }else{
  	   		for (k in from:to){
  	    		image(Data[, k], cex.main = cexval)
  	    	}
  			mtext(paste("2D virtual image of raw data intensities",
  			   ifelse(nPages>1,paste(l,"/",nPages),"")),side = 3, outer = TRUE,
  			   font = 2, cex = 1.1)
  	   	  }

  	    dev.off()
        }
      }

      if(Resid) splotFun("resids")
      if(ResSign) splotFun("sign.resids")
      if(Raw) splotFun("raw")
      if(Weight) splotFun("weights")
  }

  #################
  ## array.image ##
  #################

  # alternative function to plot virtual arrays when PLM cannot be run (when > 6 arrays)
  #-------------------------------------------------------------------------------------
  array.image<-function(Data, pcut=NULL, relative=TRUE, symm=relative,
     balance=relative,quantitative=relative,col.mod=1,postfix="",arrays=NULL,
     WIDTH=1000, HEIGHT=1414, POINTSIZE=24){

    # example calls
    #  array.image(rawData,postfix="_balanced")
    #  array.image(rawData,quantitative=FALSE,postfix="_updown")
    #  array.image(rawData,balance=FALSE)
    #  array.image(rawData,relative=FALSE,postfix="_Intensity")
    # or: per experiment group
    #  for (i in levels(experimentFactor)) {
    #    #remove bad characters because will become part of the file name
    #    array.image(rawData[, experimentFactor == i],
    #    postfix = gsub("[:/|\\?<>*\"]","_",i), balance=FALSE)
    #  }

    gc()
    if(!exists("Data"))
      stop("ERROR on spatial graph: data parameter must be specified")
    if(class(Data)!="AffyBatch")
      stop("ERROR on spatial graph: data object must be of class 'AffyBatch'")
    if(!relative & quantitative)
      warning("When relative is set to FALSE, setting quantitative to TRUE has no effect")
    if(is.null(pcut)) {
      if(relative) {
        pcut<-0.001
      } else {
  	  pcut<-c(0,0.1)
      }
    }
    if(length(pcut)==1) pcut <- c(pcut,pcut)
    if(sum(pcut<0 | pcut>.5)>0)
      stop("ERROR on spatial graph: pcut value(s) has/ve to be in the interval [0, 0.5]")
    if(is.null(arrays)) arrays <- 1:dim(exprs(Data))[2]
    cat("---preprocessing data---\n")
    data2<-Data
    exprs(data2)<-log(exprs(data2),2)
    if(relative) {
      cat("computing median array\n")
    av_array<-rowMedians(exprs(data2),na.rm=TRUE)
    cat("computing relative expression compared to median array\n")
    exprs(data2)<-exprs(data2)-av_array
    }
    #image function takes the log when operating on AffyBatch function!
    #so we visualise the log ratios as compared to the data set median
    #balance or not?
    if(balance) {
      cat("balancing for array-wide intensity differences\n")
      med <- median(exprs(data2),na.rm=TRUE)
      exprs(data2)<-apply(exprs(data2),2,function(x)
              {x-median(x,na.rm=TRUE)+med})
    }

    if(relative & !quantitative) {
      exprs(data2) <- (exprs(data2)>0)+0
    } else {
      #symmetric or not?
      cat("computing color scale limits\n")
      if(!symm) {
        low.lim <- quantile(exprs(data2),pcut[1],na.rm=TRUE)
        up.lim <- quantile(exprs(data2),1-pcut[2],na.rm=TRUE)
      } else {
        low <- quantile(exprs(data2),pcut[1],na.rm=TRUE)
        up <- quantile(exprs(data2),1-pcut[2],na.rm=TRUE)
        low.lim <- min(low,-up)
        up.lim <- max(-low,up)
      }

      cat("replacing off-scale intensities by scale limits\n")
      exprs(data2)[exprs(data2)<low.lim] <- low.lim
      exprs(data2)[exprs(data2)>up.lim] <- up.lim
    }

    cat("---creating images---")
    col.range.part <- round(((seq(0,127^col.mod,length.out=128))^(1/col.mod)))
    if(relative) {
      col.range <- c(rgb(col.range.part, col.range.part, 127, max = 127),
          rgb(127, rev(col.range.part), rev(col.range.part), max = 127))
    } else {
      #col.range <- heat.colors(170)[1:128]
  	col.range <- c(rgb(col.range.part, col.range.part, 127, max = 127),
          rgb(127, rev(col.range.part), rev(col.range.part), max = 127))
  	#col.range <- c(rgb(127, col.range.part, 0, max = 127))
    }

    exprs(data2) <- 2^exprs(data2) #image automatically takes log

   	if(length(sampleNames(Data))<=18) {nPerPage <- 6; nCol<-2; cexval = 1 } # max 3p of 6 plots
  	if(length(sampleNames(Data))>18) {nPerPage <- 12; nCol<-3; cexval = 0.8 } # max 4p of 12 plots
  	if(length(sampleNames(Data))>48) {nPerPage <- 20; nCol<-4; cexval = 0.8 }

  	nPages <- ceiling(length(sampleNames(Data))/(nPerPage+1))

  	for(l in 1:nPages){
  		from <- nPerPage*l - (nPerPage-1)
  		to <-   min(nPerPage*l, length(sampleNames(Data)))
  		png(filename = paste("RawDataArray.image",postfix,
  		    ifelse(nPages>1,paste("_",l,sep=""),""),".png",sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
  		layout(matrix(1:nPerPage,ncol=nCol,byrow=TRUE))
  		par(oma=c(0,0,3,0),mar=c(1,1,3,1),cex.main=cexval)
  		cat("\nimage",l,sep="")

  		for (k in from:to){
  		  cat(".")
  		  image(data2[,arrays][,k],col=col.range)
  		}
  		mtext(paste("2D virtual images",ifelse(nPages>1,paste(l,"/",nPages),"")),
  			side = 3, outer = TRUE, font = 2, cex = 1.1)
  		dev.off()
  	}
  	gc()
  	cat("\n")
  }


  ###############
  ## PNposPlot ##
  ###############

  PNposPlot <- function(Data,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24) {
    png(file = "RawDataPosNegPositions.png",width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(oma=c(17,0,0,0),srt=90)
    borderQC2(Data)
    mtext("COI absolute values should be < 0.5\n", side=3, cex=0.7)
    dev.off()
  }

  #############
  ## nuseFun ##
  #############

  nuseFun <- function(Data, Data.pset=NULL, experimentFactor=NULL, plotColors=NULL, legendColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
    if(is.null(plotColors)) stop("The 'plotColors' parameter must be specified")
    if(is.null(legendColors)) stop("The 'legendColors' parameter must be specified")

    if(is.null(Data.pset)) {Data.pset <- fitPLM(Data)}

    png(file = "RawDataNUSE.png",width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(oma=c(17,0,0,0), cex.axis=1)
    NUSE(Data.pset, col = plotColors,main="Normalized Unscaled Standard Errors (NUSE)", axes = FALSE)
  	if(length(sampleNames(Data))<MAXARRAY){
  		cexval <- 0.65
  	}else{
  		cexval <- 0.45
  	}
    axis(1,at=1:length(sampleNames(Data)),las=2,labels=sampleNames(Data), cex.axis=cexval)
    axis(2, cex.axis=0.7)
    if(length(levels(experimentFactor))>1){
       legend("topright", levels(experimentFactor), col = legendColors,
       fill = legendColors, cex = 0.7, bg = "white", bty = "o")
    }

    mtext("NUSE median value should be < 1.1\n", side=3, font=1, cex=0.7)
    dev.off()
  }


  ############
  ## rleFun ##
  ############

  rleFun <- function(Data, Data.pset=NULL, experimentFactor=NULL, plotColors=NULL, legendColors=NULL,
    WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

    if(is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
    if(is.null(plotColors)) stop("The 'plotColors' parameter must be specified")
    if(is.null(legendColors)) stop("The 'legendColors' parameter must be specified")

    if(is.null(Data.pset)) {Data.pset <- fitPLM(Data)}

    png(file = "RawDataRLE.png",width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    par(oma=c(17,0,0,0), cex.axis=1)
    Mbox(Data.pset, col = plotColors, main = "Relative Log Expression (RLE)", axes = FALSE)
  	if(length(sampleNames(Data))<MAXARRAY){
  		cexval <- 0.65
  	}else{
  		cexval <- 0.45
  	}
    axis(1,at=1:length(sampleNames(Data)),las=2,labels=sampleNames(Data), cex.axis=cexval)
    axis(2, cex.axis=0.7)
    if(length(levels(experimentFactor))>1){
    legend("topright", levels(experimentFactor), col = legendColors,
       fill = legendColors, cex = 0.7, bg = "white", bty = "o")
    }
    mtext("RLE distributions should be centered around 0\n", side=3, font=1, cex=0.7)
    dev.off()
  }


  ###############
  ## correlFun ##
  ###############

  correlFun <- function(Data, clusterOption1="pearson", clusterOption2="ward", normMeth="",
    experimentFactor=NULL, legendColors=NULL, WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41){

    if(is.null(experimentFactor)) stop("the 'exerimentFactor' parameter is required")
    if(is.null(legendColors)) stop("the 'legendColors' parameter is required")

    if(class(Data) == "AffyBatch") {
      Type <- "Raw"
      text1 <- "Raw data correlation plot"
    } else {
      if(normMeth == "") stop("When providing a normalized data object, the normMeth parameter is required")
      Type <- "Norm"
      text1 <- paste("Array correlation plot\nafter",normMeth,"normalization")
    }
    if(length(sampleNames(Data))<2) {
        warning("Only one array in dataset, no correlation plot made")
    } else {
  	  png(file = paste(Type,"DataArrayCorrelation.png",sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
   	  if(length(sampleNames(Data))<MAXARRAY) {
  	    par(oma=c(17,0,0,0),cex.axis=0.7,cex.main=0.8)
  	    #subval <- 10
  	  } else {
  	    par(oma=c(17,0,0,0),srt=90,las=2,cex.axis=0.5,cex.main=0.8)
  	    #subval <- 16
  	  }

      #note: for computing array correlation, euclidean would not make sense
      #only use euclidean distance to compute the similarity of the correlation vectors for the arrays
      COpt1 <- "pearson"
      if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
      crp <- cor(exprs(Data), use="complete.obs", method=COpt1)

      text1 <- paste(text1,"\ncorrelation method:",COpt1,"\ncluster method:",clusterOption2)

      switch(tolower(clusterOption1),
        "pearson" = {
          my.dist <- function(x) cor.dist(x, abs=FALSE)
        },
        "spearman" = {
          my.dist <- function(x) spearman.dist(x, abs=FALSE)
        },
        "euclidean" = {
          my.dist <- function(x) euc(x)
        }
      )

      my.hclust <- function(d) hclust(d, method=clusterOption2)

      #in order to create some space to put colored symbols as well
      #sampleNames(Data) <- paste(sampleNames(Data)," ")

      sideColors <- legendColors[as.numeric(experimentFactor)]

      heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", symm=TRUE, density.info="density",
                main=text1, dendrogram="row", ColSideColors=sideColors)

      #correlationPlot(Data)
      #axis(1,side=3,at=seq(from=0.5, to=(length(sampleNames(Data)))-0.5,by=1),
      #    labels=substr(as.character(sampleNames(Data)),1,subval),las=2)
      #par(srt=0)
  	  #plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE,
  	  #    frame.plot = FALSE, xlim = c(0, 2), ylim = c(0,2))
  	  #text(1,1,text1,cex=1)
    	dev.off()
    }
  }


  ################
  ## clusterFun ##
  ################

  clusterFun <- function(Data, experimentFactor=NULL, clusterOption1="pearson", clusterOption2="ward", normMeth="",
    plotColors=NULL, legendColors=NULL, plotSymbols=NULL, legendSymbols=NULL, WIDTH=1000, HEIGHT=1414, POINTSIZE=24, MAXARRAY=41) {

    if(is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
    if(is.null(plotSymbols)) stop("the 'plotSymbols' parameter is required")
    if(is.null(legendSymbols)) stop("the 'legendSymbols' parameter is required")

    if(class(Data)=="AffyBatch") {
      Type <- "Raw"
      main <- "Cluster dendrogram of raw data"
    } else {
      if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
      Type <- "Norm"
      main <- paste("Cluster dendrogram of",normMeth,"normalized data")
    }
    if(length(sampleNames(Data))<3) {
      warning("Only ",length(sampleNames(Data))," sample(s) in dataset, no clustering plot made")
    } else {
      switch(tolower(clusterOption1),
        "pearson" = {
          correl <- cor.dist(t(exprs(Data)),abs=FALSE)
        },
        "spearman" = {
          correl <- spearman.dist(t(exprs(Data)),abs=FALSE)
        },
        "euclidean" = {
          correl <- euc(t(exprs(Data)))
        }
      )
      clust <- hclust(correl, method = tolower(clusterOption2))
      png(file = paste(Type,"DataCluster_",clusterOption1,".png",sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
  	  if(length(sampleNames(Data))<MAXARRAY) {
  		  cexval1 <- 0.75
  		  cexval2 <- 1.23
  		  cexval3 <- 0.55
  	  } else {
    	  cexval1 <- 0.55
    	  cexval2 <- 1.6
        cexval3 <- 0.41
  	  }
      par(cex=cexval1,oma=c(14,1,0,0))
      par(cex.axis=cexval2,cex.lab=cexval2,cex.main=cexval2)
      plot(clust, hang=-1, main=main, xlab=paste("distance:",clusterOption1), sub=paste(" cluster method:",clusterOption2))
      points(1:length(clust$order),rep(0,length(clust$order)),pch=15,col="white",cex=1.5)
      points(1:length(clust$order),rep(0,length(clust$order)),pch=plotSymbols[clust$order],col=plotColors[clust$order])
      if(length(levels(experimentFactor))>1) {
        legend("topright",levels(experimentFactor),
        pch=legendSymbols,col=legendColors)
      }
      par(cex=cexval3)
      dev.off()
    }
  }


  ############
  ## pcaFun ##
  ############

  pcaFun <- function(Data, experimentFactor=NULL, normMeth="", scaled_pca=TRUE, plotColors=NULL,
     legendColors=NULL, plotSymbols=NULL, legendSymbols=NULL, namesInPlot=FALSE, WIDTH=1000, HEIGHT=1414, POINTSIZE=24){
    # Scaled PCA by default
    if(is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
    if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
    if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
    if(is.null(plotSymbols)) stop("the 'plotSymbols' parameter is required")
    if(is.null(legendSymbols)) stop("the 'legendSymbols' parameter is required")

    if(length(sampleNames(Data))<3) {
      warning("Only",length(sampleNames(Data)),"sample(s) in dataset, no PCA plot made")
    } else {

      if(class(Data) == "AffyBatch"){
        #raw data
        tmain <- "PCA analysis of Raw data"
        Type="Raw"
      } else{
        if(normMeth=="") stop("When providing a normalised data object, the normMeth parameter is required")
        tmain <- paste("PCA analysis after", normMeth, "normalization")
        Type <- "Norm"
      }

      pca1 <- NULL
      try(pca1 <- prcomp(t(exprs(Data)[apply(exprs(Data),1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=scaled_pca),TRUE)
      if(is.null(pca1) & scaled_pca) {
        try(pca1 <- prcomp(t(exprs(Data)[apply(exprs(Data),1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=FALSE),TRUE)
        if(!is.null(pca1)) warning("pca with scaling unsuccessful, successfully retried without scaling")
      }
      if(!is.null(pca1)) {
        perc_expl1 <- round(((pca1$sdev[1:3]^2)/sum(pca1$sdev^2))*100,2)

        cex.circle <- 1.5
        cex.text <- 0.7
  	    cex.legend <- 0.75
        tcol <- "#444444"

  	    png(file = paste(Type,"DataPCAanalysis.png",sep=""), width=WIDTH+200*(!namesInPlot), height=HEIGHT+283*(!namesInPlot),
  	      pointsize=POINTSIZE)

        if(!namesInPlot) {
  	     layout(rbind(c(1,1,2,2,5),c(3,3,4,4,5)))
  	    } else {
  	     layout(rbind(c(1,1,2,2),c(1,1,2,2),c(3,3,4,4),c(3,3,4,4)))
  	    }
  	    par(oma=c(20,0,5,0))
          plot(pca1$x[,1],pca1$x[,2],cex=cex.circle,pch=plotSymbols,
            col=plotColors,xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
            ylab=paste("PC2 (",perc_expl1[2],"%)",sep=""))
          if(namesInPlot) text(pca1$x[,1],pca1$x[,2], sampleNames(Data),pos=4,cex=cex.text,col=tcol)
          plot(pca1$x[,1],pca1$x[,3],cex=cex.circle,pch=plotSymbols,
          col=plotColors,xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
          ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
          if(namesInPlot) text(pca1$x[,1],pca1$x[,3], sampleNames(Data),pos=4,cex=cex.text,col=tcol)
          plot(pca1$x[,2],pca1$x[,3],cex=cex.circle,pch=plotSymbols,
          col=plotColors,xlab=paste("PC2 (",perc_expl1[2],"%)",sep=""),
          ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
          if(namesInPlot) text(pca1$x[,2],pca1$x[,3], sampleNames(Data),pos=4,cex=cex.text,col=tcol)
          barplot((100*pca1$sdev^2)/sum(pca1$sdev^2),xlab="components",ylab="% of total variance explained")

  		  if(namesInPlot) {
          if(length(levels(experimentFactor))>1){
            legend("topright",levels(experimentFactor),
            pch=legendSymbols,col=legendColors,cex=cex.legend)
          }
        } else {
  			  par(mar=c(0,0,0,0))
  			  plot(1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  			  if(length(levels(experimentFactor))>1) {
  		      legend("topleft",c(levels(experimentFactor),"",sampleNames(Data)),
   #             pch=c(rep(20,length(unique(experimentFactor))+1),plotSymbols,
  			      pch=c(legendSymbols,20,plotSymbols),
  		          col=c(legendColors,"white",plotColors),cex=(cex.legend+0.1)
   #             ,fill=c(legendColors,rep("white",length(experimentFactor)+1)),
   #             border=c(legendColors,rep("white",length(experimentFactor)+1))
  			    )
  			  } else {
  			    legend("topleft",sampleNames(Data),pch=plotSymbols,
  			      col=plotColors,cex=0.7, bty = "n")
  			  }
        }

  		  mtext(tmain, side = 3, outer = TRUE, font = 2, cex = 1.2)
        dev.off()
      } else {
  	    warning("PCA on the",Type,"data set unsuccessful, image skipped")
      }
    }
  }


  if(!exists("rawData")) { # this is the case when the script is run locally
    setwd(DATA.DIR)
    rawData <- ReadAffy()
    print("Raw data have been loaded in R")
  }
  if(!exists("libdir")) { # libdir exists only for GenePattern usage
    setwd(SCRIPT.DIR)
    setwd(WORK.DIR)
  }

  # Make sure that the CDF environment works
  rawData <- addStandardCDFenv(rawData)   # if already works, won't be changed

  # Verify the array type (PMMM or PMonly)
  aType <- getArrayType(rawData)

  # When refName does not exist, use the empty string
  if(!exists("refName")) refName <- ""

  ###############################################################################
  # Create array groups and array names                                         #
  ###############################################################################

  if(affyPARAM$arrayGroup!=""){
    # Information is available: groups will be created
    # 1- read the arrayGroup file and trim spaces
    # 2- define the array names and classes (experimentFactor)
    if(!exists("DESC.DIR")) DESC.DIR <- ""

    descfile <- paste(DESC.DIR, affyPARAM$arrayGroup, sep="")
    extension<-strsplit(descfile,"\\.")
    extension<-paste(".",extension[[1]][length(extension[[1]])],sep="")
    description = NULL;
    switch(extension,
           ".txt" = description<-trim(read.delim(descfile, fill = FALSE, as.is=TRUE)),
           ".csv" = description<-trim(read.csv(descfile, fill = FALSE, as.is=TRUE)),
           ".xls" = {library(gdata); description<-trim(read.xls(descfile, as.is=TRUE))},
           ".xlsx" = {library(gdata); description<-trim(read.xls(descfile, as.is=TRUE))}
  	)
    if(is.null(description)) stop(paste("extension",extension,"not recognised"))

   # description <- trim(read.table(paste(DESC.DIR, affyPARAM$arrayGroup, sep=""),
   # 	  header = TRUE, as.is = TRUE, sep="\t"))

    if(length(grep(".CEL",toupper(colnames(description)[1]),
      ignore.case = TRUE))>0) {
      stop(paste("The description file may not contain a header, as the first",
       	"column header seems to be a CEL file name"))
    }
    file_order <- match(description[,1],sampleNames(rawData))
    if(sum(is.na(file_order)) > 0) stop("file names in data directory and file names in description file do not match")
    if(length(unique(file_order)) < length(file_order)) stop("file names in description file are not unique")
    rawData <- rawData[,file_order]

    sampleNames(rawData)<- as.character(description[,2])
    experimentFactor <- factor(description[,3])

    # if required reorder the arrays according to group levels in order to keep
    # groups together in all plots
    if(affyPARAM$reOrder) {
      rawData <- rawData[,order(experimentFactor)]
      experimentFactor <- experimentFactor[order(experimentFactor)]
    }
  } else {
    # No information: arrays will be computed/colored independently
    sampleNames(rawData) <- as.character(sampleNames(rawData))
    experimentFactor <- factor(rep(1, length(sampleNames(rawData))))
    description <- cbind(sampleNames(rawData),sampleNames(rawData),
      experimentFactor)
    colnames(description) <- c("ArrayDataFile","SourceName","FactorValue")
  }

  # Create colorset for the array groups
  #-------------------------------------
  colList <- colorsByFactor(experimentFactor)
  plotColors <- colList$plotColors
  legendColors <- colList$legendColors
  rm(colList)

  # Create symbolset for the array groups
  #--------------------------------------
  plotSymbols <- 18-as.numeric(experimentFactor)
  legendSymbols <- sort(plotSymbols, decreasing=TRUE)

  ###############################################################################
  # Define display parameters for the images			                          #
  ###############################################################################

  WIDTH <- 1000
  HEIGHT <- 1414
  POINTSIZE <- 24
  if(!exists("maxArray")) maxArray <- 41

  ###############################################################################
  # Calculate the indicator values and begin the report                         #
  ###############################################################################

  #create a cover sheet for the report to be created later
  #and create a page indicating the naming and grouping used
  coverAndKeyPlot(description, refName,WIDTH=WIDTH,HEIGHT=HEIGHT)

  #create a table with several QC indicators
  if(affyPARAM$samplePrep || affyPARAM$ratio || affyPARAM$hybrid || affyPARAM$percPres || affyPARAM$bgPlot || affyPARAM$scaleFact) {

    # The indicators are calculated only for PM-MM arrays as the calculation
    # based on MAS5 does not work for PM-only arrays

    quality <- NULL
    try(quality <- qc(rawData),TRUE) # calculate Affymetrix quality data for PMMM
    if(is.null(quality)) {
      warning("Plots based on the simpleaffy qc function cannot be created for this chip type")
    }

    if(affyPARAM$samplePrep) {
      # find the data
      try(yack <- yaqc(rawData),TRUE)
      if(exists("yack")) {
        spnames<-rownames(yack@morespikes[grep("(lys|phe|thr|dap).*3", # only 3'
        rownames(yack@morespikes), ignore.case = TRUE),])
        sprep<-t(yack@morespikes[spnames,])
      } else {
        sprep <- NULL
        warning("Plots based on the yaqc function cannot be created for this chip type")
      }

      try({calls<-detection.p.val(rawData)$call
      lys<-calls[rownames(calls)[grep("lys.*3",rownames(calls),ignore.case=TRUE)],]
      rm(calls)},TRUE)
      if(!exists("lys")) {
        lys <- NULL
        warning("Plots based on the detection.p.val function cannot be created for this chip type")
      }else{
  		if(length(lys) > length(sampleNames(rawData))) { lys<-lys[1,] }
      }
    }

    QCtablePlot(rawData,quality,sprep,lys,samplePrep=affyPARAM$samplePrep,ratio=affyPARAM$ratio,
        hybrid=affyPARAM$hybrid,percPres=affyPARAM$percPres,bgPlot=affyPARAM$bgPlot,scaleFact=affyPARAM$scaleFact,
  	  WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
  }


  ###############################################################################
  # Raw data Quality Control graphs                                             #
  ###############################################################################

  print("Graphs ready to be computed")

  # 1.1 Sample prep controls
  #-------------------------

  if(affyPARAM$samplePrep && !is.null(sprep) && !is.null(lys)) {
    print ("   plot sample prep controls"  )
    samplePrepPlot(rawData,sprep,lys,plotColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  # 1.2 3'/5' ratio - only for PM-MM arrays
  #----------------------------------------

  if(affyPARAM$ratio && !is.null(quality)) {
    print ("   plot beta-actin & GAPDH 3'/5' ratio")
    ratioPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  # 1.3 RNA degradation plot
  #-------------------------

  if(affyPARAM$degPlot) {
    print ("   plot degradation plot"  )
    RNAdegPlot(rawData,plotColors=plotColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  ###############################################################################
  # 2.1 Spike-in controls - only for PM-MM arrays
  #----------------------------------------------

  if(affyPARAM$hybrid && !is.null(quality)) {
    print ("   plot spike-in hybridization controls"  )
    hybridPlot(rawData,quality=quality,plotColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  # 2.2 Background intensities - only for PM-MM arrays
  #---------------------------------------------------

  if(affyPARAM$bgPlot && !is.null(quality)) {
    print ("   plot background intensities"  )
    backgroundPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  # 2.3 Percent present - only for PM-MM arrays
  #---------------------------------------------

  if(affyPARAM$percPres && !is.null(quality)) {
    print ("   plot percent present"  )
    percPresPlot(rawData,quality=quality,experimentFactor,plotColors,legendColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  # 2.4 Table of PMA-calls based on the MAS5 algorithm - only for PM-MM arrays
  #---------------------------------------------------------------------------

  if(affyPARAM$PMAcalls) {
    if(affyPARAM$customCDF) {
      if(affyPARAM$species=="") {
        warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
        species <- deduceSpecies(rawData@annotation)
      }
  	if(affyPARAM$species!=""){
  		PMAtable <- computePMAtable(rawData,affyPARAM$customCDF,affyPARAM$species,affyPARAM$CDFtype)
  	}else{
  		warning("Could not define species; the CDF will not be changed")
  		PMAtable <- computePMAtable(rawData,affyPARAM$customCDF)
  	}
    } else {
      PMAtable <- computePMAtable(rawData,affyPARAM$customCDF)
    }
    if(!is.null(PMAtable)) {
      write.table(PMAtable, "PMAtable.txt", sep="\t", row.names=FALSE,
  	  col.names=TRUE, quote=FALSE)
    }
  }

  # 2.5 Pos and Neg control distribution
  #-------------------------------------

  if(affyPARAM$posnegDistrib) {
    print ("   plot pos & neg control distribution"  )
    PNdistrPlot(rawData,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
  }

  # 2.6 affx control profiles and boxplot
  #--------------------------------------

  if(affyPARAM$controlPlot) {
    print ("   plot control profiles and/or boxplots")
    controlPlots(rawData,plotColors,experimentFactor,legendColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  ###############################################################################
  # 3.1.1 Scale factor - only for PM-MM arrays
  #-------------------------------------------

  if(affyPARAM$scaleFact && !is.null(quality)) {
    print ("   plot scale factors")
    scaleFactPlot(rawData,quality=quality,experimentFactor,plotColors,
       legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
  	 MAXARRAY=maxArray)
  }

  # 3.1.2 Boxplot of raw log-intensities
  #-------------------------------------

  if(affyPARAM$boxplotRaw){
    print ("   plot boxplot for raw intensities")
    boxplotFun(Data=rawData, experimentFactor, plotColors, legendColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  # 3.1.3 Density histogram of raw log-intensities
  #-----------------------------------------------

  if(affyPARAM$densityRaw){
    print ("   plot density histogram for raw intensities")
    densityFun(Data=rawData, plotColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
    #densityFunUnsmoothed(Data=rawData, plotColors,
    #  WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  # 3.2.1 MA-plot or raw data
  #--------------------------

  if(affyPARAM$MARaw){
    print ("   MA-plots for raw intensities")
    maFun(Data=rawData, experimentFactor, perGroup=(MAOption1=="group"),
       aType=aType,WIDTH=WIDTH,HEIGHT=HEIGHT,MAXARRAY=maxArray)
  }

  # 3.3.1 Plot of the array layout
  #-------------------------------

  if(affyPARAM$layoutPlot) {
    print ("   plot array reference layout")
    plotArrayLayout(rawData,aType,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
  }

  # 3.3.2 Pos and Neg control Position
  #-----------------------------------

  if(affyPARAM$posnegCOI){
    print ("   Pos/Neg COI")
    PNposPlot(rawData,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
  }

  # 3.3.3.1 Create PLM object
  #--------------------------

  # fit a probe level model on the raw data, used by nuse and rle plot as well
    rawData.pset <- NULL
    if(affyPARAM$spatialImage || affyPARAM$PLMimage || affyPARAM$Nuse || affyPARAM$Rle) {
    print ("   Fit a probe level model (PLM) on the raw data")
      rawData.pset <- fitPLM(rawData)
    }


  # 3.3.3.2 Spatial images
  #---------------------

  if(affyPARAM$spatialImage) {
    print ("   2D virtual images")
    valtry<-try(spatialImages(rawData, Data.pset=rawData.pset, TRUE,FALSE,FALSE,FALSE,
  	          WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE),
  		      silent=TRUE)
    if(class(valtry)=="try-error") {
  	print("      Use array.image instead of spatialImages function")
  	if(length(sampleNames(rawData))>6){
  		# Usage of a median array is interesting when there are enough arrays
  		array.image(rawData,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
  	}else{
  		# Usage when few arrays in dataset (one page for 3 arrays -> max: 2 pages)
  		array.image(rawData,relative=FALSE,col.mod=4,symm=TRUE,WIDTH=WIDTH,
  		  HEIGHT=HEIGHT,POINTSIZE=POINTSIZE)
  	}
    }
  }

  # 3.3.3.3 PLM images
  #---------------------

  if(affyPARAM$PLMimage) {
    print ("   Complete set of 2D PLM images")
    valtry<-try(spatialImages(rawData, Data.pset=rawData.pset, TRUE, TRUE, TRUE, TRUE,
  	            WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray),
  				silent=TRUE)
    if(class(valtry)=="try-error") {
  	print("      Could not create the PLM images.")
    }
  }

  # 3.4.1 NUSE
  #-----------

  if(affyPARAM$Nuse){
    print ("   NUSE boxplot")
    nuseFun(rawData, Data.pset=rawData.pset, experimentFactor, plotColors,
       legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
  	 MAXARRAY=maxArray)
  }

  # 3.4.2 RLE
  #----------

  if(affyPARAM$Rle){
    print ("   RLE boxplot")
    rleFun(rawData, Data.pset=rawData.pset, experimentFactor, plotColors,
       legendColors,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
  	 MAXARRAY=maxArray)
  }

  ###############################################################################
  # 4.1 Correlation Plot  of raw data
  #----------------------------------

  if(affyPARAM$correlRaw){
    print ("   Correlation plot of raw data")
    correlFun(Data=rawData, experimentFactor=experimentFactor, legendColors=legendColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  # 4.2 PCA analysis of raw data
  #-----------------------------

  if(affyPARAM$PCARaw){
    print("   PCA analysis of raw data")
    pcaFun(Data=rawData, experimentFactor=experimentFactor,
  	plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
  	legendSymbols=legendSymbols, namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
  	(length(sampleNames(rawData))<=(maxArray/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,
  	POINTSIZE=POINTSIZE)
  }

  # 4.3 Hierarchical Clustering of raw data
  #-----------------------------------------

  if(affyPARAM$clusterRaw){
    print ("   Hierarchical clustering of raw data")
    clusterFun(Data=rawData, experimentFactor=experimentFactor,
     clusterOption1=affyPARAM$clusterOption1, clusterOption2=affyPARAM$clusterOption2,
     plotColors=plotColors, legendColors=legendColors,
     plotSymbols=plotSymbols, legendSymbols=legendSymbols,
     WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
  }

  ###############################################################################
  # Pre-processing                                                              #
  ###############################################################################

  if (aType == "PMonly") {
    if (affyPARAM$normMeth == "MAS5") {
      warning("MAS5 cannot be applied to PMonly arrays. Changed MAS5 to PLIER")
      affyPARAM$normMeth <- "PLIER"
    }
  #  if (normMeth == "GCRMA") {
  #    warning("GCRMA cannot be applied to PMonly arrays. Changed GCRMA to RMA")
  #    normMeth <- "RMA"
  #  }
  }

  if(affyPARAM$normMeth!="" && affyPARAM$normMeth!="none") {
    if(affyPARAM$customCDF) {
      if(affyPARAM$species=="") {
        warning("Species has not been set and custom cdf requested, attempting to deduce species for chip type")
        species <- deduceSpecies(rawData@annotation)
      }
  	if(affyPARAM$species!=""){
  		normData <- normalizeData(rawData,affyPARAM$normMeth,perGroup=(normOption1=="group"),
  		  experimentFactor, aType=aType, affyPARAM$customCDF, affyPARAM$species, affyPARAM$CDFtype,WIDTH=WIDTH,
  		  HEIGHT=HEIGHT)
  	}else{
  		warning("Could not define species; the CDF will not be changed")
  		normData <- normalizeData(rawData,affyPARAM$normMeth,perGroup=(normOption1=="group"),
  		  experimentFactor, aType=aType, affyPARAM$customCDF,WIDTH=WIDTH,HEIGHT=HEIGHT)
  	}

    } else {
      normData <- normalizeData(rawData,affyPARAM$normMeth,perGroup=(normOption1=="group"),
  	  experimentFactor, aType=aType, affyPARAM$customCDF,WIDTH=WIDTH,HEIGHT=HEIGHT)
    }
  }

  if((affyPARAM$boxplotNorm || affyPARAM$densityNorm || affyPARAM$MANorm || affyPARAM$correlNorm || affyPARAM$clusterNorm ||
      affyPARAM$PCANorm) && ((affyPARAM$normMeth=="") || (affyPARAM$normMeth=="none"))) {
    warning("One or more QC plots of normalized data requested, but no normalization selected, plots will be omitted")
  } else {

  ###############################################################################
  # Evaluation of the pre-processing                                            #
  ###############################################################################

  # 5.1 Make a Box-plot of the normalized data
  #-------------------------------------------

    if(affyPARAM$boxplotNorm){
      print ("   plot boxplot for normalized intensities")
      boxplotFun(Data=normData, experimentFactor, plotColors, legendColors,
  	  normMeth=affyPARAM$normMeth,WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,
  	  MAXARRAY=maxArray)
    }

  # 5.2 Make a Density histogram of the normalized data
  #----------------------------------------------------

    if(affyPARAM$densityNorm){
      print ("   plot density histogram for normalized intensities")
      densityFun(Data=normData, plotColors, normMeth=affyPARAM$normMeth,
        WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
      #densityFunUnsmoothed(Data=normData, plotColors, normMeth=normMeth,
      #  WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
    }

  # 5.3 Make separate MA-plots for each group on normalized data
  #-------------------------------------------------------------

    if(affyPARAM$MANorm){
      print ("   MA-plots for normalized intensities")
      maFun(Data=normData, experimentFactor, perGroup=(MAOption1=="group"),
  	 normMeth=affyPARAM$normMeth,WIDTH=WIDTH,HEIGHT=HEIGHT,MAXARRAY=maxArray)
    }

  # 5.4 Make correlation plots on normalized data
  #----------------------------------------------

    if(affyPARAM$correlNorm){
      print ("   Correlation plot of normalized data")
      correlFun(Data=normData, normMeth=affyPARAM$normMeth, experimentFactor=experimentFactor, legendColors=legendColors,
       WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
    }

  # 5.5 PCA analysis of normalized data
  # -----------------------------------

    if(affyPARAM$PCANorm){
      print("   PCA graph for normalized data")
      pcaFun(Data=normData, experimentFactor=experimentFactor,normMeth=affyPARAM$normMeth,
  	  plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
  	  legendSymbols=legendSymbols, namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
  	  (length(sampleNames(rawData))<=(maxArray/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,
  	  POINTSIZE=POINTSIZE)
    }

  # 5.6 Make hierarchical clustering on normalized data
  #----------------------------------------------------

    if(affyPARAM$clusterNorm){
      print ("   Hierarchical clustering of normalized data")
      clusterFun(Data=normData, experimentFactor=experimentFactor,
      clusterOption1=affyPARAM$clusterOption1, clusterOption2=affyPARAM$clusterOption2,
      normMeth=affyPARAM$normMeth, plotColors = plotColors, legendColors = legendColors,
      plotSymbols=plotSymbols, legendSymbols=legendSymbols,
      WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=maxArray)
    }
  }

  ###############################################################################
  # Prepare the output data                          				              #
  ###############################################################################

  # Export the normalized data

  if((affyPARAM$normMeth=="") || (affyPARAM$normMeth=="none")) {
    warning("No normalization selected, normalized data table not saved")
  } else {
    print("Saving normalized data table")

    normDataTable <- createNormDataTable(normData, customCDF=(sum(featureNames(normData)!=featureNames(rawData)[1:length(featureNames(normData))])>0), affyPARAM$species, affyPARAM$CDFtype)

    #output normalised expression data to file
    refName <- sub("(_\\d{4}-\\d{2}-\\d{2}_\\d{2}-\\d{2}_\\d{2})", "", refName)
    normFileName <- paste(affyPARAM$normMeth,"NormData_",refName,".txt",sep="")
    print(paste("Normalized data table:", normFileName))
    write.table(normDataTable, normFileName, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  }

  # clean R: (or quit without saving the environment...)
  # rm(list = ls())

  if(!is.null(warnings())) warnings()


  return(DATA)
}

############################################################################
################### TIMING AND REPORTING ###################################
############################################################################

cat("\nTIMES   record=do.load(data, 'condition').....\n") #DEBUG
TIMES <- addRecord(TIMES, record_name = "processName",
                   record = system.time(gcFirst = T,
                                        DATA_rna <- do.load(affyPARAM)
                                        )
                   )

#cat("\nTIMES   record=do.X(data, 'condition').....\n") #DEBUG
#TIMES <- addRecord(TIMES, record_name = "processName",
#                   record = system.time(gcFirst = T,
#                                        DATA_rna <- do.preprocess(data,"condition")
#                                        )
#                   )


#######

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
