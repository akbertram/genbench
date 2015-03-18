args <- commandArgs(trailingOnly = TRUE)
PATH <- args[1]
NGENES <- args[2]
NPATIENTS <- args[3]
GEO <- paste(PATH, '/GEO-', NGENES, '-', NPATIENTS, sep="")
GO <- paste(PATH, '/GO-', NGENES, '-', NPATIENTS, sep="")
GENES <- paste(PATH, '/GeneMetaData-', NGENES, '-', NPATIENTS, sep="")
PATIENTS <- paste(PATH, '/PatientMetaData-', NGENES, '-', NPATIENTS, sep="")

invisible(lapply(list(GEO, GO, GENES, PATIENTS), function(f)  {
	message("doing ", f)
	saveRDS(read.csv(paste(f, '.txt', sep="")), 
		file=paste(f, '.rds', sep="")) 
}))
