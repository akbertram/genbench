# human liver cohort
# synapse account
# user=synapsedata
# pwd=synapse

source('http://depot.sagebase.org/CRAN.R')
pkgInstall(c("synapseClient"))

library(synapseClient)
synapseLogin('username','password')

# Obtain a pointer and download the data

syn3275753 <- synGet(id='syn3275753')

# Load the data
syn3275753 <- synGet(id='syn3275753', load=T)

## expression data
syn88644 <- synGet(id='syn88644')

# Load the data
syn88644 <- synGet(id='syn88644', load=T)

## genotype
syn89614 <- synGet(id='syn89614')

# Load the data
syn89614 <- synGet(id='syn89614', load=T)

## phenotype
syn4629 <- synGet(id='syn4629')

# Load the data
syn4629 <- synGet(id='syn4629', load=T)

