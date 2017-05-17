#!/usr/bin/env Rscript
param = commandArgs(trailingOnly=TRUE)
setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)

print("This job is running on ")
comput <- system("/bin/hostname")

## Load libraries
.libPaths( c("/share/nordlab/libraries/R/lib") )
library(ChIPQC)
library(BiocParallel)

threads <- param[8]

mcparam <- MulticoreParam(threads)

register(mcparam, default = TRUE)
bpparam("MulticoreParam")

## Construct a ChIPQsample object and retrieving QC metrics
setwd(param[2]);

## Sample data
samples <- read.table(param[1], header = TRUE);

blacklist.file <- param[3];

if(blacklist.file == 'null') {
    sampleset <- ChIPQC(samples, annotation = param[4], chromosomes = NULL, consensus = FALSE, bCount=FALSE, mapQCth = param[5])
} else {
    blacklist.file <- read.table(param[3]);
    sampleset <- ChIPQC(samples, annotation = param[4], chromosomes = NULL, consensus = FALSE, bCount=FALSE, mapQCth = param[5], blacklist = blacklist.file)
}


## Report generation
name <- unlist(strsplit(param[1], "[.]"));
name <- name[1];
name <- unlist(strsplit(name, "[/]"));
name <- tail(name, n=1);

name

facetlist <- param[6];
facetlist <- unlist(strsplit(facetlist, "[,]"));

facetlist

ChIPQCreport(sampleset, facet=TRUE, reportName=name, reportFolder=name, facetBy=facetlist, colourBy=param[7])

save.image(file=paste0(name, ".RData"))
