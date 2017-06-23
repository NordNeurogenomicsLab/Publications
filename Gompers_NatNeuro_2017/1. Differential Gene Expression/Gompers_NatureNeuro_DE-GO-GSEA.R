#!/usr/bin/env Rscript

################################################################################

## Gompers et al. 2017 (Nature Neuroscience)
## Differential expression analysis 
## topGO
## Permutation gene set enrichment analysis

################################################################################

## EDGER
library(edgeR)

cpm.sample.cutoff <- 2
min.cpm <- 10

y <- DGEList(counts=exp.data.matrix,group=group)
keep <- rowSums(cpm(y)>min.cpm) >= cpm.sample.cutoff
y <- y[keep, , keep.lib.sizes=FALSE]

## Perform simple exact test on genotype
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
pseudo.counts <- y$pseudo.counts
sample.count.table <- y$counts

## Perform glm on genotype, with sex and batch as covariates
design <- model.matrix(~as.factor(sex)+as.factor(seq.run)+as.factor(group))
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y,design)
lrt <- glmLRT(fit)
glm.output <- topTags(lrt, n=Inf)

## TOPGO
## All DE genes

library(topGO)
library(GO.db)
BP.output.down <- list(length=length(module.names))
BP.output.up <- list(length=length(module.names))
test.gene.min <- 20
FDR.criteria <- .2

xx <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "ensembl")
background.genes <- final.set[,"Human.Ensembl.Gene.ID"]
geneUniverse <- background.genes

## Down regulated for all DE
for (module.index in 1:length(module.names)) {
  test.genes <- final.set[which(final.set[,"moduleColors"]==module.names[module.index] & final.set[,"FDR"]<FDR.criteria & final.set[,"logFC"]<0),"Human.Ensembl.Gene.ID"]
  if (length(test.genes)>test.gene.min) {
    genesOfInterest <- test.genes
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl", nodeSize=20)
    resultGO <- runTest(myGOdata, algorithm = "weight01", statistic="fisher")
    print(paste(module.names[module.index], "Down"))
    BP.output.down[[module.index]] <- GenTable(myGOdata, resultGO, topNodes=3315)
    print(BP.output.down[[module.index]][1:10,])
  } else {
    print(module.names[module.index])
    print("Too few genes")
  }
}

## Up regulated for all DE
for (module.index in 1:length(module.names)) {
  test.genes <- final.set[which(final.set[,"moduleColors"]==module.names[module.index] & final.set[,"FDR"]<FDR.criteria & final.set[,"logFC"]>0),"Human.Ensembl.Gene.ID"]
  if (length(test.genes)>test.gene.min) {
    genesOfInterest <- test.genes
    geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
    names(geneList) <- geneUniverse
    myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl", nodeSize=20)
    resultGO <- runTest(myGOdata, algorithm = "weight01", statistic="fisher")
    print(paste(module.names[module.index], "Up"))
    BP.output.up[[module.index]] <- GenTable(myGOdata, resultGO, topNodes=3315)
    print(BP.output.up[[module.index]][1:10,])
  } else {
    print(module.names[module.index])
    print("Too few genes")
  }
}


## PERMUTATION GSEA
## Function to perform general permutation test and generate histogram output
geneset.perm.test <- function(binary.score.reference, binary.score.test, iterations, plot.name) {
  binary.score.reference <- ifelse(is.na(binary.score.reference),0,binary.score.reference)
  binary.score.test <- ifelse(is.na(binary.score.test),0,binary.score.test)
  count.criteria <- vector(length=iterations)
  test.size <- length(binary.score.test)
  for (index in 1:iterations) {
    count.criteria[index] <- sum(sample(binary.score.reference, test.size, replace = F))
  }
  x.min <- min(c(count.criteria, sum(binary.score.test))) - 20
  if (x.min < 0) {x.min <- 0}
  x.max<- max(c(count.criteria, sum(binary.score.test))) + 20
  z <- abs(mean(count.criteria) - sum(binary.score.test))/sd(count.criteria)
  p <- 2*pnorm(-abs(z))
  e <- sum(binary.score.test)/mean(count.criteria)
  prop <- sum(binary.score.test)/length(binary.score.test)
  bg <- length(binary.score.reference)
  pdf(file=paste(plot.name, ".pdf", sep=""), height=6, width=10)
  hist(count.criteria, col="gray", xlab="Count", ylab="Frequency", main=paste(plot.name, " \ncount=", sum(binary.score.test), " of ",  length(binary.score.test), "; p-value=", format(p, scientific = T, digits = 2), "; FE=", format(e, scientific = F, digits = 2), "; Prop=", format(prop, scientific = F, digits = 2), "; n(bg)=", bg, sep=""), xlim=c(x.min,x.max), cex=.75)
  abline(v=sum(binary.score.test), lwd=3, col="red")
  dev.off()
  return(list(count.criteria,z,p,e,prop,bg,sum(binary.score.test)))
}