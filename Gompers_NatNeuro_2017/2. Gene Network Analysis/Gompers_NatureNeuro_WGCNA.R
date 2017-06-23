#!/usr/bin/env Rscript

################################################################################

## Gompers et al. 2017 (Nature Neuroscience)
## Weighted Gene Co-Expression Network Analysis (WGCNA) Pipeline
## Partially adapted from Parikshak et al. 2013 (Cell)

################################################################################

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

## IMPORT RPKM x GENES MATRIX 
## Low-expression genes have been removed (from DE analysis)
metData <- read.delim("Chd8.sampleInfo.txt")
rnaData <- read.delim("Chd8.logRPKM.txt")
sampleNames <- names(data.frame(rnaData))
genData <- read.delim("mm9_geneannotation.txt")
datExpr <- data.frame(rnaData)
table(dimnames(t(datExpr))[[1]]==metData$full.names)
sampleTrait <- metData$stage
sizeGrWindow(9, 5)
plotClusterTreeSamples(datExpr = t(datExpr), y=sampleTrait)

## PICK SOFT THRESHOLD
## Soft threshold used, sft=11 (R2 > 0.95) 
bsize = 13000;
SFT = pickSoftThreshold(data=t(datExpr), networkType="signed", corFnc="bicor", verbose=5, powerVector=c(seq(2, 30, by = 1)), blockSize=bsize);
print(SFT);

## CALCULATE SIGNED TOM
softPower=11; 
tmpMulti = as.data.frame(t(datExpr)); 
netData = blockwiseModules(datExpr=tmpMulti, maxBlockSize=13000, networkType="signed",
                             power=softPower, mergeCutHeight=0.15, nThreads=10,
                             saveTOMFileBase=paste("./round1-signed/data/signed-round1"),
                             saveTOMs=TRUE, corType="bicor", minModuleSize=10, pamStage=FALSE,
                             reassignThreshold=1e-10, verbose=3, deepSplit=2); 
geneTree = hclust(1 - TOM, method="average");
mColorh = mLabelh = colorLabels = NULL; 

## SPLIT TREE PARAMETERS
## Minimum module size = 500
## deepSplit parameter = 2
## deepSplit threshold = 0.20
## PAM = off
minModSize = 500;
dthresh = 0.20;
ds = 2;
tree = cutreeHybrid(dendro=geneTree, pamStage=FALSE,
                    minClusterSize=minModSize, cutHeight=NULL,
                    deepSplit=ds, distM=as.matrix(1-TOM)); 
merged = mergeCloseModules(exprData=t(datExpr),colors=tree$labels, cutHeight=dthresh); 

## COMPARE METADATA TO GENERATED MODULES
keep = goodGenes(t(datExpr)); 
thisExpr = datExpr[keep, ];
traitMat = cbind(metData[,"stage"], metData[,"group.all"], metData[,"sex.all"], metData[,"seq.run.all"]);
rownames(traitMat) = metData$full.names;
colnames(traitMat)[1:4] = c("Stage","Group","Sex","SeqRun");
traitMat = traitMat[rownames(traitMat) %in% colnames(datExpr),];
geneSigs = matrix(NA, nrow=4, ncol=nrow(thisExpr));
colnames(geneSigs) = rownames(thisExpr);

for (i in 1:ncol(geneSigs)) {
  exprvec = as.numeric(thisExpr[i,]);
  stager = sqrt(max(summary(lm(exprvec ~ as.factor(traitMat[,"Stage"])))$adj.r.squared, 0));
  groupr = bicor(as.factor(traitMat[,"Group"]), exprvec);
  sexr = bicor(as.factor(traitMat[,"Sex"]), exprvec);
  seqrunr = sqrt(max(summary(lm(exprvec ~ as.factor(traitMat[,"SeqRun"])))$adj.r.squared, 0));
  geneSigs[,i] = c(stager, groupr, sexr, seqrunr); 
}
geneCorr = geneSigs;
rownames(geneCorr) = colnames(traitMat);

## kME AND MODULE STATISTICS
keepDat = t(datExpr[goodGenes(t(datExpr)),]); 
MEs = moduleEigengenes(expr=keepDat, colors=moduleColors)$eigengenes; 
MEmat = t(MEs); 
colnames(MEmat) = colnames(datExpr); 
kMEtable = signedKME(datExpr=keepDat, datME=t(MEmat), corFnc="bicor");
keptSymbols = genData[match(rownames(t(keepDat)), genData[,2]), 4];
geneTable   = cbind(substr(rownames(kMEtable),1,15), keptSymbols, moduleColors, t(geneCorr), kMEtable);
moduleTraitCor = cor(MEs, traitMat, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
