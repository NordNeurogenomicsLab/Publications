# PanelDoC v 1.2 
# Authors: Alex S. Nord, Alex M. Mawla
# Copyright 2011-2015


###################################################################
###################################################################
###                                                             ###
###               Main function for analysis                    ###
###                                                             ###
###################################################################
###################################################################

run.analysis <- function(criteria) {

###################################################################
### 1: Load input data for experiment                           ###
###################################################################

### Load libraries
	suppressPackageStartupMessages(library("genomeIntervals", logical.return=FALSE, warn.conflicts=FALSE))
	suppressPackageStartupMessages(library("Biostrings", logical.return=FALSE, warn.conflicts=FALSE))
	# not currently implemented in analysis
#	suppressPackageStartupMessages(library("ADaCGH2", logical.return=FALSE, warn.conflicts=FALSE)) 
	suppressPackageStartupMessages(library("mclust", logical.return=FALSE, warn.conflicts=FALSE))

### Read in experiment files
	setwd(criteria$experiment.directory)
	sample.information <- read.csv(criteria$sample.filename, header=T)
	known.cnvs <- read.csv(criteria$knowncnvs.filename, header=T, as.is=T)
	setwd(criteria$platform.directory)
	bait.information <- read.table(criteria$bait.filename, sep="\t", header=F, fill=T)
	partition.information <- read.csv(criteria$partition.filename, header=, as.is=T)
	setwd(criteria$annotation.directory)
	gene.information <- read.table(criteria$gene.filename, header=F, as.is=T)
	experiment.metadata <- list(
						sample.information=sample.information, 
						bait.information=bait.information, 
						partition.information=partition.information,
						gene.information=gene.information,
						known.cnvs=known.cnvs
						)

	bait.coverage <- sum(size(bait.union(bait.information))) 
	rm(sample.information, known.cnvs, bait.information, partition.information, gene.information) 
	
### Print run information
	print(paste("Experiment name: ", criteria$experiment.name, sep=""))
	print(paste("Number of samples in batch: ", nrow(experiment.metadata$sample.information), sep=""))
	print(paste("Number of partitions for analysis: ", nrow(experiment.metadata$partition.information), sep=""))
	print(paste("Total targeted bases to analyze: ", bait.coverage, sep=""))
								
### Generate directory structure if not already present
	setwd(criteria$experiment.directory)
	if (sum(dir()=="bedgraph")==0) {dir.create("bedgraph")}
	if (sum(dir()=="calls")==0) {dir.create("calls")}
	if (sum(dir()=="raw")==0) {dir.create("raw")}
	if (sum(dir()=="normalized")==0) {dir.create("normalized")}
	if (sum(dir()=="PDFs")==0) {dir.create("PDFs")}
	if (sum(dir()=="QC_Metrics")==0) {dir.create("QC_Metrics")}


###################################################################
### 2: If no genomic data for bed file, generate genomic data   ###
###################################################################

	if (criteria$prep.partition==TRUE) {
		print("Generating platform data")
		partition.prep(criteria, experiment.metadata)
	}
	
###################################################################
### 3: Load median data and genomic data for targeted area      ###
###################################################################

	setwd(criteria$platform.directory)
### read in genomic data for platform
	targeted.bases <- read.csv("Platform_baseData.csv", header=T, as.is=T)
### if regions correspond to genes, read in data for gene plotting
	if (criteria$gene.graphics==TRUE) {exon.data <- load.exon.data(experiment.metadata)}
### identify last targeted autosomal base for normalization (sex chromosomes not used for curve, but will be normalized)
	last.aut <- max(which(targeted.bases[,"ChrID"] != "ChrX" & targeted.bases[,"ChrID"] != "ChrY" & targeted.bases[,"ChrID"] != "ChrM"))
### Generate (if necessary) and read in median coverage data for targeted bases
	if (criteria$generate.median == "TRUE") {generate.experiment.median(targeted.bases, experiment.metadata, criteria, exon.data)}
	median.coverage <- load.median.data(criteria, experiment.metadata)
	median.coverage <- merge(targeted.bases, median.coverage)
		
###################################################################
### 4: Analyze each sample                                      ###
###################################################################

### intialize output data structures
	cnv.table <- vector("list", length=nrow(experiment.metadata$sample.information))
	qc.table <- matrix(ncol=14, nrow=nrow(experiment.metadata$sample.information))
	colnames(qc.table) <- c("Sample", "Coverage", "NormalizedCoverage", "Ratio", "SDCoverage", "SDNormalizedCoverage", "SDRatio", "Sample_SN_Raw", "Sample_SN_Normalized", "CNV_Zscore", "CNV_Count", "Percent10X", "Percent50X", "Percent100X")

### run for loop covering all samples to analyze
	for (sample.index in 1:nrow(experiment.metadata$sample.information)) {
		sample.id <- experiment.metadata$sample.information[sample.index,"Sample"]
		print(paste("Running analysis for ", sample.id, ": ", sample.index, " of ", nrow(experiment.metadata$sample.information), " Time: ", date(), sep=""))
### 4a: load sample data ###
		sample.coverage <- load.coverage.data(sample.id, criteria, experiment.metadata)
		sample.coverage <- merge(median.coverage, sample.coverage, all.x=T, by.x=c("ChrID", "Position"), by.y=c("ChrID", "Position"))
		sample.coverage[,"Coverage"] <- ifelse(is.na(sample.coverage[,"Coverage"]),0,sample.coverage[,"Coverage"])			
		sample.coverage <- order.chr(sample.coverage,1,2)
		
### 4b: remove sample known CNV positions from normalization curve (will still be normalized) ###
### only necessary for CNVs that could affect normalization due to proportion of total bases within CNV ###
		cnv.regions <- subset(experiment.metadata$known.cnvs, experiment.metadata$known.cnvs[,"SampleID"]==sample.id)
		cnv.rows <- vector(length=0)
		if (nrow(cnv.regions)>0) {
			for (cnv.calls in 1:nrow(cnv.regions)) {
				chr <- cnv.regions[cnv.calls,"Chr"]
				start.base <- as.numeric(cnv.regions[cnv.calls,"Start"])
				end.base <- as.numeric(cnv.regions[cnv.calls,"End"])
				cnv.rows <- c(cnv.rows, which(sample.coverage[,"ChrID"]==chr & as.numeric(sample.coverage[,"Position"])>=start.base & as.numeric(sample.coverage[,"Position"])<=end.base))
			}
		}


### 4c: normalize vs. median using invariant set methods ###
		normalized.coverage <- invSetNormalize(as.numeric(sample.coverage[,"Coverage"]), as.numeric(sample.coverage[,"MedianCoverage"]), last.aut, cnv.rows)
	
### 4d: normalize vs. gc content bias ###		
		gc.normalized.coverage <- run.coverage.correction(normalized.coverage, sample.coverage[,"GC100"], as.numeric(sample.coverage[,"BaitCoverage"]), as.numeric(sample.coverage[,"SelfChain"]), sample.coverage[,"MedianCoverage"], sample.coverage[,"SDCoverage"], as.numeric(last.aut), criteria$minimum.coverage, criteria$minimum.bait, criteria$minimum.zscore, max(sample.coverage[,"MedianCoverage"]), cnv.rows)
		normalized.coverage <-  normalized.coverage * gc.normalized.coverage[[3]]
		sample.coverage[,"MedianCoverage"] <- sample.coverage[,"MedianCoverage"] * gc.normalized.coverage[[3]]
		sd.diff.coverage <- sd(na.omit(normalized.coverage-as.numeric(sample.coverage[,"MedianCoverage"])))

### 4e: generate ratio relative to median coverage ###
		ratio.normalized <- normalized.coverage/sample.coverage[,"MedianCoverage"]
		sample.coverage <- cbind(sample.coverage, NormalizedCoverage=round(normalized.coverage,4), RatioNormalized=round(ratio.normalized,4))
		
### 4f: generate final coverage data and output figures ###
		setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
		#HERE 
		pdf(file=paste(sample.id, "HistogramPlot.pdf", sep="_"), width=11, height=8)
		par(mfrow = c(1, 3))
		hist(as.numeric(sample.coverage[,"Coverage"]), breaks=1000, xlim=c(0,1000), main=sample.id, xlab="Coverage")
		hist(as.numeric(sample.coverage[,"NormalizedCoverage"]), breaks=1000, xlim=c(0,1000), main=sample.id, xlab="NormalizedCoverage")
		hist(as.numeric(sample.coverage[,"RatioNormalized"]), breaks=10000, main=sample.id, xlab="RatioNormalized", xlim=c(0,2))
		dev.off()
		#HERE
		generate.ratio.plot(sample.id, sample.coverage, experiment.metadata, criteria)
		write.normalized.output(sample.id, sample.coverage, criteria, experiment.metadata)
		
### 4g: call cnvs ###
### Note: at this point normalized ratio data has been generated and a different segmentation algorithm may be used ###
		if (criteria$call.cnvs==TRUE) {
			if (criteria$sample.type == "Tumor") {
				cnv.output <- run.state.calling(sample.id, sample.coverage, criteria, experiment.metadata, exon.data)
			} else {
				cnv.output <- run.variant.calling(sample.id, sample.coverage, criteria, experiment.metadata, exon.data)
			}
			if (criteria$annotate.cnvs == TRUE) {
				cnv.output <- annotate.cnvs(cnv.output, experiment.metadata)
			}
			cnv.table[[sample.index]] <- cnv.output
			cnv.table.final <- mat.build(cnv.table)
			setwd(criteria$experiment.directory)
			write.table(cnv.table.final, paste(criteria$experiment.name, "_CNVs.csv", sep=""), sep=",", col.names=T, row.names=F, quote=F)
		}
		
### 4h: generate QC output ###
		sample.raw.sn <- median(sample.coverage[,"Coverage"], na.rm=T)/sd(sample.coverage[,"Coverage"], na.rm=T)
		sample.normalized.sn <- median(sample.coverage[,"NormalizedCoverage"], na.rm=T)/sd(sample.coverage[,"NormalizedCoverage"], na.rm=T)
		zscore.cnv <- .5/sd(sample.coverage[,"RatioNormalized"], na.rm=T)
		if (criteria$call.cnvs==TRUE) {		
			cnv.count <- nrow(cnv.table[[sample.index]])
			} else {
			cnv.count <- NULL
		}
		if (is.null(cnv.count)) {cnv.count <- 0}
		qc.table[sample.index,] <- c(
			as.character(sample.id), 
			median(sample.coverage[,"Coverage"], na.rm=T), 
			median(sample.coverage[,"NormalizedCoverage"], na.rm=T), 
			median(sample.coverage[,"RatioNormalized"], na.rm=T), 
			sd(sample.coverage[,"Coverage"], na.rm=T), 
			sd(sample.coverage[,"NormalizedCoverage"], na.rm=T), 
			sd(sample.coverage[,"RatioNormalized"], na.rm=T),
			sample.raw.sn,
			sample.normalized.sn,
			zscore.cnv,
			cnv.count,
			length(which(sample.coverage[,"Coverage"]>=10))/nrow(sample.coverage),
			length(which(sample.coverage[,"Coverage"]>=50))/nrow(sample.coverage),
			length(which(sample.coverage[,"Coverage"]>=100))/nrow(sample.coverage)			
			)
		qc.table.final <- merge(experiment.metadata$sample.information, qc.table)
		setwd(criteria$experiment.directory)
		write.table(qc.table.final, paste(criteria$experiment.name, "_QCTable.csv", sep=""), sep=",", col.names=T, row.names=F, quote=F)

	}

	# Per-gene, cross-sample ratio plots
	#setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
	color.scheme <- c("blue4", "brown2", "chartreuse4", "mediumorchid3", "coral1",
"darkcyan", "darkolivegreen1", "darkorange", "darkturquoise",
"darkslategrey", "indianred4", "goldenrod1", "lawngreen",
"lightsteelblue4", "magenta", "midnightblue", "chocolate2",
"olivedrab1", "orangered", "deepskyblue", "green", "deeppink3")
	sample.names <- experiment.metadata$sample.information[,1]
	if (length(sample.names)>length(color.scheme)) {color.scheme <- rep(color.scheme, (length(sample.names)%/%length(color.scheme)+1))}
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		partition.name <- experiment.metadata$partition.information[partition.index,1]
                png(file=paste("PDFs/",partition.name, "_RatioPlot.png", sep=""), width=11,height=8, units="in",res=150)
                sample.index <- 1
		sample.id <- sample.names[sample.index]
		ratio.norm <- read.table(paste("normalized","/",sample.id,"_",partition.name,".ratio", sep=""), header=T, colClasses=c("character","numeric","numeric"))
		plot(1:nrow(ratio.norm), ratio.norm[,3], type="p", cex=.2, pch=16, col=color.scheme[sample.index], ylim=c(0,2.2), xlab=partition.name, xaxt="n", ylab="Ratio", main="Relative Depth of Coverage Ratio")
                abline(2,0)
                abline(1.75,0)
                abline(1.5,0)
                abline(1.25,0)
                abline(1,0)
                abline(.75,0)
                abline(.5,0)
                abline(.25,0)
                abline(0,0)
                for (sample.index in 2:length(sample.names)) {
			sample.id <- sample.names[sample.index]
			ratio.norm <- read.table(paste("normalized","/",sample.id,"_",partition.name,".ratio", sep=""), header=T, colClasses=c("character","numeric","numeric"))
			points(1:nrow(ratio.norm), ratio.norm[,3], type="p", cex=.2, pch=16, col=color.scheme[sample.index])
                }
		legend("topleft",legend=sample.names,col=color.scheme, horiz=T, cex=(.2*(22/length(sample.names))), pch=15, pt.cex=1, bty="n")
                dev.off()
	}

	print(paste(criteria$experiment.name, " analysis complete", sep=""))
}


###################################################################
###################################################################
###################################################################
###                                                             ###
###                        Sub-functions                        ###
###                                                             ###
### 1.  Platform preparation                                    ###
### 2.  Coverage prep and normalization                         ###
### 3.  CNV calling                                             ###
### 4.  General                                                 ###
###                                                             ###
###################################################################
###################################################################
###################################################################



###################################################################
###################################################################
###                  Platform prep functions                    ###
###################################################################
###################################################################

##############################################################################################################

###                                                                                                        ###

### partition.prep generates a set of genomic data for the targeted regions                                ###

### this funciton reads in data downloaded from the UCSC ftp site                                          ###

### and generates a series of data for each targeted base in the .bed file                                 ###

### data generated:                                                                                        ###

### 	GC-content for varying window sizes around targeted base (25,50,100,250,500)                       ###

###		distance to nearest non-targeted base                                                              ###

###		number of self-chains overlapping targeted base                                                    ###

### this data is merged during analysis and used for CNV calling and normalization                         ###

### this function outputs a table named "Platform_basedata.csv" that is read during analysis               ###

### all other functions in this section are called by partition.prep                                       ###

###                                                                                                        ###

##############################################################################################################


partition.prep <- function(criteria, experiment.metadata) {
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
#		print(date())
		partition.name <- experiment.metadata$partition.information[partition.index,1]
		partition.chr <- experiment.metadata$partition.information[partition.index,2]
		print(paste("Generating GC-percent and bait coverage for ", partition.name, sep=""))
		setwd(criteria$genome.directory)
		partition.seq <- readDNAStringSet(paste(partition.chr, ".fa", sep=""), format="fasta", use.names=FALSE)[[1]]
	
		partition.bait <- subset(experiment.metadata$bait.information, experiment.metadata$bait.information[,1]==as.character(partition.chr) 
			& experiment.metadata$bait.information[,2]>=experiment.metadata$partition.information[partition.index,3]
			& experiment.metadata$bait.information[,3]<=experiment.metadata$partition.information[partition.index,4]
			)		
		partition.bait.intervals <- partition.bait[,2:3]
		partition.bait.intervals <- Intervals(
			partition.bait.intervals,
			closed = c(TRUE, FALSE)
			) 
		targeted.intervals <- bait.union(partition.bait)
		targeted.intervals.complement <- interval_complement(
			Intervals(
 				targeted.intervals[,1:2],
				closed = c(TRUE, FALSE)
				)
			)	
		partition.output <- interval_gc_bait_data(targeted.intervals, targeted.intervals.complement, partition.seq, partition.bait.intervals, partition.name, partition.chr)
		setwd(criteria$platform.directory)
		print(paste("length posititions: ", length(partition.output[[1]]), sep=""))
		print(paste("length GC100: ", length(partition.output[[4]]), sep=""))
		print(paste("length bait coverage: ", length(partition.output[[7]]), sep=""))
		print(paste("length bait coverage > 0: ", length(which(partition.output[[7]]>0)), sep=""))
		write(partition.output[[1]], paste(partition.name, ".positions", sep=""), ncolumns=1)
		write(partition.output[[2]], paste(partition.name, ".gcpercent500", sep=""), ncolumns=1)
		write(partition.output[[3]], paste(partition.name, ".gcpercent250", sep=""), ncolumns=1)
		write(partition.output[[4]], paste(partition.name, ".gcpercent100", sep=""), ncolumns=1)
		write(partition.output[[5]], paste(partition.name, ".gcpercent50", sep=""), ncolumns=1)
		write(partition.output[[6]], paste(partition.name, ".gcpercent25", sep=""), ncolumns=1)
		write(partition.output[[7]], paste(partition.name, ".baitcoverage", sep=""), ncolumns=1)
		#Test
		if (partition.index==1) {
			overlap.positions <- partition.output[[1]]
		} else {
		overlap.positions <- c(overlap.positions,partition.output[[1]]) 
		}
	}
	positions <- load.positions(criteria, experiment.metadata)
	gc.content <- load.gcpercent(criteria, experiment.metadata)
	bait.coverage <- load.baitcoverage(criteria, experiment.metadata)
	targeted.bases <- data.frame(
		ChrID=positions[[1]], 
		Position=positions[[2]], 
		GC100=gc.content, 
		BaitCoverage=bait.coverage	
		)	
	targeted.bases[,1] <- as.character(targeted.bases[,1])
	targeted.bases[,2] <- as.numeric(as.character(targeted.bases[,2]))
	targeted.bases[,3] <- as.numeric(as.character(targeted.bases[,3]))
	targeted.bases[,4] <- as.numeric(as.character(targeted.bases[,4]))
	targeted.bases <- unique(targeted.bases)
	targeted.bases <- generate.selfchain(criteria, experiment.metadata, targeted.bases)
	targeted.bases <- subset(targeted.bases, as.numeric(as.character(targeted.bases[,"BaitCoverage"]))>0)
	rm(positions, gc.content, bait.coverage)
	setwd(criteria$platform.directory)
	write.table(targeted.bases, "Platform_baseData.csv", sep=",", col.names=T, row.names=F, quote=F)
}			


##############################################################################################################

###                                                                                                        ###

### generate.selfchain generates data for the number of self-chains mapping to targeted bases              ###
###                                                                                                        ###

##############################################################################################################



generate.selfchain <- function(criteria, experiment.metadata, targeted.bases) {
	setwd(criteria$annotation.directory)
	selfchain.filename <- "chainSelfLink.txt"
	self.chain <- read.table(selfchain.filename, header=F, as.is=T)
	self.output <- vector(length=nrow(targeted.bases))
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		partition.name <- experiment.metadata$partition.information[partition.index,1]
		partition.chr <- experiment.metadata$partition.information[partition.index,2]
		self.chain.chr <- subset(self.chain, self.chain[,2]==partition.chr)
		setwd(criteria$platform.directory)
		partition.positions <- which(targeted.bases[,"ChrID"]==partition.chr & as.numeric(targeted.bases[,"Position"])>=experiment.metadata$partition.information[partition.index,3] & as.numeric(targeted.bases[,"Position"])<=experiment.metadata$partition.information[partition.index,4])
		positions <- as.numeric(targeted.bases[partition.positions,"Position"])
		self.output.partition <- vector(length=length(positions))
		position.count.100k <- length(positions)%/%100000
		if (length(positions)%%100000>0) {position.count.100k <- position.count.100k + 1}
		for (index in 1:position.count.100k) {
			start.position <- (100000 * (index - 1)) + 1
			end.position <- 100000 * index
			if (index == position.count.100k) {end.position <- (100000 * (index-1)) + length(positions)%%100000}
			position.temp <- positions[start.position:end.position]
			self.temp <- subset(self.chain.chr, as.numeric(self.chain.chr[,4])>=min(position.temp) & as.numeric(self.chain.chr[,3])<=max(position.temp))
			self.output.partition[start.position:end.position] <- interval_selfchain_data(partition.name, self.temp, position.temp)
		}
		write(self.output.partition, paste(partition.name, ".selfdata", sep=""), ncolumns=1)
		self.output[partition.positions] <- self.output.partition
		rm(positions, self.output.partition)
	}
	return(data.frame(targeted.bases, SelfChain=self.output))
}	

##############################################################################################################

###                                                                                                        ###

### interval_selfchain.data is a subfunction that takes an interval                                        ###

### and produces self-chain count for each position                                                        ###
###                                                                                                        ###

##############################################################################################################



interval_selfchain_data <- function(partition.name, self.temp, position.temp) {
	if (nrow(self.temp)>0) {
		selfchain.intervals <- Intervals(
			self.temp[,3:4],
			closed = c(TRUE, FALSE)
			) 
		selfchain.final <- vector(length=0)
		selfchain.final <- interval_overlap(position.temp, selfchain.intervals)
		selfchain.final <- as.numeric(as.character(lapply(selfchain.final,length)))
		selfchain.final <- selfchain.final + 1
		return(selfchain.final)
		} else {
		return(rep("1", length(position.temp)))
	}	
}	
	
##############################################################################################################

###                                                                                                        ###

### interval_gc_bait_data is a subfunction that takes an interval                                          ### 

### and produces GC and distance-from for each position                                                    ###
### for now, this analyzes 225, 50, 100, 250, 500 bp windows, but                                          ### 

### this can be changed if desired in the function                                                         ###

###                                                                                                        ###

##############################################################################################################



interval_gc_bait_data <- function(buffered.intervals, targeted.intervals.complement, partition.seq, partition.bait.intervals, partition.name, partition.chr) {
	positions.final <- vector(length=0)
	gcpercent.final.500 <- vector(length=0)
	gcpercent.final.250 <- vector(length=0)
	gcpercent.final.100 <- vector(length=0)
	gcpercent.final.50 <- vector(length=0)
	gcpercent.final.25 <- vector(length=0)
	baitcoverage.final <- vector(length=0)
	interval.regions <- vector("list", length=length(buffered.intervals[[1]])) 
	# for each interval, generate list with 3 elements: position, gc-percent, bait coverage
	for (interval.index in 1:length(buffered.intervals[[1]])) {
		start <- buffered.intervals[interval.index,1]
		end <- buffered.intervals[interval.index,2]
		positions <- seq(start, end, 1)
		bait.coverage <- distance_to_nearest(positions, targeted.intervals.complement)
		min.position <- start-249
		max.position <- end+250	
		if (interval.index==1) {
			start <- start-1
		
		}
		gc.percent.500 <- letterFrequencyInSlidingView(
			subseq(partition.seq,start=min.position,end=max.position), 
			view.width=500, 
			letters=c("GC"), 
			OR="|"
			)/500
		min.position <- start-124
		max.position <- end+125	
		gc.percent.250 <- letterFrequencyInSlidingView(
			subseq(partition.seq,start=min.position,end=max.position), 
			view.width=250, 
			letters=c("GC"), 
			OR="|"
			)/250
		min.position <- start-49
		max.position <- end+50
		gc.percent.100 <- letterFrequencyInSlidingView(
			subseq(partition.seq,start=min.position,end=max.position), 
			view.width=100, 
			letters=c("GC"), 
			OR="|"
			)/100
		min.position <- start-24
		max.position <- end+25
		gc.percent.50 <- letterFrequencyInSlidingView(
			subseq(partition.seq,start=min.position,end=max.position), 
			view.width=50, 
			letters=c("GC"), 
			OR="|"
			)/50		
		min.position <- start-12
		max.position <- end+12
		gc.percent.25 <- letterFrequencyInSlidingView(
			subseq(partition.seq,start=min.position,end=max.position), 
			view.width=25, 
			letters=c("GC"), 
			OR="|"
			)/25	
	
		positions.final <- c(positions.final, positions)
		gcpercent.final.500 <- c(gcpercent.final.500, gc.percent.500)
		gcpercent.final.250 <- c(gcpercent.final.250, gc.percent.250)
		gcpercent.final.100 <- c(gcpercent.final.100, gc.percent.100)
		gcpercent.final.50 <- c(gcpercent.final.50, gc.percent.50)
		gcpercent.final.25 <- c(gcpercent.final.25, gc.percent.25)
		baitcoverage.final <- c(baitcoverage.final, bait.coverage)
		}

	return(list(position=positions.final, gcPercent500=gcpercent.final.500, gcPercent250=gcpercent.final.250, gcPercent100=gcpercent.final.100, gcPercent50=gcpercent.final.50, gcPercent25=gcpercent.final.25, baitCoverage=baitcoverage.final))
	}


##############################################################################################################

###                                                                                                        ###

### load position data generated by partition analysis                                                     ###
###                                                                                                        ###

##############################################################################################################



load.positions <- function(criteria, experiment.metadata) {
	setwd(criteria$platform.directory)
	positions <- vector(length=0)
	chr <- vector(length=0)
	for (partition.index in 1:length(experiment.metadata$partition.information$PartitionName)) {
		partition.name <- experiment.metadata$partition.information[partition.index, "PartitionName"]
		partition.chr <- experiment.metadata$partition.information[partition.index, "PartitionChr"]
		chr.position <- read.table(paste(partition.name, ".positions", sep=""), header=F)
		positions <- c(positions, chr.position[,1])
		chr.call <- rep(as.character(partition.chr), nrow(chr.position))
		chr <- c(chr, chr.call)
	}
	return(list(chr=chr, position=positions))
}



##############################################################################################################

###                                                                                                        ###

### load GC data generated by partition analysis                                                           ###

### for now, this reads in the 100bp window data, but this can be changed if desired                       ###

###                                                                                                        ###

##############################################################################################################





load.gcpercent <- function(criteria, experiment.metadata) {
	setwd(criteria$platform.directory)
	gc.percent <- vector(length=0)
	positions <- vector(length=0)
	chr <- vector(length=0)
	for (partition.index in 1:length(experiment.metadata$partition.information$PartitionName)) {
		partition.name <- experiment.metadata$partition.information[partition.index, "PartitionName"]
		partition.gc <- read.table(paste(partition.name, ".gcpercent100", sep=""), header=T)
		partition.name <- experiment.metadata$partition.information[partition.index, "PartitionName"]
		partition.chr <- experiment.metadata$partition.information[partition.index, "PartitionChr"]
		chr.position <- read.table(paste(partition.name, ".positions", sep=""), header=F)
		positions <- c(positions, chr.position[,1])
		chr.call <- rep(as.character(partition.chr), nrow(chr.position))
		chr <- c(chr, chr.call)

		gc.percent <- c(gc.percent, partition.gc[,1])
	}
	return(gc.percent)
}

##############################################################################################################

###                                                                                                        ###

### load self-chain data generated by partition analysis                                                   ###
###                                                                                                        ###

##############################################################################################################



load.selfchain <- function(criteria, experiment.metadata) {
	setwd(criteria$platform.directory)
	self.chain <- matrix(nrow=0, ncol=1)
	for (partition.index in 1:length(experiment.metadata$partition.information$PartitionName)) {
		partition.name <- experiment.metadata$partition.information[partition.index, "PartitionName"]
		partition.sc <- read.table(paste(partition.name, ".selfdata", sep=""), header=F, skip=1, as.is=T)
		self.chain <- rbind(self.chain, partition.sc)
		}
	return(self.chain)
}
	
	
##############################################################################################################

###                                                                                                        ###

### load distance to non-targeted-base data generated by partition analysis                                ###

###                                                                                                        ###

##############################################################################################################



load.baitcoverage <- function(criteria, experiment.metadata) {
	setwd(criteria$platform.directory)
	bait.coverage <- matrix(nrow=0, ncol=1)
	for (partition.index in 1:length(experiment.metadata$partition.information$PartitionName)) {
		partition.name <- experiment.metadata$partition.information[partition.index, "PartitionName"]
		partition.bc <- read.table(paste(partition.name, ".baitcoverage", sep=""), header=F)
		bait.coverage <- rbind(bait.coverage, partition.bc)
	}
	return(bait.coverage[,1])
}


##############################################################################################################

###                                                                                                        ###

### function to generate overlap count                                                                     ###
###                                                                                                        ###

##############################################################################################################



bait.overlap.positions <- function(temp.coverage, bait.intervals) { 
	overlap <- interval_overlap(temp.coverage, bait.intervals)
	return(sapply(overlap, length))
	}


##############################################################################################################

###                                                                                                        ###

### function to write bedgraphs for the GC percent data                                                    ###
###                                                                                                        ###

##############################################################################################################



write.gcpercent.bedgraph <- function(basic.data, criteria, experiment.metadata, partition.id) {
	bedgraph.main <- basic.data[,c(1,2,2,3)]
	bedgraph.main[,3] <- as.numeric(bedgraph.main[,3]) + 1
	bedgraph.main <- order.chr(bedgraph.main,1,2)
	bedgraph.main[,2] <- format(bedgraph.main[,2], scientific=F)
	bedgraph.main[,3] <- format(bedgraph.main[,3], scientific=F)
	bedgraph.main[,4] <- format(bedgraph.main[,4], scientific=F, digits=4)
	coverage.header <- make.bedgraph.header(paste("Experiment", "GC100", sep="_"), paste("Experiment", "GC100", sep="_"))
	coverage.bedgraph.filename <- paste("Experiment", "GC100.bedgraph", sep="_")
	cat(coverage.header, file=coverage.bedgraph.filename)
	cat("\n", file=coverage.bedgraph.filename, append=T)
	write.table(bedgraph.main[,1:4], file=coverage.bedgraph.filename, col.names=F, row.names=F, quote=F, append=T)

	bedgraph.main <- basic.data[,c(1,2,2,4)]	
	bedgraph.main[,3] <- as.numeric(bedgraph.main[,3]) + 1
	begraph.main <- bedgraph.main[order(as.numeric(bedgraph.main[,2])),]
	corrected.header <- make.bedgraph.header(paste("Experiment", "BaitCoverage", sep="_"), paste("Experiment", "BaitCoverage", sep="_"))
	corrected.bedgraph.filename <- paste("Experiment", "BaitCoverage.bedgraph", sep="_")
	cat(corrected.header, file=corrected.bedgraph.filename)
	cat("\n", file=corrected.bedgraph.filename, append=T)
	write.table(bedgraph.main[,c(1:4)], file=corrected.bedgraph.filename, col.names=F, row.names=F, quote=F, append=T)
	
	### Write corrected coverage by chromsome (partition to be implemented)
	setwd(paste(criteria$platform.directory, sep="/"))
	write.table(basic.data, paste("Platform", "_", ".CorrectionData", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
}	


###################################################################
###################################################################
###        Coverage prep and normalization Functions            ###
###################################################################
###################################################################

##############################################################################################################

###                                                                                                        ###

### function to load exon data for gene graphics                                                           ###
###                                                                                                        ###

##############################################################################################################



load.exon.data <- function(experiment.metadata)	{
	exon.data <- vector("list", length=nrow(experiment.metadata$partition.information))
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		partition.name <- experiment.metadata$partition.information[partition.index,"PartitionName"]
		partition.refseq <- subset(experiment.metadata$gene.information, experiment.metadata$gene.information[,2]==experiment.metadata$partition.information[partition.index,"PartitionRefseq"])
		if (nrow(partition.refseq) > 1) {
			# Get the one on from the chromosome listed in the partition file (i.e., probably not hapmap population chrs)
			partition.refseq <- subset(partition.refseq, partition.refseq[,3]==experiment.metadata$partition.information[partition.index,2])
		}
		if (nrow(partition.refseq) > 1) {
			# If there's still more than one row, there's more than one set of coordinates for the same refseq ID - identical paralogs?
			# Get the one that matches the coordinates we used
			overlap <- apply(partition.refseq[,5:6],1, function(x,y,z) length(intersect(seq(x[1],x[2]),seq(y,z))),y=as.numeric(experiment.metadata$partition.information[partition.index,3]),z=as.numeric(experiment.metadata$partition.information[partition.index,4]))
			partition.refseq <- partition.refseq[which(overlap == max(overlap))[1],]
		}
		partition.start <- as.numeric(partition.refseq[5]) 
		partition.end <- as.numeric(partition.refseq[6]) 
### Generate coding and UTR start stop table
		partition.exon <- data.frame(start=as.numeric(strsplit(partition.refseq[,10], split=",")[[1]]), end=as.numeric(strsplit(partition.refseq[,11], split=",")[[1]]))	
		coding.start.exon <- which(partition.exon[,"start"]<=as.numeric(partition.refseq[7]) & partition.exon[,"end"]>=as.numeric(partition.refseq[7]))
		coding.end.exon <- which(partition.exon[,"start"]<=as.numeric(partition.refseq[8]) & partition.exon[,"end"]>=as.numeric(partition.refseq[8]))
		partition.untranslated <- partition.exon[c(1:coding.start.exon,coding.end.exon:nrow(partition.exon)),]
		utr.coding.start.exon <- which(partition.untranslated[,"start"]<=as.numeric(partition.refseq[7]) & partition.untranslated[,"end"]>=as.numeric(partition.refseq[7]))
		utr.coding.end.exon <- which(partition.untranslated[,"start"]<=as.numeric(partition.refseq[8]) & partition.untranslated[,"end"]>=as.numeric(partition.refseq[8]))
		partition.untranslated[utr.coding.start.exon,"end"] <- as.numeric(partition.refseq[7]) - 1
		partition.untranslated[utr.coding.end.exon,"start"] <- as.numeric(partition.refseq[8]) + 1
		partition.untranslated <- cbind(partition.untranslated, class=rep("UTR", nrow(partition.untranslated)))
		partition.exon <- partition.exon[coding.start.exon:coding.end.exon,]
		partition.exon[1,"start"] <- as.numeric(partition.refseq[7])
		partition.exon[nrow(partition.exon),"end"] <- as.numeric(partition.refseq[8])
		partition.exon <- cbind(partition.exon, class=rep("coding", nrow(partition.exon)))
		partition.table <- rbind(partition.untranslated, partition.exon)
		partition.table <- partition.table[order(as.numeric(partition.table[,1])),]
		exon.data[[partition.index]] <- partition.table
	}
	return(exon.data)
}


##############################################################################################################

###                                                                                                        ###

### assign region for base position - used for generating ratio plots                                      ###
###                                                                                                        ###

##############################################################################################################



assign.region <- function(sample, partition.data) {
	GeneRegion <- vector(length=nrow(sample))
	for (partition.index in 1:nrow(partition.data)) {
		partition.chr <- partition.data[partition.index,"PartitionChr"]
		partition.start <- partition.data[partition.index,"PartitionStart"]
		partition.end <- partition.data[partition.index,"PartitionEnd"]
		part.rows <- which(sample[,"ChrID"]==partition.chr & sample[,"Position"]>=partition.start & sample[,"Position"]<=partition.end)
		GeneRegion[part.rows] <- partition.data[partition.index,"PartitionName"]
		}
	return(GeneRegion)
	}	

##############################################################################################################

###                                                                                                        ###

### generate.ratio.plot generates .png format plots of ratio and coverage data                             ###

### colors can be changed as necessary, will recycle through if more regions than colors                   ###

### all the graphical output is controlled in this funciton for the ratio plot                             ###

###                                                                                                        ###

##############################################################################################################



generate.ratio.plot <- function(sample.id, sample.coverage, experiment.metadata, criteria) {
	gene.region <- assign.region(sample.coverage, experiment.metadata$partition.information)
	gene.names <- experiment.metadata$partition.information[,"PartitionName"]
	color.scheme <- c("blue4", "brown2", "chartreuse4", "mediumorchid3", "coral1",
"darkcyan", "darkolivegreen1", "darkorange", "darkturquoise",
"darkslategrey", "indianred4", "goldenrod1", "lawngreen",
"lightsteelblue4", "magenta", "midnightblue", "chocolate2",
"olivedrab1", "orangered", "deepskyblue", "green", "deeppink3")
	if (length(gene.names)>length(color.scheme)) {color.scheme <- rep(color.scheme, (length(gene.names)%/%length(color.scheme)+1))}
	setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
	png(file=paste(sample.id, "RatioPlot.png", sep="_"), width=11, height=8, units="in", res=150)
	par(mfrow=c(2,1))
	plot(0,0, col="white", xlim=c(0,nrow(sample.coverage)), ylim=c(0,1200), xlab="", xaxt="n", ylab="Coverage", main="Raw Coverage")
	abline(1000,0)
	abline(900,0)	
	abline(800,0)
	abline(700,0)
	abline(600,0)
	abline(500,0)	
	abline(400,0)
	abline(300,0)
	abline(200,0)	
	abline(100,0)
	abline(0,0)
	position <- 1
	for (gene.index in 1:length(gene.names)) {
		sample.rows <- which(gene.region==gene.names[gene.index])
		points(position:(position+length(sample.rows)-1), sample.coverage[sample.rows, "Coverage"], type="h", cex=.1, pch=16, col=color.scheme[gene.index], lwd=.1)
		position <- position + length(sample.rows)
		}
	legend("topleft", legend=gene.names, col=color.scheme, horiz=T, cex=(.45*(22/length(gene.names))), pch=15, pt.cex=1, bty="n", bg="white")
	plot(0,0, col="white", xlim=c(0,nrow(sample.coverage)), ylim=c(0,2.2), xlab="", xaxt="n", ylab="Ratio", main="Relative Depth of Coverage Ratio")
	abline(2,0)
	abline(1.75,0)	
	abline(1.5,0)
	abline(1.25,0)
	abline(1,0)
	abline(.75,0)	
	abline(.5,0)
	abline(.25,0)
	abline(0,0)	
	position <- 1
	for (gene.index in 1:length(gene.names)) {
		sample.rows <- which(gene.region==gene.names[gene.index])
		points(position:(position+length(sample.rows)-1), sample.coverage[sample.rows, "RatioNormalized"], type="p", cex=.1, pch=16, col=color.scheme[gene.index])
		position <- position + length(sample.rows)
		}
	legend("topleft", legend=gene.names, col=color.scheme, horiz=T, cex=(.45*(22/length(gene.names))), pch=15, pt.cex=1, bty="n", bg="white")
	dev.off()
}


##############################################################################################################

###                                                                                                        ###

### load.median.data reads in the median files for an experiment                                           ###
###                                                                                                        ###

##############################################################################################################



load.median.data <- function(criteria, experiment.metadata) {
	setwd(paste(criteria$experiment.directory, "raw", sep="/"))
	coverage <- matrix(nrow=0, ncol=3)
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		partition.name <- experiment.metadata$partition.information[partition.index,1]
		partition.chr <- experiment.metadata$partition.information[partition.index,2]
		partition.coverage <- read.table(paste("Median_", partition.name, ".raw", sep=""), header=T, colClasses=c("character", "numeric", "numeric", "numeric", "numeric"))
		coverage<- rbind(coverage, partition.coverage)
		}
	return(coverage)
}


##############################################################################################################

###                                                                                                        ###

### generate.experiment.median reads in raw data for all the samples from                                  ### 

### an experiment and generates median and sd data for the experiment                                      ###

### this function also plots the median coverage data, which can be changed                                ### 

### in the same manner as the generate.ratio.plot function                                                 ###
###                                                                                                        ###

##############################################################################################################



generate.experiment.median <- function(targeted.bases, experiment.metadata, criteria, exon.data) {
# run for loop covering all samples to analyze
	full.data <- vector("list", length=nrow(experiment.metadata$partition.information))
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		partition.name <- experiment.metadata$partition.information[partition.index,1]
		partition.chr <- experiment.metadata$partition.information[partition.index,2]
		partition.start <- experiment.metadata$partition.information[partition.index,3]
		partition.end <- experiment.metadata$partition.information[partition.index,4]
		
		## Fix here 
		final.coverage <- subset(targeted.bases, targeted.bases[,"ChrID"]==partition.chr & targeted.bases[,"Position"]>=partition.start & targeted.bases[,"Position"]<=partition.end)
		
		print(paste("Running analysis for ", partition.name, ": ", partition.index, " of ", nrow(experiment.metadata$partition.information), sep=""))

		for (sample.index in 1:nrow(experiment.metadata$sample.information)) {
			sample.id <- experiment.metadata$sample.information[sample.index,"Sample"]
### load sample data ###
			setwd(paste(criteria$experiment.directory, "coverage", sep="/"))
			sample.coverage <- read.table(paste(sample.id, "_", partition.name, ".Depth", sep=""), header=T, colClasses=c("character", "integer", "integer"))
			colnames(sample.coverage)[2] <- "Position"
			colnames(sample.coverage)[3] <- paste(sample.id, "_Coverage", sep="")
			sample.coverage <- subset(sample.coverage, sample.coverage[,paste(sample.id, "_Coverage", sep="")]!=0)
			sample.coverage <- sample.coverage[,-1]
			final.coverage <- merge(final.coverage, sample.coverage, by.x="Position", by.y="Position", all.x=T)
			final.coverage[,3] <- ifelse(is.na(final.coverage[,3]),0,final.coverage[,3])
			final.coverage[,paste(sample.id, "_Coverage", sep="")] <- ifelse(is.na(final.coverage[,paste(sample.id, "_Coverage", sep="")]),0,final.coverage[,paste(sample.id, "_Coverage", sep="")])	
			setwd(paste(criteria$experiment.directory, "raw", sep="/"))
			write.table(final.coverage[,c("ChrID", "Position", paste(sample.id, "_Coverage", sep=""))], paste(sample.id, "_", partition.name, ".raw", sep=""), col.names=T, row.names=F, quote=F) 
		}	
		median.coverage <- as.numeric(apply(final.coverage[,3:ncol(final.coverage)],1,median))
		sd.coverage <- as.numeric(apply(final.coverage[,3:ncol(final.coverage)],1,sd))
		sn.coverage <- median.coverage/sd.coverage
		median.coverage <- cbind(final.coverage[,1:2], MedianCoverage=median.coverage, SDCoverage=sd.coverage, SNCoverage=sn.coverage)
		full.data[[partition.index]] <- median.coverage
		setwd(paste(criteria$experiment.directory, "raw", sep="/"))
		write.table(median.coverage[,c(2,1,3:5)], paste("Median", "_", partition.name, ".raw", sep=""), col.names=T, row.names=F)
		setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
		pdf(file=paste(partition.name, "HistogramPlot.pdf", sep="_"), width=11, height=8)
		par(mfrow = c(1, 3))
		hist(median.coverage[,"MedianCoverage"], breaks=100, xlim=c(0,1000), main=partition.name, xlab="Median Coverage")
		hist(median.coverage[,"SDCoverage"], breaks=100, xlim=c(0,500), main=partition.name, xlab="SD Coverage")
		hist(median.coverage[,"SNCoverage"], breaks=100, main=partition.name, xlab="S:N", xlim=c(0,5))
		dev.off()
		png(file=paste(partition.name, "Coverage.png", sep="_"), width=11, height=8, units="in", res=150)
# load in gene data for plotting

	### Load data for region
# plot gene if gene.graphics = TRUE #
		if (criteria$gene.graphics == TRUE) {
			partition.exons <- exon.data[[partition.index]]
			par(mfrow=c(4,1))
			start.track <- 2
			track.height <- 2
			end.track <- start.track - track.height
			plot(median.coverage[c(1,nrow(median.coverage)), "Position"],c(0,2), col="white", xaxt="n", yaxt="n", main=partition.name, xlab="", ylab="", frame.plot=F)
			draw.gene(partition.exons,1,2,start.track,end.track,"lightblue","Y", "black", "black")
		} else {
			par(mfrow=c(3,1))
		}
# plot coverage #
		plot(median.coverage[c(1,nrow(median.coverage)), "Position"], median.coverage[c(1,nrow(median.coverage)), "MedianCoverage"], col="white", ylim=c(0,1000), xlab="Position", ylab="Coverage", main="Median Coverage")
		abline(1000,0, col="grey")
		abline(900,0, col="grey")	
		abline(800,0, col="grey")
		abline(700,0, col="grey")
		abline(600,0, col="grey")
		abline(500,0, col="grey")	
		abline(400,0, col="grey")
		abline(300,0, col="grey")
		abline(200,0, col="grey")	
		abline(100,0, col="grey")
		abline(0,0, col="grey")
		points(median.coverage[,"Position"], median.coverage[,"MedianCoverage"], col="black", pch=16, cex=.2)
# plot sd coverage #
		plot(median.coverage[c(1,nrow(median.coverage)), "Position"], median.coverage[c(1,nrow(median.coverage)), "SDCoverage"], col="white", ylim=c(0,500), xlab="Position", ylab="SD Coverage", main="SD Coverage")
		abline(500,0, col="grey")	
		abline(400,0, col="grey")
		abline(300,0, col="grey")
		abline(200,0, col="grey")	
		abline(100,0, col="grey")
		abline(0,0, col="grey")
		points(median.coverage[,"Position"], median.coverage[,"SDCoverage"], col="black", pch=16, cex=.2)
# plot sn coverage #
		plot(median.coverage[c(1,nrow(median.coverage)), "Position"], median.coverage[c(1,nrow(median.coverage)), "SNCoverage"], col="white", ylim=c(0,5), xlab="Position", ylab="SN Coverage", main="S:N Coverage")
		abline(5,0, col="grey")	
		abline(4,0, col="grey")
		abline(3,0, col="grey")
		abline(2,0, col="grey")	
		abline(1,0, col="grey")
		abline(0,0, col="grey")
		points(median.coverage[,"Position"], median.coverage[,"SNCoverage"], col="black", pch=16, cex=.2)
		dev.off()
	}
	median.coverage <- mat.build(full.data)
	setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
	pdf(file=paste("Experiment", "HistogramPlot.pdf", sep="_"), width=11, height=8)
	par(mfrow = c(1, 3))
	hist(as.numeric(median.coverage[,"MedianCoverage"]), breaks=1000, xlim=c(0,1000), main=criteria$experiment.name, xlab="Coverage")
	hist(as.numeric(median.coverage[,"SDCoverage"]), breaks=1000, xlim=c(0,500), main=criteria$experiment.name, xlab="SD Coverage")
	hist(as.numeric(median.coverage[,"SNCoverage"]), breaks=1000, main=criteria$experiment.name, xlab="S:N", xlim=c(0,10))
	dev.off()
	gene.region <- assign.region(median.coverage, experiment.metadata$partition.information)
	gene.names <- experiment.metadata$partition.information[,"PartitionName"]
	color.scheme <- c("blue4", "brown2", "chartreuse4", "mediumorchid3", "coral1",
"darkcyan", "darkolivegreen1", "darkorange", "darkturquoise",
"darkslategrey", "indianred4", "goldenrod1", "lawngreen",
"lightsteelblue4", "magenta", "midnightblue", "chocolate2",
"olivedrab1", "orangered", "deepskyblue", "green", "deeppink3"  )
	if (length(gene.names)>length(color.scheme)) {color.scheme <- rep(color.scheme, (length(gene.names)%/%length(color.scheme)+1))}
	setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
	png(file=paste("Experiment", "RatioPlot.png", sep="_"), width=11, height=8, units="in", res=150)
	par(mfrow=c(3,1))
	plot(0,0, col="white", xlim=c(0,nrow(median.coverage)), ylim=c(0,1200), xlab="", xaxt="n", ylab="Coverage", main="Median")
	abline(1000,0)
	abline(900,0)	
	abline(800,0)
	abline(700,0)
	abline(600,0)
	abline(500,0)	
	abline(400,0)
	abline(300,0)
	abline(200,0)	
	abline(100,0)
	abline(0,0)
	position <- 1:nrow(median.coverage)
	for (gene.index in 1:length(gene.names)) {
		sample.rows <- which(gene.region==gene.names[gene.index])
		points(position[sample.rows], median.coverage[sample.rows, "MedianCoverage"], type="h", cex=.1, pch=16, col=color.scheme[gene.index], lwd=.1)
		}
	legend("topleft", legend=experiment.metadata$partition.information[,"PartitionName"], col=color.scheme, horiz=T, cex=(.45*(22/length(gene.names))), pch=15, pt.cex=1, bty="n", bg="white")
	plot(0,0, col="white", xlim=c(0,nrow(median.coverage)), ylim=c(0,600), xlab="", xaxt="n", ylab="Coverage", main="SD")
	abline(500,0)	
	abline(400,0)
	abline(300,0)
	abline(200,0)	
	abline(100,0)
	abline(0,0)
	position <- 1:nrow(median.coverage)
	for (gene.index in 1:length(gene.names)) {
		sample.rows <- which(gene.region==gene.names[gene.index])
		points(position[sample.rows], median.coverage[sample.rows, "SDCoverage"], type="h", cex=.1, pch=16, col=color.scheme[gene.index], lwd=.1)
		}
	legend("topleft", legend=experiment.metadata$partition.information[,"PartitionName"], col=color.scheme, horiz=T, cex=(.45*(22/length(gene.names))), pch=15, pt.cex=1, bty="n", bg="white")
	plot(0,0, col="white", xlim=c(0,nrow(median.coverage)), ylim=c(0,6), xlab="", xaxt="n", ylab="S:N", main="S:N")
	abline(5,0)	
	abline(4,0)
	abline(3,0)
	abline(2,0)	
	abline(1,0)
	abline(0,0)
	position <- 1:nrow(median.coverage)
	for (gene.index in 1:length(gene.names)) {
		sample.rows <- which(gene.region==gene.names[gene.index])
		points(position[sample.rows], median.coverage[sample.rows, "SNCoverage"], type="p", cex=.1, pch=16, col=color.scheme[gene.index])
		}
	legend("topleft", legend=experiment.metadata$partition.information[,"PartitionName"], col=color.scheme, horiz=T, cex=.45, pch=15, pt.cex=1, bty="n", bg="white")
	dev.off()
	### Write genomic bedgraph for coverage and corrected coverage (not sure if this will work for exome)
	setwd(paste(criteria$experiment.directory,  "bedgraph", sep="/"))
	bedgraph.main <- median.coverage[,c("ChrID","Position","Position","MedianCoverage","SDCoverage", "SNCoverage")]
	bedgraph.main[,3] <- as.numeric(bedgraph.main[,3]) + 1
	bedgraph.main <- order.chr(bedgraph.main,1,2)
	bedgraph.main[,2] <- format(bedgraph.main[,2], scientific=F)
	bedgraph.main[,3] <- format(bedgraph.main[,3], scientific=F)
	bedgraph.main[,4] <- format(bedgraph.main[,4], scientific=F, digits=4)
	bedgraph.main[,5] <- format(bedgraph.main[,5], scientific=F, digits=4)
	bedgraph.main[,6] <- format(bedgraph.main[,6], scientific=F, digits=4)
	coverage.header <- make.bedgraph.header(paste(sample.id, "Median ReadDepth", sep="_"), paste(sample.id, "ReadDepth", sep="_"))
	sd.header <- make.bedgraph.header(paste(criteria$experiment.name, "SD Coverage", sep="_"), paste(sample.id, "SD", sep="_"))
	sn.header <- make.bedgraph.header(paste(criteria$experiment.name, "S:N Coverage", sep="_"), paste(sample.id, "SN", sep="_"))
	coverage.bedgraph.filename <- paste(criteria$experiment.name, "ReadDepth.bedgraph", sep="_")
	sd.bedgraph.filename <- paste(criteria$experiment.name, "SD.bedgraph", sep="_")
	sn.bedgraph.filename <- paste(criteria$experiment.name, "SN.bedgraph", sep="_")
	cat(coverage.header, file=coverage.bedgraph.filename)
	cat("\n", file=coverage.bedgraph.filename, append=T)
	write.table(bedgraph.main[,1:4], file=coverage.bedgraph.filename, col.names=F, row.names=F, quote=F, append=T)
	cat(sd.header, file=sd.bedgraph.filename)
	cat("\n", file=sd.bedgraph.filename, append=T)
	write.table(subset(bedgraph.main[,c(1:3,5)], !is.na(bedgraph.main[,5])), file=sd.bedgraph.filename, col.names=F, row.names=F, quote=F, append=T)
	cat(sn.header, file=sn.bedgraph.filename)
	cat("\n", file=sn.bedgraph.filename, append=T)
	write.table(subset(bedgraph.main[,c(1:3,6)], !is.na(bedgraph.main[,6])), file=sn.bedgraph.filename, col.names=F, row.names=F, quote=F, append=T)
}


##############################################################################################################

###                                                                                                        ###

### load.coverage.data loads raw data for a sample                                                         ###
###                                                                                                        ###

##############################################################################################################



load.coverage.data <- function(sample.id, criteria, experiment.metadata) {
	setwd(paste(criteria$experiment.directory, "raw", sep="/"))
	coverage <- matrix(nrow=0, ncol=3)
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		partition.name <- experiment.metadata$partition.information[partition.index,1]
		partition.chr <- experiment.metadata$partition.information[partition.index,2]
		partition.coverage <- read.table(paste(sample.id, "_", partition.name, ".raw", sep=""), header=T, colClasses=c("character", "numeric", "numeric"))
		coverage<- rbind(coverage, partition.coverage)
		}
	colnames(coverage) <- c("ChrID", "Position", "Coverage")
	return(coverage)
}	

##############################################################################################################

###                                                                                                        ###

### write.normalized.output writes output for all the raw and normalized data in multiple formats          ###
###                                                                                                        ###

##############################################################################################################



write.normalized.output <- function(sample.id, sample.coverage, criteria, experiment.metadata) {
### Write genomic bedgraph for coverage and normalized coverage (not sure if this will work for exome)
	if (criteria$output.bedgraph==TRUE) {
		setwd(paste(criteria$experiment.directory,  "bedgraph", sep="/"))
		bedgraph.main <- sample.coverage[,c("ChrID","Position","Position","Coverage","RatioNormalized")]
		bedgraph.main[,"Coverage"] <- ifelse(is.na(bedgraph.main[,"Coverage"]),0,bedgraph.main[,"Coverage"])
		bedgraph.main[,3] <- as.numeric(bedgraph.main[,3]) + 1
		bedgraph.main <- order.chr(bedgraph.main,1,2)
		bedgraph.main[,2] <- format(bedgraph.main[,2], scientific=F)
		bedgraph.main[,3] <- format(bedgraph.main[,3], scientific=F)
		bedgraph.main[,4] <- format(bedgraph.main[,4], scientific=F, digits=4)
		bedgraph.main[,5] <- format(bedgraph.main[,5], scientific=F, digits=4)
		coverage.header <- make.bedgraph.header(paste(sample.id, "ReadDepth", sep="_"), paste(sample.id, "ReadDepth", sep="_"))
		ratio.header <- make.bedgraph.header(paste(sample.id, "Ratio", sep="_"), paste(sample.id, "Ratio", sep="_"))
		coverage.bedgraph.filename <- paste(sample.id, "ReadDepth.bedgraph", sep="_")
		ratio.bedgraph.filename <- paste(sample.id, "Ratio.bedgraph", sep="_")
		cat(coverage.header, file=coverage.bedgraph.filename)
		cat("\n", file=coverage.bedgraph.filename, append=T)
		write.table(bedgraph.main[,1:4], file=coverage.bedgraph.filename, col.names=F, row.names=F, quote=F, append=T)
		cat(ratio.header, file=ratio.bedgraph.filename)
		cat("\n", file=ratio.bedgraph.filename, append=T)
		write.table(subset(bedgraph.main[,c(1:3,5)], !is.na(bedgraph.main[,5])), file=ratio.bedgraph.filename, col.names=F, row.names=F, quote=F, append=T)
	}
	
	### Write corrected coverage by chromsome (partition to be implemented)
	setwd(paste(criteria$experiment.directory,  "normalized", sep="/"))
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		partition.name <- experiment.metadata$partition.information[partition.index,1]
		partition.chr <- experiment.metadata$partition.information[partition.index,2]
		sample.partition <- subset(sample.coverage, sample.coverage[,"ChrID"]!="")
		sample.partition <- subset(sample.partition, as.character(sample.partition[,"ChrID"])==as.character(partition.chr) & sample.partition[,"Position"]>=experiment.metadata$partition.information[partition.index,"PartitionStart"] & sample.partition[,"Position"]<=experiment.metadata$partition.information[partition.index,"PartitionEnd"])
		sample.partition[,"NormalizedCoverage"] <- ifelse(is.na(sample.partition[,"NormalizedCoverage"]),0,sample.partition[,"NormalizedCoverage"])
		if (criteria$output.normalized==TRUE) {
			write.table(sample.partition[,c("ChrID", "Position", "Coverage", "NormalizedCoverage", "RatioNormalized")], paste(sample.id, "_", partition.name, ".normalized", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
		}
		if (criteria$output.ratio==TRUE) {
			write.table(sample.partition[,c("ChrID", "Position", "RatioNormalized")], paste(sample.id, "_", partition.name, ".ratio", sep=""), sep="\t", col.names=T, row.names=F, quote=F)
		}	
	}
}


##############################################################################################################

###                                                                                                        ###

### run.coverage.correction was used for GC and distance-to-untargeted correction                          ###

### this function is not currently implemented, but may be used to optimize normalization for              ###

### a specific platform                                                                                    ###

###                                                                                                        ###

##############################################################################################################



run.coverage.correction <- function(coverage, gc.data, probe.coverage, selfchain, median.coverage, sd.coverage, last.aut, min.correct.coverage, min.correct.bait, min.correct.zscore, max.correct.coverage, cnv.rows){
	count <- length(coverage[1:last.aut])
	zscore <- abs((coverage - median.coverage)/sd.coverage)
	yx <- data.frame(y=as.numeric(coverage), x=as.numeric(gc.data), z=as.numeric(probe.coverage), zscore=zscore, sc=selfchain)
	yx$y <- ifelse(yx$y==0,1,yx$y)
	subsample <- yx
	if (length(cnv.rows)>0) {subsample <- subsample[-cnv.rows,]}
	subrows <- c(1:last.aut)
	if (count > 10000) {
		subsample.step <- count%/%10000
		subrows <- seq(1,last.aut,subsample.step)
		}
	subsample <- subsample[subrows,]
	subsample <- subset(subsample, subsample$y>=min.correct.coverage & subsample$z>=min.correct.bait & subsample$zscore<=min.correct.zscore & subsample$y<=max.correct.coverage & subsample$sc==1)
	subsample$z <- log(subsample$z + 10)
	probe.coverage.converted <- log(probe.coverage + 10)
	subsample$x <- -(subsample$x - .35)^2
	gc.data.converted <- -(gc.data - .35)^2
	g <- lm(log(subsample$y) ~ subsample$x + subsample$z)
	gs <- summary(g)
	b.coef <- gs["coefficients"]
	r.sq <- gs["adj.r.squared"]
	output <- log(yx$y)-(g$coefficients[1]+g$coefficients[2]*gc.data.converted + g$coefficients[3]*probe.coverage.converted)
	output <- exp(output)
	return(list(b.coef, r.sq, output))
}


##############################################################################################################

###                                                                                                        ###

### normalize.invariantset is the normalization curve generation described by Li                           ###

###                                                                                                        ###

##############################################################################################################



normalize.invariantset <- function (data, ref, prd.td = c(0.003, 0.007)) {
    np <- length(data)
    r.ref <- rank(ref)
    r.array <- rank(data)
    prd.td.adj <- prd.td * 10
    i.set <- rep(TRUE, np)
    ns <- sum(i.set)
    ns.old <- ns + 50 + 1
    while ((ns.old - ns) > 50) {
        air <- (r.ref[i.set] + r.array[i.set])/(2 * ns)
        prd <- abs(r.ref[i.set] - r.array[i.set])/ns
        threshold <- (prd.td.adj[2] - prd.td[1]) * air + prd.td.adj[1]
        i.set[i.set] <- (prd < threshold)
        ns.old <- ns
        ns <- sum(i.set)
        if (prd.td.adj[1] > prd.td[1]) 
		prd.td.adj <- prd.td.adj * 0.9
    }
    n.curve <- smooth.spline(ref[i.set], data[i.set], tol = max(1e-6, 1e-6 * IQR(ref[i.set])))
    return(list(n.curve = n.curve, i.set = i.set))
}


##############################################################################################################

###                                                                                                        ###

### invSetNormalize is a function that calls the invariant set normalization to                            ###

### generate the curve and then applies it to normalize the data                                           ###
###                                                                                                        ###

##############################################################################################################



invSetNormalize <- function(rawIntensity,rawRef,lastAuto, cnv.rows) {
# This function gets invariant set and fitting curve from autosomes and do adjust for all probes including x.
# rawIntensity is autosome and x chromosome.
	logIntAuto <- log(rawIntensity[1:lastAuto]+10)
	refAuto <- log(rawRef[1:lastAuto]+10)
	if (length(cnv.rows)>0) {
		logIntAuto <- log(rawIntensity[-cnv.rows])
		refAuto <- log(rawRef[-cnv.rows])	
		}
# correct for 0 intensity/coverage measurements
	logIntAuto <- ifelse(logIntAuto==-Inf,0,logIntAuto)
	refAuto <- ifelse(refAuto==-Inf,0,refAuto)
	normalized <- normalize.invariantset(logIntAuto,refAuto, prd.td=c(0.003,0.007))
	intNorm <- log(rawIntensity+10) + log(rawRef+10) - predict(normalized$n.curve,log(rawRef+10))$y
	data.norm <- exp(intNorm)
}




###################################################################
###################################################################
###                    CNV calling functions                    ###
###################################################################
###################################################################

##############################################################################################################

###                                                                                                        ###

### Function to run cnv calling                                                                            ###
###                                                                                                        ###

##############################################################################################################



run.variant.calling <- function(sample.id, sample.coverage, criteria, experiment.metadata, exon.data) {
	sd.diff.coverage <- sd(as.numeric(sample.coverage[,"Coverage"])-as.numeric(sample.coverage[,"MedianCoverage"]))
	cnv.calls.partition <- vector("list", length=nrow(experiment.metadata$partition.information))
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		cnv.calls.partition[[partition.index]] <- matrix(nrow=0,ncol=16)
		}
	setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
	png(file=paste(sample.id, "RatioPlots.png", sep="_"), width=8, height=2.5*nrow(experiment.metadata$partition.information), units="in", res=150)
	par(mfrow=c(nrow(experiment.metadata$partition.information),1))
	
	#Generate Coverage QC for Original and Corrected Minimum Depth of Coverage Threshold
	partition.qc <- matrix(NA, nrow=nrow(experiment.metadata$partition.information), ncol=7)
		colnames(partition.qc) <- c("Partition Name", "Original Minimum Coverage Threshold", "Corrected Minimum Coverage Threshold", "Mean Coverage Threshold", "Median Coverage Threshold", "Original Percent Coverage Passing Threshold", "Corrected Percent Coverage Passing Threshold")
	partition.qc[,"Original Minimum Coverage Threshold"]=criteria$minimum.coverage
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		
### Run for each partition ###
		partition.name <- experiment.metadata$partition.information[partition.index, "PartitionName"]
#		print(as.character(partition.name))
		partition.chr <- experiment.metadata$partition.information[partition.index, "PartitionChr"]
		ratio.table <- subset(sample.coverage, 
			  sample.coverage[,"ChrID"]==partition.chr & 				  
			  sample.coverage[,"Position"]>=experiment.metadata$partition.information[partition.index,"PartitionStart"] & 
			  sample.coverage[,"Position"]<=experiment.metadata$partition.information[partition.index,"PartitionEnd"]
			  )
		partition.qc[partition.index,"Partition Name"] <- partition.name
		partition.qc[partition.index,"Mean Coverage Threshold"] <- round(mean(ratio.table[,"MedianCoverage"]),2)
		partition.qc[partition.index,"Median Coverage Threshold"] <- round(median(ratio.table[,"MedianCoverage"]),2)
		partition.qc[partition.index,"Original Percent Coverage Passing Threshold"] <- round(100*sum(ratio.table[,"MedianCoverage"]>=criteria$minimum.coverage)/nrow(ratio.table),2)
		
		#If Original Depth of Coverage Threshold not met, new threshold is asigned, coverage adjusted to use fixed percentage of ratio.table 
		while ( (sum(ratio.table[,"MedianCoverage"]>=criteria$minimum.coverage) / nrow(ratio.table) ) < 0.8) {
			criteria$minimum.coverage <- criteria$minimum.coverage-1;
		}
		partition.qc[partition.index,"Corrected Minimum Coverage Threshold"]=criteria$minimum.coverage
		partition.qc[partition.index,"Corrected Percent Coverage Passing Threshold"] <- round(100*sum(ratio.table[,"MedianCoverage"]>=criteria$minimum.coverage)/nrow(ratio.table),2)
		
	# Threshold too high REMOVED 
		ratio.table <- subset(ratio.table, 
		as.numeric(ratio.table[,"MedianCoverage"])>=criteria$minimum.coverage & 
		as.numeric(ratio.table[,"BaitCoverage"])>=criteria$minimum.bait & 
		as.numeric(ratio.table[,"SelfChain"])<=criteria$maximum.selfchain
		)
				
		if (nrow(ratio.table)>(criteria$window.size*2)) {
			test.data <- ratio.table[,"RatioNormalized"] # add in z-score later

## Run calling functions for 0 and 4+ copies ###
			called.probes.proband <- run.seed.call(
			   test.data, 
			   criteria$loss.none.seed.signal, 
			   criteria$gain.many.seed.signal, 
			   criteria$loss.none.extend.signal, 
			   criteria$gain.many.extend.signal,
			   criteria$window.size,
			   criteria$pass.count
			   )
			cnv.regions <- call.cnv.regions(ratio.table[,"Position"], called.probes.proband)
			cnv.regions <- event.merge(cnv.regions,5000,.1)
			cnv.region.details.1 <- get.cnv.details(
				as.character(partition.chr), 
				cnv.regions, 
				test.data, 
				ratio.table, 
				criteria$loss.none.seed.signal, 
				criteria$gain.many.seed.signal,
				sd.diff.coverage
				)
		### Run calling functions for 1 and 3 copies ###
			called.probes.proband <- run.seed.call(
			   test.data, 
			   criteria$loss.seed.signal, 
			   criteria$gain.seed.signal, 
			   criteria$loss.extend.signal, 
			   criteria$gain.extend.signal,
			   criteria$window.size,
			   criteria$pass.count
			   )
			cnv.regions <- call.cnv.regions(ratio.table[,"Position"], called.probes.proband)
			cnv.regions <- event.merge(cnv.regions,5000,.1)
			cnv.region.details.2 <- get.cnv.details(
				as.character(partition.chr), 
				cnv.regions, 
				test.data, 
				ratio.table, 
				criteria$loss.seed.signal, 
				criteria$gain.seed.signal,
				sd.diff.coverage
				)
					
### Combine tables and process calls	
			cnv.region.details <- rbind(cnv.region.details.1, cnv.region.details.2)	
			copy.number <- vector(length=nrow(cnv.region.details))
			if (nrow(cnv.region.details)>0) {
				copy.number <- ifelse(as.numeric(cnv.region.details[,"median.ratio"])<=criteria$loss.none.call.signal,0,1)
				copy.number <- ifelse(as.numeric(cnv.region.details[,"median.ratio"])>=criteria$gain.call.signal,3,copy.number)
				copy.number <- ifelse(as.numeric(cnv.region.details[,"median.ratio"])>=criteria$gain.many.call.signal,4,copy.number)
			}
			cnv.region.details <- cbind(cnv.region.details, CopyNumber=copy.number)
			
			cnv.region.details <- subset(cnv.region.details, as.numeric(cnv.region.details[,"median.ratio"])>=criteria$gain.call.signal | as.numeric(cnv.region.details[,"median.ratio"])<=criteria$loss.call.signal) 
			cnv.region.details <- subset(cnv.region.details, as.numeric(cnv.region.details[,"base.count"])>=criteria$minimum.base.window & as.numeric(cnv.region.details[,"base.count.criteria"])>=criteria$minimum.base.pass)
			
			#if (nrow(cnv.region.details)>0) {
				#setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
				diploid.positions <- rep(2,nrow(ratio.table))
				region.length <- nrow(ratio.table)
				if (nrow(cnv.region.details)>0) {
				for (call.index in 1:nrow(cnv.region.details)) {
					if (as.numeric(cnv.region.details[call.index,"median.ratio"]) < 1) {
						call.rows <- which(ratio.table[,"Position"]>=as.numeric(cnv.region.details[call.index,"start"]) & ratio.table[,"Position"]<=as.numeric(cnv.region.details[call.index,"end"]))
						diploid.positions[call.rows] <- 1
					}
					if (as.numeric(cnv.region.details[call.index,"median.ratio"]) > 1) {
						call.rows <- which(ratio.table[,"Position"]>=as.numeric(cnv.region.details[call.index,"start"]) & ratio.table[,"Position"]<=as.numeric(cnv.region.details[call.index,"end"]))
						diploid.positions[call.rows] <- 3
					}
				}
				}	
				loss.positions <- which(diploid.positions==1)
				gain.positions <- which(diploid.positions==3)
				diploid.positions <- which(diploid.positions==2)
				#HERE
				png(file=paste(sample.id, "_", partition.name, ".png", sep=""), width=11, height=4, units="in", res=1000)
				if (criteria$gene.graphics==TRUE) {
					# Add exons 
					plot(ratio.table[,"Position"], rep(1,region.length), col="white", ylim=c(0,3.6), main=partition.name, xlab="Position", ylab="Ratio")
					start.track <- 3.6
					track.height <- .5
					end.track <- start.track - track.height
					draw.gene(exon.data[[partition.index]],1,2,start.track,end.track,"lightblue","Y", "black", "black")
				} else {
					plot(ratio.table[,"Position"], rep(1,region.length), col="white", ylim=c(0,3), main=partition.name, xlab="Position", ylab="Ratio")
				}	
				fig.lines <- c(0,.5,1,1.5,2,2.5,3)
				x.min <- min(ratio.table[,"Position"])
				x.max <- max(ratio.table[,"Position"])
				s.line <- x.min
				e.line <- x.max
				coord.mat <- matrix(nrow=2, ncol=2)
				for (line.index in 1:length(fig.lines)) {
					coord.mat[1,] <- c(s.line, fig.lines[line.index])
					coord.mat[2,] <- c(e.line, fig.lines[line.index])
					lines(coord.mat,col="dark grey")
					}	
				points(ratio.table[diploid.positions,"Position"], ratio.table[diploid.positions,"RatioNormalized"], cex=.3, pch=16)
				points(ratio.table[loss.positions,"Position"], ratio.table[loss.positions,"RatioNormalized"], cex=.3, col="red", pch=16)
				points(ratio.table[gain.positions,"Position"], ratio.table[gain.positions,"RatioNormalized"], cex=.3, col="blue", pch=16)
				#HERE
				dev.off()
			#}
			

### Output calls ###
			cnv.calls.partition[[partition.index]] <- data.frame(SAMPLE_ID=rep(as.character(sample.id),nrow(cnv.region.details)), RegionName=rep(as.character(partition.name),nrow(cnv.region.details)), cnv.region.details)
		}
	}
	dev.off() # finished with figure
	cnv.calls <- mat.build(cnv.calls.partition)
	cnv.calls <- data.frame(cnv.calls, SampleCalls=rep(nrow(cnv.calls), nrow(cnv.calls)))
	cnv.calls <- unique(cnv.calls)
	##Output metrics ##			
		setwd(paste(criteria$experiment.directory, "QC_Metrics", sep="/"))
		write.table(partition.qc, paste(sample.id, "Partition_QC.csv", sep="_"), sep=",", col.names=T, row.names=F, quote=F)
			setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))

		
	setwd(paste(criteria$experiment.directory, "calls", sep="/"))
	write.table(cnv.calls, paste(sample.id, "_CNVs.csv", sep=""), sep=",", col.names=T, row.names=F, quote=F)
	return(cnv.calls)
}	


##############################################################################################################

###                                                                                                        ###

### Function to run cnv calling for tumor data                                                             ###
###                                                                                                        ###

##############################################################################################################



run.state.calling <- function(sample.id, sample.coverage, criteria, experiment.metadata, exon.data) {
	sd.diff.coverage <- sd(as.numeric(sample.coverage[,"Coverage"])-as.numeric(sample.coverage[,"MedianCoverage"]))
	cnv.calls.partition <- vector("list", length=nrow(experiment.metadata$partition.information))
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		cnv.calls.partition[[partition.index]] <- matrix(nrow=0,ncol=16)
		}
	
	for (partition.index in 1:nrow(experiment.metadata$partition.information)) {
		
### Run for each partition ###
		partition.name <- experiment.metadata$partition.information[partition.index, "PartitionName"]
		partition.chr <- experiment.metadata$partition.information[partition.index, "PartitionChr"]
		ratio.table <- subset(sample.coverage, 
			  sample.coverage[,"ChrID"]==partition.chr & 				  
			  sample.coverage[,"Position"]>=experiment.metadata$partition.information[partition.index,"PartitionStart"] & 
			  sample.coverage[,"Position"]<=experiment.metadata$partition.information[partition.index,"PartitionEnd"]
			  )

		ratio.table <- subset(ratio.table, 
			as.numeric(ratio.table[,"MedianCoverage"])>=criteria$minimum.coverage & 
			as.numeric(ratio.table[,"BaitCoverage"])>=criteria$minimum.bait & 
			as.numeric(ratio.table[,"SelfChain"])<=criteria$maximum.selfchain
			)
				
		if (nrow(ratio.table)>(criteria$window.size*2)) {
			test.data <- ratio.table[,"RatioNormalized"] 
			relative.model <- Mclust(test.data, modelNames="E", G=1:criteria$mclust.states)
			final.classification <- relative.model$classification
			if (median(relative.model$uncertainty) > criteria$max.uncertainty) {final.classification <- rep(1, length(final.classification))}
			state.table <- build.states(final.classification, ratio.table)
			cnv.region.details <- get.state.details(
				as.character(partition.chr), 
				state.table, 
				test.data, 
				ratio.table,
				sd.diff.coverage
				)

			if (nrow(cnv.region.details)>1) {
				setwd(paste(criteria$experiment.directory, "PDFs", sep="/"))
				diploid.positions <- rep(1,nrow(ratio.table))
				for (call.index in 1:nrow(cnv.region.details)) {
					if (as.numeric(cnv.region.details[call.index,"state"]) != 1) {
						call.rows <- which(ratio.table[,"Position"]>=as.numeric(cnv.region.details[call.index,"start"]) & ratio.table[,"Position"]<=as.numeric(cnv.region.details[call.index,"end"]))
						diploid.positions[call.rows] <- 2
					}
				}	
				state2.positions <- which(diploid.positions==2)
				state1.positions <- which(diploid.positions==1)
				pdf(file=paste(sample.id, "_", partition.name, ".pdf", sep=""), width=11, height=4)
				if (criteria$gene.graphics==TRUE) {
					# Add exons 
					plot(ratio.table[,"Position"], rep(1,region.length), col="white", ylim=c(0,3.6), main=partition.name, xlab="Position", ylab="Ratio")
					start.track <- 3.6
					track.height <- .5
					end.track <- start.track - track.height
					draw.gene(exon.data[[partition.index]],1,2,start.track,end.track,"lightblue","Y", "black", "black")
				} else {
					plot(ratio.table[,"Position"], rep(1,region.length), col="white", ylim=c(0,3), main=partition.name, xlab="Position", ylab="Ratio")
				}	
				fig.lines <- c(0,.5,1,1.5,2,2.5,3)
				x.min <- min(ratio.table[,"Position"])
				x.max <- max(ratio.table[,"Position"])
				s.line <- x.min
				e.line <- x.max
				coord.mat <- matrix(nrow=2, ncol=2)
				for (line.index in 1:length(fig.lines)) {
					coord.mat[1,] <- c(s.line, fig.lines[line.index])
					coord.mat[2,] <- c(e.line, fig.lines[line.index])
					lines(coord.mat,col="dark grey")
					}	
				points(ratio.table[state1.positions,"Position"], ratio.table[state1.positions,"RatioNormalized"], cex=.3, pch=16)
				points(ratio.table[state2.positions,"Position"], ratio.table[state2.positions,"RatioNormalized"], cex=.3, col="red", pch=16)
				dev.off()
				}
		
### Output calls ###
			cnv.calls.partition[[partition.index]] <- data.frame(SAMPLE_ID=rep(as.character(sample.id),nrow(cnv.region.details)), RegionName=rep(as.character(partition.name),nrow(cnv.region.details)), cnv.region.details)
		}
	}
	cnv.calls <- mat.build(cnv.calls.partition)
	cnv.calls <- data.frame(cnv.calls, SampleCalls=rep(nrow(cnv.calls), nrow(cnv.calls)))
	cnv.calls <- unique(cnv.calls)
	setwd(paste(criteria$experiment.directory, "calls", sep="/"))
	write.table(cnv.calls, paste(sample.id, "_CNVs.csv", sep=""), sep=",", col.names=T, row.names=F, quote=F)
	return(cnv.calls)
}	


##############################################################################################################

###                                                                                                        ###

### build.states takes as input the state call generated by mclust and segments into CNVs                  ###
###                                                                                                        ###

##############################################################################################################



build.states <- function(final.classification, ratio.table) {
	if (length(unique(final.classification))>1) {
		end.probe <- 0
		run.count <- 0
		run.assignment <- vector(length=length(final.classification))
		while (end.probe < length(final.classification)) {
			run.count <- run.count + 1
			start.probe <- end.probe + 1
			start.state <- final.classification[start.probe]
			first.different.probe <- min(which(final.classification[start.probe:length(final.classification)]!=start.state)) + start.probe - 1
			end.probe <- first.different.probe - 1
			if (end.probe=="Inf") {end.probe <- length(final.classification)}
			run.assignment[start.probe:end.probe] <- run.count
		}
		final.run.count <- run.count
		run.start <- vector(length=final.run.count)
		run.end <-  vector(length=final.run.count)
		run.count <-  vector(length=final.run.count)
		run.state <-  vector(length=final.run.count)
		run.mean <-  vector(length=final.run.count)
		for (run.index in 1:final.run.count) {
			run.start[run.index] <- min(which(run.assignment==run.index))
			run.end[run.index] <- max(which(run.assignment==run.index))
			run.count[run.index] <- run.end[run.index]-run.start[run.index]+1
			run.state[run.index] <- final.classification[run.start[run.index]]
		}
		run.table <- data.frame(run.start, run.end, run.count, run.state)	
		run.table <- subset(run.table, run.table[,"run.count"]>=200)
		run.index <- 1
		while (run.index < nrow(run.table)) {
			next.row <- run.index + 1
			if (run.table[next.row,"run.state"] == run.table[run.index,"run.state"]) {
				# merge with run.index + 2
				run.table[run.index,"run.end"] <- run.table[next.row,"run.end"]
				run.table[run.index,"run.count"] <- run.table[run.index,"run.end"]-run.table[run.index,"run.start"]+1
				run.table <- run.table[-next.row,]
			} else {
				run.index <- run.index + 1
			}
		}	
		colnames(run.table) <- c("Start", "End", "BaseCount", "State")	
		return(run.table)
	} else  {
		run.table <- matrix(nrow=1, ncol=4)
		run.table[1,] <- c(1, nrow(ratio.table), nrow(ratio.table), "1")
		colnames(run.table) <- c("Start", "End", "BaseCount", "State")	
		return(run.table)
	}
}


##############################################################################################################

###                                                                                                        ###

### get.state.details generates details for each CNV call                                                  ###

###                                                                                                        ###

##############################################################################################################


get.state.details <- function(chr.id, cnv.regions, test.data, ratio.table, sd.diff.coverage) {
### initialize variables ###
	cnv.regions.temp <- cnv.regions
	if (class(cnv.regions)=="character") {cnv.regions.temp <- t(as.matrix(cnv.regions.temp))}
	output.matrix <- matrix(nrow=nrow(cnv.regions.temp), ncol=19)
### get data for each CNV ###
	if (nrow(cnv.regions.temp)>0) {
		for (cnv.index in 1:nrow(cnv.regions.temp)) {
			start <- as.numeric(cnv.regions.temp[cnv.index,1])
			end <- as.numeric(cnv.regions.temp[cnv.index,2])
			start <- ratio.table[start,"Position"]
			end <- ratio.table[end, "Position"]
			size <- end - start + 1 
			class <- cnv.regions.temp[cnv.index,4]
			cnv.probes <- which(as.numeric(ratio.table[,"Position"])>=start & as.numeric(ratio.table[,"Position"])<=end)
			base.count <-length(cnv.probes)
			gc.content <- median(as.numeric(ratio.table[cnv.probes,"GC100"]))
			mean.base.spacing <- as.numeric(size)/base.count

			median.selfchain <- median(as.numeric(ratio.table[cnv.probes,"SelfChain"])) 	
			percent.unique.sequence <- length(which(as.numeric(ratio.table[cnv.probes,"SelfChain"])==1))/base.count

			median.probe.ratio.normalized <- median(as.numeric(ratio.table[cnv.probes,"RatioNormalized"]))
			sd.probe.ratio.normalized <- sd(as.numeric(ratio.table[cnv.probes,"RatioNormalized"]))
			sn.diff.coverage <- median(as.numeric(ratio.table[cnv.probes,"NormalizedCoverage"])-as.numeric(ratio.table[cnv.probes,"MedianCoverage"]))/sd.diff.coverage

			# change to mean when ready
	
			median.coverage <- median(as.numeric(ratio.table[cnv.probes,"Coverage"]))
			median.coverage.normalized <- median(as.numeric(ratio.table[cnv.probes,"NormalizedCoverage"]))
			median.experiment.coverage <- median(as.numeric(ratio.table[cnv.probes,"MedianCoverage"]))
			median.experiment.sd <- median(as.numeric(ratio.table[cnv.probes,"SDCoverage"]))
			median.experiment.sn <- median(as.numeric(ratio.table[cnv.probes,"SNCoverage"]))
			median.experiment.relativesn <- median(as.numeric(ratio.table[cnv.probes,"SNCoverage"]))/median(as.numeric(ratio.table[,"SNCoverage"]))

			output.row <- c(
				chr.id, 
				start, 
				end, 
				class, 
				size, 
				base.count, 
				gc.content, 
				median.selfchain,
				percent.unique.sequence,
				mean.base.spacing, 
				median.probe.ratio.normalized, 
				sd.probe.ratio.normalized, 
				sn.diff.coverage,
				median.coverage, 
				median.coverage.normalized, 
				median.experiment.coverage,
				median.experiment.sd,
				median.experiment.sn,
				median.experiment.relativesn
				)
			output.matrix[cnv.index,] <- output.row
		}
	}
	colnames(output.matrix) <- c(
		"chr.id", "start", "end", "state", "size", "base.count", "gc.content", "median.selfchain.count", "percent.unique.sequence", "mean.base.spacing",  
		"median.ratio", "sd.ratio", "sn.diff.coverage", "median.coverage.raw", "median.coverage.normalized", "median.experiment.coverage", "median.experiment.sd", "median.experiment.sn", "median.experiment.relativesn")
	return(output.matrix)
}



##############################################################################################################

###                                                                                                        ###

### run.seed.call executes sliding window test for                                                         ###
### probes that pass criteria                                                                              ###
### criteria currently set by del.crit and dup.crit                                                        ###
### function tests 50 probe window                                                                         ###
### for call to be made, 45 probes must pass criteria                                                      ###
### function also cleans up edges                                                                          ###
###                                                                                                        ###

##############################################################################################################



run.seed.call <- function(chr.data, del.crit, dup.crit, extend.del.crit, extend.dup.crit, window.size, pass.count) {
#	print("identifying CNV seeds")
	half.window <- round(window.size/2,0)
#################
### Deletions ###
#################
# initialize variables
	pass.window <- ifelse(chr.data <= del.crit,1,0)
	pass.extend.window <- ifelse(chr.data <= extend.del.crit,1,0)
	call.window <- rep(0,length(chr.data))
	probe.c <- length(chr.data)
# start window scanning	
	for (probe.index in half.window:(probe.c-half.window)) {if (sum(pass.window[(probe.index-(half.window-1)):(probe.index+half.window)])>=pass.count) {call.window[(probe.index-(half.window-1)):(probe.index+half.window)] <- 1}}
# fix probes that do not pass cutoff    
	call.window <- ifelse(call.window == 1 & pass.window != 1,0,call.window)
	call.window.temp <- call.window
# Catch extended regions
	for (probe.index in half.window:(probe.c-half.window)) {
		if (call.window.temp[probe.index]==1) {
			call.window[(probe.index-(half.window-1)):(probe.index+half.window)] <- ifelse(pass.extend.window[(probe.index-(half.window-1)):(probe.index+half.window)]==1,1,call.window[(probe.index-(half.window-1)):(probe.index+half.window)])
			}
		}
# summarize deletion calls
	chr.probe.loss <- rep(0,probe.c)
	chr.probe.loss <- ifelse(call.window==1,1,chr.probe.loss)
	chr.criteria.loss <- rep(0, probe.c)
	chr.criteria.loss <- ifelse(chr.probe.loss==del.crit,del.crit,chr.criteria.loss)
	chr.probe.loss <- ifelse(chr.probe.loss>0,1,0)
####################
### Duplications ###
####################
# initialize variables
	probe.c <- length(chr.data)
	pass.window <- ifelse(chr.data >= dup.crit,1,0)
	pass.extend.window <- ifelse(chr.data >= extend.dup.crit,1,0)
	call.window <- rep(0,length(chr.data))
	# start window scanning	
	for (probe.index in half.window:(probe.c-half.window)) {if (sum(pass.window[(probe.index-(half.window-1)):(probe.index+half.window)])>=pass.count) {call.window[(probe.index-(half.window-1)):(probe.index+half.window)] <- 1}}
# fix probes that do not pass cutoff    
	call.window <- ifelse(call.window == 1 & pass.window != 1,0,call.window)
	call.window.temp <- call.window
# Catch extended regions
	for (probe.index in half.window:(probe.c-half.window)) {
		if (call.window.temp[probe.index]==1) {
			call.window[(probe.index-(half.window-1)):(probe.index+half.window)] <- ifelse(pass.extend.window[(probe.index-(half.window-1)):(probe.index+half.window)]==1,1,call.window[(probe.index-(half.window-1)):(probe.index+half.window)])
			}
		}
# summarize duplication calls
	chr.probe.gain <- rep(0,probe.c)
	chr.probe.gain <- ifelse(call.window==1,1,chr.probe.gain)
	chr.criteria.gain <- rep(0, probe.c)
	chr.criteria.gain <- ifelse(chr.probe.gain==dup.crit,dup.crit,chr.criteria.gain)
	chr.probe.gain <- ifelse(chr.probe.gain>0,1,0)
#####################	
### return output ###
#####################
	return(list(chr.probe.gain, chr.criteria.gain, chr.probe.loss, chr.criteria.loss))
}


##############################################################################################################

###                                                                                                        ###
### get.cnv.details generates information for CNV                                                          ###
### call including signal distribution data                                                                ###
###                                                                                                        ###

##############################################################################################################



get.cnv.details <- function(chr.id, cnv.regions, test.data, ratio.table, del.crit, dup.crit, sd.diff.coverage) {
#	print("getting cnv details")
### initialize variables ###
	cnv.regions.temp <- cnv.regions
	if (class(cnv.regions)=="character") {cnv.regions.temp <- t(as.matrix(cnv.regions.temp))}
	output.matrix <- matrix(nrow=nrow(cnv.regions.temp), ncol=21)
### get data for each CNV ###
	if (nrow(cnv.regions.temp)>0) {
		for (cnv.index in 1:nrow(cnv.regions.temp)) {
			start <- as.numeric(cnv.regions.temp[cnv.index,1])
			end <- as.numeric(cnv.regions.temp[cnv.index,2])
			size <- end - start + 1 
			class <- cnv.regions.temp[cnv.index,3]
			merge.count <- cnv.regions.temp[cnv.index,4]
			cnv.probes <- which(as.numeric(ratio.table[,"Position"])>=start & as.numeric(ratio.table[,"Position"])<=end)
			base.count <-length(cnv.probes)
			gc.content <- median(as.numeric(ratio.table[cnv.probes,"GC100"]))
			mean.base.spacing <- as.numeric(size)/base.count

			median.selfchain <- median(as.numeric(ratio.table[cnv.probes,"SelfChain"])) 	
			percent.unique.sequence <- length(which(as.numeric(ratio.table[cnv.probes,"SelfChain"])==1))/base.count

			median.probe.ratio.normalized <- median(as.numeric(ratio.table[cnv.probes,"RatioNormalized"]))
			sd.probe.ratio.normalized <- sd(as.numeric(ratio.table[cnv.probes,"RatioNormalized"]))
			sn.diff.coverage <- median(as.numeric(ratio.table[cnv.probes,"NormalizedCoverage"])-as.numeric(ratio.table[cnv.probes,"MedianCoverage"]))/sd.diff.coverage

			# change to mean when ready
	
			median.coverage <- median(as.numeric(ratio.table[cnv.probes,"Coverage"]))
			median.coverage.normalized <- median(as.numeric(ratio.table[cnv.probes,"NormalizedCoverage"]))
			median.experiment.coverage <- median(as.numeric(ratio.table[cnv.probes,"MedianCoverage"]))
			median.experiment.sd <- median(as.numeric(ratio.table[cnv.probes,"SDCoverage"]))
			median.experiment.sn <- median(as.numeric(ratio.table[cnv.probes,"SNCoverage"]))
			median.experiment.relativesn <- median(as.numeric(ratio.table[cnv.probes,"SNCoverage"]))/median(as.numeric(ratio.table[,"SNCoverage"]))


			if (class == "GAIN") {
				base.criteria <- length(which(as.numeric(test.data[cnv.probes])>=dup.crit))
			} else {
				base.criteria <- length(which(as.numeric(test.data[cnv.probes])<=del.crit))
			}
			output.row <- c(
				chr.id, 
				start, 
				end, 
				class, 
				size, 
				base.count, 
				gc.content, 
				median.selfchain,
				percent.unique.sequence,
				mean.base.spacing, 
				base.criteria, 
				merge.count, 
				median.probe.ratio.normalized, 
				sd.probe.ratio.normalized, 
				sn.diff.coverage,
				median.coverage, 
				median.coverage.normalized, 
				median.experiment.coverage,
				median.experiment.sd,
				median.experiment.sn,
				median.experiment.relativesn
				)
#			print(output.row)
			output.matrix[cnv.index,] <- output.row
		}
	}
	colnames(output.matrix) <- c(
		"chr.id", "start", "end", "class", "size", "base.count", "gc.content", "median.selfchain.count", "percent.unique.sequence", "mean.base.spacing", "base.count.criteria", "merge.count", 
		"median.ratio", "sd.ratio", "sn.diff.coverage", "median.coverage.raw", "median.coverage.normalized", "median.experiment.coverage", "median.experiment.sd", "median.experiment.sn", "median.experiment.relativesn")
	return(output.matrix)
}


##############################################################################################################

###                                                                                                        ###

### call.cnv.regions takes the sliding window output and joins called bases into CNVs                      ###

###                                                                                                        ###

##############################################################################################################



call.cnv.regions <- function(positions, called.probes) {
#	print("extending cnv regions")
# first call duplications
	gain.probe.calls <- called.probes[[1]]
	if (length(summary(as.factor(gain.probe.calls)))>1) {
		result.matrix <- matrix(ncol=2, nrow=10000)
		start.position <- 0
		end.position <- 0
		cnv.index <- 0
		next.start.probe <- 1
		test.condition <- 0
		while (test.condition==0) {
			cnv.index <- cnv.index + 1
			start.probe <- min(which(gain.probe.calls==1 & as.numeric(positions)>as.numeric(end.position))) 
			temp.calls <- gain.probe.calls[start.probe:length(gain.probe.calls)]
			temp.positions <- positions[start.probe:length(gain.probe.calls)]
			if (length(which(temp.calls==0 & as.numeric(temp.positions)>as.numeric(end.position)))>0) {
				next.start.probe <- min(which(temp.calls==0 & as.numeric(temp.positions)>as.numeric(end.position))) + start.probe - 1
				end.position <- positions[(next.start.probe-1)]
				start.position <- positions[start.probe]
				probe.count <- (next.start.probe-1) - start.probe + 1
				if (probe.count > 1) {
					result.matrix[cnv.index,] <-c(start.position, end.position)
				}
#				print(result.matrix[cnv.index,])
				if ((length(gain.probe.calls) - next.start.probe >=0)==FALSE | (sum(gain.probe.calls[next.start.probe:length(gain.probe.calls)])>0)==FALSE) {
					test.condition <- 1
				}
			} else {
				end.position <- max(as.numeric(positions))
				start.position <- positions[start.probe]
				probe.count <- length(positions) - start.probe + 1
# test for parent overlap
				if (probe.count > 1) {result.matrix[cnv.index,] <-c(start.position, end.position)}
				test.condition <- 1
			}
		}
		gain.result.matrix <- subset(result.matrix, !is.na(result.matrix[,1])) 
		if (nrow(gain.result.matrix)>1) {gain.result.matrix <- extend.cnvs(gain.result.matrix, positions)}
		gain.result.matrix <- cbind(gain.result.matrix, rep("GAIN", nrow(gain.result.matrix)))
	} else {		
		if (sum(gain.probe.calls)==length(gain.probe.calls)) {
			gain.result.matrix <- cbind(min(positions),max(positions), "GAIN")
		} else {
			gain.result.matrix <- matrix(nrow=0, ncol=3)
		}
	}
# then call deletions
	loss.probe.calls <- called.probes[[3]]
	if (length(summary(as.factor(loss.probe.calls)))>1) {
# first test for called probes
		result.matrix <- matrix(ncol=2, nrow=10000)
		start.position <- 0
		end.position <- 0
		cnv.index <- 0
		next.start.probe <- 1
		test.condition <- 0
		while (test.condition==0) {
			cnv.index <- cnv.index + 1
			start.probe <- min(which(loss.probe.calls==1 & as.numeric(positions)>as.numeric(end.position))) 
			temp.calls <- loss.probe.calls[start.probe:length(loss.probe.calls)]
			temp.positions <- positions[start.probe:length(loss.probe.calls)]
			if (length(which(temp.calls==0 & as.numeric(temp.positions)>as.numeric(end.position)))>0) {
				next.start.probe <- min(which(temp.calls==0 & as.numeric(temp.positions)>as.numeric(end.position))) + start.probe - 1
				end.position <- positions[(next.start.probe-1)]
				start.position <- positions[start.probe]
				probe.count <- (next.start.probe-1) - start.probe + 1
				if (probe.count > 1) {
					result.matrix[cnv.index,] <-c(start.position, end.position)
				}
				
#				print(paste(p1.call.perc, p2.call.perc, probe.count))
# print(result.matrix[cnv.index,])
				if ((length(loss.probe.calls) - next.start.probe >=0)==FALSE | (sum(loss.probe.calls[next.start.probe:length(loss.probe.calls)])>0)==FALSE) {
					test.condition <- 1
				}
			} else {
				end.position <- max(as.numeric(positions))
				start.position <- positions[start.probe]
				probe.count <- length(positions) - start.probe + 1
				if (probe.count > 1) {result.matrix[cnv.index,] <-c(start.position, end.position)}
				test.condition <- 1
			}
			
		}
		loss.result.matrix <- subset(result.matrix, !is.na(result.matrix[,1])) 
		if (nrow(loss.result.matrix)>1) {loss.result.matrix <- extend.cnvs(loss.result.matrix, positions)}
		loss.result.matrix <- cbind(loss.result.matrix, rep("LOSS", nrow(loss.result.matrix)))
	} else {
		if (sum(loss.probe.calls)==length(loss.probe.calls)) {
			loss.result.matrix <- cbind(min(positions),max(positions), "LOSS")
		} else {
			loss.result.matrix <- matrix(nrow=0, ncol=3)
		}
	}
# then return data
	result.matrix <- rbind(gain.result.matrix, loss.result.matrix)
	colnames(result.matrix) <- c("start", "end", "class")
	return(result.matrix)
}




##############################################################################################################

###                                                                                                        ###

### extend.cnvs extends regions to neighboring bases                                                       ###

###                                                                                                        ###

##############################################################################################################



extend.cnvs <- function(cnv.matrix, positions) {
	new.output.matrix <- matrix(nrow=nrow(cnv.matrix), ncol=ncol(cnv.matrix))
	new.output.matrix[1,] <- cnv.matrix[1,]
	new.cnv.index <- 0
	new.row <- 1
	for (original.row in 2:nrow(cnv.matrix)) {
# input first cnv
		first.cnv.start <- as.numeric(new.output.matrix[new.row,1])
		first.cnv.end <- as.numeric(new.output.matrix[new.row,2])
# input next cnv
		second.cnv.start <- as.numeric(cnv.matrix[original.row,1])
		second.cnv.end <- as.numeric(cnv.matrix[original.row,2])
# test for overlap
		break.probe.count <- length(which(as.numeric(positions)>first.cnv.end & as.numeric(positions)<second.cnv.start))
		if (break.probe.count <= 200) {
			new.output.matrix[new.row,] <- c(first.cnv.start, second.cnv.end)
		} else {
			new.output.matrix[new.row,] <- c(first.cnv.start, first.cnv.end)
			new.output.matrix[(new.row+1),] <- c(second.cnv.start, second.cnv.end)
			new.row <- new.row + 1
		}
	}
	new.output.matrix <- subset(new.output.matrix, !is.na(new.output.matrix[,1])) 
	colnames(new.output.matrix) <- c("START", "END")
	return(new.output.matrix)
}


##############################################################################################################

###                                                                                                        ###

### event.merge merges called CNVs                                                                         ###

###                                                                                                        ###

##############################################################################################################



event.merge <- function(input.data, merge.max.spacing, event.vs.gap.max) {
	the.data <- input.data
	the.data <- cbind(the.data, size=(as.numeric(the.data[,"end"])-as.numeric(the.data[,"start"])+1))
#	print("merging events")
#	print(paste("Starting event count:", nrow(the.data)))
# for gain, loss, then for each chromosome
	output <- matrix(nrow=nrow(the.data), ncol=(ncol(the.data)+1))
	output.start.row <- 0
	types <- c("GAIN", "LOSS")
	total.row.run <- 0
	row.total.count <- 0
	row.merged.count <- 0
	merged.count <- 0
	for (type.index in 1:2) {
		type <- types[type.index]
#		print(type)
		type.data <- subset(the.data, the.data[,"class"]==type)
		total.row.run <- total.row.run + nrow(type.data)
#		print(total.row.run)
#		print(output.start.row)
		if (nrow(type.data)>1) {
			distance <- NA
			percent.1 <- NA
			percent.2 <- NA
			test.row.index <- 1
			new.row.index <- 2
			current.row.data <- type.data[1,]
			new.row.data <- type.data[2,]
			event.count <- 1
			while (test.row.index < nrow(type.data)) {
#				print(output.start.row)
#				print(test.row.index)
#				print(new.row.index)
				distance.test <- as.numeric(new.row.data["start"]) - as.numeric(current.row.data["end"])

				percent.1.test <- as.numeric(current.row.data["size"])/distance.test
				percent.2.test <- as.numeric(new.row.data["size"])/distance.test
#				print(distance.test)
#				print(percent.1.test) 
#				print(percent.2.test)
				if (distance.test<=as.numeric(merge.max.spacing) & distance.test>0 & percent.1.test >= as.numeric(event.vs.gap.max) & percent.2.test >= as.numeric(event.vs.gap.max)) {
#					print("merge.event")
					merged.count <- merged.count + 1
# update current row data
					event.count <- event.count + 1 
					current.row.data["end"] <- new.row.data["end"]
					current.row.data["size"] <- as.numeric(current.row.data["end"]) - as.numeric(current.row.data["start"]) + 1
					new.row.index <- new.row.index + 1
					if (new.row.index>nrow(type.data)) {
						output.start.row <- output.start.row + 1
						output[output.start.row,] <- c(as.character(current.row.data), 1)
#						print(c(as.character(current.row.data[1:10]), 1))
						test.row.index<- new.row.index
						} else {
							new.row.data <- type.data[new.row.index,]
						}
					} else {
#						print("no merge")
						output.start.row <- output.start.row + 1
						output[output.start.row,] <- c(as.character(current.row.data), event.count)
#						print(c(as.character(current.row.data[1:10]), event.count))
						test.row.index <- new.row.index
						if (test.row.index==nrow(type.data)) {
							current.row.data <- type.data[test.row.index,]
							output.start.row <- output.start.row + 1
							output[output.start.row,] <- c(as.character(current.row.data), 1)
#							print(c(as.character(current.row.data[1:10]), event.count))
						} else {
								new.row.index <- new.row.index + 1
								current.row.data <- type.data[test.row.index,]
								new.row.data <- type.data[new.row.index,]
								event.count <- 1
							}
					}
				}
				
			} else {
				if (nrow(type.data)>0) {
					output[(output.start.row + 1),] <- c(type.data[1,],1)
					output.start.row <- output.start.row + 1
					}
			}
	}
	output <- subset(output, !is.na(output[,1]))
	colnames(output) <- c(colnames(the.data), "merge.count")
	output <- output[,c("start", "end", "class", "merge.count")]
	return(output) 
} 

##############################################################################################################

###                                                                                                        ###

### create.summary.bed.file generates a bed file of all called regions colors can be changed if desired    ###

###                                                                                                        ###

##############################################################################################################



create.summary.bed.file <- function(experiment.data, criteria) {
	chr.column <- experiment.data[,"chr.id"]
	start.col <- experiment.data[,"start"]
	end.col <- experiment.data[,"end"]
	name.col <- paste(
					  substr(experiment.data[,"SAMPLE_ID"], 1,10), 
					  paste(as.character(experiment.data[,"chr.id"]), as.numeric(as.character(paste(experiment.data[,"start"])), as.numeric(as.character(experiment.data[,"end"])), sep=".."), sep=":"), 
					 as.character(experiment.data[,"class"]), 
					  as.numeric(as.character(experiment.data[,"base.count"])), 
					  as.numeric(as.character(round(as.numeric(experiment.data[,"median.ratio"]), digits=2))),
					  as.numeric(as.character(experiment.data[,"size"])), 
					  sep=";"
					  )
	score.col <- experiment.data[,"median.ratio"]
	strand.col <- rep("+", nrow(experiment.data))
	thickStart.col <- experiment.data[,"start"]
	thickEnd.col <- experiment.data[,"end"]
	itemRgb.col <- ifelse(experiment.data[,"class"]=="GAIN", "0,0,255", "255,0,0")
	data.table <- cbind(
						chr.column, 
						start.col,
						end.col,
						name.col, 
						score.col, 
						strand.col, 
						thickStart.col, 
						thickEnd.col,
						itemRgb.col
						)
	output.filename <- paste(criteria$experiment.name, "_cnvcalls.bed", sep="")
	header <- make.bed.header(criteria$experiment.name, "CNV calls")
	setwd(criteria$experiment.directory)
	cat(header, file=output.filename)
	cat("\n", file=output.filename, append=T)
	write.table(data.table, output.filename, append=T, col.names=F, row.names=F, quote=F, sep="\t")
}

##############################################################################################################

###                                                                                                        ###

### Function to annotate events with gene and exon data                                                    ###
###                                                                                                        ###

##############################################################################################################


annotate.cnvs <- function(cnv.output, experiment.metadata) {
	study.temp <- cnv.output
	chr.ids <- unique(cnv.output[,"chr.id"])
	output <- vector("list", length=length(chr.ids))
	if (nrow(study.temp)>0) {
		for (chr.index in 1:length(chr.ids)) {
			chr.id <- as.character(chr.ids[chr.index])
			chr.rows <- which(study.temp[,"chr.id"]==chr.id)
			study.chr <- study.temp[chr.rows,]
			gene.overlap.data <- matrix(nrow=nrow(study.chr), ncol=4)
			gene.chr <- subset(experiment.metadata$gene.information, experiment.metadata$gene.information[,3]==chr.id)
			for (locus.index in 1:nrow(study.chr)) {
				gene.overlap.data[locus.index,] <- gene.overlap.function(study.chr[locus.index,], gene.chr)
				colnames(gene.overlap.data) <- c("Gene", "Refseq", "Overlap", "Exons")
			}
			output[[chr.index]] <- data.frame(study.chr, gene.overlap.data)
		}
		final.output <- mat.build(output)	
		return(final.output)
	} else {
		gene.overlap.data <- matrix(nrow=0, ncol=4)
		colnames(gene.overlap.data) <- c("Gene", "Refseq", "Overlap", "Exons")
		final.output <- data.frame(study.temp, gene.overlap.data)
		return(final.output)
	}	
}  

##############################################################################################################

###                                                                                                        ###

### function to test for gene overlap with cnv                                                             ###     
###                                                                                                        ###

##############################################################################################################


gene.overlap.function <- function(locus.row, gene.chr) {
	range.start <- as.numeric(as.character(locus.row[,"start"]))
	range.end <- as.numeric(as.character(locus.row[,"end"]))
	overlap.table <- subset(gene.chr, as.numeric(gene.chr[,5]) <= range.end & as.numeric(gene.chr[,6]) >= range.start)
	overlap.gene <- vector(length=nrow(overlap.table))
	overlap.refseq <- vector(length=nrow(overlap.table))
	overlap.class <- vector(length=nrow(overlap.table))
	overlap.exon <- vector(length=nrow(overlap.table))
	if (nrow(overlap.table)>0) {
		for (overlap.index in 1:nrow(overlap.table)) {
			overlap.gene[overlap.index] <- overlap.table[overlap.index,1]
			overlap.refseq[overlap.index] <- overlap.table[overlap.index,2]
			overlap.class[overlap.index] <- overlap.call.function(range.start, range.end, overlap.table[overlap.index,5], overlap.table[overlap.index,6], overlap.table[overlap.index,4])
			overlap.exon[overlap.index] <- overlap.exon.function(range.start, range.end, overlap.table[overlap.index,10], overlap.table[overlap.index,11])
		}
		return(c(paste(overlap.gene,collapse=";"),paste(overlap.refseq,collapse=";"),paste(overlap.class,collapse=";"), paste(overlap.exon, collapse=";")))
	} else {
		return(c(NA,NA,NA,NA))
		}
}

##############################################################################################################

###                                                                                                        ###

### function to characterize gene overlap with cnv with regards to exons                                   ###
###                                                                                                        ###

##############################################################################################################


overlap.exon.function <- function(range.start, range.end, starts, ends) {
	exon.starts <- as.numeric(strsplit(starts,split=",")[[1]])
	exon.ends <- as.numeric(strsplit(ends,split=",")[[1]])
	exons <- data.frame(exon.starts, exon.ends, stringsAsFactors=FALSE)
	exon.count <- which(exons[,1]<=range.end & exons[,2]>=range.start)
	return(length(exon.count))
	}



##############################################################################################################

###                                                                                                        ###

### overlap.call.function determines what portion of a gene is overlapped by a cnv                         ###

###                                                                                                        ###

##############################################################################################################



overlap.call.function <- function(range.start, range.end, overlap.start, overlap.end, overlap.strand) {
	call <- NA
	if (overlap.strand=="+") {
# complete
		if (range.start <= overlap.start & range.end >= overlap.end) {call <- "complete"}
# internal
		if (range.start > overlap.start & range.end < overlap.end) {call <- "internal"}
# 5partial
		if (range.start <= overlap.start & range.end < overlap.end) {call <- "5partial"}
# 3partial
		if (range.start > overlap.start & range.end >= overlap.end) {call <- "3partial"}
	} else {
# complete
		if (range.start <= overlap.start & range.end >= overlap.end) {call <- "complete"}
# internal
		if (range.start > overlap.start & range.end < overlap.end) {call <- "internal"}
# 3partial
		if (range.start <= overlap.start & range.end < overlap.end) {call <- "3partial"}
# 5partial
		if (range.start > overlap.start & range.end >= overlap.end) {call <- "5partial"}
	}
	return(call)
} 



###################################################################
###################################################################
###                      General functions                      ###
###################################################################
###################################################################


##############################################################################################################

###                                                                                                        ###

### make.bed.header generates a header for a .bed file                                                     ###

###                                                                                                        ###

##############################################################################################################



make.bed.header <- function(array.id, suffix) {
	file.label <- paste(array.id, suffix, sep="_")
	file.label <- paste('track type=bed name="', file.label, sep="")
	file.label <- paste(file.label, '" description="', sep=" ") 
	file.label <- paste(file.label, paste(array.id, suffix, sep=" "), sep="")
	file.label <- paste(file.label, '" visibility=full itemRgb="On"', sep="")
	return(file.label)
}

##############################################################################################################

###                                                                                                        ###

### mat.build takes a list with identical matrices as the elements and binds them by rows into one table   ###

### used for breaking up by partition and then rejoining                                                   ###
###                                                                                                        ###

##############################################################################################################



mat.build <- function(x) {
	temp.data <- x
	element.names <- names(temp.data)
    for (el.index in 1:length(temp.data)) {
		if (class(temp.data[[el.index]])=="character") {temp.data[[el.index]] <- t(as.matrix(temp.data[[el.index]]))}
	}
 	# find first element with dimensions
    element.col <- 0
    test.col <- 0
    while (element.col == 0) {
    	test.col <- test.col + 1
    	number.cols <- ncol(temp.data[[test.col]])
    	if (class(number.cols) == "NULL") {element.col <- 0} else {element.col <- 1}
    	}
	number.rows <- 0
	for (el.index in 1:length(temp.data)) {
		if (!is.null(nrow(temp.data[[el.index]]))) {
			if (nrow(temp.data[[el.index]])>0) {
				start.row <- number.rows + 1
				number.rows <- number.rows + nrow(temp.data[[el.index]])
			}
		}
	} 
	output <- matrix(nrow=number.rows, ncol=number.cols)
	colnames(output) <- colnames(temp.data[[1]])
	end.row <- 0
	for (el.index in 1:length(temp.data)) {
		if (!is.null(nrow(temp.data[[el.index]]))) {
			if (nrow(temp.data[[el.index]])>0) {
				el.temp <- temp.data[[el.index]]
				start.row <- end.row + 1
				end.row <- end.row + nrow(el.temp)
				output[start.row:end.row,] <- as.matrix(el.temp)
			}
		}
	} 
	return(output)	
}

##############################################################################################################

###                                                                                                        ###

### Function to generate union of targeted bases                                                           ###
###                                                                                                        ###

##############################################################################################################



bait.union <- function(bait.information) {
	bait.intervals <- new(
		"Genome_intervals",
		as.matrix(bait.information[,2:3]),
		closed = as.matrix(cbind(
			rep(TRUE, nrow(bait.information)),
			rep(FALSE, nrow(bait.information))
			)),
		annotation = data.frame(
			seq_name = as.factor(bait.information[,1]),
			inter_base = rep(FALSE, nrow(bait.information))
			)
		)
	bait.union <- interval_union(bait.intervals)
	return(bait.union)
	}
	
	
##############################################################################################################

###                                                                                                        ###

### make.bedgraph.header generates a header for a .bedgraph file                                           ###

###                                                                                                        ###

##############################################################################################################



make.bedgraph.header <- function(file.name, file.description) {
	file.label <- paste('track type=bedGraph name="', file.name, sep="")
	file.label <- paste(file.label, '" description="', sep=" ") 
	file.label <- paste(file.label, file.description, sep="")
	file.label <- paste(file.label, '" visibility=hide', sep="")
	return(file.label)
}
	
	
	
##############################################################################################################

###                                                                                                        ###

### order.chr takes a genomic coordinate table and organizes by chromosome and position                    ###

###                                                                                                        ###

##############################################################################################################



order.chr <- function(data, chr.column, position.column) {
 	sort.order <- vector(length=nrow(data))
	chr.ids <- paste("chr", c(1:22,"X", "Y"), sep="")
	for (chr.index in 1:length(chr.ids)) {
		sort.order <- ifelse(data[,chr.column]==chr.ids[chr.index],chr.index,sort.order)
		}
	return(data[order(sort.order, data[,position.column]),]) 
	}	
		
		
##############################################################################################################

###                                                                                                        ###

### draw.gene is used to generate gene graphics by drawing the exons and linking them with a line          ###

###                                                                                                        ###

##############################################################################################################



draw.gene <- function(features, feature.start.col, feature.end.col, y.top, y.bottom, feature.color, attach, border.color, line.color) {
	if (attach=="Y") {
		for (line.index in 2:nrow(features)) {
			x.start <- as.numeric(features[(line.index-1),feature.end.col])
			x.end <- as.numeric(features[line.index,feature.start.col])
			y.val <- ((y.top-y.bottom)/2) + y.bottom
			coord.plot <- matrix(nrow=2,ncol=2)
			coord.plot[1,] <- c(x.start,y.val)
			coord.plot[2,] <- c(x.end,y.val)
			lines(coord.plot, col=line.color, lwd=3)
		}
	}
	for (feature.index in 1:nrow(features)) {
		x.start <- as.numeric(features[feature.index, feature.start.col])
		x.end <- as.numeric(features[feature.index, feature.end.col])
		if (features[feature.index,"class"]=="UTR") {
			rect(x.start, (y.bottom+((y.top-y.bottom)*.25)), x.end, (y.top-((y.top-y.bottom)*.25)), col=feature.color, border=border.color)
			}
		if (features[feature.index,"class"]=="coding") {
			rect(x.start, y.bottom, x.end, y.top, col=feature.color, border=border.color)

			}
		
		}
	}		
sex.chromosome <- function(sample.id, criteria,experiment.metadata ) {
	
	
}
