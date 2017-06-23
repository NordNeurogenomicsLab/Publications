#!/bin/bash

################################################################################

## Gompers et al. 2017 (Nature Neuroscience)
## RNA-seq Alignment Parameters

################################################################################

## DECLARE TOOL VERSIONS	
STAR=star/2.4.2a
FASTQC=fastqc/v0.11.2
SAMTOOLS=samtools/1.3
SUBREAD=subread/1.5.0-p1
PYTHON=python2.7
RSEQC=rseqc/2.6.3
KENTUTILS=kentutils/308

## GENOME PARAMETERS
GENOME=mm9
SIZES=mm9.chrom.sizes
FASTA=WholeGenomeFasta/genome.fa
REFGTF=Genes/genes.gtf
REFBED=Genes/mm9_knownGene.bed

## STAR PARAMETERS
SAINDEX=14
SJDBOVERHANG=99
RUNTHREADN=6
GENOMELOAD=LoadAndKeep
OUTREADSUM=Fastx
OUTSAMTYPE="BAM Unsorted SortedByCoordinate"
OUTWIGTYPE=bedGraph
OUTWIGSTRAND=Unstranded
OUTWIGREFPREFIX=chr
OUTWIGNORM=None
CHIMSEGMENTMIN=15
QUANTMODE="TranscriptomeSAM GeneCounts"

################################################################################

## RNA-SEQ EFFECTIVE COMMAND LINE
STAR --runMode alignReads --runThreadN 6 --genomeDir </STAR/genomes/mm9/> --genomeLoad LoadAndKeep --readFilesIn <SAMPLE/R1.fastq> --outFileNamePrefix <SAMPLENAME.> --outReadsUnmapped Fastx --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Unstranded --outWigReferencesPrefix chr --outWigNorm None --chimSegmentMin 15 --quantMode TranscriptomeSAM GeneCounts

## FEATURECOUNTS COMMAND LINE
featureCounts -t exon -g gene_id -s 0 -a $REFGTF <SAMPLE.Aligned.out.bam>
featureCounts -t exon -g gene_id -s 1 -a $REFGTF <SAMPLE.Aligned.out.bam>
featureCounts -t exon -g gene_id -s 2 -a $REFGTF <SAMPLE.Aligned.out.bam>

## FASTQC COMMAND LINE
fastqc -o <fastqc/R1_data> -f fastq --contaminants <contaminant_list.txt> <SAMPLE/R1.fastq>

## RSEQC COMMAND LINE
mismatch_profile.py -l 50 -i <SAMPLE.Aligned.sortedByCoord.out.bam>
insertion_profile.py -s \"SE\" -i <SAMPLE.Aligned.sortedByCoord.out.bam>
insertion_profile.py -s \"PE\" -i <SAMPLE.Aligned.sortedByCoord.out.bam>
deletion_profile.py -i <SAMPLE.Aligned.sortedByCoord.out.bam> -l 50
read_distribution.py -r $REFBED -i <SAMPLE.Aligned.sortedByCoord.out.bam>
read_duplication.py  -i <SAMPLE.Aligned.sortedByCoord.out.bam>
RPKM_saturation.py  -r $REFBED -i <SAMPLE.Aligned.sortedByCoord.out.bam>
bam_stat.py -i <SAMPLE.Aligned.sortedByCoord.out.bam>
clipping_profile.py -i <SAMPLE.Aligned.sortedByCoord.out.bam> -s \"SE\"
clipping_profile.py -i <SAMPLE.Aligned.sortedByCoord.out.bam> -s \"PE\"
inner_distance.py -r $REFBED -i <SAMPLE.Aligned.sortedByCoord.out.bam>
junction_annotation.py -r $REFBED -i <SAMPLE.Aligned.sortedByCoord.out.bam>
junction_saturation.py -r $REFBED -i <SAMPLE.Aligned.sortedByCoord.out.bam>
read_GC.py -i <SAMPLE.Aligned.sortedByCoord.out.bam>
read_NVC.py -i <SAMPLE.Aligned.sortedByCoord.out.bam>
read_quality.py -i <SAMPLE.Aligned.sortedByCoord.out.bam>
read_hexamer.py -r $FASTA -i <SAMPLE.Aligned.sortedByCoord.out.bam>
FPKM_count.py -r $REFBED -i <SAMPLE.Aligned.sortedByCoord.out.bam>
split_bam.py -r $REFBED -i <SAMPLE.Aligned.sortedByCoord.out.bam>
split_bam.py -r mm9_rRNA.bed -i <SAMPLE.Aligned.sortedByCoord.out.bam>
geneBody_coverage.py -r $REFBED -i <SAMPLE.Aligned.sortedByCoord.out.bam>
