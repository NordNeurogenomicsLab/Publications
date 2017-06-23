#!/bin/bash

################################################################################

## Gompers et al. 2017 (Nature Neuroscience)
## MISO Run Parameters

################################################################################

## DECLARE TOOL VERSIONS	
MISO=miso/09b2a9
ANNOTATIONS=annotations/v1
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

source myenv/bin/activate

## INDEX ANNOTATIONS
index_gff --index </mm9_ensGene.gff3> </Isoform/>
index_gff --index </mm9_miso_annotations_v1/SE.mm9.gff3> </Exon-SE/>
index_gff --index </mm9_miso_annotations_v1/A3SS.mm9.gff3> </Exon-A3SS/>
index_gff --index </mm9_miso_annotations_v1/A5SS.mm9.gff3> </Exon-A5SS/>
index_gff --index </mm9_miso_annotations_v1/MXE.mm9.gff3> </Exon-MXE/>
index_gff --index </mm9_miso_annotations_v1/TandemUTR.mm9.gff3> </Exon-TandemUTR/>
index_gff --index </mm9_miso_annotations_v1/RI.mm9.gff3> </Exon-RI/>
index_gff --index </mm9_miso_annotations_v1/AFE.mm9.gff3> </Exon-AFE/>
index_gff --index </mm9_miso_annotations_v1/ALE.mm9.gff3> </Exon-ALE/>

## MERGE BAMS
## We merged sorted individual .bams based on stage and genotype and 
## sorted/indexed the Merged.bam for MISO analysis. To check the impact of any
## one sample on the downstream MISO analysis, sample genotypes were shuffled
## and bams were merged on shuffled IDs. We also generated merged bams by 
## removing one sample at a time and rerunning analysis.
## Based on this analysis, we removed e17.5 sample Chd8-e17-5-2_S159 
## due to lower overall coverage that was impacting the genotype comparisons.

## RUN MISO ON SORTED, INDEXED MERGED BAMS
ANNOTATION={Isoform A3SS A5SS AFE ALE MXE RI SE TandemUTR}
STAGE={e12 e14 e17 P0 P56}
GENOTYPE={WT HT}

miso --run </$ANNOTATION/> <$STAGE.$GENOTYPE.Merged.bam> --output-dir </$STAGE.$GENOTYPE/$ANNOTATION/> --read-len <50/100> --settings-filename <miso_settings.txt>
summarize_miso --summarize-samples </$STAGE.$GENOTYPE/$ANNOTATION/> </$STAGE.$GENOTYPE/$ANNOTATION/>

## COMPARE MISO OUTPUTS BETWEEN GENOTYPE AND STAGE
## Genotype comparisons: WT vs HT
compare_miso --compare-samples </$STAGE.WT/$ANNOTATION/> </$STAGE.HT/$ANNOTATION/> </$STAGE.WT_vs_HT/$ANNOTATION/Comparisons/>
## Stage comparisons: pairwise across all stages
compare_miso --compare-samples </${STAGE[$i]}.WT/$ANNOTATION/> </${STAGE[$j]}/$ANNOTATION> </${STAGE[$i]}_vs_{STAGE[$j]}/$ANNOTATION/Comparisons/>

## FILTER MISO RESULTS
filter_events --filter <$STAGE.WT_vs_HT/$ANNOTATION/Comparisons/compare_miso_output.miso_bf> --num-total 100 --num-sum-inc-exc 10 --delta-psi <0 or 0.1> --bayes-factor 100 --apply-both --output-dir <$STAGE.WT_vs_HT/$ANNOTATION/Comparisons/filtered/>
filter_events --filter <${STAGE[$i]}_vs_{STAGE[$j]}/$ANNOTATION/Comparisons/compare_miso_output.miso_bf> --num-total 100 --num-sum-inc-exc 10 --delta-psi <0,0.1> --bayes-factor 100 --apply-both --output-dir <${STAGE[$i]}_vs_{STAGE[$j]}/$ANNOTATION/Comparisons/filtered/>