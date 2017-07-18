# PanelDoC v 1.2
A targeted-sequencing CNV Calling Script written in R
<br/>Authors: Alex S. Nord, Alex M. Mawla 2011-2015

## Note: We are no longer updating these scripts.  While we are happy to provide guidance and support if you are interested in using or adapting, please note that we have not performed testing to verify compatibility with updates to R and relevant packages.

This file contains information for running and interpreting results from depth of coverage based CNV calling as described in Nord et al., 2011 (Pubmed ID: 21486468).  All the R scripts described are as described in the manuscript, however, no guarantees are made as to their use on other coverage data. Any questions can be sent to asnord@ucdavis.edu if you need using these scripts. 

There are two scripts files that are necessary:

DoC_functions.R
DoC_parameters.R

The parameters file is annotated and controls the analysis.  The functions file is loaded as a source file and all the actual analysis functions are contained in it.  For custom analysis, the functions will need to be modified.  

The two scripts can be downloaded on their own, or together with a sample dataset.

The scripts perform a series of operations on raw coverage data:

1.  Generate genome parameters for targeted regions (number of self chains (see UCSC for details), GC0content in 100-base window, distance from non-targeted base)
2.  Generate median data from all samples
3.  Normalize sample data to the median
4.  Normalize sample data for GC Content Bias
5.  Call CNVs based on relative depth of coverage

Note: The manuscript also describes a secondary analysis using BLAST to confirm CNV calls.  These scripts do not perform the secondary analysis.  I recommend either custom BLAST (as was described in the manuscript) or a program that performs global split-read analysis, such as SLOPE (Pubmed ID: 20876606).  

The DoC_parameter.R file allows the user to specify which operations to perform, as well as containing necessary variables such as the directory path and file names to use.  All the parameters in this file are annotated, and as such should be comprehendible.  Please email with questions if they are not.  The DoC_paramters.R file can either be run line by line in R or as a whole through R or the command line.

For analysis, the data must be split into sample and region, where the regions are discreet genomic regions that are assumed to be independent (i.e. not contiguous) and are analyzed individually.  Regions may be a specific contiguous block or a stretch of mixed targeted and non-targeted bases (e.g. where only exons are targeted and introns are skipped).  While the data must be subdivided, all regions from a sample are considered together during normalization.  Because of this, I recommend limiting the total number of datapoints in any given analysis run.  The scripts are optimized for 1-3Mb sized datasets, and larger datasets should be broken up as described in the notes below.  The sample dataset is a subset of the data reported in the manuscript, with two genes included.  All coding and non-coding bases are targeted with the exception of repeat-masked DNA, which is skipped.   

Some basic sequence and annotation files must also be present.  These files can be downloaded from the UCSC genome browser FTP site.  The example dataset contains all the necessary files in the UCSC folder for the sample dataset, but for any new datasets full sequence data must be added for all the chromosomes (the sample dataset only includes chr10 and chr16).  This data can be found for hg19 at:
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/
All the necessary annotation files are included in the sample dataset, and these match the files of the same name available on the UCSC FTP site.



Input Files:

There are 4 required input files plus the raw coverage data for each sample which must be in the specified experiment or platform directory.  The basic sequence and annotation data must also be present (included in the "UCSC" folder in the sample dataset). Full example files are available in the sample dataset.  Below is a general description and a few sample rows from each file.  In order for the analysis to run, all input files must match exactly the format shown below, including the header.  The input filenames must be entered in the parameter file.  The platform directory contains information about the pool of targeted bases (the platform).  This allows multiple experiments to reference the same platform directory.  The experiment directory contains data about the samples and the raw coverage data.

Platform Directory Files:

Captured region bed file (CaptureRegions.bed in the sample dataset):
This file contains the tarted regions in bed format.  All targeted bases must be within at least one row of the bed file, and no non-targeted or skipped bases must be with regions listed in this file.  For situations where multiple capture probes/oligos/regions overlap, each independent one can be listed.  For more information about the .bed format, see the UCSC genome browser.

First rows:<br>
chr10	89617793	89617913	BI425136563_1869	1000	+ <br>
chr10	89618163	89618283	BI425136563_1872	1000	+ <br>
chr10	89618203	89618323	BI425136563_1873	1000	+ <br>
chr10	89618503	89618623	BI425136563_1877	1000	+ <br>
chr10	89619263	89619383	BI425136563_1896	1000	+ <br>
chr10	89619063	89619183	BI425136563_1891	1000	+ 

Partition/Region file (Platform_partitionInformation.csv):
This file contains the discreet regions for analysis.  In the example file, these are two genes, but they may be genomic regions.  For gene-based graphics, a RefSeq ID must be included in this file.  For CNV calling, each region is handled individually.  For exome analysis or large targeted regions with multiple neighboring genes may need to be analyzed together.  
Fist rows:<br>
PartitionName,PartitionChr,PartitionStart,PartitionEnd,PartitionRefseq<br>
PTEN,chr10,89613195,89738532,NM_000314<br>
CDH1,chr16,68761195,68879444,NM_004360<br>

Experiment Directory Files:

Sample information file (SampleInformation.csv in the sample dataset):
This file contains general information about the samples.  The only required field is the first, with the header "Sample."  Data in this file is merged into the output files.
First rows:<br>
"Sample","Lane","Index"<br>
"1",7,"TAGCTT"<br>
"2",2,"TGACCA"<br>
"3",4,"TGACCA"<br>

Known CNV file (knowncnvs.csv): 
This contains any known CNVs, which are then removed during the normalization process.  Generally, only large CNVs need to be included here.  While this file must be present, it is usually only the header with no data.  If analysis is being performed on a sample with a known large CNV, the coordinates can be entered to increase normalization performance.
First rows:<br>
SampleID,Chr,Start,End<br>

Coverage files:
These must be within a directory named "coverage" in the experiment directory.  There must be one file for each sample and region combination.  For example, if there are three genes targeted and five samples, there would be 15 files within the coverage directory.  Files must be named using the format "XXXSampleXXX_XXXRegionXXX.depth" where XXXSampleXXX is the sample ID in the sample file and XXXRegionXXX is the region ID listed in the partition information file described below.  Each targeted base must be included, or coverage for that base is assumed to be zero.
First rows: <br>
ChrID	Position	Coverage <br>
chr16	68761817	233<br>
chr16	68761818	243<br>
chr16	68761819	246<br>
chr16	68761820	249<br>
chr16	68761821	250<br>





Output Files:

Output files are located in the Experiment directory listed in the parameters files.  Information about main output files:

ExperimentName_QCTable.csv:  QC data for all samples<br>
Detailed field information:<br>
Sample: Sample ID	<br>
Coverage: Median coverage for sample	<br>
NormalizedCoverage: Median normalized coverage for sample<br>
Ratio: Median ratio for sample	<br>
SDCoverage: Standard deviation for raw coverage	<br>
SDNormalizedCoverage: Standard deviation for median coverage	<br>
SDRatio: Standard deviation for ratio 	<br>
Sample_SN_Raw: Signal to noise ratio for raw coverage(median/sd)	<br>
Sample_SN_Normalized: Signal to noise ratio for normalized coverage (median/sd)	<br>
CNV_Zscore: Number of standard deviations from median for ratio value of .5<br>
CNV_Count: Number of CNV calls	<br>
Percent10X: Percent of targeted bases with at least 10X raw coverage	<br>
Percent50X: Percent of targeted bases with at least 50X raw coverage<br>
Percent100X: Percent of targeted bases with at least 100X raw coverage<br>


ExperimentName_CNVs.csv: CNV calls with parameters<br>
Detailed field information:<br>
SAMPLE_ID: sample ID (from sample information input file)<br>
RegionName: partition name (from partition information input file)<br>
chr.id: chromosome<br>
start:  variant start base<br>
end: variant end base<br>
class:  gain or loss<br>
size: variant size<br>
base.count:  total number targeted bases within call<br>
gc.content:  median GC content within call (calculated in 100bp windows from UCSC genome sequence .fa files)<br>
median.selfchain.count: median number of self-chains within call (from UCSC chainSelfLink.txt)<br>
percent.unique.sequence: percent of call not overlappin self-chains<br>
mean.base.spacing:  base.count/size<br>
base.count.criteria: number of targeted bases meeting call criteria<br>
merge.count: number of calls merged to make final call<br>
median.ratio: median relative depth of coverage ratio within call<br>
sd.ratio: standard deviation for relative depth of coverage ratio within call<br>
sn.diff.coverage: (sample coverage - experiment median coverage)/experiment median standard deviation.  This is a measure of the strength of signal from the sample normalized coverage data relative to the noise across all samples in the experiment, similar to a z-score.  Values close to 0 show little difference between the signal and noise, whereas low or high values indicate strong signal.  True positive rare CNVs should have high or low values.<br>
median.coverage.raw: median raw sample coverage within call<br>
median.coverage.normalized: median normalized sample coverage within call<br>
median.experiment.coverage: median raw coverage across all samples within call<br>
median.experiment.sd: standard deviation for coverage for all samples within call<br>
median.experiment.sn: signal:noise ratio (median/sd) for coverage across all samples in experiment within call.<br>
median.experiment.relativesn: median s:n within call/median s:n for full partition.  Values close to 1 indicate that s:n within the call is similar to overall s:n.  Low numbers  indicate lower relative s:n within call. Low numbers are more likely to be false positives.<br>
CopyNumber: copy number estimate (0, 1, 3, 4).  Value of 4 indicates >3 copies - not a quantitative estimate of 4 copies.<br>
SampleCalls: number of calls for sample<br>
Gene: genes overlapped by variant (from UCSC refFlat.txt)<br>
Refseq: refseq values matching overlapped genes (from UCSC refFlat.txt)<br>
Overlap: description of overlap for each gene.  5’ partial and 3’ partial indicate overlap of the 5’ or 3’ end, respectively. <br>Internal indicates that CNV is within gene only.  Complete indicates that CNV extends beyond gene boundaries.<br>
Exons: number of exons overlapped for each gene/refseq ID<br>


Results folders:<br>
bedgraph: bedgraph format files for raw and ratio data for use with a genome browser such as UCSC.<br>
calls: individual output files with all CNV calls for each sample.<br>
normalized: normalized data for each sample and each region.  Can be used as ratio input for other segmentation algorithms.<br>
PDFs: graphical output for each subject and all CNV calls normalized.<br>
raw: raw coverage data for targeted regions generated from the raw input files.<br>





Notes on running analysis:

1.  The DoC_parameters.R file must contain correct directory paths and file names, and the input files must match the expected format exactly.  Error messages are typically due to incorrect format or missing files.  For examples of each file, see the sample dataset.

2.  Low data quality/high inter-subject variability will cause significant problems for any relative depth of coverage analysis.  Some samples may need to be removed if CNV call counts are too high.  

3.  User beware: the analysis scripts are designed for high-coverage targeted data.  Only minimal attempt was made to be memory efficient.  For data sets with very high numbers of samples or very large targeted regions (e.g. exome data), the data should be run in manageable subsets.  For example, for exome data, you could run each chromosome as an independent analysis or even further split the data into smaller genomic blocks.  For robust normalization, larger numbers of samples are preferable for a given region.  For users interested in analyzing exam data, low coverage/high noise may produce such a high false positive rate as to make the analysis problematic.  You will need to significantly increase the minimum size by updating the minimum.base.window and minimum.base.pass in the parameters file.    

4.  The sample dataset contains coverage data from a real experiment.  There are two targeted regions and data for 12 samples for those regions.  The data contains  one true positive call for each region.  The SLOPE output for the regions are also included as an example of split-read confirmation data.  The results from this data set are also available.  

5.  These scripts have also been tested on tumor data, which is expected to have a much higher level of copy number variation.  I implemented a mixture model based method, which may perform better than the expected ratio values used in the main algorithm.  For more information, please contact me.


  
