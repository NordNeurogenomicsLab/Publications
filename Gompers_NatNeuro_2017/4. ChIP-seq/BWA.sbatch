#!/bin/bash
# BWA.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=BWA_SE		   			# A single job name for the array

# Load Modules, and Run Programs

echo 
echo BWA.sbatch version 1.3.0
echo ... started running on $(date)
echo 
echo

echo
echo This job is running on
/bin/hostname
echo
echo


i=$((8 + $SLURM_ARRAY_TASK_ID));

BWA_NAME=${!i}

echo
echo bwa_name is $BWA_NAME
echo  

INPUT_BWA=$(find $3 -iname $BWA_NAME)
echo INPUT_BWA is $INPUT_BWA

module load bwa
module load samtools/1.3
module load fastx_toolkit
module load fastqc
module load cutadapt
module load bedtools2
module load java
echo

export TMPDIR=./$SLURM_JOBID

sleep $SLURM_ARRAY_TASK_ID

mkdir -p "${4}$BWA_NAME"

if [ "$6" = "yes" ]; 
then
	perl trim_galore.pl --fastqc_args "-t $1" -o $4 -q $5 --illumina --stringency 1 --length 20 $INPUT_BWA
		
	SUFF=".fastq.gz"
	suff="_trimmed.fq.gz"	
	NEW_INPUT_BWA=$4${BWA_NAME/$SUFF/$suff}
	echo  
	echo Running bwa on sample: $NEW_INPUT_BWA
	echo  
	echo   
	
	cd "${4}$BWA_NAME"
	
	bwa aln -t $1 $2 $NEW_INPUT_BWA | bwa samse $2 - $NEW_INPUT_BWA | samtools view -q $5 -Shu - | \
		samtools sort -@ $1 -o $NEW_INPUT_BWA.srt.bam -
	
	echo
		
	samtools index $NEW_INPUT_BWA.srt.bam
	echo

	echo
	
	echo " ... deleting unnecessary files from bwa"
	echo 
	rm $NEW_INPUT_BWA
	echo "$NEW_INPUT_BWA ... removed"
	
	fastqcDIR="$4../fastqc/"
	zipSUFF="_trimmed_fastqc.zip"
	mv ${NEW_INPUT_BWA/$suff/$zipSUFF} $fastqcDIR
	echo "${NEW_INPUT_BWA/$suff/$zipSUFF} moved to $fastqcDIR"
	
	htmlSUFF="_trimmed_fastqc.html"
	mv ${NEW_INPUT_BWA/$suff/$htmlSUFF} $fastqcDIR
	echo "${NEW_INPUT_BWA/$suff/$htmlSUFF} moved to $fastqcDIR"
	
	txtSUFF="_trimming_report.txt"
	mv $4$BWA_NAME$txtSUFF $fastqcDIR
	echo "$4$BWA_NAME$txtSUFF moved to $fastqcDIR"
	echo 

	echo Finished running BWA align $NEW_INPUT_BWA
else
	echo  
	echo Running bwa on sample: $INPUT_BWA
	echo  
	echo  
	
	OUTPUT_BWA=$4$BWA_NAME
	 
	cd "${4}$BWA_NAME"
	
	bwa aln -t $1 $2 $INPUT_BWA | bwa samse $2 - $INPUT_BWA | samtools view -q $5 -Shu - | \
		samtools sort -@ $1 -o $OUTPUT_BWA.srt.bam -
	echo
	
	samtools index $OUTPUT_BWA.srt.bam
	echo
	
	bamToBed -i $OUTPUT_BWA.srt.bam > $OUTPUT_BWA.bed
	echo
	
	echo Finished running BWA align $INPUT_BWA
	
	echo " ... deleting unnecessary files from bwa"
	echo 
	rm $INPUT_BWA.bam
	echo "$INPUT_BWA.bam ... removed"
	rm $INPUT_BWA
	echo "$INPUT_BWA ... removed"
	echo 

fi

echo  
echo
echo 
echo Completed execution on $(date)
echo 

sleep 1

rm -rf "${4}$BWA_NAME"

# changelog
# 1.3.0	code searches for file to bwa
# 1.2.1	included minor corrections on file output
# 1.2.0	enabled selection of raw or trimmed sequences for alignment
# 1.1.0	introduced quality filtering, adapter trimming and alignment of trimmed sequences 
#       (with name changing)
# 1.0.2	introduced piping for disk space saving
# 1.0.1	introduced soft coding
# 1.0.0	first working edition

