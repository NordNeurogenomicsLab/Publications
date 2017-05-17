#!/bin/bash
# BWA_PE.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=BWA_PE		   			# A single job name for the array

# Load Modules, and Run Programs

echo 
echo BWA.sbatch version 1.0.0
echo ... started running on $(date)
echo 
echo 

echo
echo This job is running on
/bin/hostname
echo
echo


i=$(( 7 + 3 * $SLURM_ARRAY_TASK_ID));
j=$(( 8 + 3 * $SLURM_ARRAY_TASK_ID));
m=$(( 9 + 3 * $SLURM_ARRAY_TASK_ID));

BWA_R1=${!i}
BWA_R2=${!j}
BWA_out=${!m}

echo
echo BWA_names are $BWA_R1
echo   and $BWA_R2
echo  
echo

INPUT_BWA_R1=$(find $3 -iname $BWA_R1)
echo INPUT_BWA_R1 is $INPUT_BWA_R1
echo 

INPUT_BWA_R2=$(find $3 -iname $BWA_R2)
echo INPUT_BWA_R2 is $INPUT_BWA_R2
echo 

module load bwa
module load samtools
module load fastx_toolkit
module load fastqc
module load cutadapt
module load bedtools2
echo

if [ $6 = "yes" ]; 
then
	perl trim_galore.pl --fastqc_args "-t $1" -o $4 -q $5 $INPUT_BWA_R1
	perl trim_galore.pl --fastqc_args "-t $1" -o $4 -q $5 $INPUT_BWA_R2
		
	SUFF=".fastq.gz"
	suff="_trimmed.fq.gz"	
	NEW_INPUT_BWA_R1=$4${BWA_R1/$SUFF/$suff}
	NEW_INPUT_BWA_R2=$4${BWA_R2/$SUFF/$suff}
	echo  
	echo Running bwa on samples: $NEW_INPUT_BWA_R1 and $NEW_INPUT_BWA_R2
	echo  
	echo   

	NEW_INPUT_SAI_R1="$NEW_INPUT_BWA_R1.sai"
	NEW_INPUT_SAI_R2="$NEW_INPUT_BWA_R2.sai"
	
	BWA_PE_BAM="$4$BWA_out"

# 	bwa aln -t $1 $2 $NEW_INPUT_BWA_R1 > $NEW_INPUT_SAI_R1
# 	bwa aln -t $1 $2 $NEW_INPUT_BWA_R2 > $NEW_INPUT_SAI_R2
	
	echo
	echo starting sampe on $(date)
	echo
	echo
	sleep 1

# 	bwa sampe -P $2 $NEW_INPUT_SAI_R1 $NEW_INPUT_SAI_R2 $NEW_INPUT_BWA_R1 $NEW_INPUT_BWA_R2 | samtools view -Shu - > $BWA_PE_BAM.bam
	bwa mem -t $1 $2 $NEW_INPUT_BWA_R1 $NEW_INPUT_BWA_R2 | samtools view -Shu - > $BWA_PE_BAM.bam
	
	echo
	sleep 1
	
	samtools sort -@ $1 $BWA_PE_BAM.bam $BWA_PE_BAM.srt
	echo
	sleep 1
	
	samtools index $BWA_PE_BAM.srt.bam
	echo
	sleep 1

	bamToBed -i $BWA_PE_BAM.srt.bam > $BWA_PE_BAM.bed
	echo
	sleep 1
	
	echo " ... deleting unnecessary files from bwa"
	echo 
	rm $NEW_INPUT_SAI_R1
	echo "$NEW_INPUT_SAI_R1 ... removed"
	rm $NEW_INPUT_SAI_R2
	echo "$NEW_INPUT_SAI_R2 ... removed"
# 	rm $BWA_PE_BAM
# 	echo "$BWA_PE_BAM ... removed"
	
# 	fastqcDIR="$4../fastqc/"
# 	zipSUFF="_trimmed_fastqc.zip"
# 	mv ${NEW_INPUT_BWA_R1/$suff/$zipSUFF} $fastqcDIR
# 	echo "${NEW_INPUT_BWA_R1/$suff/$zipSUFF} moved to $fastqcDIR"
# 	
# 	mv ${NEW_INPUT_BWA_R2/$suff/$zipSUFF} $fastqcDIR
# 	echo "${NEW_INPUT_BWA_R2/$suff/$zipSUFF} moved to $fastqcDIR"
# 	
# 	htmlSUFF="_trimmed_fastqc.html"
# 	mv ${NEW_INPUT_BWA_R1/$suff/$htmlSUFF} $fastqcDIR
# 	echo "${NEW_INPUT_BWA_R1/$suff/$htmlSUFF} moved to $fastqcDIR"
# 	
# 	htmlSUFF="_trimmed_fastqc.html"
# 	mv ${NEW_INPUT_BWA_R2/$suff/$htmlSUFF} $fastqcDIR
# 	echo "${NEW_INPUT_BWA_R2/$suff/$htmlSUFF} moved to $fastqcDIR"
# 	
# 	txtSUFF="_trimming_report.txt"
# 	mv $4$NEW_INPUT_BWA_R1$txtSUFF $fastqcDIR
# 	echo "$4$NEW_INPUT_BWA_R1$txtSUFF moved to $fastqcDIR"
# 
# 	mv $4$NEW_INPUT_BWA_R2$txtSUFF $fastqcDIR
# 	echo "$4$NEW_INPUT_BWA_R2$txtSUFF moved to $fastqcDIR"
# 	echo 

	echo Finished running BWA sampe $BWA_PE_BAM
else
	echo  
	echo Running bwa on samples: $INPUT_BWA_R1 and $INPUT_BWA_R2
	echo  
	echo   
	INPUT_SAI_R1="$INPUT_BWA_R1.sai"
	INPUT_SAI_R2="$INPUT_BWA_R2.sai"
	
# 	bwa aln -t $1 $2 $INPUT_BWA_R1 > $INPUT_SAI_R1
# 	bwa aln -t $1 $2 $INPUT_BWA_R2 > $INPUT_SAI_R2
	
	echo
	echo starting sampe on $(date)
	echo
	echo

# 	bwa sampe -P $2 $INPUT_SAI_R1 $INPUT_SAI_R2 $INPUT_BWA_R1 $INPUT_BWA_R2 | samtools view -Shu - > $BWA_PE_BAM.bam
	bwa mem -t $1 $2 $INPUT_BWA_R1 $INPUT_BWA_R2 | samtools view -Shu - > $BWA_PE_BAM.bam
	
	echo
	
	samtools sort -@ $1 $BWA_PE_BAM.bam $BWA_PE_BAM.srt
	echo
	
	samtools index $BWA_PE_BAM.srt.bam
	echo

	bamToBed -i $BWA_PE_BAM.srt.bam > $BWA_PE_BAM.bed
	echo
		
	echo Finished running BWA sampe $BWA_PE_BAM
	
	echo " ... deleting unnecessary files from bwa"
	echo 
	rm $INPUT_SAI_R1
	echo "$INPUT_SAI_R1 ... removed"
	rm $INPUT_SAI_R2
	echo "$INPUT_SAI_R2 ... removed"
	rm $INPUT_BWA
	echo "$BWA_PE_BAM ... removed"
	echo 

fi

echo  
echo
echo 
echo Completed execution on $(date)
echo 

sleep 1


# changelog
# 1.0.0	first working edition

