#!/bin/bash
# Version 0 by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=FSTQCpost

echo
echo 
echo FASTQC.sbatch version 1.1.0
echo ... started running on $(date)
echo 
echo 

echo
echo This job is running on
/bin/hostname
echo
echo


# Load Modules, and Run Programs
module load fastqc/v0.11.2

i=$((4 + $SLURM_ARRAY_TASK_ID));

FASTA_NAME=${!i}

OLDPREF=".fastq.gz"
ORIGPREF="_trimmed.fq.gz.srt.bam"

FASTA_NAME=${FASTA_NAME/$OLDPREF/$ORIGPREF}

FASTA=$(find $2 -iname $FASTA_NAME)

echo 
echo FASTA_NAME is $FASTA_NAME
echo 
echo FASTA is $FASTA
echo

fastqc -t $1 -o $3 $FASTA

NADA=""
INITNAME="_trimmed_fastqc"
FINALNAME="_trimmed_post_fastqc"

OLDFASTA="${FASTA_NAME/${ORIGPREF}/$INITNAME}"
NEWFASTA="${FASTA_NAME/${ORIGPREF}/$FINALNAME}"

MYONE=$(find $3 -iname "${OLDFASTA}.html")
mv $MYONE "$3${NEWFASTA}.html"

MYONE=$(find $3 -iname "${OLDFASTA}.zip")
mv $MYONE "$3${NEWFASTA}.zip"


echo  
echo
echo 
echo Completed execution on $(date)
echo 

echo

# changelog
# 1.1.0	code searches for file to fastqc
# 1.0.0	first working edition
