#!/bin/bash
#
#SBATCH --job-name=ChIPQC

echo
echo 
echo ChIPQC.sbatch version 1.0.0
echo ... started running on $(date)
echo 

factorfile=$1
chipqc_dir_popped=$2
blacklist_file=$3
genome_ver=$4
quality=$5
facetby=$6
colorby=$7
threads=$8

/usr/bin/Rscript ChIPQC.R $factorfile $chipqc_dir_popped $blacklist_file $genome_ver $quality $facetby $colorby $threads

rm $1

echo  
echo
echo 
echo Completed execution on $(date)
echo 

