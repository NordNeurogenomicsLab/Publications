#!/bin/bash
# Homer.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=Homer

echo
echo 
echo Homer.sbatch version 1.0.0
echo ... started running on $(date)
echo  
echo  
echo
echo This job is running on
/bin/hostname
echo
echo

# Load Modules, and Run Programs
module load homer/4.7
module load python-libs/2.7.6-ubuntu

t=$((16 + 5 * $SLURM_ARRAY_TASK_ID));
c=$((17 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((18 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((19 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((20 + 5 * $SLURM_ARRAY_TASK_ID));

trt=${!t}
ctl=${!c}
otp=${!o}

if [[ "$3" == ".bayesPeaks.bed" ]];
then
	homerfile=${15}
else
	homerfile=${12}
fi

homerfile="$homerfile$otp$3"
PREFIX="-"
SIZE_ARG=""
OUTDIR="${2}$otp"


echo  
echo
echo homer name is $homerfile
echo 

if [[ "$7" == "void" ]];
then
	MASK=""
else
	MASK="$7"
fi

if [[ "$6" == "void" ]];
then
	BKGRND=""
else
	BKGRND="$6"
fi

if [[ "$8" == "void" ]];
then
	RNA=""
else
	RNA="$8"
fi

if [[ "${13}" == "binomial" ]];
then
	SCORING=""
else
	SCORING="-h"
fi

if [[ "$9" == "void" ]];
then
	NORM=""
else
	NORM="$PREFIX$9"
fi

if [[ "$4" == "given" ]];
then
	SIZE_ARG="size given"
	SIZE_ARG="$PREFIX$SIZE_ARG"
else
	if [[ "$5" == "given" ]];
	then
		SIZE_ARG="size given"
		SIZE_ARG="$PREFIX$SIZE_ARG"
	else
		SIZE_ARG="size"
		SIZE_ARG="$PREFIX$SIZE_ARG -$4,$5"
	fi
fi

# Perform motif discovery on peak file
findMotifsGenome.pl $homerfile $1 $OUTDIR $SIZE_ARG $BKGRND $MASK $RNA $NORM $SCORING ${10} -p ${11}

if [[ ! "${14}" == "null" ]];
then
	homerfile=${15}
	mkdir -p "${OUTDIR}.tmp"
	cd "${OUTDIR}.tmp"
	homerfile="$homerfile$otp${14}"
	findMotifsGenome.pl $homerfile $1 $OUTDIR $SIZE_ARG $BKGRND $MASK $RNA $NORM $SCORING ${10} -p ${11}
	cd ..
	rm -rf "${OUTDIR}.tmp"
fi

echo  
echo
echo 
echo Completed execution on $(date)
echo 
echo 
echo 

# changelog
# 1.0.0	first working edition
