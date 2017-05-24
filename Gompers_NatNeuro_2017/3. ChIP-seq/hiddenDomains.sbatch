#!/bin/bash
# hiddenDomains.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=hiddenDomains

echo
echo 
echo hiddenDomains.sbatch version 1.0.0
echo ... started running on $(date)
echo  
echo
echo
echo This job is running on
/bin/hostname
echo
echo


# Load Modules, and Run Programs
module load python/2.7.6-ubuntu
module load samtools/1.3

t=$(( 6 + 5 * $SLURM_ARRAY_TASK_ID));
c=$(( 7 + 5 * $SLURM_ARRAY_TASK_ID));
o=$(( 8 + 5 * $SLURM_ARRAY_TASK_ID));
n=$(( 9 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((10 + 5 * $SLURM_ARRAY_TASK_ID));

MACS_t=${!t}
MACS_c=${!c}
MACS_o=${!o}

MACS_sple=${!n}
MACS_ctrl=${!m}

SUFF="\-vs\.Input"
SUFF2="\-vs\.NegCtl"
suff=""

EXT=$1
PVALUE=$2
CHRSIZES=$3
CODE=$4
HIDDIR=$5

cd $CODE

if [ "$MACS_c" != "null" ]
then
	"${CODE}hiddenDomains" -g $CHRSIZES -b $EXT -p $PVALUE -t $MACS_t -c $MACS_c \
						   -o "${HIDDIR}${MACS_o}.b${EXT}.p${PVALUE}"
else
	MACS_o=${MACS_o/$SUFF/$suff}
	MACS_o=${MACS_o/$SUFF2/$suff}
	"${CODE}hiddenDomains" -g $CHRSIZES -b $EXT -p $PVALUE -t $MACS_t \
						   -o "${HIDDIR}${MACS_o}.b${EXT}.p${PVALUE}"
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
