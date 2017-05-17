#!/bin/bash
# Homer.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=Homer_annot

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

t=$(( 8 + 5 * $SLURM_ARRAY_TASK_ID));
c=$(( 9 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((10 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((11 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((12 + 5 * $SLURM_ARRAY_TASK_ID));

trt=${!t}
ctl=${!c}
otp=${!o}

homerfile="$7$otp$3"

allmotiffs="homerMotifs.all.motifs"
allmotiffs="${2}$otp/$allmotiffs"

ANNOT_HOMER=".motif.counts"
ANNOT_HOMER="${2}$otp/$otp$ANNOT_HOMER"


echo  
echo
echo homer name is $homerfile
echo allmottifs is $allmotiffs
echo Annotationfile is $ANNOT_HOMER
echo 

PREFIX="-"

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

# Perform known motif analysis on all subset bed files
# (not implemented yet)


# Map individual motifs of interest to peak bedfiles
# motif counts for all motifs
mkdir -p "$7${otp}${3}.tmp"
cd "$7${otp}${3}.tmp"
annotatePeaks.pl $homerfile $1 -m $allmotiffs -nmotifs > $ANNOT_HOMER


# motif co-occurence for all motifs
MATRIX_HOMER=".cooccur"
MATRIX_HOMER="$ANNOT_HOMER$MATRIX_HOMER"

annotatePeaks.pl $homerfile $1 -m $allmotiffs -matrix $MATRIX_HOMER > $MATRIX_HOMER.out


# motif histogram for primary motifs
annotatePeaks.pl $homerfile $1 -m $allmotiffs $SIZE_ARG -hist $6 > $ANNOT_HOMER.hist.txt
annotatePeaks.pl $homerfile $1 -m $allmotiffs $SIZE_ARG -hist $6 -ghist > $ANNOT_HOMER.heatmap.hist.txt
annotatePeaks.pl $homerfile $1 -m $allmotiffs $SIZE_ARG -hist $6 -CpG > $ANNOT_HOMER.CpGhist.txt

# GO analysis

GO_dir="GO"
GO_dir="${2}$otp/$GO_dir"

mkdir -p $GO_dir

annotatePeaks.pl $homerfile $1 -m $allmotiffs -nmotifs -go $GO_dir -genomeOntology $GO_dir > $MATRIX_HOMER.GO

cd ..
rm -rf "$7${otp}${3}.tmp"
echo  
echo
echo 
echo Completed execution on $(date)
echo 
echo 
echo 

# changelog
# 1.0.0	first working edition