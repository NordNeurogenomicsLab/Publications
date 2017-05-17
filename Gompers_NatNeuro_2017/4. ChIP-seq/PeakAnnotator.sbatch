#!/bin/bash
# PeakAnnotator.sh by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=PeakAnnotator		   			# A single job name for the array

# Load Modules, and Run Programs

echo 
echo PeakAnnotator.sbatch version 1.1.0
echo ... started running on $(date)
echo
echo 
module load java
echo
echo
echo This job is running on
/bin/hostname
echo
echo

t=$(( 8 + 5 * $SLURM_ARRAY_TASK_ID));
c=$(( 9 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((10 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((11 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((12 + 5 * $SLURM_ARRAY_TASK_ID));

MACS_o=${!o}     
MACS_sple=${!n}		
MACS_ctrl=${!m}		

term=".no_model_peaks.narrowPeak"
termG=".no_model_peaks.gappedPeak"

CHIP="$3$MACS_sple$term"
CHIP_B="$6$MACS_sple$term2"
CONTROL="$3$MACS_ctrl$term"
CTRL_B="$6$MACS_ctrl$term2"

OUTDIR="$2"

echo
cd $1

if [[ "$7" == "MACS2" ]];
then
    term2=".no_model_peaks.broadPeak"
    BROAD="$6$MACS_o$term2"
    NARROW="$3$MACS_o$term"
    GAPPED="$6$MACS_o$termG"
    
    java -jar -Xmx512m PeakAnnotator.jar -u NDG -p $GAPPED -a $4 -o $OUTDIR -g all -r $5
    java -jar -Xmx512m PeakAnnotator.jar -u TSS -p $GAPPED -a $4 -o $OUTDIR -g all -r $5
    
fi

echo "... starting gene symbol annotation"

perl /share/nordlab/users/rinaldo/chipseq/code/Gene_Annotator.pl $OUTDIR "${MACS_o}.no_model"


echo
echo
echo
echo

echo
echo Completed execution on $(date)
echo

