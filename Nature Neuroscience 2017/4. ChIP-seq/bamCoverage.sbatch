#!/bin/bash
# bamCoverage.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=bamCoverage

echo
echo 
echo bamCoverage.sbatch version 1.0.0
echo ... started running on $(date)
echo  
echo

echo
echo This job is running on 
/bin/hostname
echo
echo


# Load Modules, and Run Programs
module load python-libs/2.7.6-ubuntu
module load deeptools/2.3.6
module load samtools/1.3

t=$(( 9 + 5 * $SLURM_ARRAY_TASK_ID));
c=$((10 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((11 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((12 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((13 + 5 * $SLURM_ARRAY_TASK_ID));

MACS_t=${!t}
MACS_c=${!c}
MACS_o=${!o}

MACS_sple=${!n}
MACS_ctrl=${!m}

BROAD=$1
NARROW=$2
DEEPDIR=$3
THREADS=$4
EXT=$5
IGN_DUP=$6
NORMALIZ=$7
GENOME=$8

echo
echo $BROAD
echo $NARROW
echo $DEEPDIR
echo $THREADS
echo $EXT
echo $IGN_DUP
echo $NORMALIZ
echo $GENOME
echo

echo
echo MACS_c is $MACS_c
echo MACS_t is $MACS_t
echo 

if [[ "$NORMALIZ" == "normalizeTo1x" ]]
then
	NORMAL_SUFF="normalized1x"
elif [[ "$NORMALIZ" == "normalizeUsingRPKM" ]]
then
	NORMAL_SUFF="normalizedRPKM"
fi

if [[ "$GENOME" == "mm9" ]]
then
	GENFORNORM="2150570000"
elif [[ "$GENOME" == "hg19" ]]
then
	GENFORNORM="2451960000"
elif [[ "$GENOME" == "rn5" ]]
then
	GENFORNORM="2572853723"
	# http://genomewiki.ucsc.edu/index.php/Mm10_Genome_size_statistics
else
	echo
	echo
	echo
	echo " ... Curently supported genomes are mm9, hg19 and rn5"
	echo
	echo
fi

if [[ "$NORMALIZ" == "normalizeTo1x" ]]
then
	NORMALIZ="${NORMALIZ} $GENFORNORM"
fi

INDEX='.bai'

echo
echo

if [[ ! -f "$MACS_t$INDEX" ]]
then
	samtools index $MACS_t &
else
	echo
	echo "... indexing $MACS_t not necessary"
	echo
fi

if [ "$MACS_c" != "null" ]
then
	if [[ ! -f "$MACS_c$INDEX" ]]
	then
		samtools index $MACS_c
	else
		echo
		echo "... indexing $MACS_c not necessary"
		echo
	fi
fi

wait

cd $DEEPDIR

declare -a arr=("${MACS_t}" "${MACS_c}")

if [ "${MACS_o/NegCtl}" = "$MACS_o" ]
then
	arr=("${MACS_c}" "${MACS_t}")
fi

for COMPARE in "${arr[@]}"
do
	if [ "${COMPARE}" != "null" ]
	then
		echo "... working on ${COMPARE} coverage files"
		echo 

		if [[ "$NORMALIZ" == "none" ]]
		then
			OUTFILE="${COMPARE}.coverage.bw"
			if [[ "$IGN_DUP" == "no" ]]
			then
				bamCoverage -b "${COMPARE}" -o "$OUTFILE" -of bigwig -p "$THREADS" -e $EXT
			else
				bamCoverage -b "${COMPARE}" -o "$OUTFILE" -of bigwig -p "$THREADS" -e $EXT "$IGN_DUP"
			fi
		else
			OUTFILE="${COMPARE}.${NORMAL_SUFF}.coverage.bw"
			if [[ "$IGN_DUP" == "no" ]]
			then
				bamCoverage -b "${COMPARE}" -o "$OUTFILE" -of bigwig -p "$THREADS" -e $EXT --"$NORMALIZ"
			else
				bamCoverage -b "${COMPARE}" -o "$OUTFILE" -of bigwig -p "$THREADS" -e $EXT "$IGN_DUP" --"$NORMALIZ"
			fi
		fi

		mv $OUTFILE ../deeptools/

		echo "... done with ${OUTFILE}"
		echo 
		echo 
	fi
done          


echo  
echo
echo 
echo Completed execution on $(date)
echo 
echo 
echo 

# changelog
# 1.0.0	first working edition