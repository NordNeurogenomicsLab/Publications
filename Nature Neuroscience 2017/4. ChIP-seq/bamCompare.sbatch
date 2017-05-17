#!/bin/bash
# bamCompare.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=bamCompare

echo
echo 
echo bamCompare.sbatch version 1.0.0
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

t=$(( 9 + 3 * $SLURM_ARRAY_TASK_ID));
c=$((10 + 3 * $SLURM_ARRAY_TASK_ID));
n=$((11 + 3 * $SLURM_ARRAY_TASK_ID));

PRIMARY=${!t}
SECONDARY=${!c}
NAME=${!n}

PSEUDONUM=$1
RATIOTYPE=$2
THREADS=$3
NORMALIZ=$4
IGN_DUP=$5
DEEPDIR=$6
GENOME=$7
EXTEND=$8

echo
echo $PSEUDONUM
echo $RATIOTYPE
echo $THREADS
echo $NORMALIZ
echo $IGN_DUP
echo $DEEPDIR
echo $GENOME
echo

echo
echo Primary is $PRIMARY
echo Secondary is $SECONDARY
echo Output prefix is $NAME

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

if [ "$PRIMARY" != "null" ]
then
	if [[ ! -f "$PRIMARY$INDEX" ]]
	then
		samtools index $PRIMARY &
	else
		echo
		echo "... indexing $PRIMARY not necessary"
		echo
	fi
fi

if [ "$SECONDARY" != "null" ]
then
	if [[ ! -f "$SECONDARY$INDEX" ]]
	then
		samtools index $SECONDARY
	else
		echo
		echo "... indexing $SECONDARY not necessary"
		echo
	fi
fi

wait

cd $DEEPDIR

echo
echo " ... starting bamCompare on $(date)"
echo

if [ "$PRIMARY" != "null" ] && [ "$SECONDARY" != "null" ]
then
	if [[ "$NORMALIZ" == "none" ]]
	then
		OUTFILE="${DEEPDIR}${NAME}.${RATIOTYPE}.pseudocount${PSEUDONUM}.bw"
		if [[ "$IGN_DUP" == "no" ]]
		then
			bamCompare  -b1 "$PRIMARY" -b2 "$SECONDARY" -o "$OUTFILE" -of bigwig --ratio "$RATIOTYPE" --pseudocount "$PSEUDONUM" -e $EXTEND
		else
			bamCompare  -b1 "$PRIMARY" -b2 "$SECONDARY" -o "$OUTFILE" -of bigwig --ratio "$RATIOTYPE" --pseudocount "$PSEUDONUM" "$IGN_DUP" -e $EXTEND
		fi
	else
		OUTFILE="${DEEPDIR}${NAME}.${NORMAL_SUFF}.${RATIOTYPE}.pseudocount${PSEUDONUM}.bw"
		if [[ "$IGN_DUP" == "no" ]]
		then
			bamCompare  -b1 "$PRIMARY" -b2 "$SECONDARY" -o "$OUTFILE" -of bigwig --ratio "$RATIOTYPE" --pseudocount "$PSEUDONUM" --"$NORMALIZ" -e $EXTEND
		else
			bamCompare  -b1 "$PRIMARY" -b2 "$SECONDARY" -o "$OUTFILE" -of bigwig --ratio "$RATIOTYPE" --pseudocount "$PSEUDONUM" "$IGN_DUP" --"$NORMALIZ" -e $EXTEND
		fi	
	fi

	echo Output is "$OUTFILE"
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