#!/bin/bash
# MACS_diff.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=MACS2diff

echo
echo 
echo MACS_diff.sbatch version 1.0.0
echo ... started running on $(date)
echo  
echo  

echo
echo This job is running on
/bin/hostname
echo
echo


# Load Modules, and Run Programs
module load macs2/2.1.0.20140616

echo 

BROAD_DIR=$1
EXTSIZE=$2
DIFFCALL_DIR=$3
MEMORY=$4
CUTOFF=$5
MIN_LENGTH=$6
MAX_GAP=$7
DEPTH1=$8
DEPTH2=$9
LOG_DIR=${10}
CONTROLS=${11}
NARROW_DIR=${12}


echo
echo
echo $BROAD_DIR
echo $NARROW_DIR
echo $EXTSIZE 
echo $DIFFCALL_DIR
echo $MEMORY
echo $CUTOFF
echo $MIN_LENGTH
echo $MAX_GAP
echo $DEPTH1
echo $DEPTH2
echo $LOG_DIR
echo $CONTROLS
echo
echo

shift; 
shift; 
shift; 
shift; 
shift; 
shift; 
shift; 
shift; 
shift; 
shift; 
shift; 
shift

array=( "$@" )

batches=$(( ${#array[@]} / 5 ))

if [ "$batches" -gt 0 ]
then
	batches=$(($batches - 1))
fi

echo "no. of batches is $batches"
echo

samples=()
samples_trt=()
samples_ctl=()

origterm=""
pile=".no_model_treat_pileup.bdg"
lambda=".no_model_control_lambda.bdg"
xls=".no_model_peaks.xls"
newg=$(( EXTSIZE / 2 ))
#TYPES="$BROAD_DIR $NARROW_DIR"
TYPES="$NARROW_DIR"

SUFF="\-vs\.Input"
suff=""

if [ "$CONTROLS" != "0" ]
then
	for i in `seq 0 $batches`;
	do
		test_sple="${array[$i * 5 + 2]}"
		teste="vs\.Input"
	
		samples+=( ${array[$i * 5 + 2]} )
		echo ${array[$i * 5 + 2]}
		samples_trt+=( ${array[$i * 5 + 0]} )
		echo ${array[$i * 5 + 0]}
		samples_ctl+=( ${array[$i * 5 + 1]} )
		echo ${array[$i * 5 + 1]}
		echo "..."
	done
else
	for i in `seq 0 $batches`;
	do
		test_sple="${array[$i * 5 + 0]}"
		
		samples+=$test_sple
		echo $test_sple
		samples_trt+=( ${array[$i * 5 + 0]} )
		echo ${array[$i * 5 + 0]}
		samples_ctl+=( ${array[$i * 5 + 1]} )
		echo "${array[$i * 5 + 1]} not used"
		echo "..."
	done
fi

for k in $TYPES
do
	cd $k

	for i in `seq 0  $(( ${#samples[@]} - 2 ))`;
	do
		for j in `seq $(( i + 1 ))  $(( ${#samples[@]} - 1 ))`;
		do
			hot="${samples[$i]}"
			echo "$hot"
			cold="${samples[$j]}"
			echo "$cold"
	
			hotprefix="${hot/$SUFF/$suff}"
			coldprefix="${cold/$SUFF/$suff}"
			
			if [ "$k" == "$BROAD_DIR"  ]
			then
				oprefix="${coldprefix}_vs_${hotprefix}_diff_BROAD"
			else
				oprefix="${coldprefix}_vs_${hotprefix}_diff_NARROW"
			fi
			
			echo "$oprefix"
			echo 
			echo 
	
			if [ "$CONTROLS" != "0" ]
			then
				TRT=(`egrep "tags after filtering in treatment" "$k$cold$xls"`)
				CTL=(`egrep "tags after filtering in control" "$k$cold$xls"`)
				trt=${TRT[-1]}
				ctl=${CTL[-1]}

				D1=$( [ $trt -le $ctl ] && echo "$trt" || echo "$ctl" )

				TRT=(`egrep "tags after filtering in treatment" "$k$hot$xls"`)
				CTL=(`egrep "tags after filtering in control" "$k$hot$xls"`)
				trt=${TRT[-1]}
				ctl=${CTL[-1]}

				D2=$( [ $trt -le $ctl ] && echo "$trt" || echo "$ctl" )
		
				echo "D1: $D1   D2: $D2"
				echo
				echo
		
				cd $DIFFCALL_DIR
		
				KOchecker="KO"
		
				if [ "${cold/$KOchecker}" = "$cold" ];
				then
					macs2 bdgdiff --t1 "$k$cold$pile" --c1 "$k$cold$lambda" \
								  --t2 "$k$hot$pile" --c2 "$k$hot$lambda" \
								  --d1 "$D1" --d2 "$D2" -g $newg -l "$EXTSIZE" --o-prefix "$oprefix" \
								  --outdir "$DIFFCALL_DIR" -C $CUTOFF -l $MIN_LENGTH \
								  -g $MAX_GAP  
				else
					macs2 bdgdiff --t1 "$k$hot$pile" --c1 "$k$hot$lambda" \
								  --t2 "$k$cold$pile" --c2 "$k$cold$lambda" \
								  --d1 "$D2" --d2 "$D1" -g $newg -l "$EXTSIZE" --o-prefix "$oprefix" \
								  --outdir "$DIFFCALL_DIR" -C $CUTOFF -l $MIN_LENGTH \
								  -g $MAX_GAP 
				fi
			else
				TRT=(`egrep "tags after filtering in treatment" "$k$cold$xls"`)
				D1=${TRT[-1]}

				TRT=(`egrep "tags after filtering in treatment" "$k$hot$xls"`)
				D2=${TRT[-1]}

				echo "D1: $D1   D2: $D2"
				echo
				echo
		
				cd $DIFFCALL_DIR
		
				KOchecker="KO"
		
				if [ "${cold/$KOchecker}" = "$cold" ];
				then
					macs2 bdgdiff --t1 "$k$cold$pile" --c1 "$k$cold$lambda" \
								  --t2 "$k$hot$pile" --c2 "$k$hot$lambda" \
								  --d1 "$D1" --d2 "$D2" -g $newg -l "$EXTSIZE" --o-prefix "$oprefix" \
								  --outdir "$DIFFCALL_DIR" -C $CUTOFF -l $MIN_LENGTH \
								  -g $MAX_GAP 
				else
					macs2 bdgdiff --t1 "$k$hot$pile" --c1 "$k$hot$lambda" \
								  --t2 "$k$cold$pile" --c2 "$k$cold$lambda" \
								  --d1 "$D2" --d2 "$D1" -g $newg -l "$EXTSIZE" --o-prefix "$oprefix" \
								  --outdir "$DIFFCALL_DIR" -C $CUTOFF -l $MIN_LENGTH \
								  -g $MAX_GAP 
				fi
			fi
				echo 
				echo "... finished bdgdiff on $oprefix on $(date)"
		done
	done

done

wait

echo  
echo
echo 
echo Completed execution with peaks on $(date)
echo 
echo 
echo 



# changelog
# 1.1.0 added narrow peak differential calls
# 1.0.0	first working edition