#!/bin/bash
# MACS_sple.sbatch by Rinaldo Catta-Preta @ UC Davis

#SBATCH --job-name=MACS2sple

echo
echo 
echo MACS_sple.sbatch version 1.4.0
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

t=$((12 + 5 * $SLURM_ARRAY_TASK_ID));
c=$((13 + 5 * $SLURM_ARRAY_TASK_ID));
o=$((14 + 5 * $SLURM_ARRAY_TASK_ID));
n=$((15 + 5 * $SLURM_ARRAY_TASK_ID));
m=$((16 + 5 * $SLURM_ARRAY_TASK_ID));

MACS_t=${!t}
MACS_c=${!c}
MACS_o=${!o}

MACS_sple=${!n}
MACS_ctrl=${!m}

MEM=${10}
LOG=${11}

PAIR="_PE"

if [ "${MACS_o/$PAIR}" = "$MACS_o" ];
then
	FORMAT="BAM"
else
	FORMAT="BAMPE"
fi


echo  
echo
echo macs2 name is $MACS_o
echo 

term="_model.r"
modelterm=".model"
nomodelterm=".no_model"

PVALUE=$( echo "$6 > 0" | bc -l )
EXTSIZE=$( echo "$8 > 0" | bc -l )

if [[ "$3" == 'rn' ]];
then
	GENSIZE='1.87e9'
	# source: doi: 10.1038/srep33736
else
	GENSIZE=$3
fi


# BROAD PEAK CALL
cd $1

if [[ $EXTSIZE -eq 1 ]];
then
	if [[ $PVALUE -eq 1 ]];
	then
		macs2 callpeak -t $MACS_t -n "$MACS_sple$nomodelterm" --broad -f $FORMAT \
			-g $GENSIZE --broad-cutoff $4 --bdg -p $6 $7 --nomodel --extsize=$8 --outdir $1
		
		echo

		if [[ ! -f "${MACS_ctrl}${nomodelterm}_peaks.broadPeak" ]]
		then
			if [ "$MACS_c" != "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$nomodelterm" --broad -f $FORMAT \
					-g $GENSIZE --broad-cutoff $4 --bdg -p $6 $7 --nomodel --extsize=$8 --outdir $1
				echo
			fi
		fi

	else
		macs2 callpeak -t $MACS_t -n "$MACS_sple$nomodelterm" --broad -f $FORMAT \
			-g $GENSIZE --broad-cutoff $4 --bdg -q $5 $7 --nomodel --extsize=$8 --outdir $1 
		
		echo

		if [[ ! -f "${MACS_ctrl}${nomodelterm}_peaks.broadPeak" ]]
		then
			if [ "$MACS_c" != "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$nomodelterm" --broad -f $FORMAT \
					-g $GENSIZE --broad-cutoff $4 --bdg -q $5 $7 --nomodel --extsize=$8 --outdir $1 
					
			echo
			fi
		fi
	fi
else
	if [[ $PVALUE -eq 1 ]];
	then
		macs2 callpeak -t $MACS_t -n "$MACS_sple$modelterm" --broad -f $FORMAT \
			-g $GENSIZE --broad-cutoff $4 --bdg -p $6 $7 --outdir $1
		
		echo

		if [[ ! -f "${MACS_ctrl}${modelterm}_peaks.broadPeak" ]]
		then
			if [ "$MACS_c" != "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$modelterm" --broad -f $FORMAT \
					-g $GENSIZE --broad-cutoff $4 --bdg -p $6 $7 --outdir $1
				
				echo
			fi
		fi	
	else
		macs2 callpeak -t $MACS_t -n "$MACS_sple$modelterm" --broad -f $FORMAT \
			-g $GENSIZE --broad-cutoff $4 --bdg -q $5 $7 --outdir $1 
		
		echo

		if [[ ! -f "${MACS_ctrl}${modelterm}_peaks.broadPeak" ]]
		then
			if [ "$MACS_c" != "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$modelterm" --broad -f $FORMAT \
					-g $GENSIZE --broad-cutoff $4 --bdg -q $5 $7 --outdir $1 
				
				echo
			fi
		fi
	fi
fi	
echo  
echo  ... completed broad peak call on $(date)
echo  
echo  


cd $2

if [[ $EXTSIZE -eq 1 ]];
then
	if [[ $PVALUE -eq 1 ]];
	then
		macs2 callpeak -t $MACS_t -n "$MACS_sple$nomodelterm" -f $FORMAT -g $GENSIZE \
			--call-summits --bdg -p $6 $7 --nomodel --extsize=$8 --outdir $2 
		
		echo
			
		macs2 callpeak -t $MACS_t -n "$MACS_sple$modelterm" -f $FORMAT -g $GENSIZE \
			--call-summits --bdg -p $6 --outdir $2  
		
		echo
		
		if [[ ! -f "${MACS_ctrl}${nomodelterm}_peaks.broadPeak" ]]
		then
			if [ ! "$MACS_c" = "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$nomodelterm" -f $FORMAT -g $GENSIZE \
					--call-summits --bdg -p $6 $7 --nomodel --extsize=$8 --outdir $2  
		
				echo
			fi
		fi
			
		if [[ ! -f "${MACS_ctrl}${modelterm}_peaks.broadPeak" ]]
		then
			if [ ! "$MACS_c" = "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$modelterm" -f $FORMAT -g $GENSIZE \
					--call-summits --bdg -p $6 --outdir $2  
		
				echo
			fi
		fi
	else
		macs2 callpeak -t $MACS_t -n "$MACS_sple$nomodelterm" -f $FORMAT -g $GENSIZE \
			--call-summits --bdg -q $5 $7 --nomodel --extsize=$8 --outdir $2  
		
		echo
			
		macs2 callpeak -t $MACS_t -n "$MACS_sple$modelterm" -f $FORMAT -g $GENSIZE \
			--call-summits --bdg -q $5 --outdir $2  
		
		echo
		
		if [[ ! -f "${MACS_ctrl}${nomodelterm}_peaks.broadPeak" ]]
		then
			if [ ! "$MACS_c" = "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$nomodelterm" -f $FORMAT -g $GENSIZE \
					--call-summits --bdg -q $5 $7 --nomodel --extsize=$8 --outdir $2  
		
				echo
			fi
		fi
		
		if [[ ! -f "${MACS_ctrl}${modelterm}_peaks.broadPeak" ]]
		then
			if [ ! "$MACS_c" = "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$modelterm" -f $FORMAT -g $GENSIZE \
					--call-summits --bdg -q $5 --outdir $2  
		
				echo
			fi
		fi
	fi
else
	if [[ $PVALUE -eq 1 ]];
	then
		macs2 callpeak -t $MACS_t -n "$MACS_sple$modelterm" -f $FORMAT -g $GENSIZE \
			--call-summits --bdg -p $6 $7 --outdir $2  
		
		echo

		if [[ ! -f "${MACS_ctrl}${modelterm}_peaks.broadPeak" ]]
		then
			if [ ! "$MACS_c" = "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$modelterm" -f $FORMAT -g $GENSIZE \
					--call-summits --bdg -p $6 $7 --outdir $2  
		
				echo
			fi
		fi
	else
		macs2 callpeak -t $MACS_t -n "$MACS_sple$modelterm" -f $FORMAT -g $GENSIZE \
			--call-summits --bdg -q $5 $7 --outdir $2  
		
		echo

		if [[ ! -f "${MACS_ctrl}${modelterm}_peaks.broadPeak" ]]
		then
			if [ ! "$MACS_c" = "null" ]
			then
				macs2 callpeak -t $MACS_c -n "$MACS_ctrl$modelterm" -f $FORMAT -g $GENSIZE \
					--call-summits --bdg -q $5 $7 --outdir $2  
		
				echo
			fi
		fi
	fi
fi	

wait

echo  
echo  ... completed narrow peak call on $(date)

if [ -f $2$MACS_sple$nomodelterm$term ];
then
	Rscript $2$MACS_sple$nomodelterm$term
	echo  
	echo 
	echo Completed pdf generation of broad peak files on $(date)
	echo 

fi

if [ -f $2$MACS_ctrl$nomodelterm$term ];
then
	Rscript $2$MACS_ctrl$nomodelterm$term
	echo  
	echo 
	echo Completed pdf generation of broad peak files on $(date)
	echo 

fi

if [ -f $2$MACS_sple$modelterm$term ];
then
	Rscript $2$MACS_sple$modelterm$term
	echo  
	echo 
	echo Completed pdf generation of broad peak files on $(date)
	echo 

fi
echo 

if [ -f $2$MACS_ctrl$modelterm$term ];
then
	Rscript $2$MACS_ctrl$modelterm$term
	echo  
	echo 
	echo Completed pdf generation of broad peak files on $(date)
	echo 

fi
echo 
echo 
echo Completed pdf generation of narrow peak model files on $(date)
echo 
cd $1

if [ -f $1$MACS_sple$nomodelterm$term ];
then
	Rscript $1$MACS_sple$nomodelterm$term
	echo  
	echo 
	echo Completed pdf generation of broad peak files on $(date)
	echo 

fi

if [ -f $1$MACS_ctrl$nomodelterm$term ];
then
	Rscript $1$MACS_ctrl$nomodelterm$term
	echo  
	echo 
	echo Completed pdf generation of broad peak files on $(date)
	echo 

fi

if [ -f $1$MACS_sple$modelterm$term ];
then
	Rscript $1$MACS_sple$modelterm$term
	echo  
	echo 
	echo Completed pdf generation of broad peak files on $(date)
	echo 

fi
echo 

if [ -f $1$MACS_ctrl$modelterm$term ];
then
	Rscript $1$MACS_ctrl$modelterm$term
	echo  
	echo 
	echo Completed pdf generation of broad peak files on $(date)
	echo 

fi
echo 


echo  
echo
echo 
echo Completed execution on $(date)
echo 
echo 
echo 

# changelog
# 1.4.0	introduced peak calling on pair-ended reads
# 1.3.0	introduced additional parameter settings
# 1.2.0	introduced soft coding for peak calling parameters
# 1.1.0	included single ChIP, inputDNA and NegCtrl peak calling
# 1.0.0	first working edition

