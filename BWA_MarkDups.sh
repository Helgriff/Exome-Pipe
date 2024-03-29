#! /bin/bash
#$ -cwd
#$ -j y
#$ -m e
#$ -pe smp 4

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job
#cwd=`pwd`; export cwd

#/************************************
# note: fastq file patterns are used in the perl script detectSampleLanes.pl and READ_FILE1 & READ_FILE2
# they need to be adjusted accordingly
#************************************/

echo $'\n'"["`date`"]: Job started."

SAMPLE_ID=$1
SAMPLE_PATH=$2
SCRIPTS_DIR=$3
REF_DIR=$4
SCRATCH_DIR=$5
TARGETS=$6
SEQ_PLATFORM=$7 # Valid values are: ILLUMINA, SOLID, LS454, HELICOS and PACBIO
Library_ID=$8 #Better to use {PI name}.{Batch Number} to identify the library
COV_DIR_NAME=$9
GATK_OUT_DIR_NAME=${10}
DUP_FREE_BAM_DIR_NAME=${11}
WRKGDIR_NAME=${12}
JAVA_TMP_DIR_NAME=${13}
FREE_DIR_NAME=${14}
OUT_PATH=${15}
REF_FILE2=${16}

REF_FILE="${REF_DIR}/${REF_FILE2}"
# GATK base quality recalibration covariates
COVARIATES="-cov ReadGroupCovariate -cov QualityScoreCovariate -cov ContextCovariate -cov CycleCovariate"

#/************************ fmscluster add modules
module load apps/samtools/0.1.19
module load apps/bwa/0.7.12
module load apps/picard/1.130
module load apps/perl/5.18.2
module add apps/java/jre-1.8.0_25
export _JAVA_OPTIONS="-XX:-UseLargePages -Xms8g"
#PICARDDIR="$PICARD_PATH"
PICARDDIR="/opt/software/bsu/versions/picard-tools-1.119"
GATKDIR="/opt/software/bsu/versions/GenomeAnalysisTK-3.2-2"
INDIR="${SAMPLE_PATH}/${SAMPLE_ID}"
#INDIR="${SAMPLE_PATH}/Sample_${SAMPLE_ID}"
#INDIR="${SAMPLE_PATH}/PFC_${SAMPLE_ID}"
OUTDIR="${OUT_PATH}/${SAMPLE_ID}"
#OUTDIR="${OUT_PATH}/Sample_${SAMPLE_ID}"
#OUTDIR="${OUT_PATH}/PFC_${SAMPLE_ID}"
WRKGDIR="${SCRATCH_DIR}/${SAMPLE_ID}${WRKGDIR_NAME}"
JAVA_TMP_DIR="${SCRATCH_DIR}/${SAMPLE_ID}.Java.tmp.dir.${JAVA_TMP_DIR_NAME}"
DUP_FREE_BAM_DIR="${OUTDIR}/${DUP_FREE_BAM_DIR_NAME}"

Picard_sort="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/SortSam.jar VALIDATION_STRINGENCY=LENIENT"
Picard_nodups="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT"
Picard_CleanSam="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/CleanSam.jar VALIDATION_STRINGENCY=LENIENT"

# prepare folders
echo $'\n'"mkdir $JAVA_TMP_DIR"
if [ ! -d $JAVA_TMP_DIR ]; then
	mkdir $JAVA_TMP_DIR
else
	echo "$JAVA_TMP_DIR exists"
fi

echo $'\n' "mkdir $OUTDIR"
if [ ! -d $OUTDIR ]; then
	mkdir $OUTDIR
else 
	echo "$OUTDIR exists"	
fi

echo $'\n' "mkdir $WRKGDIR"
if [ ! -d $WRKGDIR ]; then
	mkdir $WRKGDIR
else 
	echo "$WRKGDIR exists"	
fi

echo $'\n' "mkdir $DUP_FREE_BAM_DIR"
if [ ! -d $DUP_FREE_BAM_DIR ]; then
	mkdir ${DUP_FREE_BAM_DIR}
else
        echo "${DUP_FREE_BAM_DIR} exists"   
fi

PICARD_TEMP="$WRKGDIR/Picard_Temp"
echo $'\n'"mkdir $PICARD_TEMP"
if [ ! -d $WRKGDIR ]; then
	mkdir $PICARD_TEMP
else
        echo "$PICARD_TEMP exists"   
fi
PICARD_LOG="$WRKGDIR/${SAMPLE_ID}_picard.log"

################## BWA Bit ####################
LANES_STRING=`perl ${SCRIPTS_DIR}/detectSampleLanes.pl $SAMPLE_PATH $SAMPLE_ID`
LANES=($LANES_STRING)

#echo "$LANES ...are the lanes detected!!"
echo "${LANES[@]} ...are the lanes detected!!"

if [ ${#LANES[@]} -lt 1 ]; then
	exit "Error: No lane info detected, check file name patterns and paths"
fi

BAM_FILE_LIST=""
echo $'\n'"started to loop through lanes"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
for LANE in "${LANES[@]}"
do
	echo "---------------- hua li de fen jie xian ----------------"
	echo "process lane $LANE for sample ${SAMPLE_ID}"
		
	READ_FILE1=( $INDIR/${SAMPLE_ID}*L00${LANE}_R1_001.fast* )
	READ_FILE2=( $INDIR/${SAMPLE_ID}*L00${LANE}_R2_001.fast* )
	
	READ_FILE1=${READ_FILE1[0]}; READ_FILE1=${READ_FILE1%'.gz'}
	READ_FILE2=${READ_FILE2[0]}; READ_FILE2=${READ_FILE2%'.gz'}
	
	echo $READ_FILE1; echo $READ_FILE2

	SAM_FILE1="$WRKGDIR/${SAMPLE_ID}_${LANE}.sam"
	
	#bwa mapping using bwa mem
	# echo $'\n'"["`date`"]: bwa starts to align the lane ${LANE}."
	# echo bwa mem -R "@RG\tID:FlowCell.${LANE}.${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:${SEQ_PLATFORM}\tLB:${Library_ID}.${SAMPLE_ID}" -t 4 -M $REF_FILE $READ_FILE1 $READ_FILE2 "> $SAM_FILE1"
	# bwa mem -R "@RG\tID:FlowCell.${LANE}.${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:${SEQ_PLATFORM}\tLB:${Library_ID}.${SAMPLE_ID}" -t 4 -M $REF_FILE $READ_FILE1 $READ_FILE2 > $SAM_FILE1 # "-M" is for Picard compatibility (mark shorter split hits as secondary)
	
	#bwa mapping using original 2-step aln and sampe commands
	SAI_FILE1_1="$WRKGDIR/${SAMPLE_ID}_L00${LANE}_1.sai"
	SAI_FILE1_2="$WRKGDIR/${SAMPLE_ID}_L00${LANE}_2.sai"  
	echo $'\n'"["`date`"]: bwa starts to align the lane ${LANE}."
	echo "bwa aln -t 4 $REF_FILE $READ_FILE1 > $SAI_FILE1_1"
	bwa aln -t 4 $REF_FILE $READ_FILE1 > $SAI_FILE1_1
	echo "bwa aln -t 4 $REF_FILE $READ_FILE2 > $SAI_FILE1_2"
	bwa aln -t 4 $REF_FILE $READ_FILE2 > $SAI_FILE1_2
	echo "bwa sampe -r "@RG\tID:FlowCell.${LANE}.${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:${SEQ_PLATFORM}\tLB:${Library_ID}.${SAMPLE_ID}" $REF_FILE $SAI_FILE1_1 $SAI_FILE1_2 $READ_FILE1 $READ_FILE2 > $SAM_FILE1"
	bwa sampe -r "@RG\tID:FlowCell.${LANE}.${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:${SEQ_PLATFORM}\tLB:${Library_ID}.${SAMPLE_ID}" $REF_FILE $SAI_FILE1_1 $SAI_FILE1_2 $READ_FILE1 $READ_FILE2 > $SAM_FILE1

	rm ${SAI_FILE1_1}
	rm ${SAI_FILE1_2}

	#Sam to bam, Sort and Remove Duplicates
	RAW_BAM="$SAM_FILE1.bam"
	SORTED_BAM_FILE_NODUPS="$SAM_FILE1_nodups.bam"
	DUP_FREE_BAM_FILE_LIST="$DUP_FREE_BAM_FILE_LIST-I ${SORTED_BAM_FILE_NODUPS} "

	echo $'\n'"["`date`"]: PICARD to sort the sam file ${SAM_FILE1}"
	echo "$Picard_sort INPUT=$SAM_FILE1 OUTPUT=$RAW_BAM SORT_ORDER=coordinate TMP_DIR=$PICARD_TEMP"
	$Picard_sort INPUT=$SAM_FILE1 OUTPUT=$RAW_BAM SORT_ORDER=coordinate
	echo $'\n'"samtools index $RAW_BAM"
    samtools index $RAW_BAM
	echo $'\n'"["`date`"]: PICARD to remove duplicates from the bam file ${RAW_BAM}"
	echo "$Picard_nodups INPUT=$RAW_BAM OUTPUT=$SORTED_BAM_FILE_NODUPS METRICS_FILE=$PICARD_LOG REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$PICARD_TEMP"
	$Picard_nodups INPUT=$RAW_BAM OUTPUT=$SORTED_BAM_FILE_NODUPS METRICS_FILE=$PICARD_LOG REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=$PICARD_TEMP
	echo $'\n'"samtools index $SORTED_BAM_FILE_NODUPS"
    samtools index $SORTED_BAM_FILE_NODUPS
    
    rm "${RAW_BAM}"
    rm "${RAW_BAM}.bai"
    
done
echo $'\n'"------------------------------"
echo "all lanes have been processed!"
echo $'\n'"------------------------------"

echo $'\n'"Merge alignment and move files"
# Merge Multi-lane bams
DUP_FREE_BAM="${DUP_FREE_BAM_DIR}/${SAMPLE_ID}_nodups.bam"
DUP_FREE_BAM_IDX="${DUP_FREE_BAM_DIR}/${SAMPLE_ID}_nodups.bam.bai"
if [ ${#LANES[@]} -eq 1 ]
then
	echo "mv $SORTED_BAM_FILE_NODUPS $DUP_FREE_BAM"
	mv $SORTED_BAM_FILE_NODUPS $DUP_FREE_BAM
	mv "${SORTED_BAM_FILE_NODUPS}.bai" "${DUP_FREE_BAM_IDX}"
else # merge with GATK printReads
	echo $'\n'java -Xmx8g -Djava.io.tmpdir=${JAVA_TMP_DIR} -jar $GATKDIR/GenomeAnalysisTK.jar -T PrintReads -nct 4 -R $REF_FILE $DUP_FREE_BAM_FILE_LIST -o $DUP_FREE_BAM
    java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $GATKDIR/GenomeAnalysisTK.jar -T PrintReads -nct 4 -R $REF_FILE $DUP_FREE_BAM_FILE_LIST -o $DUP_FREE_BAM
    echo "samtools index $DUP_FREE_BAM"
    samtools index $DUP_FREE_BAM
fi

# cleaning bit #
echo rm -r $WRKGDIR
rm -r $WRKGDIR
echo rm -r $JAVA_TMP_DIR
rm -r $JAVA_TMP_DIR

echo $'\n'"["`date`"]: The whole job is DONE!!"

# runtime calculation
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
echo "----------------------------- hua li li de fen jie xian -----------------------------"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
echo "exit status $?"
