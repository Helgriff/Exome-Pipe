#! /bin/bash
#$ -cwd
#$ -j y
#$ -m e

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

#/************ less frequently changed options
REF_FILE="${REF_DIR}/${REF_FILE2}"
# GATK base quality recalibration covariates
COVARIATES="-cov ReadGroupCovariate -cov QualityScoreCovariate -cov ContextCovariate -cov CycleCovariate"
#************/

#/************************ fmscluster add modules
FastUniq="fastuniq"
module load apps/fastuniq/1.1
module load apps/samtools/0.1.19
module load apps/bwa/0.7.12
module load apps/picard/1.130
module load apps/bedtools/2.19.0
module load apps/perl/5.18.2
module add apps/java/jre-1.8.0_25
#module load apps/gatk/3.2-protected
#module load apps/gatk-queue/3.2-protected
export _JAVA_OPTIONS="-XX:-UseLargePages -Xms8g"
PICARDDIR="$PICARD_PATH"
GATKDIR="/sharedlustre/IGM/GATK/GATK3"
INDIR="${SAMPLE_PATH}/Sample_${SAMPLE_ID}"
OUTDIR="${OUT_PATH}/Sample_${SAMPLE_ID}"
WRKGDIR="${SCRATCH_DIR}/${SAMPLE_ID}${WRKGDIR_NAME}"
JAVA_TMP_DIR="${SCRATCH_DIR}/${SAMPLE_ID}.Java.tmp.dir.${JAVA_TMP_DIR_NAME}"
Picard_sort="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/SortSam.jar VALIDATION_STRINGENCY=LENIENT"
Picard_nodups="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT"
Picard_CleanSam="java -Djava.io.tmpdir=$JAVA_TMP_DIR -Xmx8g -jar $PICARDDIR/CleanSam.jar VALIDATION_STRINGENCY=LENIENT"
#*************************/

##Output folders' names
COV_DIR="${OUTDIR}/${COV_DIR_NAME}"
DUP_FREE_BAM_DIR="${OUTDIR}/${DUP_FREE_BAM_DIR_NAME}"

# prepare folders
echo $'\n'mkdir $OUTDIR
if [ ! -d $OUTDIR ]; then
	mkdir $OUTDIR
else 
	echo "$OUTDIR exists"	
fi


echo $'\n'mkdir $WRKGDIR
if [ ! -d $WRKGDIR ]; then
	mkdir $WRKGDIR
else 
	echo "$WRKGDIR exists"	
fi

echo $'\n'mkdir $JAVA_TMP_DIR
if [ ! -d $JAVA_TMP_DIR ]; then
	mkdir $JAVA_TMP_DIR
else
	echo "$JAVA_TMP_DIR"
fi

################## BWA Bit ####################
LANES_STRING=`perl ${SCRIPTS_DIR}/detectSampleLanes.pl $SAMPLE_PATH $SAMPLE_ID`
LANES=($LANES_STRING)

echo "$LANES ...are the lanes detected!!"

if [ ${#LANES[@]} -lt 1 ]; then
	exit "Error: No lane info detected, check file name patterns and paths"
fi

#BAM_FILE_LIST=""
echo $'\n'"started to loop through lanes"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
echo "----------------------------- hua li li de fen jie xian -----------------------------"
for LANE in "${LANES[@]}"
do
	echo "---------------- hua li de fen jie xian ----------------"
	echo "process lane $LANE for sample ${SAMPLE_ID}"

	READ_FILE1=( $INDIR/${SAMPLE_ID}*_L00${LANE}_R1_001.fast* )
	READ_FILE2=( $INDIR/${SAMPLE_ID}*_L00${LANE}_R2_001.fast* )
	
	READ_FILE1=${READ_FILE1[0]}; READ_FILE1=${READ_FILE1%'.gz'}
	READ_FILE2=${READ_FILE2[0]}; READ_FILE2=${READ_FILE2%'.gz'}
	
	echo $READ_FILE1; echo $READ_FILE2
	# can handle zipped file now, will do unzip automatically
	if [ -f "$READ_FILE1.gz" ]; then
		echo "find gzipped fastq file $READ_FILE1.gz, unzip now.."
		gunzip -c "$READ_FILE1.gz" > "$READ_FILE1"
	else
		if [ -f $READ_FILE1 ]; then
			echo "find file $READ_FILE1"
		else
			exit "fastq file name error, can not find file $READ_FILE1"
		fi
	fi

	if [ -f "$READ_FILE2.gz" ]; then
		echo "find gzipped fastq file $READ_FILE2.gz, unzip now.."
        gunzip -c "$READ_FILE2.gz" > "$READ_FILE2"
	else
		if [ -f $READ_FILE2 ]; then
                        echo "find file $READ_FILE2"
                else
                        exit "fastq file name error, can not find file $READ_FILE2"
                fi
        fi
	
	READ_FILE1_nodup="$INDIR/${SAMPLE_ID}_L00${LANE}_R1_001.nodup.fastq"
    READ_FILE2_nodup="$INDIR/${SAMPLE_ID}_L00${LANE}_R2_001.nodup.fastq"
	FASTQ_LIST_FILE="$INDIR/${SAMPLE_ID}_L00${LANE}.temp.qlist"
	
	# remove dups with FastUniq
	echo -e "$READ_FILE1\n$READ_FILE2" > $FASTQ_LIST_FILE
	echo $FastUniq -i "${FASTQ_LIST_FILE}" -t q -o "${READ_FILE1_nodup}" -p "${READ_FILE2_nodup}"
	$FastUniq -i "${FASTQ_LIST_FILE}" -t q -o "${READ_FILE1_nodup}" -p "${READ_FILE2_nodup}"
	
done

echo $'\n'"["`date`"]: Fastuniq part is DONE!!"

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
printf "Fastuniq runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
echo "exit status $?"
