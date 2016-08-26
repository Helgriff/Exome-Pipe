#! /bin/bash

## Script to launch Exome Pipeline (FastUniq - BWA - Picard/GATK/Samtools) for multiple samples
## INALL is an array of directory names, with one directory per sample
## SAMPLE_PATH is the pathway to the sample directories ${SAMPLE_PATH}/${INALL[$SAMPLE_NUMBER]}

#/************************* Paras need to be adjusted for diff samples
SAMPLE_PATH="/home/path"
OUT_PATH="/home/path"

SCRIPT_PATH="/home/Scripts"
TARGETS="/home/CCDS_20160421.txt"
Library_ID="Test_plus.26Aug2016" #e.g. {PI name}.{Batch Number} to identify the library
REFDIR="/home/hg38"

REF_FILE="hg38.fa"
SCRATCH_DIR="/home/pipe.temp"
SEQ_PLATFORM="ILLUMINA" # Valid values are: ILLUMINA, SOLID, LS454, HELICOS and PACBIO
COV_OUT_DIR_NAME="coverage" #coverage file output directory name. Just a name, it will be placed under the sample directory
DUP_FREE_BAM_DIR_NAME="dup_free_bam" # folder to ouput bam file without any GATK process, this can be input of freebayes
GATK_OUT_DIR_NAME="gatk" #gatk tuned alignment file output directory name
WRKGDIR_NAME="_no_QC" # tail of temp working directory under scratch dir, in case same sample running at the same using the same temp directory
JAVA_TMP_DIR_NAME="_no_QC" #tail of temp java working directory under scratch dir
#**************************/

if [ ! -d ${SCRATCH_DIR} ]; then
	mkdir ${SCRATCH_DIR}
else 
	echo "${SCRATCH_DIR} exists"	
fi

###############################

## Submitting jobs ##
INALL=( "${SAMPLE_PATH}"/SampleID_prefix )
for SampleFolder in "${INALL[@]}"; do
	SAMPLE_ID="${SampleFolder##*/}"
	       	JOB_ID2="bwa_${SAMPLE_ID}"
	       	JOB_IDc="cov_${SAMPLE_ID}"
		arr=("${SAMPLE_ID}" "${SAMPLE_PATH}" "${SCRIPT_PATH}" "$REFDIR" "${SCRATCH_DIR}" "$TARGETS" "${SEQ_PLATFORM}" "${Library_ID}" "${COV_OUT_DIR_NAME}" "${GATK_OUT_DIR_NAME}" "${DUP_FREE_BAM_DIR_NAME}" "${WRKGDIR_NAME}" "${JAVA_TMP_DIR_NAME}" "${FREE_DIR_NAME}" "${OUT_PATH}" "${REF_FILE}")
        arr2=("${SAMPLE_ID}" "${SAMPLE_PATH}" "${SCRIPT_PATH}" "$REFDIR" "$TARGETS" "${COV_OUT_DIR_NAME}" "${DUP_FREE_BAM_DIR_NAME}" "${OUT_PATH}" "${REF_FILE}" "${Library_ID}")
        echo "${SAMPLE_ID}"
		##Add separate submission for fastuniq.sh with larger memory allocation!
		qsub -N "${JOB_ID2}" -l "h_vmem=60G" -M "helen.griffin@ncl.ac.uk" ${SCRIPT_PATH}/BWA_MarkDups.sh "${arr[@]}"
		qsub -hold_jid "${JOB_ID2}" -N "${JOB_IDc}" -l "h_vmem=40G" -M "helen.griffin@ncl.ac.uk" ${SCRIPT_PATH}/Coverage.sh "${arr2[@]}"
done
