#! /bin/bash

#*********************************************************************************************************************************************************##
## Script to Launch Variant Calling accross multiple bam files, creating per Chr Multi-Sample VCF Files, Concatenate per Chr files and Annotation by Annovar																##
#*********************************************************************************************************************************************************##

##N.B. Bam files moved (mv) into directory (${SAMPLE_PATH}/${SAMPLE_ID}/) and $BAMLIST created (ls ${SAMPLE_PATH}/${SAMPLE_ID}/*${REF_BUILD}.bam > $BAMLIST)
##Bam lists restricted to groups of < approx. 800 samples ($GROUP_$BAMLIST) ... otherwise jobs now fail on cluster, memory restrictions!?!

SAMPLE_ID="BAMs"

INALL=(chrM chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY) #Array containing group names

SAMPLE_PATH="/home/Aligned_reads"

BAMLIST="${SAMPLE_PATH}/${SAMPLE_ID}_bams.txt"
ls ${SAMPLE_PATH}/${SAMPLE_ID}/*.bam > ${BAMLIST}

SCRIPT_PATH="/home/Scripts"
REF_BUILD="hg38"
REF_FILE="/home/hg38/hg38.fa"
Q=25
MaxD=200
IntGENES="/home/hg38_MitochondrialGenes.txt"
INHOUSE="/home/InHouse_MAFs_378exomes_bcftools13.txt"
PATHS="/home/CPDB_pathways_genes.tab" ##tab delim file of pathways with genelists downloaded from cpdb site
BMQ="all.q,bigmem.q"
LMEM="h_vmem=60G"
EMAIL="email@place.com"
SMT13="apps/samtools/1.3"
PRL="apps/perl/5.18.2"
VTL="apps/vcftools/0.1.12b"
ANV="/home/annovar_jan2016"
OMIM="/home/Ensembl_OMIM_AllGeneNames_withChrPos_hg38.txt"

#1# Call Variants using Samtools-bcftools in batches of < ~800 sampleIDs (now fails with larger sample numbers!!)
JOBHOLD="NONE"
for GROUP in "${INALL[@]}"; do
JOB_ID1="v_${GROUP}_${SAMPLE_ID}"
arr13=("${SAMPLE_ID}" "${SAMPLE_PATH}" "$Q" "$MaxD" "${REF_BUILD}" "${REF_FILE}" "${BAMLIST}" "$SMT13" "$GROUP")
qsub -N "${JOB_ID1}"  -l "$LMEM" -M "$EMAIL" ${SCRIPT_PATH}/Samtools1-3Bcftools.sh "${arr13[@]}"
JOBHOLD="${JOBHOLD},${JOB_ID1}"
done

#2# bgzip/tabix each vcf
JOBHOLD2="NONE"
for GROUP in "${INALL[@]}"; do
JOB_ID2="tb_${GROUP}_${SAMPLE_ID}"
arr2=("${SAMPLE_ID}" "${SAMPLE_PATH}" "$Q" "${REF_BUILD}" "$SMT13" "$GROUP")
qsub -hold_jid "${JOBHOLD}" -N "${JOB_ID2}"  -l "$LMEM" -M "$EMAIL" ${SCRIPT_PATH}/Bgzip_tabix.sh "${arr2[@]}"
JOBHOLD2="${JOBHOLD2},${JOB_ID2}"
done

#3# Concatenate all group vcfs
JOB_ID3="cc_${SAMPLE_ID}"
arr3=("${SAMPLE_ID}" "${SAMPLE_PATH}" "$Q" "${REF_BUILD}" "$SMT13")
qsub -hold_jid "${JOBHOLD},${JOBHOLD2}" -N "${JOB_ID3}"  -l "$LMEM" -M "$EMAIL" ${SCRIPT_PATH}/Concat_vcfs.sh "${arr3[@]}"

#5# Annotate Variant Calls for one (concatenated, multisample) vcf file
JOB_ID5="Annovar_${SAMPLE_ID}"
arr5=("${SAMPLE_ID}" "${SAMPLE_PATH}" "$Q" "$ANV" "${REF_BUILD}" "$PRL" "$VTL" "$SMT")
qsub -hold_jid "${JOBHOLD},${JOBHOLD2},${JOB_ID2},${JOB_ID2b},${JOB_ID3},${JOB_ID3b},${JOB_ID4}" -N "${JOB_ID5}" -q "$BMQ" -l "$LMEM" -M "$EMAIL" ${SCRIPT_PATH}/Annovar.sh "${arr5[@]}"

#6# Get per Gene, per Pathway counts of variants
JOB_ID6="Count_${SAMPLE_ID}"
FILE="${SAMPLE_ID}_Concat_SourceBio_Q${Q}_Bcftools13_${REF_BUILD}.${REF_BUILD}_multianno.vcf"
arr6=("${SAMPLE_ID}" "${SAMPLE_PATH}" "$Q" "${REF_BUILD}" "${IntGENES}" "$PRL" "$PATHS" "${SCRIPT_PATH}" "${FILE}" "${INHOUSE}")
qsub -hold_jid "${JOBHOLD},${JOBHOLD2},${JOB_ID2},${JOB_ID3},${JOB_ID4},${JOB_ID5}" -N "${JOB_ID6}" -q "$BMQ" -l "$LMEM" -M "$EMAIL" ${SCRIPT_PATH}/Get_Counts.sh "${arr6[@]}"

#7# Use VCF_to_Excel.pl to produce tab delim file with desired headers for exonic/splicing but not synonymous variants
JOB_ID7="ConV_${SAMPLE_ID}"
FILE="${SAMPLE_ID}_Q${Q}_Bcftools13_${REF_BUILD}_ProteinAlteringVariants.vcf"
arr7=("${SAMPLE_ID}" "${SAMPLE_PATH}" "${IntGENES}" "$PRL" "${SCRIPT_PATH}" "${FILE}" "${INHOUSE}" "${OMIM}")
qsub -hold_jid "${JOBHOLD},${JOBHOLD2},${JOB_ID2},${JOB_ID3},${JOB_ID4},${JOB_ID5},${JOB_ID6}" -N "${JOB_ID7}" -q "$BMQ" -l "$LMEM" -M "$EMAIL" ${SCRIPT_PATH}/VCF_to_Excel_run.sh "${arr7[@]}"
