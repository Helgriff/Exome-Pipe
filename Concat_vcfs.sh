#! /bin/bash
#$ -cwd -V
#$ -j y
#$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job

echo $'\n'"["`date`"]: Job started."

## Add Modules ########################
module load $5
#######################################

SAMPLE_ID=$1
SAMPLE_PATH=$2
INDIR="${SAMPLE_PATH}/${SAMPLE_ID}"
Q=$3
REF_BUILD=$4

## Create list file of vcf files
VCFLIST="${INDIR}/${SAMPLE_ID}_vcfs.txt"
ls "${INDIR}/${SAMPLE_ID}_SourceBio_Q${Q}_Bcftools13_${REF_BUILD}_chr"*".vcf.gz" > "${VCFLIST}"

## Merge vcfs with bcftools
OUTVCF="${INDIR}/${SAMPLE_ID}_Concat_SourceBio_Q${Q}_Bcftools13_${REF_BUILD}.vcf"
bcftools concat -o "${OUTVCF}" -O v -f "${VCFLIST}"

echo $'\n'"["`date`"]: Merging of multiple vcfs is Complete!!"

# runtime calculation
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds
echo "exit status $?"
