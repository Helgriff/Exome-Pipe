#! /bin/bash
#$ -cwd -V
#$ -j y
#$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job

echo $'\n'"["`date`"]: Job started."

## Add Modules ########################
module load $8 $9 ${10}
#######################################

SAMPLE_ID=$1
SAMPLE_PATH=$2
Q=$3
ANNOVAR_PATH=$4
REF_BUILD=$5
IntGENES=$6
INHOUSE_MAF=$7

INDIR="${SAMPLE_PATH}/${SAMPLE_ID}"
SAM13bVCF="${INDIR}/${SAMPLE_ID}_Concat_SourceBio_Q${Q}_Bcftools13_${REF_BUILD}.vcf"
S13bOUT="${INDIR}/${SAMPLE_ID}_Concat_SourceBio_Q${Q}_Bcftools13_${REF_BUILD}"

##Annovar
${ANNOVAR_PATH}/table_annovar.pl "${SAM13bVCF}" --outfile "${S13bOUT}" "${ANNOVAR_PATH}/humandb" \
--vcfinput \
--buildver ${REF_BUILD} \
--protocol knownGene,refGene,exac03,esp6500siv2_all,avsnp144,cosmic70,dbnsfp30a,dbnsfp31a_interpro,dbscsnv11,hrcr1,kaviar_20150923,nci60,mitimpact24 \
--operation g,g,f,f,f,f,f,f,f,f,f,f,f \
--nastring . \
--otherinfo

SAMPLE_S13bVCF="${S13bOUT}.${REF_BUILD}_multianno.vcf"

## Final Cleaning
rm ${S13bOUT}.avinput ${S13bOUT}.knownGene* ${S13bOUT}.refGene* ${S13bOUT}.log ${S13bOUT}.${REF_BUILD}_*_filtered ${S13bOUT}.${REF_BUILD}_*_dropped

echo $'\n'"["`date`"]: Annotation of ${SAMPLE_ID} variants is Complete!!"

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
