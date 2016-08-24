#! /bin/bash
#$ -cwd -V
#$ -j y
#$ -m e

set -e
res1=$(date +%s.%N) # use to calculate whole run time of the job

echo $'\n'"["`date`"]: Job started."

## Add Modules ########################
module load $4
#######################################

SAMPLE_ID=$1
SAMPLE_PATH=$2
INTGENES=$3
SCRIPTPATH=$5
FILE=$6
INHOUSE=$7
OMIM=$8

##Convert VCF to tab delimited text file using perl script
perl "${SCRIPTPATH}/VCF_to_Excel.pl" ${SAMPLE_PATH} ${FILE} ${SAMPLE_ID} ${INTGENES} ${INHOUSE} ${OMIM}

echo $'\n'"["`date`"]: VCF to Excel is complete!!"

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
