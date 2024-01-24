#!/bin/bash
#SBATCH --job-name=job1.1
#SBATCH --output=job1.1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=6-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


# https://www.ebi.ac.uk/ena/browser/view/PRJNA545408
# 'Show Column Selection' --> add/remove columns --> Download report: TSV


ARR=( $(awk -F'\t' '$10 ~ /^Lon_Virome_/ {print $9}' PRJNA545408_Shkoporov_2019.tsv) )

cd raw_reads_Shkoporov_2019

for x in "${ARR[@]}"
do

    f1=$(echo ${x} | cut -d';' -f1)
    f2=$(echo ${x} | cut -d';' -f2)

    wget --no-verbose ${f1}
    wget --no-verbose ${f2}

done
