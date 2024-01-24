#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# ------------------------------ build dotplots ------------------------------ #
ARR1=( $(ls *fasta | sed -E 's/\.fasta$//') )
ARR2=( $(ls *hhm | sed -E 's/\.hhm$//') )

module purge; module load Anaconda3; module list
source activate HHsuite; conda list

for b1 in "${ARR1[@]}"
do

    for b2 in "${ARR1[@]}" "${ARR2[@]}"
    do

        if [[ -f ${b2}.hhm ]]
        then
            f1=HHpred_output/hhpred_full_${b1}.a3m
            f2=${b2}.hhm
        else
            f1=${b1}.fasta
            f2=${b2}.fasta
        fi


        hhalign \
            -i ${f1} \
            -t ${f2} \
            -alt 100 \
            -norealign \
            -o ${b1}___${b2}.hhr \
            -atab ${b1}___${b2}.atab


        grep -A1 '^>' ${b1}___${b2}.hhr |
        grep '^Probab=' |
        cut -d' ' -f1 |
        sed 's/^Probab=//' \
        > ${b1}___${b2}.probab


        grep -v '^>' ${b1}___${b2}.atab |
        grep -v '^missing dssp' |
        sed -E 's/^ +//' |
        sed -E 's/ +/\t/g' \
        > ${b1}___${b2}.atab1

    done

done

mkdir HHalign_output; chmod 750 HHalign_output
chmod 440 *___*; mv *___* HHalign_output/




# ------------------------------ plot dotplots ------------------------------ #
module purge; module load R/4.2.2-foss-2022b; module list

Rscript dotplots.R

chmod 440 DGR_targets_dotplots.pdf
