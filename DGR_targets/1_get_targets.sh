#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




module purge; module load seqtk; module list

for p in 'Tulp' 'Roos' 'Pioen' 'Lelie' 'CP2' 'CP11'
do

    if [[ "${p}" == "CP"* ]]
    then
        f=../CPs_proteomes/${p}_proteome.fasta
    else
        f=../Phages_proteomes/${p}_proteome.fasta
    fi


    TARGETS=( $(sed '1d' ../DGR_detection/${p}_DGR_repeats.txt | cut -d$'\t' -f9) )

    for t in "${TARGETS[@]}"
    do

        echo ${t} > target.id

        seqtk subseq \
            -l 80 \
            ${f} \
            target.id \
            > ${t}.fasta

        rm target.id

        chmod 440 ${t}.fasta

    done

done
