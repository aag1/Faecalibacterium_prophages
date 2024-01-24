#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load Anaconda3; module list
conda activate mgcod; conda list

for phage in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie'
do

    if [[ "${phage}" == "Tulp" ]]
    then

        mgcod.py \
            --path_to_genome ../Induction_assembly/${phage}.fasta \
            --isoforms \
            --path_to_mgm_predictions ${phage}_mgcod \
            --path_to_output ${phage}_mgcod \
            --logfile ${phage}_mgcod/${phage}_mgcod.log \
            --amino_acids \
            --delete \
            --verbose

    else

        mgcod.py \
            --path_to_genome ../Induction_assembly/${phage}.fasta \
            --isoforms \
            --circular \
            --path_to_mgm_predictions ${phage}_mgcod \
            --path_to_output ${phage}_mgcod \
            --logfile ${phage}_mgcod/${phage}_mgcod.log \
            --amino_acids \
            --delete \
            --verbose

    fi


    awk \
        -F'\t' \
        -v OFS='\t' \
        '($1 !~ /^#/) && ($1 !~ /^$/) && ($6 != ".") {print $9,$4,$5,$7}' \
       ${phage}_mgcod/${phage}.gff |
        sed -E "s/^gene_id ([0-9\.]+);[^\t]+(\t.+)$/${phage}_\1\2/" \
        > ${phage}_proteome.txt


    sed \
        -E "s/^>.+ gene_id ([0-9\.]+); .+$/>${phage}_\1/" \
        ${phage}_mgcod/proteins_aa_${phage}.fasta \
        > ${phage}_proteome.fasta

done



grep 'genetic code ' *_mgcod/*log



mkdir Mgcod_OUT; chmod 750 Mgcod_OUT; mv *_mgcod Mgcod_OUT
chmod 440 *_proteome.txt *_proteome.fasta
