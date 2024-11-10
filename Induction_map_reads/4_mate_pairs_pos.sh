#!/bin/bash
#SBATCH --job-name=job4
#SBATCH --output=job4_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load SAMtools; module list

while read -r phage_name strain contig coo1 coo2
do

    if [[ "${phage_name}" == "phage_name" ]] || [[ "${phage_name}" == "Tulp" ]] || [[ "${contig}" == "NODE_24_length_43119_cov_24.736044" ]]
    then
        continue
    fi


    pp_5p=${contig}':'${coo1}'-'$((${coo1} + 99))

    ARR=( $(awk -F'\t' -v strain=${strain} '($2 == strain) && ($7 == "Yes") {print $1}' ../Induction_raw_reads/Novogene_summary_table.txt) )

    for sample in "${ARR[@]}"
    do

        samtools view \
            MAPPED_READS/${sample}_mapped_to_${strain}.sorted.bam \
            ${pp_5p} \
            > ${phage_name}_5p_reads_${sample}.sam

    done

done < ../Induction_assembly/integr_sites_coo.txt



module purge; module load R; module list

Rscript mate_pairs_pos.R



chmod 440 *sam mate_pairs_pos.txt
