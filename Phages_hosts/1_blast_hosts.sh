#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-220
#SBATCH --export=NONE
#SBATCH --get-user-env=L




N=$(awk -F'\t' '$2 == "Isolate"' ../UHGG_database/genomes-all_metadata.tsv | wc -l)

from=$(((${SLURM_ARRAY_TASK_ID} - 1) * 50 + 1))

to=$((${SLURM_ARRAY_TASK_ID} * 50))

if [[ "${to}" -gt "${N}" ]]; then to=${N}; fi

ARR=( $(awk -F'\t' '$2 == "Isolate" {print $1}' ../UHGG_database/genomes-all_metadata.tsv | sed -n "${from},${to}p") )




cd BLASTN_out

module purge; module load BLAST+; module list

for b in "${ARR[@]}"
do

    echo 'Working with '${b}' genome...'


    makeblastdb \
        -in ../../UHGG_database/fna_files/${b}.fna \
        -dbtype nucl \
        -out ${b}


    for p in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie'
    do

        blastn \
            -task 'blastn' \
            -db ${b} \
            -query ../../Induction_assembly/${p}.fasta \
            -evalue 0.001 \
            -num_threads ${SLURM_CPUS_PER_TASK} \
            -out ${p}_vs_${b}_blastn.out \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'


        awk -F'\t' '$10 >= 95' ${p}_vs_${b}_blastn.out > ${p}_vs_${b}_blastn.qcovs95.out


        if [[ -s ${p}_vs_${b}_blastn.qcovs95.out ]]; then chmod 440 ${p}_vs_${b}_blastn.qcovs95.out; else rm ${p}_vs_${b}_blastn.qcovs95.out; fi

        rm ${p}_vs_${b}_blastn.out

    done


    rm ${b}.n*

done
