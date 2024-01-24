#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L



base='crispr_spacers_uhgg_isolates'

cat ../UHGG_database/CRISPRidentify_out/*_spacer_dataset.fasta > ${base}.fasta



cd CRISPR_out

module purge; module load BLAST+; module list

for p in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie'
do

    makeblastdb \
        -in ../../Induction_assembly/${p}.fasta \
        -dbtype nucl \
        -out ${p}


    blastn \
        -task 'blastn-short' \
        -db ${p} \
        -query ../${base}.fasta \
        -evalue 1 \
        -max_target_seqs 1000000 \
        -num_threads ${SLURM_CPUS_PER_TASK} \
        -out ${base}_vs_${p}.out \
        -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'


    awk -F'\t' '$17 / $9 >= 0.95' ${base}_vs_${p}.out > ${base}_vs_${p}.spacer95match.txt


    if [[ -s ${base}_vs_${p}.spacer95match.txt ]]; then chmod 440 ${base}_vs_${p}.spacer95match.txt; else rm ${base}_vs_${p}.spacer95match.txt; fi

    rm ${base}_vs_${p}.out

    rm ${p}.n*

done



rm ../${base}.fasta
