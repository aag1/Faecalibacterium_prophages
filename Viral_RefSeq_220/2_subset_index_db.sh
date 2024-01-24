#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



base='viral_refseq_complete_above10kb_noRiboviria'



awk \
    -F'\t' \
    '(NR > 1) && ($2 ~ /complete genome/) && ($3 > 10000) && ($4 !~ /Riboviria/) {print $1}' \
    viral_refseq_taxo.txt > ${base}.ids



module purge; module load seqtk; module list

seqtk subseq \
    -l 80 \
    viral.1.1.genomic.fna \
    ${base}.ids > ${base}.fasta



chmod 440 ${base}*



module purge; module load BLAST+; module list

makeblastdb \
    -in ${base}.fasta \
    -dbtype nucl \
    -out ${base}
