#!/bin/bash
#SBATCH --job-name=job7
#SBATCH --output=job7_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load seqtk; module list

seqtk trimfq -b 24 -e 22 PhageTerm_OUT/Tulp_contig.fasta | seqtk seq -l 80 | sed -E 's/^>.+$/>Tulp/' > Tulp.fasta

chmod 440 Tulp.fasta



for phage in 'Roos' 'Pioen' 'Aster' 'Lelie'
do

    seqtk seq -l 80 PhageTerm_OUT/${phage}_sequence.fasta | sed -E 's/^(>[^ ]+) .+$/\1/' > ${phage}.fasta

    chmod 440 ${phage}.fasta

done
