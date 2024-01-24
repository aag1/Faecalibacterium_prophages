#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=40gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L



awk \
    -F'\t' \
    '(($1 == $14) && ($15 ~ /g__Faecalibacterium/)) {print "../UHGG_database/fna_files/"$1".fna"}' \
    ../UHGG_database/genomes-all_metadata.tsv > genomes_to_dereplicate.txt

ls *_contigs.fasta >> genomes_to_dereplicate.txt



module purge; module load dRep; module list

dRep dereplicate \
    dRep_out \
    --genomes genomes_to_dereplicate.txt \
    -pa 0.9 -sa 0.95 -nc 0.30 -cm larger \
    --processors ${SLURM_CPUS_PER_TASK} \
    --skip_plots



module purge; module load R; module list

Rscript assign_taxo.R



chmod 440 genomes_to_dereplicate.txt Fprau_strains_taxo.txt
