#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



wget --no-verbose ftp://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER

wget \
    --no-verbose \
    --recursive \
    --no-directories \
    ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/



for f in *gz; do gzip -d ${f}; done



module purge; module load Perl BioPerl; module list

perl get_viral_refseq_taxo.pl viral.1.genomic.gbff > viral_refseq_taxo.txt

perl get_viral_refseq_aa_nt_ids.pl viral.1.protein.gpff > viral_refseq_aa_nt_ids.txt

module purge; module load R; module list

Rscript fix_ids.R



chmod 440 RELEASE_NUMBER viral*
