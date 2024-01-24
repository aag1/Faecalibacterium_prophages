#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=10gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE
#SBATCH --get-user-env=L



### download the latest version of the human genome
# https://github.com/biobakery/kneaddata ->
# "NCBI project page" ->
# "GRCh38.p14 (latest minor release) FTP" ->
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz



### build index
module purge; module load Bowtie2; module list

bowtie2-build \
    GCA_000001405.29_GRCh38.p14_genomic.fna.gz \
    GRCh38p14 \
    --threads ${SLURM_CPUS_PER_TASK}
