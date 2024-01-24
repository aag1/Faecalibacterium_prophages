#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




module purge; module load R/4.2.2-foss-2022b; module list

for p in 'Tulp' 'Roos' 'Pioen' 'Lelie' 'CP2' 'CP11'
do

    Rscript get_repeats_seq.R ${p}

done




module purge; module load MAFFT; module list

for p in 'Tulp' 'Roos' 'Pioen' 'Lelie' 'CP2' 'CP11'
do

    mafft \
        --nuc \
        --maxiterate 1000 \
        ${p}_DGR_repeats_seq.fasta > ${p}_DGR_repeats_ali.fasta

done




module purge; module load R/4.2.2-foss-2022b; module list

Rscript plot_DGR_repeats.R




if [[ -d ALI_out ]]; then rm -rf ALI_out; fi
mkdir ALI_out; chmod 750 ALI_out
mv *_DGR_repeats_*.fasta *_DGR_repeats_summary.txt ALI_out/
chmod 440 ALI_out/* DGR_repeats_ali.pdf
