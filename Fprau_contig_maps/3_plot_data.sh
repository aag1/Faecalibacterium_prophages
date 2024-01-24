#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load R/4.2.2-foss-2022b; module list

ARR=( $(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt) )

for STRAIN_ID in "${ARR[@]}"
do
    Rscript plot_all.R ${STRAIN_ID}
done

mkdir PDF_PER_STRAIN; chmod 750 PDF_PER_STRAIN
chmod 440 *_contig_maps.pdf; mv *_contig_maps.pdf PDF_PER_STRAIN/



Rscript plot_selected.R

chmod 440 contigs_with_candidate_prophages.pdf



Rscript plot_example.R

chmod 440 CP11_example.pdf
