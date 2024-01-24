#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load R; module list

for x in 'VanEspen_2021' 'Shkoporov_2019'
do
    md5sum raw_reads_${x}/* | sed -E 's/ +/\t/' > md5sum_reads_${x}.txt

    Rscript check_md5.R ${x}
done

chmod 440 md5sum_reads_*.txt
