#!/bin/bash
#SBATCH --job-name=job7
#SBATCH --output=job7_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load Anaconda3; module list
source activate MetaPhlAn4; conda list

merge_metaphlan_tables.py \
    IBD_metaphlan/*_metaphlan.txt \
    > IBD_metaphlan.txt

merge_metaphlan_tables.py \
    LLD_metaphlan/*_metaphlan.txt \
    > LLD_metaphlan.txt

chmod 440 *_metaphlan.txt
