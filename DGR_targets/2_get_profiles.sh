#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# ------------------------------ pfam profiles ------------------------------ #
wget --no-verbose http://ftp.tuebingen.mpg.de/pub/protevo/toolkit/databases/hhsuite_dbs/pfama_v35.tar.gz

tar -xvzf pfama_v35.tar.gz


module purge; module load Anaconda3; module list
source activate HHsuite; conda list

ffindex_get Pfam35/pfama_hhm.ffdata Pfam35/pfama_hhm.ffindex PF19789.2 > DUF6273.hhm

ffindex_get Pfam35/pfama_hhm.ffdata Pfam35/pfama_hhm.ffindex PF13750.9 > Big_3_3.hhm


rm -rf pfama_v35.tar.gz Pfam35




# ------------------------------ pdb profiles ------------------------------ #
wget --no-verbose http://ftp.tuebingen.mpg.de/pub/protevo/toolkit/databases/hhsuite_dbs/pdb70_from_mmcif_2023-12-24.tar.gz

mkdir pdb70
tar -xvzf pdb70_from_mmcif_2023-12-24.tar.gz --directory pdb70


for x in '1YU0_A' '6HHK_A'
do
    ffindex_get pdb70/pdb70_hhm.ffdata pdb70/pdb70_hhm.ffindex ${x} > ${x}.hhm
done


chmod 440 *hhm
rm -rf pdb70_from_mmcif_2023-12-24.tar.gz pdb70
