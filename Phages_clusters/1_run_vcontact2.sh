#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# ------------------------------ collect proteomes ------------------------------ #
if [[ -f gene_to_genome.csv ]]; then rm -f gene_to_genome.csv; fi

if [[ -f proteomes.fasta ]]; then rm -f proteomes.fasta; fi

echo 'protein_id,contig_id,keywords' > gene_to_genome.csv



### induced phages
for phage in 'Tulp' 'Roos' 'Pioen' 'Aster' 'Lelie'
do

    awk -F'\t' -v phage="${phage}" '{print $1","phage",None_provided"}' ../Phages_proteomes/${phage}_proteome.txt >> gene_to_genome.csv

    cat ../Phages_proteomes/${phage}_proteome.fasta >> proteomes.fasta

done



### non-induced CPs
for i in 2 4 5 6 7 11
do

    CP='CP'${i}

    awk -F'\t' -v CP="${CP}" '{print $1","CP",None_provided"}' ../CPs_proteomes/${CP}_proteome.txt >> gene_to_genome.csv

    cat ../CPs_proteomes/${CP}_proteome.fasta >> proteomes.fasta

done




# ------------------------------ run vConTACT2 ------------------------------ #
module purge; module load Miniconda3; module list
source activate vConTACT2; conda list

vcontact2 \
	--raw-proteins proteomes.fasta \
	--proteins-fp gene_to_genome.csv \
	--db 'ProkaryoticViralRefSeq211-Merged' \
	--output-dir vcontact2_output \
	--c1-bin ${HOME}/SOFTWARE/cluster_one-1.0.jar \
	--threads ${SLURM_CPUS_PER_TASK}

conda deactivate

for f in 'c1.ntw' 'genome_by_genome_overview.csv'; do mv vcontact2_output/${f} .; chmod 440 ${f}; done
rm -rf gene_to_genome.csv proteomes.fasta vcontact2_output
