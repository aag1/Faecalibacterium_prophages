#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-23
#SBATCH --export=NONE
#SBATCH --get-user-env=L



STRAIN_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../Fprau_raw_reads/Fprau_strains.txt | cut -d$'\t' -f2)



module purge; module load SPAdes; module list

spades.py \
    --isolate \
    -1 ../Fprau_clean_reads/${STRAIN_ID}_kneaddata_paired_1.fastq.gz \
    -2 ../Fprau_clean_reads/${STRAIN_ID}_kneaddata_paired_2.fastq.gz \
    -o ${STRAIN_ID}_spades_out \
    --threads ${SLURM_CPUS_PER_TASK}



for x in 'contigs.fasta' 'scaffolds.fasta' 'assembly_graph.fastg' 'contigs.paths' 'scaffolds.paths' 'spades.log'
do
    mv ${STRAIN_ID}_spades_out/${x} ${STRAIN_ID}_${x}
done

rm -rf ${STRAIN_ID}_spades_out
chmod 440 ${STRAIN_ID}_*
