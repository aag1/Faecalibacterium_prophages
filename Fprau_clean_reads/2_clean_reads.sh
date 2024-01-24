#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --array=1-23
#SBATCH --export=NONE
#SBATCH --get-user-env=L



STRAIN_ID=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../Fprau_raw_reads/Fprau_strains.txt | cut -d$'\t' -f2)



module purge; module load Anaconda3; module list
source activate KneadData; conda list

trimmomatic_path=${HOME}/SOFTWARE/Trimmomatic-0.33

kneaddata \
    --input ../Fprau_raw_reads/${STRAIN_ID}_L001_R1_001.fastq.gz \
    --input ../Fprau_raw_reads/${STRAIN_ID}_L001_R2_001.fastq.gz \
    --reference-db ../human_genome_ref/ \
    --output-prefix ${STRAIN_ID}_kneaddata \
    --output ${STRAIN_ID}_kneaddata_out \
    --trimmomatic ${trimmomatic_path} \
    --trimmomatic-options "ILLUMINACLIP:${trimmomatic_path}/adapters/NexteraPE-PE.fa:2:30:10:1:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
    --sequencer-source NexteraPE \
    --bypass-trf \
    --reorder \
    --threads ${SLURM_CPUS_PER_TASK}

conda deactivate



mv ${STRAIN_ID}_kneaddata_out/${STRAIN_ID}_kneaddata_paired_1.fastq .
mv ${STRAIN_ID}_kneaddata_out/${STRAIN_ID}_kneaddata_paired_2.fastq .
mv ${STRAIN_ID}_kneaddata_out/${STRAIN_ID}_kneaddata.log .

chmod 440 ${STRAIN_ID}_kneaddata.log
rm -rf ${STRAIN_ID}_kneaddata_out

gzip ${STRAIN_ID}_kneaddata_paired_1.fastq
gzip ${STRAIN_ID}_kneaddata_paired_2.fastq
