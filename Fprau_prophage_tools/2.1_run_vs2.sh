#!/bin/bash
#SBATCH --job-name=job2.1
#SBATCH --output=job2.1_%A_%a.out
#SBATCH --mem=20gb
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-23
#SBATCH --export=NONE
#SBATCH --get-user-env=L



STRAIN_ID=$(cut -d$'\t' -f2 ../Fprau_raw_reads/Fprau_strains.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo '# ------------------------------ WORKING WITH '${STRAIN_ID}' STRAIN ------------------------------ #'

cd vs2_results



module purge; module load Anaconda3; module list
source activate vs2; conda list

virsorter run \
    --working-dir ${STRAIN_ID} \
    --seqfile ../../Fprau_assemblies/${STRAIN_ID}_contigs.fasta \
    --include-groups 'dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae' \
    --jobs ${SLURM_CPUS_PER_TASK} \
    all



for x in 'final-viral-combined.fa' 'final-viral-score.tsv' 'final-viral-boundary.tsv'
do
    mv ${STRAIN_ID}/${x} ${STRAIN_ID}_${x}
    chmod 440 ${STRAIN_ID}_${x}
done

rm -rf ${STRAIN_ID}/
