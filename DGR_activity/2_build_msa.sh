#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A_%a.out
#SBATCH --mem=4gb
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=1-3%1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]; then p='Roos'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 2 ]]; then p='CP2'; fi
if [[ ${SLURM_ARRAY_TASK_ID} -eq 3 ]]; then p='CP11'; fi




if [[ "${p}" == "CP"* ]]
then
    f=../CPs_refine_coo/${p}_refined.fna
else
    f=../Induction_assembly/${p}.fasta
fi




# ------------------------------ cognate genomes ------------------------------ #
module purge; module load seqtk; module list

if [[ -f ${p}_cognate_genomes.fasta ]]; then rm -f ${p}_cognate_genomes.fasta; fi

for db in 'viral_refseq' 'IMGVR'
do

    awk -F'\t' '$10 >= 95 {print $2}' ${p}_vs_${db}.txt | sort | uniq > ${p}_cognate_genomes_${db}.ids

    if [[ ! -s ${p}_cognate_genomes_${db}.ids ]]; then rm ${p}_cognate_genomes_${db}.ids; continue; fi


    if [[ "${db}" == "viral_refseq" ]]
    then
        dbF=../Viral_RefSeq_220/viral_refseq_complete_above10kb_noRiboviria.fasta
    else
        dbF=../IMG_VR_2022-12-19_7.1/IMGVR_above10kb_DTR_noRefSeq.fasta
    fi


    seqtk subseq \
        -l 80 \
        ${dbF} \
        ${p}_cognate_genomes_${db}.ids \
        >> ${p}_cognate_genomes.fasta

done




# ------------------------------ cognate genomes DTRs ------------------------------ #
perl ../Induction_assembly/identify_circular.pl \
    ${p}_cognate_genomes.fasta \
    ${p}_cognate_genomes_DTR.txt \
    10




# ------------------------------ reference genome 5'-end ------------------------------ #
echo -e ${p}'\t0\t500' > ${p}_5p.bed

seqtk subseq -l 80 ${f} ${p}_5p.bed > ${p}_5p.fasta

rm ${p}_5p.bed




# ------------------------------ find new 5'-ends ------------------------------ #
module purge; module load BLAST+; module list

makeblastdb \
    -in ${p}_cognate_genomes.fasta \
    -dbtype nucl \
    -out ${p}_cognate_genomes

blastn \
    -task 'blastn' \
    -query ${p}_5p.fasta \
    -db ${p}_cognate_genomes \
    -num_threads ${SLURM_CPUS_PER_TASK} \
    -out ${p}_5p_vs_cognate_genomes.txt \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend qlen qcovs sstart send slen sstrand evalue bitscore nident'

rm ${p}_5p.fasta ${p}_cognate_genomes.n*




# ------------------------------ apply new 5'-ends ------------------------------ #
module purge; module load R/4.2.2-foss-2022b; module list

Rscript new_5prime.R \
    ${f} \
    ${p}_cognate_genomes.fasta \
    ${p}_cognate_genomes_DTR.txt \
    ${p}_5p_vs_cognate_genomes.txt \
    ${p}_cognate_genomes_for_msa.fasta




# ------------------------------ align genomes ------------------------------ #
module purge; module load MAFFT; module list

mafft \
    --nuc \
    --maxiterate 1000 \
    --thread ${SLURM_CPUS_PER_TASK} \
    ${p}_cognate_genomes_for_msa.fasta \
    > ${p}_cognate_genomes_msa.fasta

rm ${p}_cognate_genomes_for_msa.fasta
chmod 440 ${p}_cognate_genomes* ${p}_5p_vs_cognate_genomes.txt
