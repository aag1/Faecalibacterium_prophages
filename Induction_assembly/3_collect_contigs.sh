#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



module purge; module load seqtk; module list

if [[ -f induced_phages.fasta ]]; then rm -f induced_phages.fasta; fi
if [[ -f induced_phages.txt ]]; then rm -f induced_phages.txt; fi
echo -e 'contig_name\tcontig_len\tcontig_ends\tsample_name\tstrain\ttreatment' > induced_phages.txt

while IFS=$'\t' read -r Sample_Name Strain DNA_type Treatment Concentration Libr_QC Sequencing
do

    if [[ "${DNA_type}" == "bacterial DNA" ]] || [[ "${Sequencing}" != "Yes" ]]; then continue; fi


    file=geNomad_OUT/${Sample_Name}/${Sample_Name}_contigs_summary/${Sample_Name}_contigs_virus_summary.tsv

    if [[ ! -f ${file} ]]; then continue; fi


    while IFS=$'\t' read -r seq_name length topology coordinates n_genes genetic_code virus_score fdr n_hallmarks marker_enrichment taxonomy
    do

        if [[ "${seq_name}" == "seq_name" ]]; then continue; fi

        echo -e ${Sample_Name}'_'${seq_name}'\t'${length}'\t'${topology}'\t'${Sample_Name}'\t'${Strain}'\t'${Treatment} >> induced_phages.txt

        echo ${seq_name} > sele.id

        seqtk subseq \
            -l 80 \
            ${Sample_Name}_contigs.fasta \
            sele.id |
        sed -E "s/^>(.+)$/>${Sample_Name}_\1/" >> induced_phages.fasta

    done < ${file}

done < ../Induction_raw_reads/Novogene_summary_table.txt



perl identify_circular.pl induced_phages.fasta induced_phages_DTR.txt 10



rm sele.id
chmod 440 induced_phages*
