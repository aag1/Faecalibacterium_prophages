#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L



ls *.fastq.gz |
sed -E 's/^(.+)(_S[0-9]+)_L001_R[12]_001\.fastq\.gz$/\1\t\1\2/' |
sort | uniq > Fprau_strains.txt

chmod 440 Fprau_strains.txt



echo -e 'strain_id\tfile_R1\tmd5_R1\tfile_R2\tmd5_R2' > Fprau_raw_reads_md5.txt

while read x1 x2
do
    f1=${x2}_L001_R1_001.fastq.gz
    f2=${x2}_L001_R2_001.fastq.gz

    m1=$(md5sum ${f1} | cut -d' ' -f1)
    m2=$(md5sum ${f2} | cut -d' ' -f1)

    echo -e ${x1}'\t'${f1}'\t'${m1}'\t'${f2}'\t'${m2} >> Fprau_raw_reads_md5.txt
done < Fprau_strains.txt

chmod 440 Fprau_raw_reads_md5.txt
