sessionInfo()



tab <- read.table('viral_refseq_aa_nt_ids.txt', header = T, sep = '\t', stringsAsFactors = F)



idx <- which(tab$genome_id %in% tab$protein_id)   # polyprotein instead of genome

for (i in idx) {

    j <- which(tab$protein_id == tab$genome_id[i])

    tab$genome_id[i] <- tab$genome_id[j]

}



write.table(tab, row.names = F, quote = F, sep = '\t', file = 'viral_refseq_aa_nt_ids.2.txt')
