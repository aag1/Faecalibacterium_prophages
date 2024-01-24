sessionInfo()



tab <- read.table('../Fprau_raw_reads/Fprau_strains.txt', sep = '\t', header = F, stringsAsFactors = F)

V <- tab[, 2]



DF <- data.frame(NULL)

for (s in V) {

    file <- paste0('geNomad_results/', s, '/', s, '_contigs_summary/', s, '_contigs_virus_summary.tsv')

    t <- read.table(file, sep = '\t', header = T, stringsAsFactors = F)

    if (nrow(t) == 0) { next }


    df <- data.frame(
        strain_id = s,
        contig_id = sub('\\|.+$', '', t$seq_name),
        from = sapply(1:nrow(t), function (i) ifelse(is.na(t$coordinates[i]), 1, strsplit(t$coordinates[i], '-')[[1]][1])),
        to = sapply(1:nrow(t), function (i) ifelse(is.na(t$coordinates[i]), t$length[i], strsplit(t$coordinates[i], '-')[[1]][2])),
        stringsAsFactors = F
    )

    DF <- rbind(DF, df)

}



write.table(DF, sep = '\t', row.names = F, quote = F, file = 'Fprau_prophages_geNomad.txt')
