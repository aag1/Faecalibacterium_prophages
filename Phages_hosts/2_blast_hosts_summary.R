sessionInfo()



tab <- read.table('../UHGG_database/genomes-all_metadata.tsv', sep = '\t', header = T, stringsAsFactors = F)



DF <- data.frame(NULL)

for (p in c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie')) {

    files <- list.files(path = 'BLASTN_out/', pattern = paste0('^', p, '_vs_.+_blastn.qcovs95.out'))

    if (length(files) == 0) { next }

    for (f in files) {

        b <- strsplit(f, '_')[[1]][3]

        x <- tab$Lineage[tab$Genome == b]
        x <- strsplit(x, ';')[[1]]
        x <- rev(x)[1]

        df <- data.frame(
            phage = p,
            host_id = b,
            host_taxo = x,
            stringsAsFactors = F
        )

        DF <- rbind(DF, df)

    }

}



write.table(DF, sep = '\t', quote = F, row.names = F, file = 'blast_hosts_summary.txt')
