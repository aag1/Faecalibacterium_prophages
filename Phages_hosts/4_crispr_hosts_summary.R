sessionInfo()



tab <- read.table('../UHGG_database/genomes-all_metadata.tsv', sep = '\t', header = T, stringsAsFactors = F)



DF <- data.frame(NULL)

for (p in c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie')) {

    f <- paste0('CRISPR_out/crispr_spacers_uhgg_isolates_vs_', p, '.spacer95match.txt')

    if (!file.exists(f)) { next }



    t <- read.table(f, sep = '\t', header = F, stringsAsFactors = F)

    B <- unlist(lapply(strsplit(t[, 1], '_'), function (v) v[1]))

    B <- unique(B)



    X <- sapply(B, function (b) {

        x <- tab$Lineage[tab$Genome == b]
        x <- strsplit(x, ';')[[1]]
        x <- rev(x)[1]
        return(x)

    })



    df <- data.frame(
        phage = p,
        host_id = B,
        host_taxo = X,
        stringsAsFactors = F
    )

    DF <- rbind(DF, df)

}



write.table(DF, sep = '\t', quote = F, row.names = F, file = 'crispr_hosts_summary.txt')
