sessionInfo()



tab <- read.table('../Fprau_raw_reads/Fprau_strains.txt', sep = '\t', header = F, stringsAsFactors = F)

V <- tab[, 2]



DF <- data.frame(NULL)

for (s in V) {

    file <- paste0('vs2_results/', s, '_final-viral-boundary.tsv')

    t <- read.table(file, sep = '\t', header = T, stringsAsFactors = F)

    if (nrow(t) == 0) { next }


    df <- data.frame(
        strain_id = s,
        contig_id = t$seqname,
        from = t$trim_bp_start,
        to = t$trim_bp_end,
        stringsAsFactors = F
    )

    DF <- rbind(DF, df)

}



write.table(DF, sep = '\t', row.names = F, quote = F, file = 'Fprau_prophages_vs2.txt')
