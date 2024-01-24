sessionInfo()



tab <- read.table('../Fprau_raw_reads/Fprau_strains.txt', sep = '\t', header = F, stringsAsFactors = F)

V <- tab[, 2]



DF <- data.frame(NULL)

for (s in V) {

    x <- gsub('-', '_', s)


    file <- paste0('ct2_results/', x, '_PRUNING_INFO_TABLE.tsv')

    if (!file.exists(file)) { next }

    t <- read.table(file, sep = '\t', header = T, stringsAsFactors = F)


    df <- data.frame(
        strain_id = s,
        contig_id = t$INPUT_PARENT_NAME,
        from = sapply(1:nrow(t), function (i) ifelse(t$CHROM_REMOVED[i] == 'True', t$LEFT_JUNCTION[i], 1)),
        to = sapply(1:nrow(t), function (i) ifelse(t$CHROM_REMOVED[i] == 'True', t$RIGHT_JUNCTION[i], t$PARENT_LENGTH[i])),
        stringsAsFactors = F
    )

    DF <- rbind(DF, df)

}



write.table(DF, sep = '\t', row.names = F, quote = F, file = 'Fprau_prophages_ct2.txt')
