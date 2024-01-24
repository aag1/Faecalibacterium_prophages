sessionInfo()



for (K in c('LLD', 'IBD')) {

    file <- ifelse(K == 'LLD', '../../LLD_reads/LLD_raw_reads_number_sele.txt', '../../IBD_reads/selected_ibd_samples.txt')

    samples <- read.table(file, sep = '\t', header = F, stringsAsFactors = F)[, 1]


    DF <- data.frame(NULL)

    for (s in samples) {

        t <- read.table(paste0(K, '_mapped_reads/', s, '_coverage.txt'), sep = '\t', header = F, stringsAsFactors = F)

        v <- setNames(t[, 7], t[, 1])

        df <- as.data.frame(t(v), stringsAsFactors = F)
        rownames(df) <- s

        DF <- rbind(DF, df)

    }


    write.table(DF, quote = F, sep = '\t', file = paste0(K, '_five_phages_cov.txt'))

}
