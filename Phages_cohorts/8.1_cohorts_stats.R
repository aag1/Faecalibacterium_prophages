sessionInfo()



DF <- data.frame(NULL)

for (K in c('LLD', 'IBD')) {

    ### phages
    f <- paste0(K, '_five_phages_cov.txt')
    t <- read.table(f, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

    df <- as.data.frame(t >= 0.75, stringsAsFactors = F)
    df$Cohort <- K



    ### host
    q <- read.table(paste0(K, '_metaphlan.txt'), sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

    idx <- grep('\\|g__Faecalibacterium$', rownames(q))

    df$HostMetaphlan <- sapply(rownames(df), function (x) {

        s <- paste0(ifelse(K == 'LLD', 'X', ''), x, '_metaphlan')
        q[idx, s]

    })



    ### reads
    f <- paste0('../../', K, '_reads/', K, '_clean_reads_multiqc_report_data/mqc_fastqc_sequence_counts_plot_1.txt')
    t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)

    df$NumCleanReads <- sapply(rownames(df), function (x) {

        idx <- which(t$Sample %in% paste0(x, '_kneaddata_paired_', 1:2))
        sum(t[idx, 2:3])

    })



    DF <- rbind(DF, df[, c(6, 8, 1:5, 7)])

}



write.table(DF, sep = '\t', quote = F, file = 'LLD_IBD_reads_phages_host.txt')
