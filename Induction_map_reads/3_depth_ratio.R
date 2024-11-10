sessionInfo()



P <- read.table('../CPs_refine_coo/prophages_refined_coo.txt', sep = '\t', header = T, stringsAsFactors = F)

P <- P[P$candidate_prophage %in% c('CP3', 'CP8', 'CP9', 'CP10'), c(1:3, 7:8)]

p <- data.frame(
    candidate_prophage = 'Roos_fragm',
    strain_id = 'FM8_S2',
    contig_id = 'NODE_1_length_167466_cov_21.860236',
    new_from = 157269,
    new_to = 167466,
    stringsAsFactors = F
)

P <- rbind(P, p)



Q <- read.table('../Induction_raw_reads/Novogene_summary_table.txt', sep = '\t', header = T, stringsAsFactors = F)

Q <- Q[Q$Sequencing == 'Yes', ]



DF <- data.frame(NULL)

for (i in 1:nrow(P)) {

    idx <- which(Q$Strain == P$strain_id[i])


    df <- P[rep(i, length(idx)), ]

    df <- cbind(df, Q[idx, c(1, 3:4)])


    f <- paste0('../Fprau_map_reads/BED_FILES/', P$strain_id[i], '_contigs.bed')

    b <- read.table(f, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)


    df$depth_ratio <- NA

    for (j in 1:nrow(df)) {

        f <- paste0('MAPPED_READS/', df$Sample.Name[j], '_mapped_to_', df$strain_id[j], '.depth.txt')

        t <- read.table(f, sep = '\t', header = F, stringsAsFactors = F)


        k <- which(t[, 1] == df$contig_id[j])
        coo <- t[k, 2]
        cov <- t[k, 3]

        len <- b[df$contig_id[j], 2]

        COV <- rep(0, len)
        COV[coo] <- cov


        pp <- df$new_from[j]:df$new_to[j]

        df$depth_ratio[j] <- mean( COV[ pp ] ) / mean( COV[ -pp ] )

    }


    DF <- rbind(DF, df)

}



DF$depth_ratio <- round(DF$depth_ratio)

DF <- DF[, c('candidate_prophage', 'DNA_type', 'Treatment', 'depth_ratio')]

write.table(DF, sep = '\t', row.names = F, quote = F, file = 'depth_ratio.txt')
