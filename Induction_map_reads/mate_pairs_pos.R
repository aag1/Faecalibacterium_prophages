sessionInfo()



S <- read.table('../Induction_assembly/integr_sites_coo.txt', sep = '\t', header = T, stringsAsFactors = F)

Q <- read.table('../Induction_raw_reads/Novogene_summary_table.txt', sep = '\t', header = T, stringsAsFactors = F)

Q <- Q[Q$Sequencing == 'Yes', ]



DF <- data.frame(NULL)

for (i in 3:6) {

    idx <- which(Q$Strain == S$strain[i])


    df <- S[rep(i, length(idx)), ]

    df <- cbind(df, Q[idx, c(1, 3:4)])


    df$mates_up_1000 <- NA
    df$mates_down_1000 <- NA
    df$mates_3p_1000 <- NA
    df$mates_other <- NA

    for (j in 1:nrow(df)) {

        f <- paste0(df$phage_name[j], '_5p_reads_', df$Sample.Name[j], '.sam')

        max_len <- max(count.fields(f, sep = '\t', comment.char = ''))

        t <- read.table(f, sep = '\t', comment.char = '', fill = T, header = F, col.names = 1:max_len, stringsAsFactors = F)


        sele <- which((t[, 7] == '=') & (t[, 8] >= (df$coo1[j] - 1000)) & (t[, 8] < df$coo1[j]))
    
        df$mates_up_1000[j] <- length(sele)


        sele <- which((t[, 7] == '=') & (t[, 8] >= df$coo1[j]) & (t[, 8] < df$coo1[j] + 1000))
    
        df$mates_down_1000[j] <- length(sele)


        if (i == 3) {

            sele <- which((t[, 7] == 'NODE_24_length_43119_cov_24.736044') & (t[, 8] > S$coo2[2] - 1000) & (t[, 8] <= S$coo2[2]))

        } else {

            sele <- which((t[, 7] == '=') & (t[, 8] > df$coo2[j] - 1000) & (t[, 8] <= df$coo2[j]))

        }

        df$mates_3p_1000[j] <- length(sele)


        df$mates_other[j] <- nrow(t) - df$mates_up_1000[j] - df$mates_down_1000[j] - df$mates_3p_1000[j]


        df[j, 9:12] <- round(df[j, 9:12] / nrow(t) * 100)

    }


    DF <- rbind(DF, df)

}



DF <- DF[, c(1, 7:12)]

write.table(DF, sep = '\t', row.names = F, quote = F, file = 'mate_pairs_pos.txt')
