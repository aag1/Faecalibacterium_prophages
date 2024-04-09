sessionInfo()



t <- read.table(
    'PhIFM7_102A3_split_info.txt',
    sep = '\t',
    header = F,
    comment.char = '',
    fill = T,
    col.names = paste0('V', 1:18),  # https://stackoverflow.com/a/18922914
    stringsAsFactors = F
)



DF <- data.frame(NULL)

for (i in 1:nrow(t)) {

    v <- strsplit(t[i, 17], ',')[[1]]

    df <- data.frame(
        contig1   = t[i, 3],
        pos1      = t[i, 4],
        cigar1    = t[i, 6],
        pos2      = as.numeric(v[2]),
        strand2   = v[3],
        cigar2    = v[4],
        read_flag = t[i, 2],
        read_seq  = t[i, 10],
        stringsAsFactors = F
    )

    DF <- rbind(DF, df)

}



DF$cigar1_type <- gsub('[0-9]', '', DF$cigar1)
table(DF$cigar1_type)

DF$cigar2_type <- gsub('[0-9]', '', DF$cigar2)
table(DF$cigar2_type)



CIGAR <- c('MH', 'HM', 'MS', 'SM')

sele <- which((DF$cigar1_type %in% CIGAR) & (DF$cigar2_type %in% CIGAR))

DF <- DF[sele, ]



DF$pos1_adj <- sapply(1:nrow(DF), function (i) {

    x <- DF$pos1[i]

    if (DF$cigar1_type[i] %in% c('MH', 'MS')) {

        x <- x + as.numeric(sub('M[0-9]+[HS]$', '', DF$cigar1[i])) - 1

    } else {

        x <- x - 1

    }   # nucleotide preceeding the transposition site

    return(x)

})



DF$pos2_adj <- sapply(1:nrow(DF), function (i) {

    x <- DF$pos2[i]

    if (DF$cigar2_type[i] == 'MS') {

        x <- x + as.numeric(sub('M[0-9]+S$', '', DF$cigar2[i])) - 1

    }   # nucleotide fused to a bacterial genome fragment

    return(x)

})



write.table(DF, sep = '\t', row.names = F, quote = F, file = 'PhIFM7_102A3_split_coo.txt')
