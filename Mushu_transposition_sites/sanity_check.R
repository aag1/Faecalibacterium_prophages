.libPaths(paste0(Sys.getenv('HOME'), '/SOFTWARE/R_LIB'))
library(seqinr)
library(SamSeq)
sessionInfo()



t <- read.table('PhIFM7_102A3_split_coo.txt', sep = '\t', header = T, stringsAsFactors = F)

S <- read.fasta('../Fprau_assemblies/FM7_S1_contigs.fasta', seqtype = 'DNA', forceDNAtolower = F)



t$strand1 <- sapply(t$read_flag, function (x) ifelse(samFlags(x)['READ_REVERSE_STRAND'], '-', '+'))



t$test1 <- -1
t$test2 <- -1

for (i in 1:nrow(t)) {

    read <- strsplit(t$read_seq[i], '')[[1]]



    if (t$cigar1_type[i] %in% c('MS', 'MH')) {

        l <- as.numeric(sub('M[0-9]+[SH]$', '', t$cigar1[i]))

        seq1 <- read[1:l]

    } else {

        l <- as.numeric(sub('^[0-9]+[SH]([0-9]+)M$', '\\1', t$cigar1[i]))

        if (t$cigar1_type[i] == 'SM') { seq1 <- rev(rev(read)[1:l]) }
        if (t$cigar1_type[i] == 'HM') { seq1 <- read }

    }

    seq2 <- S[[ t$contig[i] ]][ t$pos1[i] + 0:(l-1) ]

    if (length(seq1) != length(seq2)) { stop('Unequal lengths seq1 & seq2!') }

    t$test1[i] <- sum(seq1 != seq2)



    if (t$cigar1_type[i] %in% c('HM', 'MH')) { next }   # hard-clipped part of the read is absent from the bam file

    if (t$strand1[i] == '-') { read <- rev(comp(read, forceToLower = F, ambiguous = T)) }
    if (t$strand2[i] == '-') { read <- rev(comp(read, forceToLower = F, ambiguous = T)) }

    if (t$cigar2_type[i] == 'MS') {

        l <- as.numeric(sub('M[0-9]+S$', '', t$cigar2[i]))

        SEQ1 <- read[1:l]

    } else {

        l <- as.numeric(sub('^[0-9]+S([0-9]+)M$', '\\1', t$cigar2[i]))

        SEQ1 <- rev(rev(read)[1:l])

    }

    SEQ2 <- S[[ 'NODE_11_length_112199_cov_46.673290' ]][ t$pos2[i] + 0:(l-1) ]

    if (length(SEQ1) != length(SEQ2)) { stop('Unequal lengths SEQ1 & SEQ2!') }

    t$test2[i] <- sum(SEQ1 != SEQ2)

}



table(t$test1)
table(t$test2)
