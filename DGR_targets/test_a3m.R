library(seqinr)
sessionInfo()



for (f in list.files(pattern = 'fasta')) {

    b <- sub('\\.fasta', '', f)


    S1 <- read.fasta(f, seqtype = 'AA', set.attributes = F)[[1]]
    S1 <- toupper(S1)

    S2 <- read.fasta(paste0('HHpred_output/hhpred_full_', b, '.a3m'), seqtype = 'AA', set.attributes = F)[[3]]


    if (identical(S1, S2)) {
        cat(b, 'FASTA and A3M sequences are identical!\n')
    } else {
        cat(b, 'FASTA and A3M sequences are NOT identical!\n')
    }

}
