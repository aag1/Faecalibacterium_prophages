library(seqinr)
sessionInfo()




args <- commandArgs(trailingOnly = T)

p <- args[1]


f <- paste0(p, '_DGR_repeats.txt')

t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F, colClasses = setNames('character', 'TR_nt'))


f <- ifelse(
    grepl('^CP', p),
    paste0('../CPs_refine_coo/', p, '_refined.fna'),
    paste0('../Induction_assembly/', p, '.fasta')
)

S <- read.fasta(f, seqtype = 'DNA', forceDNAtolower = F)




if (nrow(t) == 1) {
        
    TR_seq <- S[[1]][ t$TR_coo1[1]:t$TR_coo2[1] ]
    VR_seq <- S[[1]][ t$VR_coo1[1]:t$VR_coo2[1] ]


    write.fasta(
        names = paste0(p, c('_TR', '_VR')),
        sequences = list(TR_seq, VR_seq),
        file.out = paste0(p, '_DGR_repeats_seq.fasta'),
        nbchar = 80
    )


    df <- data.frame(
        genome = p,
        TR_coo = paste0(t$TR_coo1, '-', t$TR_coo2),
        VR_coo = paste0(t$VR_coo1, '-', t$VR_coo2),
        stringsAsFactors = F
    )

}




if (nrow(t) == 2) {

    if (p == 'Lelie') {

        t[1, c('TR_coo2', 'VR_coo2')] <- t[1, c('TR_coo2', 'VR_coo2')] + 3
        t[2, c('TR_coo1', 'VR_coo1')] <- t[2, c('TR_coo1', 'VR_coo1')] - 21

    }

    if (p == 'CP2') {

        t[1, c('TR_coo1', 'VR_coo1')] <- t[1, c('TR_coo1', 'VR_coo1')] - 10

    }

    if (p == 'CP11') {

        t[1, c('TR_coo2', 'VR_coo2')] <- t[1, c('TR_coo2', 'VR_coo2')] + 1

    }

    if ((t$TR_coo1[1] != t$TR_coo1[2]) | (t$TR_coo2[1] != t$TR_coo2[2])) {

        stop(p, ' TRs do not match!')

    }


    TR_seq <- S[[1]][ t$TR_coo1[1]:t$TR_coo2[1] ]
    VR1_seq <- S[[1]][ t$VR_coo1[1]:t$VR_coo2[1] ]
    VR2_seq <- S[[1]][ t$VR_coo1[2]:t$VR_coo2[2] ]


     write.fasta(
        names = paste0(p, c('_TR', '_VR1', '_VR2')),
        sequences = list(TR_seq, VR1_seq, VR2_seq),
        file.out = paste0(p, '_DGR_repeats_seq.fasta'),
        nbchar = 80
     )


    df <- data.frame(
        genome = p,
        TR_coo = sapply(1:2, function (i) paste0(t$TR_coo1[i], '-', t$TR_coo2[i])),
        VR_coo = sapply(1:2, function (i) paste0(t$VR_coo1[i], '-', t$VR_coo2[i])),
        stringsAsFactors = F
    )

}

write.table(df, sep = '\t', quote = F, row.names = F, file = paste0(p, '_DGR_repeats_summary.txt'))
