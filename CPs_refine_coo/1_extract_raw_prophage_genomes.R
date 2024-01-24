library(seqinr)
sessionInfo()



DF <- read.table('../Fprau_contig_maps/prophages_raw_coo.txt', sep = '\t', header = T, stringsAsFactors = F)



for (i in 1:nrow(DF)) {

    S <- read.fasta(
        paste0('../Fprau_assemblies/', DF$strain_id[i], '_contigs.fasta'),
        seqtype = 'DNA',
        forceDNAtolower = F
    )


    p <- S[[ DF$contig_id[i] ]][ DF$from[i]:DF$to[i] ]


    write.fasta(
        names = DF$candidate_prophage[i],
        sequences = p,
        nbchar = 80,
        file.out = paste0(DF$candidate_prophage[i], '_raw.fna')
    )

}
