library(seqinr)
sessionInfo()



DF <- read.table('prophages_refined_coo.txt', sep = '\t', header = T, stringsAsFactors = F)



for (i in 1:nrow(DF)) {

    S <- read.fasta(
        paste0('../Fprau_assemblies/', DF$strain_id[i], '_contigs.fasta'),
        seqtype = 'DNA',
        forceDNAtolower = F
    )


    p <- S[[ DF$contig_id[i] ]][ DF$new_from[i]:DF$new_to[i] ]


    if (DF$strand[i] == 'r') {

        p <- rev(comp(p, forceToLower = F, ambiguous = T))

    }


    write.fasta(
        names = DF$candidate_prophage[i],
        sequences = p,
        nbchar = 80,
        file.out = paste0(DF$candidate_prophage[i], '_refined.fna')
    )

}
