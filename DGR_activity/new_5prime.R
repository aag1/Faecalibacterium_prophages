library(seqinr)
sessionInfo()




args <- commandArgs(trailingOnly = T)
ref_fasta <- args[1]
seq_fasta <- args[2]
seq_dtr   <- args[3]
blast_res <- args[4]
out_fasta <- args[5]




# ------------------------------ read data ------------------------------ #
S <- read.fasta(
    ref_fasta,
    seqtype = 'DNA',
    forceDNAtolower = F
)


L <- read.fasta(
    seq_fasta,
    seqtype = 'DNA',
    forceDNAtolower = F
)


df <- read.table(
        seq_dtr,
        sep = '\t',
        header = F,
        stringsAsFactors = F
)

endsOverlap <- setNames(df[,2], df[,1])


tab <- read.table(
    blast_res,
    sep = '\t',
    header = F,
    stringsAsFactors = F
)

colnames(tab) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'qlen', 'qcovs', 'sstart', 'send', 'slen', 'sstrand', 'evalue', 'bitscore', 'nident')

tab <- tab[tab$length > 400, ]

if (!all(names(L) %in% tab$sseqid)) { stop('No new 5\'-termini for some genomes!') }

if (nrow(tab) > length(L)) { stop('Multiple new 5\'-termini for some genomes!') }




# ------------------------------ process data ------------------------------ #
for (i in 1:nrow(tab)) {

    # genome sequence
    ACC <- tab$sseqid[i]
    SEQ <- L[[ ACC ]]


    # new start coordinate
    START <- tab$sstart[i]


    # reverse complement sequence
    if (tab$sstrand[i] == 'minus') {

        SEQ <- rev(comp(SEQ, forceToLower = F, ambiguous = T))

        START <- length(SEQ) - START + 1
    }


    # remove 3'-end nucleotides matching the 5'-end
    if (ACC %in% names(endsOverlap)) {

        SEQ <- SEQ[ 1:(length(SEQ) - endsOverlap[ACC]) ]

        if (length(SEQ) < START) { stop('After removing the 3\'-end nucleotides matching the 5\'-end, the new start coordinate is gone!') }

    }


    # shift the start coordinate
    SEQ <- SEQ[ c(START:length(SEQ), 1:(START-1)) ]


    # update genome
    L[[ ACC ]] <- SEQ

}


L <- c(S, L)


write.fasta(
    sequences = L,
    names = names(L),
    nbchar = 80,
    file.out = out_fasta
)
