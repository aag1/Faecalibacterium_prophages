library(igraph)
sessionInfo()




# ------------------------------ read data ------------------------------ #
DF <- read.table('genome_by_genome_overview.csv', sep = ',', header = T, row.names = 1, stringsAsFactors = F)

tab <- read.table('../Viral_RefSeq_220/viral_refseq_taxo.txt', sep = '\t', header = T, quote = '', fill = T, stringsAsFactors = F)



PHAGES <- c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie')

CPS <- paste0('CP', c(2, 4, 5, 6, 7, 11))

REFSEQ <- rownames(DF)[ !(rownames(DF) %in% c(PHAGES, CPS)) ]



E <- read.table('c1.ntw', header = F, stringsAsFactors = F)

colnames(E) <- c('from', 'to', 'weight')

G <- graph_from_data_frame(E, directed = F, vertices = c(PHAGES, CPS, REFSEQ))




# ------------------------------ membership of the phages & CPs in the connected components of the graph ------------------------------ #
cat('\nMembership of the phages & CPs in the connected components of the graph:\n\n')

K <- clusters(G)

M <- K$membership[ c(PHAGES, CPS) ]
M; cat('\n')

for (i in unique(M)) {

    cat('\nTaxonomy of the RefSeq viruses belonging to the connected component', i, ':\n')

    V <- names(K$membership)[ K$membership == i ]

    V <- V[ V %in% REFSEQ ]

    df <- DF[V, ]

    print(table(df$Class)); cat('\n')

    idx <- which(df$Class != 'Caudoviricetes')

    if (length(idx) > 0) { print(df[idx, 'Class', drop = F]); cat('\n') }

}




# ------------------------------ membership of the phages & CPs in the VCs ------------------------------ #
cat('\nMembership of the phages & CPs in the VCs:\n\n')

IDX <- c()

for (x in c(PHAGES, CPS)) {

    if (DF[x, 'VC.Status'] != 'Clustered') { next }

    idx <- which(DF$VC == DF[x, 'VC'])

    IDX <- c(IDX, idx)

}

df <- DF[unique(IDX), 'VC', drop = F]
df; cat('\n')




# ------------------------------ pairs for dotplots ------------------------------ #
d <- data.frame(NULL)

for (x in c(PHAGES, CPS)) {

    if (DF[x, 'VC.Status'] != 'Clustered') { next }

    idx <- which((E$from == x) & (E$to %in% REFSEQ) & (DF[E$to, 'VC'] == DF[x, 'VC']))

    d <- rbind(d, E[idx, 1:2])

}

d$to <- sapply(d$to, function (x) { tab$genome_id[ grepl(paste0('^', gsub('~', ' ', x)), tab$genome_desc) ] })

write.table(d, row.names = F, col.names = F, quote = F, sep = '\t', file = 'pairs_for_dotplots.txt')
