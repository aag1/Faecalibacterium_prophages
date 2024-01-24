sessionInfo()




t1 <- read.table('../Fprau_assemblies/Fprau_strains_taxo.txt', sep = '\t', header = T, stringsAsFactors = F)

t2 <- read.table('../Fprau_contig_maps/prophages_raw_coo.txt', sep = '\t', header = T, stringsAsFactors = F)

t3 <- read.table('new_coo.txt', sep = '\t', header = T, stringsAsFactors = F)




DF <- t2

for (i in 1:nrow(DF)) {

    j <- which(t3$prophage == DF$candidate_prophage[i])

    DF$new_from[i] <- DF$from[i] - 1 + t3$from[j]

    DF$new_to[i]   <- DF$from[i] - 1 + t3$to[j]

}

write.table(DF, sep = '\t', row.names = F, quote = F, file = 'prophages_refined_coo.txt')




DF$contig_id_short <- sub('^(NODE_[0-9]+)_.+$', '\\1', DF$contig_id)

DF$new_coo <- paste(DF$new_from, '-', DF$new_to)

DF$strain_taxo <- sapply(DF$strain_id, function (x) {

    sp <- rev(strsplit(t1$taxonomy[ t1$strain_id == x ], ';')[[1]])[1]
    sp <- sub('^s__Faecalibacterium', 'F.', sp)
    sp <- gsub('_', ' ', sp)
    return(sp)

})

DF <- DF[, c('candidate_prophage', 'strain_id', 'strain_taxo', 'contig_id_short', 'new_coo')]

write.table(DF, sep = '\t', row.names = F, quote = F, file = 'prophages_summary.txt')
