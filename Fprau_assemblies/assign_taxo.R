sessionInfo()



tab <- read.table('dRep_out/data_tables/Cdb.csv', sep = ',', header = T, stringsAsFactors = F)
tab$genome <- sub('_contigs\\.fasta$', '', tab$genome)
tab$genome <- sub('\\.fna$', '', tab$genome)

t1 <- read.table('../UHGG_database/genomes-all_metadata.tsv', sep = '\t', header = T, stringsAsFactors = F)



DF <- data.frame()

for (i in 1:nrow(tab)) {

    if (grepl('^MGY', tab$genome[i])) { next }


    rep <- tab$genome[ (tab$secondary_cluster == tab$secondary_cluster[i]) & grepl('^MGY', tab$genome) ]


    if (length(rep) == 0) {

        taxonomy <- tab$secondary_cluster[i]

    } else {

        if (length(rep) > 1) { stop('Multiple UHGG genomes in ', tab$secondary_cluster[i], ' cluster!') }

        taxonomy <- t1$Lineage[ t1$Genome == rep ]

    }


    df <- data.frame(strain_id = tab$genome[i], taxonomy, stringsAsFactors = F)
    DF <- rbind(DF, df)

}

DF <- DF[order(DF$taxonomy), ]



write.table(DF, sep = '\t', row.names = F, quote = F, file = 'Fprau_strains_taxo.txt')
