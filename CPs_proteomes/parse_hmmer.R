.libPaths(paste0(Sys.getenv('HOME'), '/SOFTWARE/R_LIB'))
library(rhmmer)
library(IRanges)
sessionInfo()




args <- commandArgs(trailingOnly = T)
genome_id <- args[1]
pfam_file <- args[2]




# ------------------------------ input files ------------------------------ #
file1 <- paste0(genome_id, '_vs_Pfam.txt')

t <- read_domtblout(file1)

t <- as.data.frame(t, stringsAsFactors = F)

p <- read.table(pfam_file, sep = '\t', header = T, quote = '', stringsAsFactors = F)




# ------------------------------ collect data per protein-profile pair ------------------------------ #
DF <- data.frame(NULL)

for (x in unique(t$domain_name)) {

    for (y in unique(t$query_accession)) {

        idx <- which((t$domain_name == x) & (t$query_accession == y))
        if (length(idx) == 0) { next }

        ir <- IRanges(start = t$env_from[idx], end = t$env_to[idx])
        cl <- reduce(ir)
        cl <- as.data.frame(cl)
        if (sum(cl$width) < 100) { next }

        df <- data.frame(
            protein_id = x,
            profile_id = y,
            profile_name = t$query_name[idx][1],
            profile_desc = p$DESC[p$ACC == y],
            protein_env_cov = sum(cl$width),
            protein_env_coo = paste(paste0(cl$start, '-', cl$end), collapse = ';'),
            stringsAsFactors = F
        )

        DF <- rbind(DF, df)

    }

}




# ------------------------------ for each protein select profile providing maximal coverage ------------------------------ #
sele <- c()

for (x in unique(DF$protein_id)) {

    idx <- which(DF$protein_id == x)
    sele <- c(sele, idx[ DF$protein_env_cov[idx] == max(DF$protein_env_cov[idx]) ][1])

}

DF <- DF[sele, ]




# ------------------------------ output file ------------------------------ #
file2 <- paste0(genome_id, '_annot.txt')

write.table(DF, sep = '\t', row.names = F, quote = F, file = file2)
