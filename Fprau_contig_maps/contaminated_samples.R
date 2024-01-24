sessionInfo()




# ------------------------------ list gut virome samples ------------------------------ #
t1 <- read.table('../VLP_reads/PRJNA722819_VanEspen_2021.tsv', sep = '\t', header = T, stringsAsFactors = F)
t2 <- read.table('../VLP_reads/PRJNA723467_VanEspen_2021.tsv', sep = '\t', header = T, stringsAsFactors = F)
t3 <- read.table('../VLP_reads/PRJNA545408_Shkoporov_2019.tsv', sep = '\t', header = T, stringsAsFactors = F)

idx <- grep('^Lon_Virome_', t3$sample_alias)

DF <- data.frame(
    study_name = c(rep('VanEspen_2021', nrow(t1) + nrow(t2)), rep('Shkoporov_2019', length(idx))),
    study_acc  = c(t1$study_accession, t2$study_accession, t3$study_accession[idx]),
    sample_id  = c(t1$run_accession, t2$run_accession, t3$run_accession[idx]),
    stringsAsFactors = F
)

rownames(DF) <- DF$sample_id




# ------------------------------ % bacterial genome nucleotides covered by reads ------------------------------ #
tab <- read.table('../Fprau_raw_reads/Fprau_strains.txt', sep = '\t', header = F, stringsAsFactors = F)
V <- tab[, 2]


for (STRAIN_ID in V) {

    cat('Working with strain', STRAIN_ID, '...\n')

    DF[, STRAIN_ID] <- NA

    for (STUDY_NAME in unique(DF$study_name)) {

        for (SAMPLE_ID in rownames(DF)[DF$study_name == STUDY_NAME]) {

            f <- paste0('../Fprau_map_reads/mapped_', STUDY_NAME, '_reads/', STRAIN_ID, '_', STUDY_NAME, '_', SAMPLE_ID, '_reads.coverage.txt')
            t <- read.table(f, sep = '\t', header = F, stringsAsFactors = F)

            DF[SAMPLE_ID, STRAIN_ID] <- sum(t[, 5]) / sum(t[, 6]) * 100

        }

    }

}




# ------------------------------ summary ------------------------------ #
DF$contaminated <- apply(DF[, 3 + seq_along(V)], 1, function (v) any(v > 25))


write.table(DF, sep = '\t', row.names = F, quote = F, file = 'contaminated_samples.txt')


aggregate(DF$sample_id, by = list(DF$study_name, DF$study_acc, DF$contaminated), FUN = length)
