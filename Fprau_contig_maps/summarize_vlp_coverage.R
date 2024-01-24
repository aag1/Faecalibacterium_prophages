sessionInfo()



STRAIN_ID <- commandArgs(trailingOnly = T)[1]


f1 <- paste0('../Fprau_map_reads/BED_FILES/', STRAIN_ID, '_contigs.bed')
t1 <- read.table(f1, sep = '\t', row.names = 1, header = F, stringsAsFactors = F)


t2 <- read.table('contaminated_samples.txt', sep = '\t', header = T, stringsAsFactors = F)


W <- list()
D <- list()
B <- list()



# for each VLP study...
for (STUDY_NAME in c('VanEspen_2021', 'Shkoporov_2019')) {


    # for each VLP sample...
    for (SAMPLE_ID in t2$sample_id[(t2$study_name == STUDY_NAME) & (!t2$contaminated)]) {


        # coverage depth file
        f3 <- paste0('../Fprau_map_reads/mapped_', STUDY_NAME, '_reads/', STRAIN_ID, '_', STUDY_NAME, '_', SAMPLE_ID, '_reads.depth.txt')

        t3 <- tryCatch(
            read.table(f3, sep = '\t', header = F, stringsAsFactors = F),
            error = function (e) { return(data.frame(NULL)) }
        )

        if (nrow(t3) == 0) {
            cat('Depth file corresponding to the', STUDY_NAME, 'VLP sample', SAMPLE_ID, 'is empty! Skipping this sample...\n')
            next
        }


        # for each F. prau contig > 50 kb...
        for (CONTIG_ID in rownames(t1)[t1[, 2] > 50000]) {


            if (!(CONTIG_ID %in% t3[, 1])) { next }


            len <- t1[CONTIG_ID, 2]


            i <- which(t3[, 1] == CONTIG_ID)
            coo <- t3[i, 2]
            cov <- t3[i, 3]


            COVERAGE <- rep(0, len)
            COVERAGE[coo] <- cov


            w_center <- seq(from = 1501, to = len - 1500, by = 500)
            w_from <- w_center - 1500
            w_to <- w_center + 1500
            w_to[ length(w_to) ] <- len
            w_size <- w_to - w_from + 1


            # sliding window centers
            if (!(CONTIG_ID %in% names(W))) { W[[ CONTIG_ID ]] <- w_center }


            # depth of coverage per window
            if (!(CONTIG_ID %in% names(D))) { D[[ CONTIG_ID ]] <- list() }
            if (!(STUDY_NAME %in% names(D[[ CONTIG_ID ]]))) { D[[ CONTIG_ID ]][[ STUDY_NAME ]] <- list() }

            D[[ CONTIG_ID ]][[ STUDY_NAME ]][[ SAMPLE_ID ]] <- sapply(seq_along(w_center), function (i) {

                mean( COVERAGE[ w_from[i]:w_to[i] ] )

            })


            # breadth of coverage per window
            if (!(CONTIG_ID %in% names(B))) { B[[ CONTIG_ID ]] <- list() }
            if (!(STUDY_NAME %in% names(B[[ CONTIG_ID ]]))) { B[[ CONTIG_ID ]][[ STUDY_NAME ]] <- list() }

            B[[ CONTIG_ID ]][[ STUDY_NAME ]][[ SAMPLE_ID ]] <- sapply(seq_along(w_center), function (i) {

                sum( COVERAGE[ w_from[i]:w_to[i] ] > 0 ) / w_size[i] * 100

            })

        }

    }

}



save(W, D, B, file = paste0(STRAIN_ID, '_vlp_coverage_summary.RData'))
