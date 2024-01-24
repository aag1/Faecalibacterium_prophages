sessionInfo()



x <- commandArgs(trailingOnly = T)[1]



Q <- c()

files <- list.files(pattern = paste0(x, '.tsv'))

for (f in files) {

    t <- read.table(f, sep = '\t', header = T, stringsAsFactors = F)

    for (i in 1:nrow(t)) {

        v1 <- strsplit(t$fastq_md5[i], ';')[[1]]

        v2 <- strsplit(t$fastq_ftp[i], ';')[[1]]
        v2 <- sub('^.+\\/([^\\/]+)$', '\\1', v2)

        q <- setNames(v1, v2)
        Q <- c(Q, q)

    }

}



f <- paste0('md5sum_reads_', x, '.txt')

t <- read.table(f, sep = '\t', header = F, stringsAsFactors = F)

v1 <- t[, 1]

v2 <- t[, 2]
v2 <- sub('^.+\\/([^\\/]+)$', '\\1', v2)

W <- setNames(v1, v2)



cat('\n', x, '\n')
if (!all(names(W) %in% names(Q))) { stop('Unexpected files in the folder!') }
cat(length(W), 'files downloaded.\n')
cat('All md5sum values as expected:', identical(W, Q[names(W)]), '\n\n')
