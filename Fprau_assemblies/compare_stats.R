sessionInfo()




# ------------------------------ raw data ------------------------------ #
s <- read.table('stats_scaffolds.txt', sep = '\t', header = T, stringsAsFactors = F)
rownames(s) <- sub('^([^_]+)_.+$', '\\1', basename(s$filename))


l <- read.table('stats_Lei.txt', sep = '\t', header = T, stringsAsFactors = F)
rownames(l) <- sub('\\.fna$', '', basename(l$filename))


rownames(s)[!(rownames(s) %in% rownames(l))]    # FM4
rownames(l)[!(rownames(l) %in% rownames(s))]    # FM5

sele <- rownames(s)[rownames(s) %in% rownames(l)]
s <- s[sele, ]
l <- l[sele, ]




# ------------------------------ scaffolds: my vs Lei ------------------------------ #
pdf('scaffolds_my_vs_Lei.pdf')

layout(matrix(1:4, nrow = 2, byrow = T))


plot(x = s$scaf_N50, y = l$scaf_N50, asp = 1, pch = 20, xlab = 'my', ylab = 'Lei', main = 'scaffolds N50')
abline(coef = c(0,1))


plot(x = s$scaf_L50 / 1000, y = l$scaf_L50 / 1000, asp = 1, pch = 20, xlab = 'my', ylab = 'Lei', main = 'scaffolds L50 (kbp)')
abline(coef = c(0,1))


plot(x = s$scaf_bp / 1000, y = l$scaf_bp / 1000, asp = 1, pch = 20, xlab = 'my', ylab = 'Lei', main = 'scaffolds kbp')
abline(coef = c(0,1))


plot(x = s$scaf_n_gt50K, y = l$scaf_n_gt50K, asp = 1, pch = 20, xlab = 'my', ylab = 'Lei', main = '# scaffolds > 50 kbp')
abline(coef = c(0,1))


dev.off()
