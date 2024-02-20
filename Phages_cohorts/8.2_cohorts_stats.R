sessionInfo()




# ------------------------------ input data ------------------------------ #
DF <- read.table('LLD_IBD_reads_phages_host.txt', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)

DF$Host <- (DF$HostMetaphlan > 0)

sc <- min(DF$HostMetaphlan[DF$HostMetaphlan > 0]) / 2
DF$HostMetaphlanLogNorm <- log10(DF$HostMetaphlan + sc)



path <- 'XXX/'

key1 <- read.table(paste0(path, 'LLD_LLD2_Info/LLD_GTMicrob.txt'), sep = '\t', row.names = 1, header = F, stringsAsFactors = F)

path <- 'YYY/'

key2 <- read.table(paste0(path, 'IBD_Info/rename_IBD.txt'), sep = '\t', header = T, stringsAsFactors = F)

pheno <- read.table(paste0(path, 'IBD_Info/data_from_Renate/LLD_IBD_meta_201020.txt'), quote = '', sep = '\t', row.names = 1, header = T, stringsAsFactors = F)

idx <- which(rownames(pheno) %in% rownames(key1))
rownames(pheno)[idx] <- key1[rownames(pheno)[idx], 1]

idx <- which(rownames(pheno) %in% key2$Classic[key2$Classic != 'QQQ'])
rownames(pheno)[idx] <- sapply(rownames(pheno)[idx], function (x) key2[key2$Classic == x, 'old'])
rownames(pheno)[rownames(pheno) == 'QQQ'] <- 'NNN'



all(rownames(DF) %in% rownames(pheno))

DF$AgeAtFecalSampling <- pheno[rownames(DF), 'AgeAtFecalSampling']

DF$Sex <- pheno[rownames(DF), 'Sex']




# ------------------------------ test w/o adjustment for host ------------------------------ #
TAB <- data.frame(NULL)

for (x in c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie', 'Host')) {

    glm.sum <- summary(glm(
        DF[, x] ~ DF$Cohort + DF$AgeAtFecalSampling + DF$Sex + DF$NumCleanReads,
        family = binomial(link = 'logit')
    ))

    tab <- data.frame(
        LLD_pos = sum(DF[DF$Cohort == 'LLD', x]) / sum(DF$Cohort == 'LLD'),
        IBD_pos = sum(DF[DF$Cohort == 'IBD', x]) / sum(DF$Cohort == 'IBD'),
        Beta = glm.sum$coef[2, 1],
        SE = glm.sum$coef[2, 2],
        Z = glm.sum$coef[2, 3],
        P = glm.sum$coef[2, 4],
        stringsAsFactors = F
    )
    rownames(tab) <- x

    TAB <- rbind(TAB, tab)

}

TAB$P_adj <- p.adjust(TAB$P, method = 'BH')




# ------------------------------ test with adjustment for host ------------------------------ #
TAB$Beta.2 <- NA
TAB$SE.2 <- NA
TAB$Z.2 <- NA
TAB$P.2 <- NA
TAB$P_adj.2 <- NA

for (x in c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie')) {

    glm.sum <- summary(glm(
        DF[, x] ~ DF$Cohort + DF$AgeAtFecalSampling + DF$Sex + DF$NumCleanReads + DF$HostMetaphlanLogNorm,
        family = binomial(link = 'logit')
    ))

    TAB[x, 'Beta.2'] <- glm.sum$coef[2, 1]
    TAB[x, 'SE.2'] <- glm.sum$coef[2, 2]
    TAB[x, 'Z.2'] <- glm.sum$coef[2, 3]
    TAB[x, 'P.2'] <- glm.sum$coef[2, 4]

}

idx <- which(!is.na(TAB$P.2))
TAB[idx, 'P_adj.2'] <- p.adjust(TAB[idx, 'P.2'], method = 'BH')

write.table(TAB, sep = '\t', quote = F, file = 'LLD_IBD_phages_logistic_regression.txt')
