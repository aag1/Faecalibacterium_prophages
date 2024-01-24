library(TeachingDemos)
source('function_plot_cp_map.R')
sessionInfo()



DF <- read.table('../CPs_refine_coo/prophages_refined_coo.txt', sep = '\t', header = T, stringsAsFactors = F)
DF$len <- DF$new_to - DF$new_from + 1

STR <- read.table('detected_domains_info.txt', sep = '\t', header = T, quote = '', stringsAsFactors = F)
STR <- STR$profile_id[ STR$virion == 1 ]



pdf('CP_genome_maps.pdf', height = 8.75, width = 7.5)


par(mar = c(1, 4.5, 0.5, 0.5))

Xmax <- ceiling(max(DF$len) / 1000) * 1000
Ymax <- 15 * nrow(DF)

plot(
    NA,
    xlim = c(0, Xmax),
    ylim = c(0, Ymax),
    xaxs = 'i',
    yaxs = 'i',
    axes = F,
    ann = F
)


YB <- Ymax - 15

for (i in 1:nrow(DF)) {

    tab <- read.table(paste0(DF$candidate_prophage[i], '_proteome.txt'), sep = '\t', header = F, stringsAsFactors = F)

    dom <- read.table(paste0(DF$candidate_prophage[i], '_annot.txt'), sep = '\t', header = T, quote = '', stringsAsFactors = F)
    dom <- dom[dom$profile_id %in% STR, ]
    dom$color <- 'deepskyblue'

    plot_cp_map(
        lab = DF$candidate_prophage[i],
        len = DF$len[i],
        tab,
        dom,
        yBottom = YB
    )

    YB <- YB - 15

}


axis(
    side = 1,
    pos = Ymax - 9,
    at = max(DF$len) - c(10000, 0),
    labels = c('0', '10 kb'),
    cex.axis = 0.5,
    mgp = c(3, 0.25, 0),
    tck = -0.01
)


dev.off()
