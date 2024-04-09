library(TeachingDemos)
source('function_plot_genome_map.R')
source('function_legend_genome_map.R')
sessionInfo()



COL <- setNames(
    c('deepskyblue', 'darkorange2'),
    c('virion', 'RT')
)

D <- read.table('../Phages_proteomes/detected_domains_info.txt', sep = '\t', header = T, quote = '', stringsAsFactors = F)



pdf('phages_genome_maps.pdf', height = 6.5, width = 7.5)

par(mar = c(0.1, 6.5, 0.1, 1))

Xmax <- 55000
Ymax <- 16 * 7

plot(
    NA,
    xlim = c(0, Xmax),
    ylim = c(0, Ymax),
    xaxs = 'i',
    yaxs = 'i',
    axes = F,
    ann = F
)



# ------------------------------ genome maps ------------------------------ #
YB <- Ymax - 10

for (p in c('Tulp', 'Roos', 'Pioen', 'Aster', 'Lelie', 'CP2', 'CP11')) {

    lab <- p

    len <- read.table(paste0('../Phages_map_reads/BED_FILES/', p, '_genome.bed'), sep = '\t', header = F, stringsAsFactors = F)[, 3]

    if (grepl('^CP', p)) { dir <- '../CPs_proteomes/' } else { dir <- '../Phages_proteomes/' }
    tab <- read.table(paste0(dir, p, '_proteome.txt'), sep = '\t', header = F, stringsAsFactors = F)
    dom <- read.table(paste0(dir, p, '_annot.txt'), sep = '\t', header = T, quote = '', stringsAsFactors = F)

    dom$color <- NA
    for (x in names(COL)) {

        sele <- D$profile_id[ D[, x] == 1 ]

        dom$color[ dom$profile_id %in% sele ] <- COL[x]

    }
    dom <- dom[!(is.na(dom$color)), ]

    if (p == 'Aster') { dgr <- NULL } else { dgr <- read.table(paste0('../DGR_detection/', p, '_DGR_repeats.txt'), sep = '\t', header = T, stringsAsFactors = F) }

    save(lab, len, tab, dom, dgr, COL, file = paste0(p, '_genome_map.RData'))


    if (lab == 'Tulp') { lab <- 'Mushu' }
    if (lab == 'CP2') { lab <- 'Lagaffe_CP' }

    plot_genome_map(
        lab,
        len,
        tab,
        dom,
        dgr,
        dom_lab = F,
        dgr_col = COL['RT'],
        cex_lab = 0.35,
        yBottom = YB
    )


    YB <- YB - 16

}



# ------------------------------ legend ------------------------------ #
XL <- 44314
YB <- 48


axis(
    side = 1,
    pos = 108,
    at = c(XL, XL + 10000),
    labels = c('0', '10 kb'),
    cex.axis = 0.8,
    mgp = c(3, 0.4, 0),
    tck = -0.01
)


legend_genome_map(XL, YB, COL, 0.8)


dev.off()
