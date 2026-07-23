library('CrossTalkeR')
library('dplyr')
library("stringr")
require(circlize)
require(igraph)
library(RColorBrewer)
suppressPackageStartupMessages({require(igraph)})
suppressPackageStartupMessages({require(ggraph)})
suppressPackageStartupMessages({require(ggplot2)})


fibCore=c('Spp1+Mac', 'MK', 'TCells', 'OsteoCAR', 'Ery')
data_sub=read.csv('../res_data/revision/data_sub_coloc_glycoprots_sankey.csv', row.names=1)

pdf('../figures/revision/F4/F4A.pdf', 8, 5)

plot_sankey(lrobj_tbl = data_sub, threshold=30, ligand_cluster='OsteoCAR', receptor_cluster=fibCore)

dev.off()

#pdf('../figures/revision/sankey_PIC_GPs_FibCore_30ints_n.pdf', 8, 5)

#plot_sankey(lrobj_tbl = data_sub, threshold=30, ligand_cluster=fibCore, receptor_cluster=fibCore)

#dev.off()

pdf('../figures/revision/F5/F5A.pdf', 8, 5)

plot_sankey(lrobj_tbl = data_sub, threshold=25, target=c("Spp1|L"), ligand_cluster=c('Spp1+Mac'))

dev.off()
