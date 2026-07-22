#library(enrichR)
library(ggplot2)
require(EnhancedVolcano)

res<-read.csv('../res_data/revision/volcano/doubletsDEA_dc9.csv', row.names=1)

pdf('../figures/revision/c9_volcano_n.pdf')
update_geom_defaults("text", aes(family = "sans"))
p1=EnhancedVolcano(res,
    lab = res$gene_name,
    selectLab=c('Spp1'),
    x = 'logFC',
    y = 'p.value',
    title = '9_Fibro-OsteoCAR',
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 4.0) 
options(repr.plot.width=9,repr.plot.height=7)
p1
#p1*theme(text=element_text(family = 'sans'))
dev.off()

pdf('../figures/revision/c9_volcano_allgenes_n.pdf')
update_geom_defaults("text", aes(family = "sans"))
p1=EnhancedVolcano(res,
    lab = res$gene_name,
    #selectLab=c('Spp1'),
    x = 'logFC',
    y = 'p.value',
    title = '9_Fibro-OsteoCAR',
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 4.0) 
options(repr.plot.width=9,repr.plot.height=7)
p1
#p1*theme(text=element_text(family = 'sans'))
dev.off()

## OsteoCAR

res<-read.csv('../res_data/revision/volcano/singletsDEA_OsteoCAR.csv', row.names=1)

pdf('../figures/revision/OsteoCAR_volcano_n.pdf')
update_geom_defaults("text", aes(family = "sans"))
p1=EnhancedVolcano(res,
    lab = res$gene_name,
    selectLab=c('Spp1'),
    x = 'logFC',
    y = 'p.value',
    title = 'OsteoCAR',
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 4.0) 
options(repr.plot.width=9,repr.plot.height=7)
p1
#p1*theme(text=element_text(family = 'sans'))
dev.off()

pdf('../figures/revision/OsteoCAR_volcano_allgenes_n.pdf')
update_geom_defaults("text", aes(family = "sans"))
p1=EnhancedVolcano(res,
    lab = res$gene_name,
    #selectLab=c('Spp1'),
    x = 'logFC',
    y = 'p.value',
    title = 'OsteoCAR',
    pCutoff = 0.05,
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 4.0) 
options(repr.plot.width=9,repr.plot.height=7)
p1
#p1*theme(text=element_text(family = 'sans'))
dev.off()

