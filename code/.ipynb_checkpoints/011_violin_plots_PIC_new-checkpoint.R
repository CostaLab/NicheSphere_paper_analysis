library(ggplot2)
library(dplyr)
library(patchwork)

violin_p<-function(d, gene, ymax=5){

    d.arranged <- d %>% 
        group_by(annotation) %>% 
        mutate(gene_imp.mean = mean(gene_imp, na.rm=T)) %>% 
        arrange(desc(gene_imp.mean))

    annotation.ordered <- unique(d.arranged$annotation)

    d.reordered <- d.arranged %>% 
        mutate(annotation = factor(annotation, levels = annotation.ordered))

    ggplot(d.reordered, aes(x=annotation, y=gene_imp, fill=stage)) + geom_violin(scale='width', trim=F, na.rm=T) + 
        ylim(0, ymax)+ 
        theme(axis.text.x = element_text(angle = 75, hjust=1, size=12), axis.title=element_text(size=14,face="bold"),
            panel.background = element_rect(fill = "white"),
            panel.border = element_blank(),  # Remove panel border
            panel.grid.major = element_blank(), # Remove major grid lines
            panel.grid.minor = element_blank(), # Remove minor grid lines
            axis.line = element_line(color = "black")) +
        labs(y=gene)+scale_fill_manual(values=c("darkblue", "darkred"))

}


dat<-read.csv('../res_data/revision/violin/Thpo_sing_violinData_imp_f.csv', row.names=1)

violin_p(dat, gene='Thpo', ymax=0.35)
ggsave(file="../figures/revision/Thpo_sing_violin_imp_f.pdf", width=7, height=4, dpi=300)

###

dat<-read.csv('../res_data/revision/violin/Mpl_sing_violinData_imp_f.csv', row.names=1)

violin_p(dat, gene='Mpl', ymax=0.35)
ggsave(file="../figures/revision/Mpl_sing_violin_imp_f.pdf", width=7, height=4, dpi=300)

###

violin_p<-function(d, gene, min_ratio, ymax=5){
    #filter by min ratio
    rts<-sapply(unique(d$annotation), function(x){mean(d$gene_imp[d$stage=='TPO' & d$annotation==x])/mean(d$gene_imp[d$stage=='EV' &        d$annotation==x])})
    d<-d[d$annotation %in% names(rts)[rts>=min_ratio & !is.nan(rts)],]

    d.arranged <- d %>% 
        group_by(annotation) %>% 
        mutate(gene_imp.mean = mean(gene_imp, na.rm=T)) %>% 
        arrange(desc(gene_imp.mean))

    annotation.ordered <- unique(d.arranged$annotation)

    d.reordered <- d.arranged %>% 
        mutate(annotation = factor(annotation, levels = annotation.ordered))

    ggplot(d.reordered, aes(x=annotation, y=gene_imp, fill=stage)) + geom_violin(scale='width', trim=F, na.rm=T) + 
        ylim(0, ymax)+ 
        theme(axis.text.x = element_text(angle = 75, hjust=1, size=12), axis.title=element_text(size=14,face="bold"),
            panel.background = element_rect(fill = "white"),
            panel.border = element_blank(),  # Remove panel border
            panel.grid.major = element_blank(), # Remove major grid lines
            panel.grid.minor = element_blank(), # Remove minor grid lines
            axis.line = element_line(color = "black")) +
        labs(y=gene)+scale_fill_manual(values=c("darkblue", "darkred"))

}


###

dat<-read.csv('../res_data/revision/violin/Il1b_sing_violinData_imp_f.csv', row.names=1)

p1=violin_p(dat, gene='Il1b', min_ratio = 1, ymax=1.5)

###

dat<-read.csv('../res_data/revision/violin/Il1b_doub_violinData_imp_f.csv', row.names=1)

p2=violin_p(dat, gene='Il1b', min_ratio = 1, ymax=1.5)
p1/p2
ggsave(file="../figures/revision/Il1b_sing_doub_violin_imp_f.pdf", width=4, height=6, dpi=300)

###

###

dat<-read.csv('../res_data/revision/violin/S100a9_sing_violinData_imp_f.csv', row.names=1)

p1=violin_p(dat, gene='S100a9', min_ratio = 0, ymax=10)

###

dat<-read.csv('../res_data/revision/violin/S100a9_doub_violinData_imp_f.csv', row.names=1)

p2=violin_p(dat, gene='S100a9', min_ratio = 0, ymax=10)
p1/p2
ggsave(file="../figures/revision/S100a9_sing_doub_violin_imp_f.pdf", width=4, height=6, dpi=300)

###

dat1<-read.csv('../res_data/revision/violin/Pf4_sing_violinData_imp_f.csv', row.names=1)
dat2<-read.csv('../res_data/revision/violin/Pf4_doub_violinData_imp_f.csv', row.names=1)

p1=violin_p(dat1, gene="Pf4", min_ratio = 1, ymax=1.5)
p2=violin_p(dat2, gene="Pf4", min_ratio = 1, ymax=1.5)

p1/p2
ggsave(file="../figures/revision/Pf4_sing_doub_violin_imp_f.pdf", width=4, height=6, dpi=300)

###

violin_p<-function(d, gene, min_ratio, ymax=5){
    #filter by min ratio
    rts<-sapply(unique(d$annotation), function(x){mean(d$gene[d$stage=='TPO' & d$annotation==x])/mean(d$gene[d$stage=='EV' &        d$annotation==x])})

    d.arranged <- d %>% 
        group_by(annotation) %>% 
        mutate(gene.mean = mean(gene, na.rm=T)) %>% 
        arrange(desc(gene.mean))

    annotation.ordered <- unique(d.arranged$annotation)

    d.reordered <- d.arranged %>% 
        mutate(annotation = factor(annotation, levels = annotation.ordered))

    ggplot(d.reordered, aes(x=annotation, y=gene, fill=stage)) + geom_violin(scale='width', trim=F, na.rm=T) + 
        ylim(0, ymax)+ 
        theme(axis.text.x = element_text(angle = 75, hjust=1, size=12), axis.title=element_text(size=14,face="bold"),
            panel.background = element_rect(fill = "white"),
            panel.border = element_blank(),  # Remove panel border
            panel.grid.major = element_blank(), # Remove major grid lines
            panel.grid.minor = element_blank(), # Remove minor grid lines
            axis.line = element_line(color = "black")) +
        labs(y=gene)+scale_fill_manual(values=c("darkblue", "darkred"))

}


dat1<-read.csv('../res_data/revision/violin/Spp1_sing_violinData_imp_f.csv', row.names=1)
dat1<-dat1[dat1$annotation %in% c('OsteoCAR', 'Spp1+Mac', 'ResMac', 'InfiltMac', 'Ly6ClowMac'),]
dat2<-read.csv('../res_data/revision/violin/Spp1_doub_violinData_imp_f.csv', row.names=1)
dat2<-dat2[dat2$annotation %in% c('9_Fibro-OsteoCAR', '10_OsteoCAR-MonoPro', '6_Spp1+Mac-MonoPro', '3_Ery-MK'),]
dat3<-read.csv('../res_data/revision/violin/Itgav_sing_violinData_imp_f.csv', row.names=1)
dat3<-dat3[dat3$annotation %in% c('OsteoCAR', 'Spp1+Mac', 'ResMac', 'InfiltMac', 'Ly6ClowMac'),]
dat4<-read.csv('../res_data/revision/violin/Itgav_doub_violinData_imp_f.csv', row.names=1)
dat4<-dat4[dat4$annotation %in% c('9_Fibro-OsteoCAR', '10_OsteoCAR-MonoPro', '6_Spp1+Mac-MonoPro', '3_Ery-MK'),]
###

p1=violin_p(dat1, gene="Spp1", min_ratio = 0, ymax=10)
p2=violin_p(dat2, gene="Spp1", min_ratio = 0, ymax=10)
p3=violin_p(dat3, gene="Itgav", min_ratio = 0, ymax=5)
p4=violin_p(dat4, gene="Itgav", min_ratio = 0, ymax=5)

p1+p3+p2+p4
ggsave(file="../figures/revision/Spp1_Itgav_sing_doub_violin_f.pdf", width=8, height=6, dpi=300)
