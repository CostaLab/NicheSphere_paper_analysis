## conda activate r_trajectory
require(Seurat)
require(rstatix)
require(tidyverse)
get_aggr <- function(obj, gene.set){
  gene.set = intersect(gene.set, rownames(obj))
  raw_counts <- GetAssayData(obj, slot="counts", assay="RNA")
  gene.set.val <- log10((Matrix::colSums(raw_counts[gene.set, ])/Matrix::colSums(raw_counts))*10000 + 1)
  return(gene.set.val)
}

options(repr.plot.width=10,repr.plot.height=20)

### Data
#mesen.evthpo<-readRDS('../rds/scrna_phase_clustering_doublets.Rds')
mesen.evthpo<-readRDS('../data/PICseq/final_data/adata_doublets_final.rds')

### DB
library(msigdbr)
m_df = msigdbr(species = "Mus musculus", category = "C2")
naba <- m_df[grep("NABA", m_df$gs_name), ]


naba.l <- list()
for(i in unique(naba$gs_name)){
  sub <- naba[naba$gs_name == i, ]
  sub <- as.data.frame(sub)
  naba.l[[i]] <- as.character(sub$gene_symbol)
}

### get aggr. gene expression for each geneset
df.s <- data.frame(
  #Cluster = mesen.evthpo@meta.data$harmony_inte_clusters,
  Cluster = mesen.evthpo@meta.data$max_pair,
  condition = mesen.evthpo@meta.data$stage
)
for(i in names(naba.l)){
  df.s[, ncol(df.s) + 1] <- get_aggr(mesen.evthpo, naba.l[[i]])
  names(df.s)[ncol(df.s)] <- paste(i)
}
library(reshape2)
df.s <- melt(df.s)
df.s[df.s == -Inf] <- 0

### Test

signftbl <- df.s %>%
      mutate(condition=factor(condition)) %>%
      mutate(id=paste0(Cluster,":",variable)) %>%
      group_by(id) %>%
      wilcox_test(value~condition,ref.group = "TPO",alternative = 'greater', p.adjust.method = "fdr",
)%>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj") 

signftbl <- signftbl %>%
            column_to_rownames("id")

#### Ridge plots

library(ggridges)
for(i in c("NABA_ECM_REGULATORS", "NABA_ECM_REGULATORS")){
aux <- df.s[df.s$variable == i, ] 
aux$signf <- signftbl[paste0(aux$Cluster,":",aux$variable),]$p.adj.signif
aux$signf[aux$condition=="EV"] <- NA
aux[duplicated(paste0(aux$Cluster,aux$condition,aux$variable)),"signf"]<-NA

aux<-aux %>% mutate(Cluster = fct_reorder(.f = Cluster, .x = value, .fun = mean))

ggplot(aux,
       aes(x=value,
           y=Cluster,
           fill = factor(condition)
           )
       ) +
  geom_density_ridges(
    jittered_points=FALSE,
    scale = .95,
    rel_min_height = .01,
    quantile_lines = TRUE,
    quantiles = 2,alpha=0.5
  ) +
  scale_fill_manual(name="Stim", values = c("EV"="#b2b2b2", "TPO"="#f02a2a"))+
  scale_y_discrete(expand = c(.01, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "Expression") +
  scale_discrete_manual("point_color") +
  geom_text(aes(x=max(density(aux$value)$x),label=signf,color='text'))+
  scale_color_manual(values = c('text'="#000000"))+
  guides(fill = guide_legend(
  override.aes = list(
  color = NA, point_color = NA))
  ) +
 theme_ridges(center = TRUE) +
  xlab("Log10 of Average Expression") +ylab("") + ggtitle(i) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(filename = paste("../figures/revision/F1/F1G",".pdf"), device = "pdf",width = 9,height = 10)
}

####################

### Data
mesen.evthpo<-readRDS('../data/PICseq/intermediate_steps/rds/scrna_phase_comparing_singlet_final_noSMCs_Xist.Rds')


### DB
library(msigdbr)
m_df = msigdbr(species = "Mus musculus", category = "C2")
naba <- m_df[grep("NABA", m_df$gs_name), ]


naba.l <- list()
for(i in unique(naba$gs_name)){
  sub <- naba[naba$gs_name == i, ]
  sub <- as.data.frame(sub)
  naba.l[[i]] <- as.character(sub$gene_symbol)
}

### get aggr. gene expression for each geneset
df.s <- data.frame(
  Cluster = mesen.evthpo@meta.data$annotation,
  condition = mesen.evthpo@meta.data$stage
)
for(i in names(naba.l)){
  df.s[, ncol(df.s) + 1] <- get_aggr(mesen.evthpo, naba.l[[i]])
  names(df.s)[ncol(df.s)] <- paste(i)
}
library(reshape2)
df.s <- melt(df.s)
df.s[df.s == -Inf] <- 0

### Test

signftbl <- df.s %>%
      mutate(condition=factor(condition)) %>%
      mutate(id=paste0(Cluster,":",variable)) %>%
      group_by(id) %>%
      wilcox_test(value~condition,ref.group = "TPO",alternative = 'greater', p.adjust.method = "fdr",
)%>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj") 

signftbl <- signftbl %>%
            column_to_rownames("id")

#### Ridge plots

library(ggridges)
for(i in c("NABA_ECM_REGULATORS", "NABA_ECM_REGULATORS")){
aux <- df.s[df.s$variable == i, ] 
aux$signf <- signftbl[paste0(aux$Cluster,":",aux$variable),]$p.adj.signif
aux$signf[aux$condition=="EV"] <- NA
aux[duplicated(paste0(aux$Cluster,aux$condition,aux$variable)),"signf"]<-NA

aux<-aux %>% mutate(Cluster = fct_reorder(.f = Cluster, .x = value, .fun = mean))

ggplot(aux,
       aes(x=value,
           y=Cluster,
           fill = factor(condition)
           )
       ) +
  geom_density_ridges(
    jittered_points=FALSE,
    scale = .95,
    rel_min_height = .01,
    quantile_lines = TRUE,
    quantiles = 2,alpha=0.5
  ) +
  scale_fill_manual(name="Stim", values = c("EV"="#b2b2b2", "TPO"="#f02a2a"))+
  scale_y_discrete(expand = c(.01, 0)) +
  scale_x_continuous(expand = c(0, 0), name = "Expression") +
  scale_discrete_manual("point_color") +
  geom_text(aes(x=max(density(aux$value)$x),label=signf,color='text'))+
  scale_color_manual(values = c('text'="#000000"))+
  guides(fill = guide_legend(
  override.aes = list(
  color = NA, point_color = NA))
  ) +
 theme_ridges(center = TRUE) +
  xlab("Log10 of Average Expression") +ylab("") + ggtitle(i) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(filename = paste("../figures/revision/F1/F1D",".pdf"), device = "pdf",width = 9,height = 10)
}




