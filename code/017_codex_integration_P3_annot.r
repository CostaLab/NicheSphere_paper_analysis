library(Seurat)
library(sctransform)
library(harmony)
library(tidyverse)
library(hdf5r)
library(patchwork)

library(cowplot)
library(RColorBrewer)
library(viridis)
library(Matrix)

options(future.globals.maxSize = 5000 * 1024^2)

## Folders
fig_dir <- "../figures/revision/FS3/"
out_data_dir <- "../data/CODEX/intermediate_steps/"

## resolution to use
chosen_res='SCT_snn_res.0.3'

## cluster color palette
cluster_palette="Paired"

#######################

cur_seurat<-readRDS(paste0(out_data_dir, 'seurat_object_harmony_clusts_filt_bg_id.rds'))
cur_seurat@meta.data$seurat_clusters<-cur_seurat@meta.data[,chosen_res]	
Idents(object = cur_seurat) <- "seurat_clusters"
# Classic palette BuPu, with 4 colors
newpal <- brewer.pal(12, cluster_palette) 
# Add more colors to this palette :
newpal <- colorRampPalette(newpal)(20)

#### Annotation
#other_clusters<-rownames(cur_seurat@meta.data)[which(!cur_seurat@meta.data$seurat_clusters %in% c(6,13))]
#cur_seurat_sub <- subset(x = cur_seurat, cells = other_clusters)
cur_seurat@meta.data$clusters_named<-0
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==0]<-'Immune_cells'
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==1]<-'Osteoblasts'
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==2]<-'Sinusoidal_endothelial_cells'
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==3]<-'Megakaryocytes'
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==4]<-'Erythroid_lineage'
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==5]<-'Myofibroblasts'
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==6]<-'MSCs'
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==7]<-'Vascular_endothelial_cells'
#cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==8]<-'Perivascular_pericytes'

cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==0]<-'CD45+_hematopoietic_cells'
cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==1]<-'Perivascular_myofibroblasts'
cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==2]<-'Sinusoidal_endothelial_cells'
cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==3]<-'Megakaryocytes'
cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==4]<-'Erythroid_lineage'
cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==5]<-'Myofibroblasts'
cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==6]<-'Activated_MSCs_CAR_cells'
cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==7]<-'Vascular_endothelial_cells'
cur_seurat@meta.data$clusters_named[cur_seurat@meta.data$seurat_clusters==8]<-'Pericytes'

Idents(object = cur_seurat) <- "clusters_named"

#clustMarkers <- FindAllMarkers(cur_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#write.csv(clustMarkers, paste0('../res_data/revision/CODEX/clusterMarkers_filt_bg.csv'))

#pdf(paste0(fig_dir, 'posClustMarkers_filt_bg.pdf'), height=4, width=9)
#DotPlot(cur_seurat, features = unique(as.character(unlist(tapply(clustMarkers$gene, clustMarkers$cluster, function(x){head(x, 2)}))))) + RotatedAxis()
#dev.off()

#pdf(paste0(fig_dir, 'allMarkers_bg.pdf'), height=4, width=9)
#DotPlot(cur_seurat, features = rownames(cur_seurat)) + RotatedAxis()
#dev.off()

## Spp1+ clusters
qthr=0.9
cur_seurat$Osteopontin_expr <- FetchData(cur_seurat, vars = "Osteopontin.166")[, 1]
#thr1=quantile(cur_seurat[,'Osteopontin.166'],na.rm = T,probs = qthr)
#thr1 <- quantile(cur_seurat$Osteopontin_expr, na.rm = TRUE, probs = 0.9)
thr1=0.5376757654610957
# 1. Get the cell names (barcodes) that meet the expression threshold
#high_expr_cells <- WhichCells(cur_seurat, expression = 'Osteopontin.166' > thr1)
high_expr_cells<-names(cur_seurat$Osteopontin_expr[cur_seurat$Osteopontin_expr>thr1])

# 2. Find which of those specific cells are also currently labeled "Immune_cells"
target_cells <- high_expr_cells[cur_seurat$clusters_named[high_expr_cells] == "CD45+_hematopoietic_cells"]
# 3. Update the metadata safely using the cell names as indices
cur_seurat$clusters_named[target_cells] <- "Spp1+_Immune"

# 2. Find which of those specific cells are also currently labeled "Immune_cells"
target_cells <- high_expr_cells[cur_seurat$clusters_named[high_expr_cells] == "Perivascular_myofibroblasts"]
# 3. Update the metadata safely using the cell names as indices
cur_seurat$clusters_named[target_cells] <- "Spp1+_Osteoblasts"

# 2. Find which of those specific cells are also currently labeled "Immune_cells"
target_cells <- high_expr_cells[cur_seurat$clusters_named[high_expr_cells] == "Sinusoidal_endothelial_cells"]
# 3. Update the metadata safely using the cell names as indices
cur_seurat$clusters_named[target_cells] <- "Spp1+_Sinusoidal_pericytes"

#cur_seurat@meta.data$clusters_named[(which(cur_seurat[,'Osteopontin.166']>thr1))&(cur_seurat@meta.data$clusters_named=='Immune_cells')]='Spp1+_Immune'

Idents(object = cur_seurat) <- "clusters_named"

#####

#cluster markers with Wilcoxon test
clustMarkers <- FindAllMarkers(cur_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(clustMarkers, paste0('../res_data/revision/CODEX/clusterMarkers_filt_bg_final.csv'))

#pdf(paste0(fig_dir, 'posClustMarkers_filt_bg_Spp1_final.pdf'), height=4, width=9)
#DotPlot(cur_seurat, features = unique(as.character(unlist(tapply(clustMarkers$gene, clustMarkers$cluster, function(x){head(x, 2)}))))) + RotatedAxis()
#dev.off()

pdf(paste0(fig_dir, 'FS3B.pdf'), height=4, width=9)
DotPlot(cur_seurat, features = rownames(cur_seurat)) + RotatedAxis()
dev.off()

#write.csv(cur_seurat@meta.data, paste0(out_data_dir,'seurat_object_harmony_clusts_filt_meta_bg.csv'))
#saveRDS(cur_seurat, paste0(out_data_dir,'seurat_object_harmony_clusts_filt_bg_id.rds'))
#counts_matrix <- GetAssayData(object = cur_seurat, assay='RNA')
#writeMM(counts_matrix, paste0(out_data_dir,'seurat_object_harmony_clusts_filt_bg_id_counts.mtx'))
#write.csv(Features(cur_seurat), paste0(out_data_dir,'seurat_object_harmony_clusts_filt_bg_id_mk_names.csv'))

print('END!!!')
